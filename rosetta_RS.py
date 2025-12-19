import os
import argparse
import subprocess

import warnings
warnings.filterwarnings("ignore")

import multiprocessing as mp
import numpy as np
import pandas as pd

from random import randint
from itertools import repeat

#micromamba install python-blosc

from pyrosetta import *  # type: ignore
from pyrosetta.rosetta import *  # type: ignore

pyrosetta.init('-use_input_sc -flip_HNQ -ex1 -ex2 -relax:cartesian -ignore_unrecognized_res -mute all -ignore_zero_occupancy false')



parser = argparse.ArgumentParser()

parser.add_argument("-f", "--file", 
                    required=True,
                    help="Path to PDB file with complex")

parser.add_argument("-m", "--mode", 
                    required=True,
                    choices=['FR', 'FastRelax', 'RS', 'Residue_Scanning', 'DM', 'Double_Mut_Searching'],
                    help="Mode for a task. FastRelax - preparing input structure. Residue Scanning - calculate ddG for mutants. Double Mut Searching - experimental. Calculate ddG only for combinations of single muts agree with condition")

parser.add_argument("-r", "--receptor",
                    required=True,
                    default="HL",
                    help='Receptor chains name. If receptor consists several chains dont separate them. Default is HL')

parser.add_argument("-l", "--ligand",
                    required=True,
                    default="A",
                    help='Ligand chains name. If ligand consists several chains dont separate them. Default is HL. ddG_interface will be calculated for interaction between receptor and ligan chains')

parser.add_argument("--cpu", 
                    required=True,
                    type=int,
                    help="Count of CPU for a task")

parser.add_argument("-o", "--output",
                    required=False,
                    default='results',
                    help='Path to output files. Default: folder "results" will be created')

parser.add_argument("--not_relax",
                    required=False,
                    action='store_true',
                    help='FastRelax Minimization of every structure is included by default. Use it if u dont need this (or u wanna make fast mut analysis)')

parser.add_argument("--replics",
                    required=False,
                    default=5,
                    help='Max count of replics for FastRelax sampling. Default is 5')

parser.add_argument("--radius",
                    required=False,
                    default=8,
                    help='Radius around mutated residue which will be relaxed. Default is 8 angstrom')

parser.add_argument("--condition",
                    required=False,
                    default='ddG_complex < 0.5 and ddG_interface < 0.5',
                    help='Condition for single muts selection to generating double mut combinations. Default is "ddG_complex < 0.5 and ddG_interface < 0.5"')

parser.add_argument("--debug", 
                    required=False,
                    action='store_true',
                    help="Debug mode")

args = parser.parse_args()



input_file = os.path.abspath(args.file)

if not os.path.isabs(args.output):
    output_folder = os.path.join(os.path.dirname(input_file), args.output)
else:
    output_folder = args.output

if not os.path.exists(output_folder): os.mkdir(output_folder)
os.chdir(output_folder)

subprocess.run(f"sed -E 's/HI(D|E|P)/HIS/g' {input_file} | sed -E 's/ASH/ASP/g' | sed -E 's/GLH/GLU/g' > {output_folder}/for_pyrosetta.pdb", shell=True)

pyrosetta.toolbox.cleanATOM(f'{output_folder}/for_pyrosetta.pdb', f'{output_folder}/clear.pdb') # type: ignore
pdb = f'{output_folder}/clear.pdb'



#---------- Defaults -----------

RADIUS = args.radius # angstrom
REPLICS = args.replics
AA = 'ADEFGHIKLMNQRSTYVW' # excluded cysteine and proline

receptor_chains = args.receptor
ligand_chains = args.ligand

num_of_processes = args.cpu # count of cpus

debug = args.debug
not_relax = args.not_relax

three_to_one = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
}

one_to_three = {}

for key, value in three_to_one.items():
    one_to_three[value] = key

ddG_condition = args.condition



#---------- Functions ----------

def random_seed_pyrosetta(debug=debug):
    seed_int = randint(0, 9999999)
    pyrosetta.rosetta.basic.random.init_random_generators(seed_int, 'mt19937')
    if debug: print(f'Seed is {seed_int}')



def wildtyper(mut):
    mut = '_'.join([ f'{i[:-1]}{i[0]}' for i in mut.replace(',', '_').split('_') ])
    
    return mut



def score_fxn_checker(not_relax=not_relax):
    if not_relax:
        sfxn = pyrosetta.rosetta.core.scoring.ScoreFunctionFactory.create_score_function('ref2015')
    else:
        sfxn = pyrosetta.rosetta.core.scoring.ScoreFunctionFactory.create_score_function('ref2015_cart')
    
    return sfxn



def mutate(pose, mutation, debug=debug):
    if debug: print(f"Step mutate: {pose.pdb_info().name()} {mutation}")

    chain = mutation[1]
    from_aa = mutation[0]
    to_aa = mutation[-1]
    
    if mutation[-2].isalpha():
        insertion_code = mutation[-2]
        position = int(mutation[2:-2])
        rosetta_resid = pose.pdb_info().pdb2pose(chain, position, insertion_code)
    else:
        position = int(mutation[2:-1])
        rosetta_resid = pose.pdb_info().pdb2pose(chain, position)

    if from_aa != pose.sequence()[rosetta_resid-1]: 
        print(f'{mutation} -- Incorrect Mutation. AA in the Mutation Name is not Equal AA in the Pose')
        return KeyError, rosetta_resid

    mutate_residue = pyrosetta.rosetta.protocols.simple_moves.MutateResidue()

    residue_selector = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector()
    residue_selector.set_index(rosetta_resid)

    mutate_residue.set_selector(residue_selector)
    mutate_residue.set_res_name(one_to_three[to_aa])

    mutate_residue.apply(pose)

    return pose, rosetta_resid



def interface_analyzer(pose, receptor_chains=receptor_chains, ligand_chains=ligand_chains):
    iam = pyrosetta.rosetta.protocols.analysis.InterfaceAnalyzerMover()
    sfxn = score_fxn_checker()

    iam.set_interface(f"{''.join(ligand_chains)}_{''.join(receptor_chains)}")
    iam.set_compute_packstat(False)
    iam.set_pack_input(False)
    iam.set_pack_separated(True)
    iam.set_scorefunction(sfxn)
    iam.apply(pose)

    return iam.get_interface_dG()



def repack_and_minimize(pose, replica, mutation=False, 
                        RADIUS=RADIUS, adjacent_aa_range=1, debug=debug):
    
    if debug: print(f"Step repack_and_minimize: {pose.pdb_info().name()}")
    sfxn = score_fxn_checker()

    all_residue_selector = pyrosetta.rosetta.core.select.residue_selector.TrueResidueSelector()

    not_residue_selector = pyrosetta.rosetta.core.select.residue_selector.NotResidueSelector()
    not_residue_selector.set_residue_selector(all_residue_selector)

    tf = pyrosetta.rosetta.core.pack.task.TaskFactory()

    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.InitializeFromCommandline())
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.NoRepackDisulfides())

    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(
        pyrosetta.rosetta.core.pack.task.operation.PreventRepackingRLT(),
        not_residue_selector
        ))

    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(
        pyrosetta.rosetta.core.pack.task.operation.RestrictToRepackingRLT(),
        all_residue_selector
        ))
    
    mmf = pyrosetta.rosetta.core.select.movemap.MoveMapFactory()
    mmf.set_cartesian(True)
    mmf.all_bb(False)
    mmf.all_bondangles(True)
    mmf.all_bondlengths(True)
    mmf.all_chi(True)
    mmf.all_jumps(False)

    if mutation != False:
        residue_selector = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector()

        for i in mutation.split('_'): # type: ignore
            pose, rosetta_resid = mutate(pose=pose, mutation=i) # e.g. mutation = GH122AL
            residue_selector.append_index(rosetta_resid)
        
        adjacent_selector = pyrosetta.rosetta.core.select.residue_selector.PrimarySequenceNeighborhoodSelector()
        adjacent_selector.set_selector(residue_selector)
        adjacent_selector.set_lower_residues(adjacent_aa_range)
        adjacent_selector.set_upper_residues(adjacent_aa_range)

        mmf.add_bb_action(pyrosetta.rosetta.core.select.movemap.mm_enable, adjacent_selector)

        neighbor_selector = pyrosetta.rosetta.core.select.residue_selector.NeighborhoodResidueSelector()
        neighbor_selector.set_distance(RADIUS)
        neighbor_selector.set_focus_selector(residue_selector)
        neighbor_selector.set_include_focus_in_subset(True)

        tf.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(
            pyrosetta.rosetta.core.pack.task.operation.PreventRepackingRLT(),
            neighbor_selector, True
            ))

    if not_relax or debug:
        pack_mover = pyrosetta.rosetta.protocols.minimization_packing.PackRotamersMover()

        pack_mover.score_function(sfxn)
        pack_mover.task_factory(tf)
        pack_mover.apply(pose)
    else:
        fr = pyrosetta.rosetta.protocols.relax.FastRelax()

        fr.set_scorefxn(sfxn)
        fr.set_task_factory(tf)
        fr.set_movemap_factory(mmf)
        fr.max_iter(200)
        fr.cartesian(True)
        fr.apply(pose)

    if mutation != False:
        pose.dump_pdb(f'{pose.pdb_info().name().replace(".pdb", "")}_{mutation}_{replica}.pdb') # type: ignore
    else:
        pose.dump_pdb(f'{pose.pdb_info().name().replace(".pdb", "")}_{replica}.pdb') # type: ignore

    dG_interface = interface_analyzer(pose=pose)
    return sfxn.score(pose), dG_interface



def mutate_pose_FR(original_pose, replica, mutation):
    pose = original_pose.clone()
    
    if debug: print(f"Step mutate_pose_FR: {pose.pdb_info().name()} {mutation}")

    sfxn = score_fxn_checker()

    sfxn.score(pose)
    random_seed_pyrosetta()
    score, dG_interface = repack_and_minimize(pose=pose, replica=replica, mutation=mutation)
    
    return score, dG_interface



def prepare_mut_df(pose):
    chains = {}

    for i in range(1, pose.num_chains() + 1):
        chain_start, chain_end = pose.chain_begin(i), pose.chain_end(i)
        chain_name = pose.pdb_info().chain(chain_start)
        chains[chain_name] = [i, chain_start, chain_end]

    if debug: print(chains)

    residue_selector = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector()

    for chain in ligand_chains:
        chain_start, chain_end = chains[chain][1], chains[chain][2]
        for position in range(chain_start, chain_end+1):
            residue_selector.append_index(position)

    neigbor_selector = pyrosetta.rosetta.core.select.residue_selector.NeighborhoodResidueSelector()
    neigbor_selector.set_focus_selector(residue_selector)
    neigbor_selector.set_include_focus_in_subset(False)

    receptor_interface = neigbor_selector.get_residues(pose)

    screen_mut = []

    for position in receptor_interface:
        from_aa = pose.residue(position).name1()
        chain = pose.pdb_info().pose2pdb(position).split()[1]
        pdb_position = pose.pdb_info().pose2pdb(position).split()[0]
        insertion_code = pose.pdb_info().icode(position).replace(' ', '')
        for to_aa in AA:
            mut_name = [f"{from_aa}{chain}{pdb_position}{insertion_code}{to_aa}"]
            mut_name.extend(['NaN' for _ in range(REPLICS*2)])
            screen_mut.append(mut_name)
            
    if debug: print(len(screen_mut))

    col_names = ["Name"]
    col_names.extend([f'Replica_{i+1}_Complex' for i in range(REPLICS)])
    col_names.extend([f'Replica_{i+1}_dG' for i in range(REPLICS)])

    df = pd.DataFrame(screen_mut, columns=col_names)

    return df



def FastRelax_replics(pose, replics, mutation, not_relax=not_relax):

    if debug: print(f"Step FastRelax_replics: {pose.pdb_info().name()}")

    score_mut = [ 'NaN' for _ in range(replics*2) ]

    for replica in range(1, replics+1):
        score_mut[replica-1], score_mut[replica-1+5] = mutate_pose_FR(original_pose=pose, replica=replica, mutation=mutation)
        if replica == 2:
            if not_relax:
                continue
            elif abs(score_mut[0] - score_mut[1]) < 1: # type: ignore
                break
    if debug: print(mutation, score_mut)
    
    return score_mut



def df_FastRelax(df, pose, replics=REPLICS):

    if debug: print(f"Step df_FastRelax: {pose.pdb_info().name()}")

    df[list(df.columns[1:])] = df.apply(lambda x : FastRelax_replics(pose=pose, replics=replics, mutation=x.Name), axis=1).to_list()
    
    return df



def ddG_calculation(x, df):

    wt_REU = df.loc[ df['Name'] == f'{wildtyper(x.Name)}' , ['REU_min'] ].iat[0,0]
    wt_dG_interface = df.loc[ df['Name'] == f'{wildtyper(x.Name)}' , ['Min_dG_Interface'] ].iat[0,0]
    
    return x.REU_min - wt_REU, x.Min_dG_Interface - wt_dG_interface



def df_ddG_postprocessing(df_concatenated, REPLICS=REPLICS, postfix=''):

    df1 = df_concatenated.apply(lambda x : x.replace('NaN', 0))
    df1.loc[:, 'REU_min'] = df1.apply(lambda x : min(x[1:REPLICS+1]), axis=1)
    df1.loc[:, 'Min_structure'] = df1.apply(lambda x : list(x[1:REPLICS+1]).index( min( x[1:REPLICS+1] )) + 1, axis=1)
    df1.loc[:, 'Min_dG_Interface'] = df1.apply(lambda x : x[f'Replica_{x.Min_structure}_dG'], axis=1)

    df2 = df1.loc[:, ['Name', 'REU_min', 'Min_dG_Interface']]
    df2.loc[:, 'Is_WT'] = df2.Name.apply(lambda x : x == wildtyper(x)) # type: ignore

    if not os.path.exists(f"all_REU_pdb{postfix}"): os.mkdir(f"all_REU_pdb{postfix}")
    os.system(f'mv {pose.pdb_info().name().replace(".pdb", "")}_* all_REU_pdb{postfix}')

    if not os.path.exists(f"min_REU_pdb{postfix}"): os.mkdir(f"min_REU_pdb{postfix}")
    df1.apply(lambda x : os.system(f'cp all_REU_pdb{postfix}/{pose.pdb_info().name().replace(".pdb", "")}_{x.Name}_{x.Min_structure}.pdb min_REU_pdb{postfix}'), axis=1)

    df2.loc[:, ['ddG_complex', 'ddG_interface']] = df2.apply(lambda x: ddG_calculation(x, df2), axis=1).to_list()
    df3 = df2.query('Is_WT != True').loc[:, ['Name', 'REU_min', 'ddG_complex', 'ddG_interface']].reset_index(drop=True)

    df3 = df3.sort_values("ddG_interface")
    df3.to_csv(f"Rosetta_ddG_mut{postfix}.csv", index=False)

    return df3



#------------- Structure Preparing ---------------

pose = pyrosetta.pose_from_pdb(pdb)
original_pose = pose.clone()

sfxn = score_fxn_checker()

sfxn.score(pose)
original_pose_score = sfxn.score(original_pose)


if args.mode in ['FR', "FastRelax", 'RS', "Residue_scanning", 'DM', 'Double_Mut_Searching']:

    num_of_processes = REPLICS
    relaxed_pose_scores = {}
    lowest_energy_pose = ''
    lowest_energy = 0

    with mp.Pool(num_of_processes) as pool:
        pool.starmap(mutate_pose_FR, [ (pose, i+1, False) for i in range(REPLICS) ]) # type: ignore

    for i in range(REPLICS):
        pose_name = f'{pdb.replace(".pdb", "")}_{i+1}.pdb'
        relaxed_pose = pyrosetta.pose_from_pdb(pose_name)
        relaxed_pose_scores[pose_name] = sfxn.score(relaxed_pose)

    for relaxed_pose, energy in relaxed_pose_scores.items():
        if energy < lowest_energy:
            lowest_energy = energy
            lowest_energy_pose = relaxed_pose

    print(f'\n{os.path.basename(lowest_energy_pose)} has the smallest REU')
    print(f'Original pose -- {round(sfxn.score(original_pose), 2)} REU')
    print(f'{os.path.basename(lowest_energy_pose)} pose -- {round(lowest_energy, 2)} REU')

    pose = pyrosetta.pose_from_pdb(lowest_energy_pose)
    original_pose = pose.clone()
    sfxn.score(pose)

    os.system(f'mkdir {os.path.dirname(pdb)}/FastRelax && mv {pdb.replace(".pdb", "")}_*.pdb {os.path.dirname(pdb)}/FastRelax/.')

    pose.dump_pdb('FastRelax.pdb')



#--------------- Residue Scanning ----------------

if args.mode in ['RS', "Residue_scanning", 'DM', 'Double_Mut_Searching']:

    pose = pyrosetta.pose_from_pdb('FastRelax.pdb')
    original_pose = pose.clone()

    print(f'\n{pose.pdb_info().name()} pose -- {round(sfxn.score(pose), 2)} REU')
    df = prepare_mut_df(pose)

    df = df.loc[ ( -( df.loc[:, 'Name'].str.startswith('P') + df.loc[:, 'Name'].str.startswith('C') )), :] # type: ignore
    print(f'\nThere are {len(df)} mutations and {len(df)*5} structures with 5 replics\n')

    if debug: df = df.head( len(AA) * 4 )

    num_of_processes = args.cpu

    if len(df) < num_of_processes: num_of_processes = len(df)

    with mp.Pool(num_of_processes) as pool: 
        df_splitted = np.array_split(df, num_of_processes) 
        df_concatenated = pd.concat(pool.starmap(df_FastRelax, zip(df_splitted, repeat(pose))))

    df_concatenated.to_csv("Rosetta_results_REU.csv", index=False)

    df3 = df_ddG_postprocessing(df_concatenated)



#------------- Double Mut Searching --------------

if args.mode in ['DM', 'Double_Mut_Searching']:

    if debug: 
        selected_positions = df3.Name.apply(lambda x : x[1:-1]).unique().tolist() # type: ignore
    else:
        selected_positions = df3.query(ddG_condition).Name.apply(lambda x : x[1:-1]).unique().tolist() # type: ignore

    neighbors_dict = {}

    for i in selected_positions:
        residue_selector = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector()
        neighbor_selector = pyrosetta.rosetta.core.select.residue_selector.NeighborhoodResidueSelector()

        chain = i[0]

        if i[-1].isalpha():
            position = int(i[1:-1])
            insertion_code = i[-1]
            pyrosetta_position = pose.pdb_info().pdb2pose(chain, position, insertion_code)
        else:
            position = int(i[1:])
            pyrosetta_position = pose.pdb_info().pdb2pose(chain, position)
        
        residue_selector.set_index(pyrosetta_position)
        neighbor_selector.set_focus_selector(residue_selector)
        neighbor_selector.set_include_focus_in_subset(False)

        for neigbor in neighbor_selector.selection_positions(pose):
            neigbor_position, neigbor_chain = pose.pdb_info().pose2pdb(neigbor).split()
            if neigbor_chain in receptor_chains:
                neighbors_dict[i] = neighbors_dict.get(i, [])
                neighbors_dict[i].append(f'{neigbor_chain}{neigbor_position}')

    updated_neighbors_dict = {}

    for key, values in neighbors_dict.items():
        for value in values:
            if value in neighbors_dict.keys():
                if key not in updated_neighbors_dict.get(value, ''):
                    updated_neighbors_dict[key] = updated_neighbors_dict.get(key, [])
                    updated_neighbors_dict[key].append(value)

    if debug: 
        df4 = df3.copy() # type: ignore
    else:
        df4 = df3.query(ddG_condition).copy() # type: ignore
    
    df4.loc[:, 'Position'] = df4.apply(lambda x : x.Name[1:-1], axis=1)
    df4.loc[:, 'Mutate_from'] = df4.apply(lambda x : x.Name[0], axis=1)
    df4.loc[:, 'Mutate_to'] = df4.apply(lambda x : x.Name[-1], axis=1)

    df5 = df4.groupby("Position").agg(lambda x : (x).tolist()).loc[:, ['Mutate_to', 'Mutate_from']].reset_index()

    df_line = []
    mut_combinations = []

    for mut1, values in updated_neighbors_dict.items():
        for mut2 in values:
            for mut1_to in df5.loc[ df5['Position'] == mut1, ['Mutate_to'] ].iat[0,0]: # type: ignore
                mut1_from = df5.loc[ df5['Position'] == mut1, ['Mutate_from'] ].iat[0,0][0] # type: ignore
                for mut2_to in df5.loc[ df5['Position'] == mut2, ['Mutate_to'] ].iat[0,0]: # type: ignore
                    mut2_from = df5.loc[ df5['Position'] == mut2, ['Mutate_from'] ].iat[0,0][0] # type: ignore
                    df_line = [f'{mut1_from}{mut1}{mut1_to}_{mut2_from}{mut2}{mut2_to}']
                    df_line.extend( [f'NaN' for i in range(REPLICS * 2)] )
                    mut_combinations.append(df_line)
            df_line = [f'{mut1_from}{mut1}{mut1_from}_{mut2_from}{mut2}{mut2_from}'] # type: ignore
            df_line.extend( [f'NaN' for i in range(REPLICS * 2)] )
            mut_combinations.append(df_line)

    print(f'There are {len(mut_combinations)} double mutations\n')

    col_names = ["Name"]
    col_names.extend([f'Replica_{i+1}_Complex' for i in range(REPLICS)])
    col_names.extend([f'Replica_{i+1}_dG' for i in range(REPLICS)])

    df6 = pd.DataFrame(mut_combinations, columns=col_names)

    num_of_processes = args.cpu

    if len(df6) < num_of_processes: num_of_processes = len(df6)

    with mp.Pool(num_of_processes) as pool: 
        df_splitted = np.array_split(df6, num_of_processes) 
        df_concatenated = pd.concat(pool.starmap(df_FastRelax, zip(df_splitted, repeat(pose))))

    df_concatenated.to_csv("Rosetta_results_REU_double_mut.csv", index=False)

    df_ddG_postprocessing(df_concatenated, postfix='_double')

print('\nDONE 🐁')