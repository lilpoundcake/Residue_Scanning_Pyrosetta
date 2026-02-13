"""PDB preprocessing: residue name normalization and structure cleaning."""

from __future__ import annotations

from pathlib import Path

# Map of non-standard / protonated residue names to canonical PDB names.
# Keys are checked against 4-char names first, then 3-char names.
NONSTANDARD_RESIDUES: dict[str, str] = {
    # Histidine protonation variants (AMBER, CHARMM, GROMACS)
    "HID": "HIS",
    "HIE": "HIS",
    "HIP": "HIS",
    "HSD": "HIS",
    "HSE": "HIS",
    "HSP": "HIS",
    "HISB": "HIS",
    "HISD": "HIS",
    "HISE": "HIS",
    "HISH": "HIS",
    # Protonated aspartate
    "ASPH": "ASP",
    "ASPP": "ASP",
    "ASH": "ASP",
    # Protonated glutamate
    "GLUH": "GLU",
    "GLUP": "GLU",
    "GLH": "GLU",
    # Protonated / modified lysine
    "LYSH": "LYS",
    "LYP": "LYS",
    "LYN": "LYS",
    # Deprotonated arginine
    "ARGN": "ARG",
    # Cysteine variants
    "CYSH": "CYS",
    "CYX": "CYS",
    "CYM": "CYS",
    # Deprotonated tyrosine
    "TYRM": "TYR",
    # Protonated glutamine / asparagine (rare)
    "GLNH": "GLN",
    "ASNH": "ASN",
}


def normalize_residue_names(pdb_path: Path) -> str:
    """Read a PDB file and replace non-standard residue names in ATOM/HETATM lines.

    Checks 4-character names (cols 18-21) before 3-character names (cols 18-20)
    so that names like ASPH are caught before ASP.
    """
    lines: list[str] = []
    for line in pdb_path.read_text().splitlines():
        if line.startswith(("ATOM", "HETATM")) and len(line) >= 21:
            resname4 = line[17:21].strip()
            resname3 = line[17:20].strip()

            if resname4 in NONSTANDARD_RESIDUES:
                canonical = NONSTANDARD_RESIDUES[resname4].ljust(4)
                line = line[:17] + canonical + line[21:]
            elif resname3 in NONSTANDARD_RESIDUES:
                canonical = NONSTANDARD_RESIDUES[resname3].ljust(3)
                line = line[:17] + canonical + line[20:]

        lines.append(line)

    return "\n".join(lines) + "\n"


def prepare_pdb(input_path: Path, output_dir: Path) -> Path:
    """Normalize residue names and strip non-ATOM records.

    Writes intermediate files into *output_dir* and returns the path to the
    cleaned PDB ready for PyRosetta.
    """
    import pyrosetta  # imported here so the module can be loaded before init()

    normalized_path = output_dir / "for_pyrosetta.pdb"
    clean_path = output_dir / "clear.pdb"

    normalized_path.write_text(normalize_residue_names(input_path))
    pyrosetta.toolbox.cleanATOM(str(normalized_path), str(clean_path))  # type: ignore[attr-defined]

    return clean_path
