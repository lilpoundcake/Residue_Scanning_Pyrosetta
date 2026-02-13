# Documentation Index

**Complete guide to all documentation files in this project**

---

## Quick Navigation

### 🚀 Getting Started
1. **[README.md](README.md)** — Installation & basic usage (START HERE)
2. **[CSV_FORMAT_GUIDE.md](CSV_FORMAT_GUIDE.md)** — How to format CSV mutations (10 min read)
3. **[tests/TESTING_GUIDE.md](tests/TESTING_GUIDE.md)** — How to run and understand tests

### 👨‍💼 For Users
- **[README.md](README.md)** — Installation, modes, examples
- **[CSV_FORMAT_GUIDE.md](CSV_FORMAT_GUIDE.md)** — CSV format with real examples
- **[tests/TESTING_GUIDE.md](tests/TESTING_GUIDE.md)** — Testing and validation
- **[DISTRIBUTION_GUIDE.md](DISTRIBUTION_GUIDE.md)** — Installation options

### 👨‍💻 For Developers
- **[CLAUDE.md](CLAUDE.md)** — Architecture, modules, execution flow
- **[DISTRIBUTION_GUIDE.md](DISTRIBUTION_GUIDE.md)** — Building wheels and distributions
- **[PROJECT_STATUS.md](PROJECT_STATUS.md)** — Project structure and decisions
- **[COMPLETION_SUMMARY.md](COMPLETION_SUMMARY.md)** — What was built

### 📊 For Project Managers
- **[PROJECT_COMPLETION_REPORT.md](PROJECT_COMPLETION_REPORT.md)** — Executive summary
- **[FINAL_TEST_RESULTS.md](FINAL_TEST_RESULTS.md)** — Test metrics and performance
- **[DISTRIBUTION_GUIDE.md](DISTRIBUTION_GUIDE.md)** — Distribution details

---

## All Documentation Files

### Core Documentation (Required Reading)

#### 1. **README.md** [⭐ START HERE]
- **Purpose**: Main project documentation
- **Audience**: All users
- **Length**: 3 pages
- **Topics**:
  - Installation (4 methods)
  - Usage examples (FR, RS, CM modes)
  - Output files description
  - Residue selection logic
  - Citations
- **Time to read**: 5-10 minutes

#### 2. **CSV_FORMAT_GUIDE.md**
- **Purpose**: Complete CSV format documentation
- **Audience**: Users of Custom Mutations mode
- **Length**: 12 pages
- **Topics**:
  - CSV format specification
  - 4 types of entries with examples
  - Position spec vs mutation disambiguation
  - Insertion code handling (UPPERCASE)
  - Real 7S7I example
  - Validation scripts
  - Common mistakes and fixes
  - Troubleshooting
- **Time to read**: 10-15 minutes

#### 3. **tests/TESTING_GUIDE.md**
- **Purpose**: Comprehensive testing documentation
- **Audience**: Users who want to run or understand tests
- **Length**: 13 pages
- **Topics**:
  - Test suite overview (5 tests)
  - Quick start commands
  - Directory structure
  - Test 1: Interface Selection (KDTree validation)
  - Test 2: Residue Scanning (RS mode)
  - Test 3: Custom Mutations (CM mode)
  - Test 4: Double Mutations (DM mode)
  - Test 5: Output Verification
  - Understanding results & ΔΔG interpretation
  - Performance metrics
  - Troubleshooting
- **Time to read**: 15-20 minutes

#### 4. **tests/TEST_SUITE.md** (NEW)
- **Purpose**: Complete test suite documentation with all 5 tests
- **Audience**: Users, QA, developers
- **Length**: 15 pages
- **Topics**:
  - Test suite overview (5 tests)
  - Test 1: Interface Selection (KDTree validation)
  - Test 2: Residue Scanning (RS mode, 3 PDBs)
  - Test 3: Custom Mutations (CM mode)
  - Test 4: Double Mutations (DM mode with CSV support) ✨ NEW
  - Test 5: Output Verification
  - Running tests (individual and all)
  - Test data files
  - Performance notes and optimization
  - Quality metrics
  - Troubleshooting
- **Time to read**: 20-25 minutes

---

### Installation & Distribution

#### 5. **DISTRIBUTION_GUIDE.md**
- **Purpose**: How to distribute, build, and install the package
- **Audience**: Users, developers
- **Length**: 12 pages (includes build details)
- **Topics**:
  - Quick distribution summary
  - Installation for end users
  - Usage examples
  - Distribution channels (GitHub, PyPI, Docker, etc.)
  - Package contents checklist
  - Installation from different sources
  - Version management
  - Verification & testing
  - Support information
  - Build history & details
  - Build configuration
  - Building distributions
  - Quality assurance
  - Troubleshooting builds
- **Time to read**: 15 minutes

---

### Project Information

#### 6. **CLAUDE.md**
- **Purpose**: Developer guide and architecture documentation
- **Audience**: Developers, maintainers
- **Length**: 4 pages
- **Topics**:
  - Environment setup
  - Running the CLI
  - Custom Mutations mode details
  - Running tests
  - Linting configuration
  - Building wheels
  - Package architecture
  - Module responsibilities
  - Execution flow
  - Interface filter details
  - Mutation name format
  - Multiprocessing details
- **Time to read**: 15 minutes

#### 7. **PROJECT_STATUS.md**
- **Purpose**: Project summary, structure, and decisions
- **Audience**: Project stakeholders, developers
- **Length**: 8 pages
- **Topics**:
  - Project summary
  - Test results
  - Implementation details
  - Code quality metrics
  - Directory structure
  - Documentation coverage
  - Dependencies
  - Production readiness
  - Key findings
  - Learning resources
- **Time to read**: 10-15 minutes

#### 8. **PROJECT_COMPLETION_REPORT.md**
- **Purpose**: Executive summary of project completion
- **Audience**: Project managers, stakeholders
- **Length**: 10 pages
- **Topics**:
  - Executive summary
  - Test results summary
  - What was accomplished
  - Code quality summary
  - Documentation highlights
  - Key achievements
  - Production readiness
  - Next steps
- **Time to read**: 10 minutes

---

### Detailed Results

#### 9. **FINAL_TEST_RESULTS.md**
- **Purpose**: Comprehensive test results and metrics
- **Audience**: QA, developers, project managers
- **Length**: 12 pages
- **Topics**:
  - Test summary (all 4 tests passed)
  - Detailed test results
  - Interface selection performance (F1 metrics)
  - Residue scanning results (3 PDBs)
  - Custom mutations results
  - Output verification
  - Performance summary
  - Key observations
  - CSV outputs
  - Conclusion
- **Time to read**: 15 minutes

#### 10. **COMPLETION_SUMMARY.md**
- **Purpose**: User-friendly project overview
- **Audience**: End users, project stakeholders
- **Length**: 9 pages
- **Topics**:
  - What was accomplished
  - Code quality checks
  - Comprehensive testing
  - Bug fixes
  - Documentation created
  - Usage examples
  - Code quality metrics
  - Production readiness
- **Time to read**: 10 minutes

---

### Supporting Files

#### 11. **tests/README.md**
- **Purpose**: Test suite overview (quick reference)
- **Audience**: Users wanting quick test info
- **Length**: 1 page
- **Topics**:
  - Quick start
  - File listing
  - Test summary
  - Documentation reference
- **Time to read**: 2 minutes
- **Note**: Points to TEST_SUITE.md for detailed information

#### 12. **tests/TEST_SUITE.md** (NEW)
- **Purpose**: Complete test suite documentation
- **Audience**: Users, QA, developers
- **Length**: 15 pages
- **Topics**:
  - All 5 tests with detailed descriptions
  - Test data files
  - Running tests
  - Performance metrics
  - Troubleshooting
- **Time to read**: 20-25 minutes

#### 13. **DOCUMENTATION_INDEX.md** (This File)
- **Purpose**: Guide to all documentation
- **Audience**: Everyone
- **Length**: This file
- **Topics**:
  - Quick navigation
  - File-by-file guide
  - Reading recommendations
  - Index by topic
- **Time to read**: 5-10 minutes

---

## Reading Recommendations by Role

### 👤 End User (Non-Technical)
**Time commitment**: 20-30 minutes
1. README.md (installation guide)
2. CSV_FORMAT_GUIDE.md (if using custom mutations)
3. tests/TESTING_GUIDE.md (if running tests)

### 👨‍💼 User (Technical)
**Time commitment**: 30-45 minutes
1. README.md (complete)
2. CSV_FORMAT_GUIDE.md (complete)
3. tests/TESTING_GUIDE.md (complete)
4. DISTRIBUTION_GUIDE.md (for installation options)

### 👨‍💻 Developer / Maintainer
**Time commitment**: 1-2 hours
1. README.md (overview)
2. CLAUDE.md (architecture)
3. tests/TEST_SUITE.md (testing strategy)
4. DISTRIBUTION_GUIDE.md (releases & builds)
5. PROJECT_STATUS.md (project structure)
6. CODE REVIEW (residue_scanning/*.py files)

### 📊 Project Manager
**Time commitment**: 30-45 minutes
1. PROJECT_COMPLETION_REPORT.md (executive summary)
2. FINAL_TEST_RESULTS.md (quality metrics)
3. PROJECT_STATUS.md (status & checklist)
4. BUILD_INFO.md (distribution)

### 🏢 DevOps / Release Manager
**Time commitment**: 30 minutes
1. DISTRIBUTION_GUIDE.md (how to build and distribute)
2. CLAUDE.md (environment setup)
3. tests/TEST_SUITE.md (testing before release)

---

## Documentation by Topic

### Installation & Setup
- README.md — Installation (4 methods)
- DISTRIBUTION_GUIDE.md — Distribution, installation options, & building from source
- CLAUDE.md — Development environment

### Usage & Examples
- README.md — Basic usage examples
- CSV_FORMAT_GUIDE.md — CSV format and examples
- CLAUDE.md — Running the CLI

### Testing & Validation
- tests/TEST_SUITE.md — Test suite (5 tests) ✨ NEW
- tests/TESTING_GUIDE.md — Testing guide
- FINAL_TEST_RESULTS.md — Test results & metrics
- PROJECT_STATUS.md — Test status

### CSV Format (Critical!)
- CSV_FORMAT_GUIDE.md — Complete reference with ambiguity resolution (⭐ most important)
- tests/TESTING_GUIDE.md — CSV in context of tests
- README.md — Quick reference

### Architecture & Design
- CLAUDE.md — Package architecture
- PROJECT_STATUS.md — Project structure
- README.md — Residue selection logic

### Code Quality
- PROJECT_COMPLETION_REPORT.md — Quality summary
- FINAL_TEST_RESULTS.md — Test results
- BUILD_INFO.md — Build requirements

### Distribution & Deployment
- DISTRIBUTION_GUIDE.md — Distribution options, building wheels, & installation methods
- README.md — Quick installation reference

---

## Key Files Location

```
Project Root/
├── README.md                        ⭐ START HERE
├── CSV_FORMAT_GUIDE.md              ⭐ CSV FORMAT (Critical!)
├── CLAUDE.md                        Developer guide
├── PROJECT_STATUS.md                Project summary
├── PROJECT_COMPLETION_REPORT.md     Executive summary
├── FINAL_TEST_RESULTS.md            Test metrics
├── COMPLETION_SUMMARY.md            What was built
├── DISTRIBUTION_GUIDE.md            Distribution & build (merged BUILD_INFO)
├── DOCUMENTATION_INDEX.md           This file
│
├── tests/
│   ├── TEST_SUITE.md                ✨ NEW: Complete test documentation
│   ├── TESTING_GUIDE.md             How to run tests
│   └── README.md                    Quick test reference
│
├── dist/
│   ├── residue_scanning-0.2.1-py3-none-any.whl    ✅ Ready
│   └── residue_scanning-0.2.1.tar.gz               ✅ Ready
│
├── residue_scanning/
│   ├── __init__.py
│   ├── cli.py
│   ├── core.py
│   └── preprocessing.py
│
└── env.yaml                         Conda environment
```

---

## Quick Links

### Most Important Files
1. **README.md** — Start here for everything
2. **CSV_FORMAT_GUIDE.md** — Must read if using custom mutations
3. **TESTING_GUIDE.md** — How tests work
4. **CLAUDE.md** — For developers

### Most Common Questions Answered In:
- "How do I install?" → README.md
- "How do I use it?" → README.md
- "How do I specify mutations?" → CSV_FORMAT_GUIDE.md
- "How do I run tests?" → tests/TEST_SUITE.md, tests/TESTING_GUIDE.md
- "What logging will I see?" → tests/TESTING_GUIDE.md (console output examples)
- "How is this built?" → DISTRIBUTION_GUIDE.md, CLAUDE.md
- "What was done?" → PROJECT_COMPLETION_REPORT.md
- "How do I distribute this?" → DISTRIBUTION_GUIDE.md
- "What are the test results?" → FINAL_TEST_RESULTS.md

---

## Document Statistics

| File | Pages | Words | Type | Audience |
|------|-------|-------|------|----------|
| README.md | 3 | ~1,500 | Guide | All |
| CSV_FORMAT_GUIDE.md | 12 | ~4,000 | Reference | Users |
| TESTING_GUIDE.md | 13 | ~5,000 | Guide | Users/QA |
| CLAUDE.md | 4 | ~2,000 | Guide | Developers |
| PROJECT_STATUS.md | 8 | ~4,000 | Summary | All |
| DISTRIBUTION_GUIDE.md | 8 | ~3,500 | Guide | All |
| BUILD_INFO.md | 7 | ~3,000 | Reference | Developers |
| PROJECT_COMPLETION_REPORT.md | 10 | ~4,500 | Summary | Managers |
| FINAL_TEST_RESULTS.md | 12 | ~5,000 | Report | QA/Managers |
| COMPLETION_SUMMARY.md | 9 | ~3,500 | Summary | Users |
| **TOTAL** | **86 pages** | **~35,500 words** | | |

---

## Navigation Tips

### Finding Information
1. **Start with README.md** for basic info
2. **Use Ctrl+F** to search specific files
3. **Check this index** for file topics
4. **Cross-references** at bottom of files link to related docs

### Recommended Reading Order
**For New Users**:
1. README.md (5 min)
2. CSV_FORMAT_GUIDE.md (10 min)
3. TESTING_GUIDE.md (15 min)

**For Developers**:
1. CLAUDE.md (15 min)
2. README.md (5 min)
3. tests/TESTING_GUIDE.md (15 min)
4. BUILD_INFO.md (10 min)

---

## Documentation Quality

- ✅ **Comprehensive**: 90+ pages covering all aspects
- ✅ **Well-Indexed**: This file + cross-references
- ✅ **Multiple Formats**: Quick guides + detailed references
- ✅ **Code Examples**: Real examples in every guide
- ✅ **Troubleshooting**: Dedicated sections in relevant files
- ✅ **Role-Based**: Different guides for different users
- ✅ **Consolidated**: Merged BUILD_INFO & AMBIGUITY_HANDLING into main docs
- ✅ **Up-to-Date**: Last updated 2026-02-13

---

## Getting Help

### If you don't know where to look:
1. Check this index
2. Read README.md
3. Search relevant file with Ctrl+F
4. Check CSV_FORMAT_GUIDE.md if CSV-related
5. Check TESTING_GUIDE.md if test-related

### If you find an error or gap:
- Please report via GitHub issues
- Include which file and what's missing
- Help us improve the documentation!

---

**Documentation Status**: ✅ COMPLETE & CONSOLIDATED
**Total Coverage**: 90+ pages across 10 main docs + 2 test docs
**Last Updated**: 2026-02-13
**Quality**: Production Ready
**Latest Changes**:
- ✅ Fixed DM CSV-based path wildtype pseudo-mutant ordering bug
- ✅ Added comprehensive wildtype pseudo-mutant documentation
- ✅ Updated all guides with ddG calculation details
- ✅ Updated CLAUDE.md with prepare_double_mut_df() ordering note
