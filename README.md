# mappingff

Molecular force field parameterization pipeline for generating LAMMPS data files from segment databases.

## Features

- **Graph-based atom typing**: Classifies atoms using molecular subgraph hashing at four levels of environment detail (hop0 → hop1 → hop2 → hop3)
- **Hierarchical fallback**: Uses progressively more general environment matching when exact matches aren't available
- **Bonded parameter support**: Handles bonds, angles, dihedrals, and impropers with symmetric canonical keys
- **SQLite database**: Efficient storage and lookup of force field parameters
- **Charge adjustment**: Automatically adjusts total system charge by redistributing to non-hydrogen atoms
- **LAMMPS output**: Generates complete LAMMPS data files ready for simulation

## Installation

### Requirements

- Python 3.7+
- RDKit (`rdkit-pypi>=2022.9.5`)

### Install from source

```bash
pip install -e .
```

## Quick Start

### 1. Build a parameter database

Prepare a directory structure with sample molecules:

```
samples/
├── segment1/
│   ├── segment1.pdb       # Molecular structure
│   └── segment1.lammps.lmp # Reference LAMMPS parameters
├── segment2/
│   ├── segment2.pdb
│   └── segment2.lammps.lmp
└── ...
```

Build the database:

```bash
mappingff build-db samples/ --db-dir ./database
```

### 2. Parameterize a target molecule

```bash
mappingff parameterize target.mol --db-path ./database/db.pkl --out target_param.lmp
```

Or with charge adjustment:

```bash
mappingff parameterize target.mol --db-path ./database/db.pkl -c 0.0
```

## CLI Commands

### `build-db`

Build a parameter database from sample molecules.

```bash
mappingff build-db <samples_dir> [options]

Options:
  --db-dir PATH       Output directory for database (default: ./database)
  -v, --verbose       Print detailed progress
```

### `parameterize`

Parameterize a target molecule and generate LAMMPS file.

```bash
mappingff parameterize <mol_file> [options]

Options:
  --db-path PATH       Path to database file (default: ./database/db.pkl)
  --out PATH           Output LAMMPS file (default: <mol_file>_param.lmp)
  -c, --charge FLOAT  Target total charge (default: 0)
  -v, --verbose        Print detailed progress
```

## How It Works

### Graph-Based Environment Encoding

mappingff classifies atoms based on their molecular environment. The environment is described as a subgraph centered on the target atom, expanded to increasing radii:

| Level | Description |
|-------|-------------|
| hop0  | Center atom only |
| hop1  | Center + first neighbors |
| hop2  | + Second neighbors |
| hop3  | + Third neighbors (finest classification) |

Each environment is canonicalized and hashed to a SHA-256 key, enabling efficient database lookup.

### Four-Level Fallback

When parameterizing a target molecule, mappingff tries to match atoms at progressively more general levels:

1. **hop3 exact match** - Most specific, based on full 4-hop environment
2. **hop2 fallback** - Based on 3-hop environment
3. **hop1 fallback** - Based on 2-hop environment
4. **hop0 fallback** - Based on immediate neighbors only

This ensures every atom gets assigned parameters even if its exact environment wasn't in the training data.

### Canonical Keys for Bonded Parameters

Bonded parameters (bonds, angles, dihedrals) use canonical keys that account for molecular symmetry:

- **Bonds**: (keyA, keyB) with keyA ≤ keyB
- **Angles**: (keyA, keyB, keyC) with keyA ≤ keyC (center B is fixed)
- **Dihedrals**: (keyA, keyB, keyC, keyD) where (A,B,C,D) and (D,C,B,A) are equivalent

### Charge Adjustment

When `--charge` is specified, the system charge is adjusted by distributing the delta evenly across all non-hydrogen atoms:

```
adjustment_per_atom = (target_charge - current_charge) / num_non_hydrogen_atoms
```

## Programmatic Usage

mappingff can be imported and used as a Python library:

```python
from mappingff import buildDb, parameterize

# Build database
result = buildDb(
    samplesDir="./samples",
    dbPath="./database/db.pkl",
    verbose=True
)
print(f"Built database with {result['atoms_processed']} atoms")

# Parameterize molecule
result = parameterize(
    molPath="./target.mol",
    dbPath="./database/db.pkl",
    outPath="./target_param.lmp",
    total_charge=0.0,
    verbose=True
)
print(f"Assigned {result['unique_types']} atom types")
```

## Database Schema

The SQLite database contains:

| Table | Description |
|-------|-------------|
| `atom_types` | Atom type definitions with hop3/hop2/hop1/hop0 keys |
| `hop2_keymap` | hop2-level fallback mappings |
| `hop1_keymap` | hop1-level fallback mappings |
| `hop0_keymap` | hop0-level fallback mappings |
| `bond_params` | Bond force constants and equilibrium distances |
| `angle_params` | Angle force constants and equilibrium angles |
| `dihedral_params` | Dihedral OPLS coefficients |
| `improper_params` | Improper coefficients |
| `meta` | Database metadata (build time, sample count) |

## File Formats

### Input: Sample Molecules

Each sample directory should contain:
- **Structure file**: `.mol` (MDL Molfile) or `.pdb` (PDB format)
- **Reference LAMMPS file**: `.lmp` or `.lammps.lmp` with known parameters

### Output: LAMMPS Data Files

Standard LAMMPS data file format with:
- Header with atom/bond/angle/dihedral/improper counts
- Box dimensions
- Masses, Pair Coeffs, Bond Coeffs, Angle Coeffs, Dihedral Coeffs, Improper Coeffs
- Atoms, Bonds, Angles, Dihedrals, Impropers sections

## License

MIT License
