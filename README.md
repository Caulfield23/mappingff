# mappingff

Molecular force field parameterization pipeline for generating LAMMPS data files from segment databases.

## Features

- **Graph-based atom typing**: Classifies atoms using molecular subgraph hashing at four levels of environment detail (hop0 → hop1 → hop2 → hop3)
- **Hierarchical fallback**: Uses progressively more general environment matching when exact matches aren't available
- **Bonded parameter support**: Handles bonds, angles, dihedrals, and impropers with symmetric canonical keys
- **SQLite database**: Efficient storage and lookup of force field parameters
- **Two-step charge adjustment**: Automatically adjusts total system charge using a weighted two-step approach: first redistributes to atoms with multiple charge options, then evenly across all atoms if residual remains
- **LAMMPS output**: Generates complete LAMMPS data files ready for simulation
- **Flexible sample input**: Supports both `.mol`/`.pdb` + `.lmp` pairs, as well as standalone `.lmp` files (atom types inferred from mass via built-in element mass table)

## Installation

### Requirements

- Python 3.10+
- RDKit (`rdkit>=2022.9.1`, via conda-forge)

### Installation

```bash
git clone https://github.com/Caulfield23/mappingff.git
cd mappingff
conda install -c conda-forge rdkit
pip install -e .
```

## Quick Start

### 1. Build a parameter database

Prepare a directory structure with sample molecules. Each sample can be either a `.mol`/`.pdb` + `.lmp` pair, or a standalone `.lmp` file:

```
workdir/
├── samples/
│   ├── segment1/
│   │   ├── segment1.mol       # Molecular structure
│   │   └── segment1.lmp      # Reference LAMMPS parameters
│   ├── segment2/
│   │   ├── segment2.mol
│   │   └── segment2.lmp
│   ├── segment3.lmp          # Standalone LAMMPS file (no mol/pdb)
│   └── ...
└── target.mol               # Target molecule to parameterize
```

For standalone `.lmp` files, atom element types are inferred from the mass table in `mass.txt` (user-editable). Bond orders are automatically determined from 3D coordinates via RDKit.

Build the database:

```bash
mappingff build-db samples/ -d samples.db
```

### 2. Parameterize a target molecule

```bash
mappingff parameterize target.mol -d samples.db
```

If `-d` is not specified, mappingff automatically searches for a `.db` file in the same directory as the target molecule.

Or with charge adjustment:

```bash
mappingff parameterize target.mol -d samples.db -c 0.0
```

### Segment Coverage Requirements

Use this section as a preparation checklist before running `build-db`.

**What must be covered**
Each segment set must cover all chemical environments that appear in the target. In particular, environments at segment junctions must be preserved carefully, because end-capping can alter local connectivity and change the inferred environment.

**Boundary case: a dihedral term**
To determine bonded terms, the hop0 environment (including neighbor signatures) must be available for all involved atoms. For a dihedral A-B-C-D, matching requires a six-atom local context: X-A-B-C-D-Y.
- A-B-C-D are the dihedral atoms.
- X and Y are the outer neighbors of A and D.
- For X and Y, only the element type and the bond order to A or D are required.

If an original chain `R1-X-A-B-C-D-Y-R2` is split in the middle into `R1-X-A-B-end_capping` and `end_capping-C-D-Y-R2`, the two capped segments must be extended to at least `R1-X-A-B-C-D-Y-end_capping` and `end_capping-A-B-C-D-Y-R2` to cover the original environment with dihedral term `A-B-C-D` (if this dihedral environment is not fully represented inside a single segment).

Under this condition, adjacent segments must share at least three atoms with matching hop0 environments across the junction.


**Checklist before database build**
1. Confirm all unique target environments are present in your segment set.
2. Extend each boundary to include the dihedral context around the cut.
3. If a dihedral environment is represented only across a junction, verify that neighboring segments share at least three atoms with matching hop0 environments across that junction.

If these conditions are not met, some atoms or bonded terms can remain unmatched (`no_match`) during parameterization.

## CLI Commands

### `build-db`

Build a parameter database from sample molecules.

```bash
mappingff build-db <samples_dir> [options]

Options:
  -d, --db PATH        Output database file path (default: samples.db)
  -a, --append         Append to existing database instead of replacing
```

### `parameterize`

Parameterize a target molecule and generate LAMMPS file.

```bash
mappingff parameterize <mol_file> [options]

Options:
  -d, --db PATH        Path to database file (auto-searched if omitted)
  -o, --out PATH       Output LAMMPS file (default: <mol_file>.lmp)
  -c, --charge FLOAT  Target total charge for the system (default: no adjustment)
  -v, --verbose       Print detailed progress
```

## How It Works

### Graph-Based Environment Encoding

mappingff classifies atoms based on their molecular environment. The environment is described as a subgraph centered on the target atom, expanded to increasing radii:

| Level | Description |
|-------|-------------|
| hop0  | Center atom + neighbor signatures + bond kinds (coarse-grained) |
| hop1  | Center + first neighbors (first fallback) |
| hop2  | + Second neighbors (second fallback) |
| hop3  | + Third neighbors (finest classification) |

Each environment is canonicalized and hashed to a SHA-256 key, enabling efficient database lookup.

### Four-Level Fallback

When parameterizing a target molecule, mappingff tries to match atoms at progressively more general levels:

1. **hop3 matching** - Based on 3-hop environment (finest)
2. **hop2 fallback** - Based on 2-hop environment
3. **hop1 fallback** - Based on 1-hop environment
4. **hop0 fallback** - Based on center + immediate neighbors

mappingff makes a best-effort attempt to ensure every atom is assigned a parameter at the most specific level possible. Atoms that cannot be matched at any level are logged as `no_match` in the output.

### Bonded Parameters

Bonded parameters (bonds, angles, dihedrals, impropers) are classified using hop0 keys, which include neighbor signatures (`element:bond_type:formal_charge`) and bond kinds for precise matching.

### OPLS Improper Filtering

For impropers, mappingff follows OPLS force field conventions: only C or N center atoms with exactly 3 neighbors are considered for improper parameter assignment. This avoids generating spurious improper terms for non-trigonal centers. Other combinations are skipped silently during parameterization.

### Canonical Keys for Bonded Parameters

Bonded parameters (bonds, angles, dihedrals) use canonical keys that account for molecular symmetry:

- **Bonds**: (keyA, keyB) ordered lexicographically
- **Angles**: (keyA, keyB, keyC) with keyA ≤ keyC (center B is fixed, outer atoms swappable)
- **Dihedrals**: (keyA, keyB, keyC, keyD) where (A,B,C,D) and (D,C,B,A) are equivalent

### Charge Adjustment

When `--charge` is specified, the system charge is adjusted using a two-step weighted approach:

**Step 1 - Multi-entry types**: Atoms whose atom type has multiple charge options (`charge_list` length > 1) receive charge adjustments weighted by their current absolute charge, bounded by the available range for that type.

**Step 2 - All atoms**: Any remaining charge delta is distributed evenly across all atoms in the system.

When `--charge` is not specified, no charge adjustment is performed.

## Programmatic Usage

mappingff can be imported and used as a Python library:

```python
from mappingff import build_db, parameterize
from pathlib import Path

# Build database
build_db(
    samples_dir=Path("samples"),
    db_path=Path("database/samples.db"),
    append=False,
)
print("Database built successfully")

# Parameterize molecule
result = parameterize(
    topo_path=Path("target.mol"),
    db_path=Path("database/samples.db"),
    out_path=Path("target.lmp"),
    total_charge=None,
)
print(f"Assigned {result['unique_types']} atom types")
print(f"  hop3={result['hop3_matches']}, hop2={result['hop2_matches']}, "
      f"hop1={result['hop1_matches']}, hop0={result['hop0_matches']}, "
      f"no_match={result['no_match']}")
```

## Database Schema

The SQLite database contains:

| Table | Description |
|-------|-------------|
| `atom_types` | Atom type definitions with hop3/hop2/hop1/hop0 keys, mass, sigma, epsilon, charge |
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

Each sample can be provided in two ways:

**Paired mode**: A `.mol`/`.pdb` structure file and a `.lmp` reference file in the same subdirectory:
- **Structure file**: `.mol` (recommended) or `.pdb` (PDB format)
- **Reference LAMMPS file**: `.lmp` with known parameters

**Standalone mode**: A single `.lmp` file directly under `samples/` (no mol/pdb). Element types are inferred from atomic masses using the built-in mass table (`mass.txt`), and bond orders are automatically determined from 3D coordinates.

### Element Mass Table

The element-to-mass mapping is stored in `mass.txt` (editable):

```
H 1.008
C 12.011
N 14.007
O 15.999
...
```

When using standalone `.lmp` mode, each atom type's mass is matched against this table (with a tolerance of 0.1 by default) to determine the element symbol. Edit `mass.txt` to customize or extend the mapping.

### Output: LAMMPS Data Files

Standard LAMMPS data file format with:
- Header with atom/bond/angle/dihedral/improper counts
- Box dimensions
- Masses, Pair Coeffs, Bond Coeffs, Angle Coeffs, Dihedral Coeffs, Improper Coeffs
- Atoms, Bonds, Angles, Dihedrals, Impropers sections

## License

MIT License
