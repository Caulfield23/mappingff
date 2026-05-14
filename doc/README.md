# mappingff

`mappingff` is a molecular force-field parameter mapping tool for generating LAMMPS data files from a user-built segment parameter database.

It is designed for workflows where you already have parameterized molecular fragments or segments and want to transfer those parameters onto a target molecule with matching local chemical environments. The tool builds an SQLite database from reference LAMMPS data files, assigns atom types to a target structure using graph-based molecular environment matching, looks up bonded terms, and writes a complete LAMMPS data file.

`mappingff` does **not** fit or optimize new force-field parameters. Its output is only as reliable as the chemical coverage and consistency of the segment/reference data used to build the database.

## Features

- **Segment-based parameter database**: builds a reusable SQLite database from reference `.lmp` files and associated molecular structures.
- **Graph-based atom environment encoding**: computes canonical atom-environment hashes at four levels: `hop0`, `hop1`, `hop2`, and `hop3`.
- **Hierarchical atom-type fallback**: assigns atom types by trying `hop3` first, then falling back to `hop2`, `hop1`, and `hop0` when needed.
- **Bonded parameter transfer**: maps bond, angle, dihedral, and improper coefficients using canonical keys based on local atom environments.
- **Duplicate-parameter averaging**: merges repeated observations of the same atom or bonded environment when saving the database.
- **Standalone LAMMPS segment support**: can build database entries from root-level `.lmp` files without separate `.mol`/`.pdb` topology files by inferring elements from masses and bond orders from 3D geometry.
- **Optional total-charge adjustment**: redistributes partial charges to match a requested total system charge.
- **LAMMPS data output**: writes masses, pair coefficients, bonded coefficients, topology records, coordinates, and simulation box bounds in LAMMPS data-file format.

## Installation

### Requirements

- Python `>=3.10`
- RDKit `>=2022.9.1`
- NumPy `>=1.21`

RDKit is usually most reliable when installed from `conda-forge`:

```bash
git clone https://github.com/Caulfield23/mappingff.git
cd mappingff

conda create -n mappingff python=3.10 rdkit numpy -c condaforge
conda activate mappingff

python -m pip install -e .
```

## Quick start

### 1. Prepare a segment library

A typical project layout is:

```text
project/
├── samples/
│   ├── segment1/
│   │   ├── segment1.mol      # or segment1.pdb
│   │   └── segment1.lmp      # reference LAMMPS data with known parameters
│   ├── segment2/
│   │   ├── segment2.pdb
│   │   └── segment2.lmp
│   ├── segment3.lmp          # optional standalone LAMMPS-only segment       
└── target.mol                # or target.pdb
```

If standalone LAMMPS mode (single `.lmp` file) gives poor results, try pairing it with a `.mol` or `.pdb` topology file in the same subdirectory.

### 2. Build the parameter database

```bash
mappingff build-db samples/ -d parameters.db
```

This command creates an SQLite database containing atom types, fallback keys, pair parameters, bonded parameters, and metadata. It also writes a log file named `build-db.log` next to the database path.

To append new segments to an existing database instead of replacing it:

```bash
mappingff build-db samples/ -d parameters.db --append
```

### 3. Parameterize a target molecule

```bash
mappingff parameterize target.mol -d parameters.db -o target.lmp -v
```

If `-o/--out` is omitted, the output path defaults to the target filename with the `.lmp` suffix. If `-d/--db` is omitted, `mappingff` looks for the first `.db` file in the same directory as the target molecule.

To adjust the total charge of the generated system:

```bash
mappingff parameterize target.mol -d parameters.db -o target.lmp --charge 0.0
```

The command writes a LAMMPS data file and a `parameterize.log` file next to the output file. Always inspect the log before running production simulations.

## Command reference

### `mappingff build-db`

Build a parameter database from a segment directory.

```bash
mappingff build-db <samples_dir> [options]
```

Options:

| Option | Description | Default |
|---|---|---|
| `-d, --db PATH` | Output SQLite database path | `samples.db` |
| `-a, --append` | Append to an existing database instead of replacing it | disabled |

### `mappingff parameterize`

Assign parameters to a target molecule and write a LAMMPS data file.

```bash
mappingff parameterize <mol_file> [options]
```

Options:

| Option | Description | Default |
|---|---|---|
| `-d, --db PATH` | Parameter database path | first `.db` in target directory |
| `-o, --out PATH` | Output LAMMPS data file | `<mol_file>.lmp` |
| `-c, --charge FLOAT` | Target total charge; enables charge redistribution | no adjustment |
| `-v, --verbose` | Print detailed matching information | disabled |

## Input file expectations

### Target structures

Target molecules should be provided as `.mol` or `.pdb` files.

- `.mol` files are read with RDKit while preserving explicit hydrogens. `OBJ3D` blocks in V3000-style files are removed before parsing because some RDKit versions reject them.
- `.pdb` files are read with RDKit, and bond orders are inferred from geometry.

### Reference LAMMPS data files

Reference `.lmp` files are expected to use a LAMMPS data layout compatible with `atom_style full` atom records:

```text
atom_id molecule_tag atom_type charge x y z
```

The parser recognizes the following sections:

```text
Masses
Pair Coeffs
Bond Coeffs
Angle Coeffs
Dihedral Coeffs
Improper Coeffs
Atoms
Bonds
Angles
Dihedrals
Impropers
```

Coefficient layouts used by the current implementation are:

| Section | Expected fields after type id | Meaning |
|---|---:|---|
| `Pair Coeffs` | `epsilon sigma` | Lennard-Jones pair parameters |
| `Bond Coeffs` | `K r0` | harmonic bond coefficient and equilibrium distance |
| `Angle Coeffs` | `K theta0` | harmonic angle coefficient and equilibrium angle |
| `Dihedral Coeffs` | `c1 c2 c3 c4` | four-term dihedral coefficient record, typically used with OPLS-style inputs |
| `Improper Coeffs` | `K d n` | CVFF-style improper coefficient record |

The generated LAMMPS data file contains coefficients only; your LAMMPS input script must define compatible `pair_style`, `bond_style`, `angle_style`, `dihedral_style`, and `improper_style` settings.

## How parameter mapping works

### Atom environment keys

For every atom, `mappingff` computes four canonical graph descriptors:

| Level | Description | Role |
|---|---|---|
| `hop0` | center atom properties, neighbor signatures, and bond kinds | coarse local identity; also used for bonded parameter keys |
| `hop1` | `hop0` plus first-neighbor shell | fallback atom-type matching |
| `hop2` | `hop1` plus second-neighbor shell | fallback atom-type matching |
| `hop3` | `hop2` plus third-neighbor shell | most specific atom-type matching |

Each descriptor is canonicalized and encoded as a SHA-256 hash. During parameterization, atom types are resolved in this order:

```text
hop3 exact match → hop2 fallback → hop1 fallback → hop0 fallback → no_match
```

A `hop3` match is the most specific assignment. `hop2`, `hop1`, and `hop0` matches are progressively less specific and should be reviewed when they appear in the log.

### Bonded terms

Bond, angle, dihedral, and improper parameters are looked up using the resolved `hop0` environments of the atoms in each term.

Canonicalization handles common symmetries:

- bonds: `A-B` is equivalent to `B-A`
- angles: `A-B-C` is equivalent to `C-B-A`, with `B` fixed as the center
- dihedrals: `A-B-C-D` is equivalent to `D-C-B-A`
- impropers: the first atom is treated as the center; the remaining three atoms are sorted

If a bonded term cannot be matched, the tool still writes the output file, but the missing term receives zero/default coefficients and a warning is written to `parameterize.log`. Do not treat such an output as simulation-ready until those warnings are resolved.

### Duplicate observations

When the same environment appears multiple times in the segment library, `mappingff` stores all observations and computes averaged parameters when saving the database:

- atom `sigma`, `epsilon`, and `charge` are averaged over duplicate `hop3` environments;
- bond, angle, and dihedral coefficient arrays are averaged over duplicate keys;
- improper coefficients use the most common integer `(d, n)` pair and average the force coefficient over records with that pair.

## Segment coverage requirements

The segment database must cover the chemical environments present in the target molecule. This is especially important near segment boundaries, where end-capping can change local connectivity and therefore change the environment hash.

For a dihedral term `A-B-C-D`, the bonded-parameter lookup depends on the `hop0` environments of all four atoms. In practice, the surrounding context often needs to preserve the six-atom pattern:

```text
X-A-B-C-D-Y
```

where `X` is an outer neighbor of `A` and `Y` is an outer neighbor of `D`. If a polymer chain or large molecule is split across a boundary, include enough overlap or extended capped fragments so that the original local environments are represented in at least one segment.

Before building the database, check that:

1. all unique atom environments in the target occur in the segment set;
2. junction and end-cap environments do not unintentionally replace the target environments you want to transfer;
3. bonded terms spanning segment boundaries, especially dihedrals, are represented with sufficient local context;
4. `parameterize.log` reports acceptable `hop3/hop2/hop1/hop0/no_match` statistics.

## Charge adjustment

By default, `mappingff` preserves the charges transferred from the database.

When `--charge` is supplied, the tool adjusts atom charges in two stages:

1. It first redistributes the charge difference over atoms whose matched database atom type has multiple observed charge values. The adjustment is weighted by the absolute value of the current charge and bounded by the observed charge range for that atom type.
2. Any residual charge difference is then distributed evenly over all atoms.

Example:

```bash
mappingff parameterize target.mol -d parameters.db -o target.lmp --charge 0.0
```

The log reports the charge before adjustment, after the first stage, after the final stage, and the requested target charge.

## Output files

A successful parameterization produces:

```text
target.lmp
parameterize.log
```

The LAMMPS data file contains:

- header counts and type counts;
- orthogonal box bounds generated from the target coordinates with padding;
- `Masses` and `Pair Coeffs`;
- `Bond Coeffs`, `Angle Coeffs`, `Dihedral Coeffs`, and `Improper Coeffs`;
- `Atoms`, `Bonds`, `Angles`, `Dihedrals`, and `Impropers` sections.

The generated box is a convenience bounding box, not an equilibrated simulation box. Adjust it as needed for your simulation setup.

## Interpreting matching quality

The CLI prints a summary like:

```text
Parameterize complete: 3313 atoms, 3486 bonds, 6245 angles, 9652 dihedrals, 1014 impropers, 73 types, hop3=..., hop2=..., hop1=..., hop0=..., no_match=...
```

Use the summary and `parameterize.log` as a quality-control report:

- high `hop3` counts indicate exact high-specificity environment matches;
- `hop2`, `hop1`, and `hop0` matches indicate fallback assignments and should be chemically reviewed;
- any atom `no_match` means the atom was not represented by the database;
- any missing bonded parameter warning means a bond, angle, dihedral, or improper term was assigned zero/default coefficients.

For production use, aim for no atom or bonded-term `no_match` warnings unless you have explicitly reviewed and corrected the resulting data file.

## Python API

```python
from pathlib import Path
from mappingff import build_db, parameterize

build_db(
    samples_dir=Path("samples"),
    db_path=Path("parameters.db"),
    append=False,
)

stats = parameterize(
    topo_path=Path("target.mol"),
    db_path=Path("parameters.db"),
    out_path=Path("target.lmp"),
    total_charge=None,
)

print(stats)
```

The `parameterize` function returns a dictionary with atom/topology counts and matching statistics.

## Troubleshooting

| Symptom | Likely cause | What to check |
|---|---|---|
| `Database not found` | `--db` was omitted and no `.db` file exists next to the target | pass `-d parameters.db` explicitly |
| `Unsupported file type` | input suffix is not supported by `MolReader` | use `.mol`, `.pdb`, or supported internal `.lmp` workflows |
| `No element found for mass ...` | standalone `.lmp` atom mass is not in `mass.txt` within tolerance | edit `src/mappingff/mass.txt` or use paired `.mol/.pdb + .lmp` mode |
| many `hop1`, `hop0`, or `no_match` assignments | segment database does not cover target environments, or bond orders differ | add/extend segments, check capping, verify RDKit bond perception |
| missing bond/angle/dihedral/improper warnings | bonded environment key not present in database | add reference segments covering those local bonded contexts |
| unexpected formal charges or bond orders | RDKit inferred bond orders from geometry and assumed charge may not match your chemistry | inspect the parsed structure and prefer high-quality `.mol` files when possible |

## Additional documentation

See [`docs/USAGE.md`](docs/USAGE.md) for a deeper description of input contracts, database schema, environment matching, and validation workflow.

## License

MIT License.
