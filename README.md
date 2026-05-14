# mappingff

`mappingff` is a molecular force-field parameter mapping tool for generating LAMMPS data files from a user-built segment parameter database.

It is designed for workflows where you already have parameterized simple molecular fragments or segments and want to transfer those parameters onto a **larger, more complex target molecule** with matching local chemical environments. Directly parameterizing large and architecturally involved molecules (block copolymers, dendrimers, supramolecular assemblies, etc.) is often a difficult work — `mappingff` automates the transfer by building an SQLite database from reference LAMMPS data files, assigning atom types to a target structure using graph-based molecular environment matching, looking up bonded terms, and writing a complete LAMMPS data file.

## What mappingff does

`mappingff` performs two main tasks:

1. **Build a parameter database** from reference segment structures and their known LAMMPS data files.
2. **Parameterize a target molecule** by matching each atom and bonded term against that database, then writing a complete LAMMPS data file.

You can obtain the reference LAMMPS data files for small molecules from other general parameterization tools (e.g., **ligpargen** for OPLS-style parameters), literature-reported values, or your own builds. **Please verify that the reference parameters are chemically reasonable and consistent with your intended force field** before building the database. The key requirement is that each atom in the target molecule must have its local chemical environment represented in the reference segments — see [Segment coverage requirements](#segment-coverage-requirements) for details on what counts as sufficient coverage. `mappingff` then transfers those parameters by matching environments, enabling parameterization of arbitrarily large or architecturally complex molecules (block copolymers, dendrimers, supramolecular assemblies, etc.) as long as their atomic environments are represented in the database.

## What mappingff does NOT do

`mappingff` does **not** fit, optimize, or derive new force-field parameters. Its output is only as reliable as the chemical coverage and consistency of the segment/reference data used to build the database. The quality of the generated LAMMPS file depends directly on the coverage and consistency of your segment database.

`mappingff` does **not**:
- Access or use standard force-field parameter tables (OPLS, AMBER, CHARMM, etc.) directly
- Perform quantum mechanical calculations to derive parameters
- Fit parameters to experimental data or molecular dynamics simulations
- Validate force-field physics or compatibility
- Generate LAMMPS input scripts (pair_style, bond_style, etc.)
- Create production-ready periodic simulation boxes

**Important:** Treat generated parameters as a reproducible starting point for simulation setup, and validate them before scientific production runs. `mappingff` is currently marked as beta software.

## Supported LAMMPS styles

`mappingff` has been tested and is stable for **OPLS force field** simulations using the following LAMMPS styles:

| Section | Style | Functional form |
|---------|-------|----------------|
| `bond_style` | `harmonic` | `E = K·(r-r₀)²` |
| `angle_style` | `harmonic` | `E = K·(θ-θ₀)²` |
| `dihedral_style` | `opls` | `E = ½K₁[1+cos(φ)] + ½K₂[1−cos(2φ)] + ½K₃[1+cos(3φ)] + ½K₄[1−cos(4φ)]` |
| `improper_style` | `cvff` | `E = K·[1 + d·cos(n·φ)]` |

**Parameter averaging** (for the rare case of duplicate environments): compatible with all of the above styles. The averaging logic is purely mathematical and makes no assumption about physical meaning. In normal usage with a well-designed segment library, genuine duplicate observations are extremely rare because `hop0` atom typing for bonded terms is specific enough.

**Improper cvff coefficients `d` (+1 or -1) and `n` (periodicity, integer ∈ {0,1,2,3,4,6}) are treated as a categorical pair.** The most frequently observed `(d, n)` pair is used, and the force constant `K` is averaged only over records sharing that pair. **This improper handling is only applicable to cvff-style impropers. It is NOT compatible with CHARMM (which uses harmonic impropers), GROMOS, or other force fields with different improper functional forms.**

Be careful when using `mappingff` with other force fields or LAMMPS styles, This tool has not validated their stablility.

## Features

- **Segment-based parameter database**: Builds a reusable SQLite database from reference `.lmp` files and associated molecular structures. The database stores atom types with four-level environment keys, fallback mappings, pair parameters, bonded parameters, and build metadata.

- **Graph-based atom environment encoding**: Computes canonical atom-environment hashes at four levels (hop0, hop1, hop2, hop3). Each descriptor includes center atom properties (atomic number, formal charge, aromaticity, degree, hydrogen count, ring membership, ring count, hybridization), neighbor signatures (element:bond_type:formal_charge), and bond kinds. All descriptors are canonicalized before hashing to ensure consistent keys regardless of atom ordering.

- **Hierarchical atom-type fallback**: Assigns atom types by trying hop3 first (most specific), then falling back to hop2, hop1, and hop0 (coarsest) when an exact match is not available. This allows partial matching when the exact environment is not in the database, though such fallbacks should be reviewed.

- **Bonded parameter transfer**: Maps bond, angle, dihedral, and improper coefficients using canonical keys based on local atom hop0 environments. For impropers, the first atom in each record is treated as the center atom; the remaining three substituents are sorted to produce a canonical key — this matches the convention expected by LAMMPS `improper_style cvff`.

- **Duplicate-parameter averaging**: When the same environment appears multiple times in the segment library, mappingff stores all observations and computes averaged parameters during database finalization. Bonded coefficient arrays (K, r0, theta0, c1–c4) and pair parameters (sigma, epsilon) are highly stable and rarely produce genuine duplicates — averaging is effectively a safety net for these. **Charge is the exception**: partial charges can vary across different parameter sources or protonation states, so the averaged charge may differ noticeably from any single observation. Inspect your segment library for charge consistency and avoid mixing segments with large charge discrepancies for the same environment.**

- **Optional total-charge adjustment**: When `--charge` is supplied, redistributes partial charges to match a requested total system charge using a two-stage algorithm: first over atoms with multiple observed charge values (weighted by absolute charge, bounded by observed range), then evenly over all atoms for any residual.

- **LAMMPS data output**: Writes masses, pair coefficients, bonded coefficients, topology records, coordinates, and simulation box bounds in standard LAMMPS data-file format. The generated file contains all sections needed for a LAMMPS simulation.

- **Flexible sample input**: Supports paired `.mol`/`.pdb` + `.lmp` segment directories (recommended), plus top-level standalone `.lmp` samples for cases where only LAMMPS data is available.

## Installation

### Requirements

- Python `>=3.10`
- RDKit `>=2022.9.1`
- NumPy `>=1.21`

### Recommended installation

RDKit is most reliably installed from `conda-forge`:

```bash
git clone https://github.com/Caulfield23/mappingff.git
cd mappingff

conda create -n mappingff python=3.10 rdkit numpy -c conda-forge
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
│   └── segment3.lmp          # sigle LAMMPS data file
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

To adjust the total charge of the generated system (e.g., for a neutral polymer):

```bash
mappingff parameterize target.mol -d parameters.db -o target.lmp --charge 0.0
```

Always inspect `parameterize.log` before using the generated data file in production simulations.

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

Notes:

- Paired samples are discovered from immediate subdirectories of `<samples_dir>`.
- Top-level `.lmp` files directly inside `<samples_dir>` are treated as standalone samples (not processed if inside subdirectories).
- Existing database files are replaced unless `--append` is used.

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
| `-v, --verbose` | Print detailed matching information and write verbose log | disabled |

## Input file expectations

### Target structures

Target molecules should be provided as `.mol` or `.pdb` files. For chemically complex, charged, or unusual systems, prefer `.mol` files with explicit bond orders — RDKit bond-order inference from geometry may not be reliable for such cases.

## Charge adjustment

By default, `mappingff` preserves the charges transferred from the database without modification.

When `--charge` is supplied, the tool adjusts atom charges in two stages to reach the requested total:

**Stage 1 - Bounded weighted adjustment**: The charge difference (`target_charge - current_total_charge`) is distributed over atoms whose matched database atom type has more than one observed charge value (`charge_list` length > 1). The adjustment for each eligible atom is weighted by `abs(current_charge)` and bounded by the minimum and maximum observed charges for that atom type. Atoms that hit their bounds are fixed, and the remaining residual is redistributed to still-active atoms.

**Stage 2 - Uniform residual adjustment**: Any remaining charge difference is distributed evenly across all atoms in the system.

The log reports the charge before adjustment, after stage 1, after stage 2, and the requested target charge. This two-stage method guarantees the final written total approaches the requested target, but it can slightly perturb all charges in stage 2. Use charge adjustment as a final consistency operation, not as a substitute for validating the underlying charge model.

Example:

```bash
mappingff parameterize target.mol -d parameters.db -o target.lmp --charge 0.0
```

## Segment coverage requirements

The segment database must cover the chemical environments present in the target molecule. This is the most important step in the workflow and the most common source of quality issues.

### Atom coverage

Each chemically distinct target atom should appear in a reference segment with a compatible `hop3` environment whenever possible. (See the [Atom environment keys](#atom-environment-keys) section for what `hop3` includes.) If the exact `hop3` environment is not present, the resolver falls back to `hop2`, `hop1`, or `hop0`, but these are progressively less specific and may assign incorrect parameters for atoms in unusual environments.

### Bonded-term coverage

Bonded terms require the `hop0` environments of all participating atoms to be available in the database (see the [Atom environment keys](#atom-environment-keys) section for what `hop0` includes). For a dihedral `A-B-C-D`, matching depends on the `hop0` keys of A, B, C, and D.

Because `hop0` includes neighbor signatures (element:bond_type:formal_charge), the useful local context extends beyond the four dihedral atoms:

```text
X-A-B-C-D-Y
```

where X is an outer neighbor of A and Y is an outer neighbor of D. If a polymer chain is cut across this region and capped too aggressively, the `hop0` environments of A or D can change, and the dihedral may not match or may match incorrectly.

### Boundary design rule

When a bonded environment is only represented across a segment boundary, include enough overlapping atoms so that each bonded term appears inside at least one reference segment with the same local connectivity as in the target. For segment junctions, prefer fragments that overlap by at least three atoms around the cut and preserve the immediate neighbors required for dihedral lookup. Simple end-capping (e.g., adding hydrogens) can change the neighbor signatures and invalidate the hop0 environment.

## Output files

A successful parameterization produces:

```text
target.lmp          # LAMMPS data file
parameterize.log    # detailed matching log
```

The LAMMPS data file contains:
- Header with atom, bond, angle, dihedral, and improper counts and type counts
- Orthogonal box bounds generated from target coordinates with 5 Å padding in each direction (this is a convenience bounding box, not an equilibrated simulation box)
- `Masses` and `Pair Coeffs`
- `Bond Coeffs`, `Angle Coeffs`, `Dihedral Coeffs`, and `Improper Coeffs`
- `Atoms`, `Bonds`, `Angles`, `Dihedrals`, and `Impropers` sections

The generated file is a LAMMPS **data** file. 

### Missing parameters

When parameters are missing:
- Unmatched atoms receive an `unknown` atom type with zero mass, zero pair coefficients, and zero charge.
- Missing bond and angle parameters are written as zero coefficients.
- Missing dihedral parameters are written as four zero coefficients.
- Missing improper parameters are written with the default zeroed improper form used internally.
- Warnings are written to `parameterize.log`.

`mappingff` always produces a complete LAMMPS data file and writes all bonded terms, even when some parameters are missing — those terms are written with zero/default coefficients. Before using the output in production, inspect `parameterize.log` and confirm `no_match == 0`. A file with unresolved parameters is a diagnostic intermediate, not a production-ready input.

## Interpreting matching quality

The CLI prints a summary at the end of parameterization:

```text
Parameterize complete: 3313 atoms, 3486 bonds, 6245 angles, 9652 dihedrals, 1014 impropers, 73 types, hop3=2177, hop2=889, hop1=239, hop0=2, no_match=6
```

Use this summary and `parameterize.log` as a quality-control report:

| Statistic | What it means | Quality bar |
|---|---|---|
| `hop3=...` | Exact high-specificity environment matches | Higher is better; aim for most atoms |
| `hop2=...` | Second-level fallback; acceptable for uncommon environments | Review the specific atoms in log |
| `hop1=...` | First-level fallback; indicates missing hop2 coverage | Investigate and add segments if frequent |
| `hop0=...` | Coarse fallback; indicates significant environment gaps | Add segments covering these atoms |
| `no_match=...` | Atoms with no database representation at any level | Must be zero for production use; add segments |
| Missing bonded warnings | Bonded terms not found in database | Must be resolved; add reference segments |

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
# {'atoms': 1552, 'bonds': 1601, 'angles': 2900, 'dihedrals': 4291,
#  'impropers': 350, 'unique_types': 45, 'hop3_matches': 1184,
#  'hop2_matches': 268, 'hop1_matches': 100, 'hop0_matches': 0, 'no_match': 0}
```

The `parameterize` function returns a dictionary containing topology counts and matching statistics. Use this for automated quality checks or reporting.

## License

MIT License.