# mappingff

`mappingff` is a molecular force-field parameter mapping tool for generating LAMMPS data files from a user-built segment parameter database.

It is designed for workflows where you already have parameterized simple molecular fragments or segments and want to transfer those parameters onto a **larger, more complex target molecule** with matching local chemical environments. Directly parameterizing large and architecturally involved molecules (block copolymers, dendrimers, supramolecular assemblies, etc.) is often a difficult work — `mappingff` automates the transfer by building an SQLite database from reference LAMMPS data files, assigning atom types to a target structure using graph-based molecular environment matching, looking up bonded terms, and writing a complete LAMMPS data file.

## What mappingff does

`mappingff` performs two main tasks:

1. **Build a segment parameter database** from reference molecular structures and their known LAMMPS data files.
2. **Parameterize a target molecule** by matching each atom and bonded term against that database, then writing a complete LAMMPS data file.

The reference LAMMPS files may come from external parameterization tools(e.g., **ligpargen** for OPLS-style parameters), literature values, or your own validated builds. **Please verify that the reference parameters are chemically reasonable and consistent with your intended force field** before building the database. The key requirement is that the target molecule's local chemical environments are represented in the segment database. — see [Segment coverage requirements](#segment-coverage-requirements) for details on what counts as sufficient coverage.

## What mappingff does not do

`mappingff` does not:

- Fit, optimize, or derive new force-field parameters.
- Run quantum-mechanical calculations.
- Validate the physical correctness of a force field.
- Have built-in standard parameter tables such as OPLS, AMBER, CHARMM, or GAFF for general molecules directly.
- Generate LAMMPS input scripts such as `pair_style`, `bond_style`, or `run` commands.
- Build production-ready periodic simulation boxes.

Treat the generated LAMMPS data file as a reproducible starting point for simulation setup. Inspect the logs, check all fallbacks and missing parameters, and validate the system before production simulations.

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

- **Segment-based parameter database**: Builds a reusable SQLite database from reference `LAMMPS data` files and associated molecular structures. The database stores atom types with four-level environment keys, fallback mappings, pair parameters, bonded parameters, and build metadata.
- **Rooted graph-based atom environment encoding**: computes four atom-environment keys (`hop0`, `hop1`, `hop2`, `hop3`) using canonicalized local molecular environments.
- **Rooted induced subgraph matching for hop1-hop3**: the center atom is marked as the root; all atoms and bonds within the selected graph radius are represented without using original RDKit atom indices in the serialized fingerprint.
- **Hierarchical atom-type fallback**: target atoms are matched by trying `hop3` first, then `hop2`, `hop1`, and `hop0`.
- **Bonded parameter transfer**: bonds, angles, dihedrals, and impropers are looked up using canonical combinations of participating atoms' `hop0` keys.
- **Duplicate-observation aggregation**: repeated observations of the same environment are reduced by averaging numeric coefficients during database finalization.

- **Optional total-charge adjustment**: When `--charge` redistributes charges to match a requested total system charge using a two-stage algorithm: first over atoms with multiple observed charge values (weighted by absolute charge, bounded by observed range), then evenly over all atoms for any residual.

- **LAMMPS data output**: writes topology, coordinates, masses, pair coefficients, bonded coefficients, and box bounds.

- **Flexible sample input**: supports paired `.mol`/`.pdb` + `LAMMPS data` segment directories and top-level standalone `LAMMPS data` samples.

## Installation

### Requirements

- Python `>=3.10`
- RDKit `>=2022.9.1`
- NumPy `>=1.21`

RDKit is usually most reliable from `conda-forge`.

```bash
git clone https://github.com/Caulfield23/mappingff.git
cd mappingff

conda create -n mappingff python=3.10 rdkit numpy -c conda-forge
conda activate mappingff

pip install -e .
```

## Quick start

### 1. Prepare a segment library

A typical project layout is:

```text
project/
├── samples/
│   ├── segment1/
│   │   ├── segment1.mol      # or segment1.pdb
│   │   └── segment1.data      # reference LAMMPS data with known parameters
│   ├── segment2/
│   │   ├── segment2.pdb
│   │   └── segment2.data
│   └── segment3.data          # sigle LAMMPS data file
└── target.mol                # or target.pdb
```

If standalone LAMMPS mode (single `.data` file) gives poor results, try pairing it with a `.mol` or `.pdb` topology file in the same subdirectory.

### 2. Build the parameter database

```bash
mappingff build samples/ -d parameters.db
```

This creates an SQLite database and writes `build.log` next to the database path.

To append new samples to an existing database

```bash
mappingff build samples/ -d parameters.db --append
```

### 3. Parameterize a target molecule

```bash
mappingff par target.mol -d parameters.db -o target.data -v
```

If `-o/--out` is omitted, the output path defaults to the target filename with a `.data`, `.lammps` or `.lmp` suffix. If `-d/--db` is omitted, `mappingff` uses the first `.db` file found in the target molecule's directory.

To request a specific total charge:

```bash
mappingff par target.mol -d parameters.db -o target.data -c 0.0 -v
```

Always inspect `par.log` before using the generated LAMMPS data file.

## Command reference

### `mappingff build`

```bash
mappingff build <samples_dir> [options]
```

| Option | Description | Default |
|---|---|---|
| `-d, --db PATH` | Output SQLite database path | `samples.db` |
| `-a, --append` | Append to an existing database instead of replacing it | disabled |

Notes:

- Paired samples are discovered from immediate subdirectories of `<samples_dir>`.
- Top-level `LAMMPS data` files directly inside `<samples_dir>` are treated as standalone samples.
- Existing database files are replaced unless `--append` is used.

### `mappingff par`

```bash
mappingff par <mol_file> [options]
```

| Option | Description | Default |
|---|---|---|
| `-d, --db PATH` | Parameter database path | first `.db` in target directory |
| `-o, --out PATH` | Output LAMMPS data file | `<mol_file>.data` |
| `-c, --charge FLOAT` | Requested total system charge | no adjustment |
| `-v, --verbose` | Print detailed progress and write verbose log | disabled |

## Input expectations

### Target structures

Target molecules should be provided as `.mol` or `.pdb` files. For charged, aromatic, conjugated, or otherwise chemically sensitive systems, prefer `.mol` files with explicit bond orders. PDB-derived bond orders depend on RDKit inference and may not be reliable for all systems.

## Segment coverage requirements

The segment database must cover the chemical environments present in the target molecule.

### Atom coverage

Each chemically distinct target atom should appear in a reference segment with a compatible `hop3` rooted environment whenever possible. Segment boundaries, end caps, aromaticity perception, protonation, hydrogen representation, and bond-order inference can all change environment keys.

### Bonded-term coverage

Bonded terms are matched by the `hop0` keys of participating atoms.

For a dihedral:

```text
A-B-C-D
```

matching depends on the `hop0` keys of A, B, C, and D. Because `hop0` includes first-neighbor signatures, the practical context can extend to:

```text
X-A-B-C-D-Y
```

where X is an outer neighbor of A and Y is an outer neighbor of D.

If a polymer chain is cut and capped too aggressively, A or D may receive different neighbor signatures from the target molecule. For segment junctions, include enough overlap so that each bonded term appears inside at least one reference segment with the same local connectivity as in the target.

## Charge adjustment

By default, transferred charges are written unchanged.

When `--charge` is supplied, `mappingff` adjusts the total system charge in two stages:

1. **Bounded weighted adjustment** over atoms whose matched database type has multiple observed charge values. The adjustment is weighted by `abs(current_charge)` and bounded by that atom type's observed charge range.
2. **Uniform residual adjustment** over all atoms if any charge difference remains.

This makes the final written charge approach the requested total, but it can perturb all charges in stage 2. Use this as a consistency operation, not as a substitute for validating the charge model.

## Output files

A successful parameterization produces:

```text
target.data          # LAMMPS data file
par.log    # matching and warning log
```

The LAMMPS data file contains:

- Header counts and type counts.
- Orthogonal box bounds.
- Masses.
- Pair coefficients.
- Bond, angle, dihedral, and improper coefficients.
- Atoms, bonds, angles, dihedrals, and impropers.

Missing bonded parameters are written as zero/default coefficients and reported in the log. Fix these by adding better-covered segments rather than treating the output as production ready.

## Recommended validation workflow

Before production use:

1. Confirm `no_match == 0` whenever possible.
2. Review all hop2/hop1/hop0 fallbacks.
3. Check warnings for missing bonded parameters.
4. Inspect total charge, especially when `--charge` was used.
5. Verify type counts and coefficient values.
6. Run a short LAMMPS minimization or short MD check.
7. Compare against reference systems, experiments, or higher-level calculations for important applications.

## Programmatic usage

```python
from pathlib import Path
from mappingff import build_db, parameterize

build_db(
    samples_dir=Path("samples"),
    db_path=Path("parameters.db"),
    append=False,
)

result = parameterize(
    topo_path=Path("target.mol"),
    db_path=Path("parameters.db"),
    out_path=Path("target.data"),
    total_charge=None,
)

print(result)
```

## Database overview

The SQLite database contains:

| Table | Description |
|---|---|
| `atom_types` | Atom type definitions with hop3/hop2/hop1/hop0 keys, parameters, source metadata, and stored hop0/hop3 graphs |
| `hop2_keymap` | hop2-level fallback mapping for inspection |
| `hop1_keymap` | hop1-level fallback mapping for inspection |
| `hop0_keymap` | hop0-level fallback mapping for inspection |
| `bond_params` | Bond force constants and equilibrium distances |
| `angle_params` | Angle force constants and equilibrium angles |
| `dihedral_params` | OPLS dihedral coefficient arrays |
| `improper_params` | Improper coefficients |
| `meta` | Database metadata |

## License

MIT License.
