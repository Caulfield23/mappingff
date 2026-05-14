# mappingff technical overview

This document provides an architectural overview and implementation details for `mappingff`. It is intended for users who need to understand how the tool works internally, why a parameter assignment succeeds or fails, and how to diagnose quality issues.

## Source-level architecture

The package is organized around a small set of modules with focused responsibilities:

| Module | Responsibility |
|---|---|
| `mappingff.cli` | Command-line parsing, argument validation, and log setup. Handles both `build-db` and `parameterize` entry points. |
| `mappingff.workflow` | High-level `build_db()` and `parameterize()` functions. Orchestrates database construction and target parameterization workflows. |
| `mappingff.mol` | RDKit-based molecular parsing (`.mol`, `.pdb`, and internal `.lmp` conversion via `lmp2rdkitmol`), topology enumeration (atoms, bonds, angles, dihedrals, impropers), and hop-key computation. |
| `mappingff.encode` | Local graph construction, canonicalization, and SHA-256 key generation for hop0/hop1/hop2/hop3 environments. |
| `mappingff.db` | SQLite database schema, insertion, aggregation (averaging of duplicate environments), and lookup operations. |
| `mappingff.fallback` | Atom type resolution from hop3 to hop0 during target parameterization. Queries the database in fallback order and returns the first match. |
| `mappingff.lmp` | LAMMPS data file parsing (`parse_lammps`), writing (`generate_lammps`), and charge adjustment algorithm (`adjust_total_charge`). Also defines the `LammpsData` dataclass. |
| `mappingff.lmp2rdkitmol` | Standalone `.lmp` to RDKit molecule conversion. Maps atom-type masses to elements using `mass.txt`, reconstructs connectivity from the Bonds section, and uses RDKit to infer bond orders from 3D coordinates. |
| `mappingff.fallback` | Atom type resolution from hop3 to hop0. |

The public Python API is:

```python
from mappingff import build_db, parameterize
```

All other modules are private and subject to change.

## Database construction workflow

`build_db(samples_dir, db_path, append=False)` performs the following steps:

### Step 1: Sample collection

Two types of samples are collected:

1. **Paired samples** from immediate subdirectories of `samples_dir`.
   Each subdirectory is scanned for the first `.mol` or `.pdb` file and the first `.lmp` file found. To avoid ambiguity, keep only one topology file and one LAMMPS file per segment directory.

2. **Standalone LAMMPS samples** from `.lmp` files placed directly inside `samples_dir` (not in subdirectories).
   These are converted to RDKit molecules by mapping atom-type masses to elements using `mass.txt` and inferring bond orders from 3D coordinates.

### Step 2: Per-sample processing

For each sample:

1. Parse the LAMMPS reference data into a `LammpsData` object.
2. Parse or reconstruct the molecular graph with RDKit.
3. Enumerate atoms, bonds, angles, dihedrals, and impropers from the RDKit graph.
4. Compute hop3, hop2, hop1, and hop0 keys for every atom.
5. Insert atom parameters (hop keys, mass, sigma, epsilon, charge, source) into the `atom_types` table.
6. Insert bonded parameters (bond, angle, dihedral, improper) into their respective tables using canonicalized hop0 keys.
7. Build hop keymap entries for external validation.

### Step 3: Database finalization

During `save()`, accumulated duplicate observations are reduced to scalar parameters by averaging:

| Parameter | Aggregation method | Rounding |
|---|---|---|
| Atom `sigma` | Arithmetic mean | 7 decimals |
| Atom `epsilon` | Arithmetic mean | 3 decimals |
| Atom `charge` | Arithmetic mean | 6 decimals |
| Bond coefficients `K, r0` | Arithmetic mean | 4 decimals |
| Angle coefficients `K, theta0` | Arithmetic mean | 3 decimals |
| Dihedral coefficients `c1, c2, c3, c4` | Arithmetic mean per term | 3 decimals |
| Improper `K` | Arithmetic mean after choosing the modal `(d, n)` pair | 3 decimals |
| Improper `d, n` | Most frequent pair (categorical, not averaged) | none |

**Practical note on averaging**: Atom typing at `hop3` resolution and bonded parameters keyed by `hop0` are already highly specific. It is rare for the same environment combination to produce genuinely different reference values requiring averaging. In normal usage, this aggregation step is effectively a safety net — not a routine occurrence. Users who want to verify or manually adjust any coefficient can inspect the database directly at any time.

## Environment key encoding

Each atom receives four hashed graph descriptors representing increasingly broad local chemical environments.

### hop0 (coarsest)

`hop0` is a compact local descriptor containing:

- **Center atom properties**: atomic number, formal charge, aromaticity, hybridization, total degree, total hydrogens, ring membership, ring count.
- **Neighbor signatures**: For each neighbor, `neighbor_atomic_number : bond_type_code : neighbor_formal_charge`.
- **Bond kinds**: Sorted list of bond type codes (S, D, T, A) to neighbors.

Bond type codes:

| RDKit bond type | Code |
|---|---|
| single | `S` |
| double | `D` |
| triple | `T` |
| aromatic | `A` |
| unknown/unsupported | `U` |

### hop1, hop2, hop3 (expanding shells)

Higher-level descriptors expand the graph outward:

- `hop1`: hop0 plus first-neighbor shell (immediate neighbors of the center atom).
- `hop2`: hop1 plus second-neighbor shell (neighbors of the first-neighbor atoms).
- `hop3`: hop2 plus third-neighbor shell (neighbors of the second-neighbor atoms).

Each level includes:
- Atom properties for atoms in the shell.
- Parent-bond information (bond type code to the parent atom in the previous level).
- Cross-level bonds (bonds connecting atoms within the shell).
- Intra-shell bonds (bonds between atoms both within the same shell).

### Canonicalization

All graphs are canonicalized before hashing:
1. Sort atoms by their properties.
2. Remap atom indices according to the sorted order.
3. Sort bond lists.
4. Serialize as canonical JSON.
5. Hash with SHA-256.

This ensures the same chemical environment always produces the identical key regardless of atom ordering in the input file.

`hop0` is the coarsest fallback and is also used as the key basis for all bonded parameter lookups. `hop3` is the most specific atom environment and is used for primary atom type assignment.

## Atom type resolution

During target parameterization, the resolver queries the `atom_types` table in the following order:

```text
hop3_key → hop2_key → hop1_key → hop0_key
```

The first match is used. The returned database row provides:
- Internal LAMMPS type ID from the database.
- Matched `hop0_key`.
- Element symbol.
- Mass, pair coefficients (sigma, epsilon), and charge (averaged from duplicate environments).
- Source metadata.

If no match is found at any level, the atom is assigned a generated `unknown` output type with zero mass, zero pair coefficients, and zero charge, reported as `no_match` in the output.

The `hop0_keymap`, `hop1_keymap`, and `hop2_keymap` tables are maintained for validation and external inspection, but the fallback resolver currently queries indexed columns in the `atom_types` table directly.

## Bonded parameter lookup

Bonded parameters are keyed by atom `hop0` environments rather than by output atom type IDs. This allows bonded parameters to be independent of atom type numbering.

### Bonds

Bond keys are order-independent. The stored key is:

```text
(min(keyA, keyB), max(keyA, keyB))
```

### Angles

The center atom is fixed; the two outer atoms are swappable:

```text
min((keyA, keyB, keyC), (keyC, keyB, keyA))
```

### Dihedrals

Forward and reverse paths are treated as equivalent:

```text
min((keyA, keyB, keyC, keyD), (keyD, keyC, keyB, keyA))
```

### Impropers

The first atom is treated as the center atom. The three substituents are sorted:

```text
(center_key, sorted(substituent_key1, substituent_key2, substituent_key3))
```

During target parameterization, impropers are filtered before lookup: only C or N center atoms with exactly three neighbors are considered. This matches OPLS-style improper-generation conventions and avoids generating impropers for non-trigonal centers.

### Force field compatibility

`mappingff` has been tested with OPLS-style force field parameters using these LAMMPS styles:

| Section | Style | Functional form | Averaging compatible? |
|---------|-------|-----------------|----------------------|
| `bond_style` | `harmonic` | `E = K·(r-r₀)²` | ✅ 安全 |
| `angle_style` | `harmonic` | `E = K·(θ-θ₀)²` | ✅ 安全 |
| `dihedral_style` | `opls` | `E = ½K₁[1+cos(φ)] + ½K₂[1−cos(2φ)] + ½K₃[1+cos(3φ)] + ½K₄[1−cos(4φ)]` | ✅ 安全（数学平均） |
| `improper_style` | `cvff` | `E = K·[1 + d·cos(n·φ)]` | ✅ 对 CVFF 安全；❌ 其他 style 不兼容 |

**Improper averaging is specific to CVFF style.** The `d` and `n` fields are treated as a categorical pair. The most common pair is selected, and `K` is averaged only over records sharing that pair. This matches CVFF/OPLS conventions where `d` and `n` encode discrete phase and periodicity values. **This approach is NOT compatible with CHARMM-style harmonic impropers, GROMOS impropers, or any style where improper coefficients are all purely numerical.**

Using segment data with different force fields or LAMMPS styles in the same database may produce incorrect averaged parameters. The tool does not validate style consistency.

## Charge adjustment algorithm

Charge adjustment is applied only when `total_charge` or `--charge` is explicitly provided.

Let `delta = target_charge - current_total_charge`.

### Stage 1: Bounded weighted adjustment

Atoms whose matched database atom type has more than one observed charge value (`charge_list` length > 1) are eligible for adjustment.

The adjustment is distributed in proportion to `abs(current_charge)` for each eligible atom, and bounded by the minimum and maximum observed charge values for that database atom type. Atoms that hit their bounds are fixed, and the remaining residual is redistributed to still-active atoms.

### Stage 2: Uniform residual adjustment

Any remaining charge difference is distributed evenly across all atoms in the system.

This guarantees that the final written total approaches the requested target, but it can slightly perturb all charges in stage 2. Always inspect the log to see before/after charge values.

## Output construction

The generated `LammpsData` object is written as a standard LAMMPS data file.

Implementation details:

- Atom coordinates are taken directly from the target structure file.
- The simulation box bounds are set to coordinate min/max plus 5 Å padding in each direction (orthogonal box only).
- All atoms use molecule tag `1`.
- Output atom, bond, angle, dihedral, and improper type IDs are renumbered consecutively starting from 1.
- Bonded terms are always written to the output, even when parameters are missing; missing values are written as zero/default coefficients with warnings in the log.
- Coefficient type deduplication uses numerical tolerances (bond: K < 0.01, r0 < 0.001; angle: K < 0.01, theta0 < 0.1; dihedral/improper: per-element tolerance of 0.01).

## Recommended validation workflow

After generating a LAMMPS file:

1. Read `parameterize.log` and confirm `no_match == 0` whenever possible for production use.
2. Check all missing bonded-parameter warnings and add segments for any missing environments.
3. Inspect total charge, especially if `--charge` was used.
4. Verify the generated type counts and coefficient values are chemically reasonable.
5. Run a short LAMMPS energy minimization or short MD run to check for crashes or obviously incorrect geometries.
6. Compare energies, bonded distributions, or radii of gyration against known reference systems when available.
7. For critical applications, validate against high-level calculations or experimental data for small model compounds.

## Fragment design guidelines

### Atom environments

Each chemically distinct target atom environment should appear in a reference segment with a compatible `hop3` environment whenever possible. The more environments covered at hop3, the fewer fallbacks and the more specific the parameterization.

### Bonded environments

For a dihedral `A-B-C-D`, matching depends on the `hop0` environment of all four atoms. The practical context is:

```text
X-A-B-C-D-Y
```

where X is an outer neighbor of A and Y is an outer neighbor of D. If a polymer chain is cut across this region and capped, the `hop0` environments of A or D can change, and the dihedral may not match.

When a bonded environment is only represented across a segment boundary, include enough overlapping atoms so that each bonded term appears inside at least one reference segment with the same local connectivity as in the target. For junctions, prefer fragments that overlap by at least three atoms around the cut.

### Common failure modes

- **Different bond orders**: If the target uses double bonds but the segment uses aromatic bonds (or vice versa) for the same connectivity, the hop keys will differ and matching will fail or fall back.
- **Protonation differences**: If the target has a deprotonated carboxyl group and the segment has a neutral carboxyl, the neighbor signatures differ.
- **Aromaticity perception**: RDKit may perceive aromaticity differently between a standalone `.lmp` and a `.mol` file for the same structure.
- **Ring membership**: Atoms on ring boundaries can have different ring counts depending on how the ring was cut in the segment.
- **Hydrogen handling**: Explicit vs. implicit hydrogens produce different hop keys. Prefer segments with explicit hydrogens that match the target's hydrogen representation.