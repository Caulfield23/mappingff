# mappingff technical notes

This document describes how `mappingff` works internally. It is intended for users and developers who need to understand why a parameter assignment succeeds or fails, how environment keys are generated, and how to diagnose segment coverage problems.

## Source-level architecture

| Module | Responsibility |
|---|---|
| `mappingff.cli` | Command-line parsing, logging setup, and dispatch for `build-db` and `parameterize`. |
| `mappingff.workflow` | High-level `build_db()` and `parameterize()` workflows. |
| `mappingff.mol` | RDKit-based parsing, topology enumeration, improper detection, and hop-key computation. |
| `mappingff.encode` | Local atom-environment encoding: `hop0`, rooted induced subgraph construction, canonicalization, and SHA-256 key generation. |
| `mappingff.db` | SQLite schema, insertion, duplicate aggregation, lookup, and metadata management. |
| `mappingff.fallback` | Atom type resolution in the order `hop3 → hop2 → hop1 → hop0`. |
| `mappingff.lmp` | LAMMPS data parsing, LAMMPS data writing, and total-charge adjustment. |
| `mappingff.lmp2rdkitmol` | Conversion of standalone `.lmp` samples into RDKit molecules using mass-to-element mapping and bond-order inference. |

The public Python API is:

```python
from mappingff import build_db, parameterize
```

All other modules should be treated as implementation details unless explicitly documented.

## Database construction workflow

`build_db(samples_dir, db_path, append=False)` performs the following high-level steps.

### 1. Sample discovery

Two sample modes are supported:

1. **Paired samples**: immediate subdirectories of `samples_dir`, each containing one `.mol` or `.pdb` topology file and one `.lmp` reference LAMMPS data file.
2. **Standalone LAMMPS samples**: `.lmp` files placed directly under `samples_dir`.

Standalone `.lmp` files inside subdirectories are not treated as standalone samples. Keep sample directories unambiguous.

### 2. Per-sample processing

For each sample:

1. Parse the reference LAMMPS data into an internal `LammpsData` object.
2. Parse or reconstruct the RDKit molecular graph.
3. Enumerate atoms, bonds, angles, dihedrals, and impropers from the molecular graph.
4. Compute `hop3`, `hop2`, `hop1`, and `hop0` keys for every atom.
5. Insert atom parameters into `atom_types`.
6. Insert bonded parameters using canonicalized combinations of participating atoms' `hop0` keys.
7. Insert fallback keymap entries for inspection and validation.

### 3. Database finalization

During `save()`, duplicate observations are reduced to scalar values.

| Parameter | Aggregation |
|---|---|
| Atom `sigma` | arithmetic mean |
| Atom `epsilon` | arithmetic mean |
| Atom `charge` | arithmetic mean |
| Bond coefficients | arithmetic mean per coefficient |
| Angle coefficients | arithmetic mean per coefficient |
| Dihedral coefficients | arithmetic mean per coefficient |
| Improper `K` | arithmetic mean within the modal `(d, n)` pair |
| Improper `d, n` | most frequent pair, treated categorically |

Duplicate observations should be rare in a well-designed segment library. When they appear, inspect whether they represent genuine repeated observations or inconsistent reference data.

## Environment key encoding

Each atom receives four hashed descriptors:

```text
hop3_key, hop2_key, hop1_key, hop0_key
```

The fallback resolver uses them from most specific to most general:

```text
hop3 → hop2 → hop1 → hop0
```

The current rooted-ego encoder version is:

```text
mappingff-rooted-ego-v2
```

Any database built with a different encoder version should be rebuilt.

## Public encode API

The current encode module preserves only the API needed by the rest of the package:

```python
get_hop0_subgraph(mol, atom) -> dict
get_hop3_subgraph(mol, atom) -> dict
compute_graph_hop_keys(mol, atom) -> tuple[str, str, str, str]
```

`compute_graph_hop_keys()` returns:

```python
(hop3_key, hop2_key, hop1_key, hop0_key)
```

Private helper functions should not be considered stable.

## hop0 descriptor

`hop0` is a compact descriptor centered on one atom. It is intentionally coarser than the rooted subgraph descriptors.

It contains:

- Center atom atomic number.
- Formal charge.
- Aromaticity flag.
- Hybridization string.
- Total degree.
- Total hydrogens.
- Ring membership.
- Ring count.
- Sorted first-neighbor signatures.
- Sorted bond kind list.

A neighbor signature has the form:

```text
neighbor_atomic_number : bond_type_code : neighbor_formal_charge
```

Bond type codes are:

| RDKit bond type | Code |
|---|---|
| single | `S` |
| double | `D` |
| triple | `T` |
| aromatic | `A` |
| unsupported or unknown | `U` |

`hop0` is used both as the weakest atom-type fallback and as the key basis for bonded-parameter lookup.

## hop1-hop3 rooted induced subgraphs

`hop1`, `hop2`, and `hop3` are generated from rooted induced molecular subgraphs.

For a center atom `c` and radius `r`:

1. Compute graph distances from `c`.
2. Select all atoms with distance `<= r`.
3. Keep every bond whose two endpoints are both selected.
4. Label the center atom as the root.
5. Store each selected atom's distance from the root.
6. Canonicalize the graph without serializing original atom indices.
7. Serialize the canonical representation and hash it with SHA-256.

This differs from a parent-tree shell expansion. In a parent-tree representation, each atom in hop2 or hop3 is attached to one arbitrary parent atom. That is fragile for rings and symmetric environments because the same atom may be reachable through multiple equivalent paths. The rooted induced graph instead preserves all local bonds inside the selected radius.

## Node and edge labels

### Node labels

Rooted subgraph node labels include:

| Field | Meaning |
|---|---|
| `z` | atomic number |
| `fc` | formal charge |
| `ar` | aromaticity flag |
| `deg` | total degree |
| `h` | total hydrogens |
| `ring` | ring-membership flag |
| `ring_count` | number of rings containing the atom |
| `dist` | shortest graph distance from the root |
| `root` | `1` for the center atom, `0` otherwise |

The original RDKit atom index is never included in the serialized graph or hash.

### Edge labels

Edge labels include:

| Field | Meaning |
|---|---|
| `bt` | bond type code: `S`, `D`, `T`, `A`, or `U` |
| `ar` | aromaticity flag |
| `conj` | conjugation flag |
| `ring` | ring-bond flag |

## Canonicalization strategy

The rooted-ego encoder uses deterministic color refinement to produce a stable, index-free serialized graph.

At a high level:

1. Build radius-3 local graph data once for the center atom.
2. For each requested radius, filter the same local graph down to radius 1, 2, or 3.
3. Initialize node colors from node labels.
4. Iteratively refine node colors using sorted neighboring color and edge-label signatures.
5. Group equivalent nodes into node classes.
6. Preserve multiplicity using a `count` field.
7. Represent edges between node classes with multiplicities.
8. Serialize the canonical graph as deterministic JSON.
9. Hash the JSON using SHA-256.

The canonical graph has this conceptual shape:

```json
{
  "version": "mappingff-rooted-ego-v2",
  "kind": "rooted_induced_subgraph",
  "radius": 3,
  "root": "n0",
  "node_classes": [
    {
      "id": "n0",
      "color": "...",
      "label": {
        "z": 6,
        "fc": 0,
        "ar": 0,
        "deg": 4,
        "h": 1,
        "ring": 0,
        "ring_count": 0,
        "dist": 0,
        "root": 1
      },
      "count": 1
    }
  ],
  "edges": [
    {
      "u": "n0",
      "v": "n1",
      "label": {
        "bt": "S",
        "ar": 0,
        "conj": 0,
        "ring": 0
      },
      "count": 1
    }
  ],
  "stats": {
    "nodes": 10,
    "bonds": 10,
    "color_rounds": 3
  }
}
```

The exact hash is an implementation detail. Users should interpret `hop3`, `hop2`, `hop1`, and `hop0` as environment-match levels, not as chemically meaningful names.

## Performance behavior

The optimized path is `compute_graph_hop_keys()`.

Instead of independently building hop1, hop2, and hop3 for each atom, it:

1. Runs one BFS to radius 3.
2. Collects only local bonds by inspecting neighbors of atoms inside that local environment.
3. Reuses the resulting `_EgoData` object for radius 1, 2, and 3.
4. Uses frozen tuple signatures internally to avoid repeated JSON serialization during refinement.
5. Serializes to canonical JSON only for the final hash.

This matters for large polymers because scanning all molecular bonds for every atom is unnecessarily expensive.

## Canonicalization limitations

The current encoder is deterministic and removes dependency on original atom indices, but it is not a full nauty/bliss-style graph canonical labeling backend.

Practical implications:

- It is designed to be stable for typical local organic molecular environments.
- It handles rings, same-shell bonds, cross-level bonds, and symmetric node classes more robustly than a single-parent hop tree.
- Extremely pathological graph pairs that defeat Weisfeiler-Lehman-style color refinement are theoretically possible.
- If strict graph-isomorphism canonical labeling becomes necessary, the encoder can be extended with an optional backend such as RDKit-assisted ranking, igraph/bliss, or nauty/Traces.

For the intended use case, the most important practical improvement is that chemically equivalent local environments should not fail `hop3` matching merely because they appear at different atom indices or in differently sized molecules.

## Atom type resolution

During target parameterization, the resolver queries the database in this order:

```text
hop3_key → hop2_key → hop1_key → hop0_key
```

The first match is used. The returned database row provides:

- Internal database LAMMPS type.
- Element symbol.
- Mass.
- Pair coefficients.
- Charge.
- Source metadata.
- Stored hop0/hop3 graph information.

A `hop3` match means the most specific stored environment was found. A `hop2`, `hop1`, or `hop0` match means a fallback was used and the assignment should be reviewed.

If no match is found, the atom is marked as `no_match` and assigned a generated unknown output type with zero/default values where necessary.

## Bonded parameter lookup

Bonded parameters are keyed by atom `hop0` environments rather than output atom type IDs. This makes lookup independent of output type renumbering.

### Bonds

Bond keys are order-independent:

```text
(min(keyA, keyB), max(keyA, keyB))
```

### Angles

The center atom is fixed; the two outer atoms are swappable:

```text
min((keyA, keyB, keyC), (keyC, keyB, keyA))
```

### Dihedrals

Forward and reverse paths are equivalent:

```text
min((keyA, keyB, keyC, keyD), (keyD, keyC, keyB, keyA))
```

### Impropers

The first atom is treated as the center atom. The three substituents are sorted:

```text
(center_key, sorted(substituent_key1, substituent_key2, substituent_key3))
```

During target parameterization, impropers are filtered before lookup: only C or N center atoms with exactly three neighbors are considered.

## Force-field compatibility

`mappingff` is tested around OPLS-style usage patterns. The tool can store and transfer numeric coefficients, but it does not know the physical meaning of arbitrary force-field styles.

For CVFF-style impropers, the final two coefficients `d` and `n` are treated as a categorical pair. The most frequent `(d, n)` pair is selected, and `K` is averaged only across records sharing that pair.

This is not compatible with all improper forms. In particular, CHARMM-style harmonic impropers and other styles with different coefficient semantics should not be mixed without code changes and validation.

## Charge adjustment algorithm

Charge adjustment is applied only when `total_charge` or `--charge` is explicitly supplied.

Let:

```text
delta = target_charge - current_total_charge
```

### Stage 1: bounded weighted adjustment

Eligible atoms are those whose matched database atom type has more than one observed charge value in `charge_list`.

The adjustment is distributed in proportion to `abs(current_charge)` and bounded by the minimum and maximum observed charge values for that atom type. Atoms that hit bounds are fixed, and the remaining residual is redistributed among still-active atoms.

### Stage 2: uniform residual adjustment

Any remaining charge difference is distributed evenly across all atoms.

This makes the final written total approach the requested target, but stage 2 can slightly perturb every atom charge. Users should inspect the log and validate the charge model.

## Output construction

The generated `LammpsData` object is written as a standard LAMMPS data file.

Implementation details:

- Atom coordinates are taken from the target structure file.
- The simulation box is orthogonal and set to coordinate min/max plus padding.
- All atoms use molecule tag `1`.
- Output atom, bond, angle, dihedral, and improper type IDs are renumbered consecutively from 1.
- Bonded terms are written even when coefficients are missing; missing values are written as zero/default coefficients with warnings.
- Coefficient type deduplication uses numerical tolerances.

## Fragment design guidelines

### Atom environments

Each chemically distinct target atom should appear in a reference segment with a compatible `hop3` environment whenever possible.

Because hop1-hop3 are rooted induced subgraphs, the database should preserve not only shell membership but also local bonds inside each radius. Ring cuts, end caps, protonation changes, aromaticity changes, or altered hydrogen representation can all change the key.

### Bonded environments

Bonded environments depend on `hop0` keys.

For a dihedral:

```text
A-B-C-D
```

the practical context is often:

```text
X-A-B-C-D-Y
```

where X and Y are first neighbors of the outer dihedral atoms A and D.

If a chain is cut and capped such that A or D gains a different neighbor signature, the dihedral key may not match. When a bonded environment is represented only across a junction, include enough overlap so that the bonded term appears inside at least one reference segment with the same local connectivity as the target.

## Common failure modes

- **Old database after encoder changes**: rebuild the database whenever the environment encoding logic changes.
- **Different bond orders**: single/double/aromatic differences change keys.
- **Different aromaticity perception**: standalone `.lmp` reconstruction and `.mol` input may not perceive aromaticity identically.
- **Different hydrogen representation**: explicit and implicit hydrogens affect atom labels.
- **Protonation or charge differences**: formal charges and neighbor signatures affect keys.
- **Segment boundary artifacts**: end caps can change `hop0` and rooted subgraph environments.
- **Insufficient hop3 coverage**: frequent hop2/hop1/hop0 fallbacks indicate missing or incomplete segment coverage.
- **Missing bonded coefficients**: add reference segments that contain the missing bonded environment with matching `hop0` context.

## Recommended validation workflow

After parameterization:

1. Confirm `no_match == 0` whenever possible.
2. Review all fallback counts: `hop2`, `hop1`, and `hop0`.
3. Inspect missing bonded-parameter warnings.
4. Check total charge and per-type charge consistency.
5. Verify masses, pair coefficients, and bonded coefficients.
6. Run a short LAMMPS minimization or short MD test.
7. Compare against a reference system or model compound when possible.
