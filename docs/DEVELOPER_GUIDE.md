# MacroMapFF Developer Guide

This document is for developers and maintainers. It describes internal architecture, implementation boundaries, and module responsibilities.

## 1. Internal Rules

- Function-first design: prefer plain functions unless persistent mutable state is required.
- Boundary clarity:
  - `cli`: argument parsing and command dispatch only.
  - `pipeline`: orchestration and workflow sequencing.
  - `domain`: pure business logic and matching rules.
  - `io`: file/database parsing, serialization, and logging.
- No pseudo-configurability:
  - internal matching strategy, fallback order, and constants stay in code.
  - only user-facing runtime inputs are exposed as CLI arguments.
- No pass-through parameter chains:
  - pass values only where consumed.

## 2. Package Layout

- `src/macromapff/cli.py`: CLI entrypoint.
- `src/macromapff/pipeline/`: workflow orchestration.
- `src/macromapff/domain/`: matching and merge logic.
- `src/macromapff/io/`: parsers, writers, logs.

## 3. Workflow Architecture

### `build-db`
1. Discover sample pairs (`structure`, `.lmp`).
2. Build sample atom maps (`*_AtomMap.csv`).
3. Build sample bonded observations (`*_BondedTerms.csv`).
4. Merge to global atom map (`Global_AtomMap.csv`).
5. Build hop fallback databases (`hop2/hop1/hop0_KeyMap.csv`).
6. Merge global bonded database (`Global_BondedTerms.csv`).

### `add-samples`
- Same sequence as `build-db` after adding/overwriting sample env folders.

### `parameterize`
1. Load hop databases + global bonded database.
2. Atom matching with fallback (`hop2 -> hop1 -> hop0`).
3. Enumerate and match bonded terms.
4. Write final LAMMPS data file + logs.

## 4. Module Responsibilities

### 4.1 CLI
- `macromapff.cli`
  - `main`: command parsing and dispatch.

### 4.2 Pipeline
- `pipeline.workflow`
  - `discover_samples`: find sample inputs.
  - `build_sample_envs`: create per-sample env outputs.
  - `discover_sample_env_records`: collect `_env` records.
  - `merge_database`: rebuild global db artifacts.
  - `build_db`, `add_samples`, `parameterize`: top-level workflows.
- `pipeline.atommap_sample`
  - `build_sample_atommap`: generate one sample atom map.
- `pipeline.bondedterms_sample`
  - `extract_bondedterms_sample`: generate one sample bonded observation table.
- `pipeline.keymap_hop`
  - `build_keymap`: merge global atom map.
  - `build_hop_map`, `build_hop_databases`: fallback db construction.
- `pipeline.global_bonded`
  - `build_global_bonded`: merge global bonded table.
- `pipeline.parameterize`
  - `build_atom_match_with_logs`: atom matching wrapper.
  - `assign_bonded_params_with_logs`: bonded matching wrapper.
  - `parameterize_lammps`: end-to-end parameterization.

### 4.3 Domain
- `domain.env_key_codec`: env key canonicalization and split-column codec.
- `domain.env_key_match`: atom environment feature extraction.
- `domain.atom_match`: atom-level matching and fallback.
- `domain.keymap_merge`: keymap merge and stats.
- `domain.term_enumeration`: bonds/angles/dihedrals/impropers enumeration.
- `domain.bonded_observed`: sample bonded observation mapping.
- `domain.bonded_global_merge`: global bonded merge.
- `domain.bonded_match`: bonded-term coefficient matching.

### 4.4 IO
- `io.input`: structure/lammps/db parsers and loaders.
- `io.output`: CSV/LAMMPS writers.
- `io.log`: build, keymap, matching, conflict logs.

## 5. Key Internal Constants

- `USER_DEFAULT_DB_DIR = "database"` in workflow boundary.
- Atom fallback order in parameterization: `(1, 0)` after hop2 primary matching.
- Hop depth for env feature generation: internal fixed depth (`2`).

## 6. Development Validation

Run after code changes:

```bash
pytest -q
```

For standard fixture validation:

```bash
pytest tests/test_standard_validation.py -q
```
