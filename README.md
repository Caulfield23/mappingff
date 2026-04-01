# MacroMapFF

MacroMapFF is an open-source molecular parameterization workflow that builds atom/multi-atom parameter databases from segment-level data and generates full LAMMPS data files from complete molecular structures.

## Project status

- Version: v0.1.0 (initial open-source layout)
- Current scope: working pipeline migration and repository standardization
- Algorithm logic: unchanged from validated internal scripts

## Why this structure

This repository uses a hybrid pattern:

- Python package (`src/macromapff`): for versioning, testing, extension, and reuse
- CLI (`MacroMapFF`): for reproducible production workflows
- Legacy bundle in `examples/ps_odms7poss_legacy`: fully isolated historical workflow and assets

## Quick start

### 1. Install

```bash
python -m pip install -e .
```

### 2. Build database from samples

```bash
MacroMapFF build-db examples/ps_odms7poss_legacy
```

### 3. Append new samples (optional)

```bash
MacroMapFF add-samples /path/to/new_samples_folder
```

`add-samples` writes sample `*_env` folders directly under `--db-dir` (same level as existing sample env folders), overwrites on module-name conflict, and then auto-rebuilds global databases by scanning all `*_env` folders.

### 4. Parameterize a new molecule

```bash
MacroMapFF parameterize /path/to/new_molecule.mol
```

### 5. Legacy example (reference only)

The complete original runnable workflow is isolated in:

- `examples/ps_odms7poss_legacy/`

Run from its own scripts folder if you need the historical workflow:

```bash
cd examples/ps_odms7poss_legacy/scripts
./rebuild_databases.sh
./generate_from_current_db.sh
```

## Repository layout

See the full documentation:

- docs/OPEN_SOURCE_GUIDE_zh.md

## Example assets

All historical assets are now bundled under one isolated example:

- `examples/ps_odms7poss_legacy/PS-oDMS7POSS.mol`
- `examples/ps_odms7poss_legacy/segment1/`
- `examples/ps_odms7poss_legacy/segment2/`
- `examples/ps_odms7poss_legacy/segment3/`
- `examples/ps_odms7poss_legacy/segment4/`
- `examples/ps_odms7poss_legacy/outputs/`
- `examples/ps_odms7poss_legacy/scripts/`

This keeps the production source code (`src/`) and the legacy runnable snapshot completely separated.

## Development workflow

- Branch model: `main` + feature branches
- CI on push/PR: install, CLI smoke test, pytest
- Release trigger: git tag `v*` (build wheel/sdist artifacts)

## Standard regression validation

This project includes a fixed end-to-end regression test with local fixtures:

- segment dataset: `tests/fixtures/standard/segdata`
- target molecule: `tests/fixtures/standard/target/PS-oDMS7POSS.mol`

Run:

```bash
pytest tests/test_standard_validation.py -q
```

This validates both core user workflows:

1. Build merged mapping databases from all segment `.lammps.lmp` samples.
2. Parameterize the target full molecule using the built database.

All build outputs are written to a fixed repository-local path for inspection:

- `tests/artifacts/standard/`

This artifacts folder is ignored by Git, so results are easy to inspect locally and never uploaded to GitHub.

## License

MIT. See LICENSE.
