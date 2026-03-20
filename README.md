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

### 2. Explore CLI

```bash
MacroMapFF --help
```

### 3. Run isolated legacy example

The complete original runnable workflow is isolated in:

- `examples/ps_odms7poss_legacy/`

Run from its own scripts folder:

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

## License

MIT. See LICENSE.
