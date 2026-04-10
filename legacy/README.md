# MacroMapFF

MacroMapFF is a local CLI tool that builds parameter databases from sample segments and generates parameterized LAMMPS data files for new molecules.

## Installation

```bash
python -m pip install -e .
```

## CLI Commands

All commands follow:

```bash
MacroMapFF <command> [arguments] [options]
```

### 1. `build-db`

Build a full database from a sample folder.

```bash
MacroMapFF build-db <samples_dir> [--db-dir <db_dir>]
```

Arguments:
- `samples_dir` (required): folder containing sample data.

Options:
- `--db-dir` (optional, default `./database`): output database directory.

Examples:

```bash
MacroMapFF build-db examples/ps_odms7poss_legacy
MacroMapFF build-db /data/samples --db-dir /data/macro_db
```

### 2. `add-samples`

Append new samples and rebuild merged global databases.

```bash
MacroMapFF add-samples <samples_dir> [--db-dir <db_dir>]
```

Arguments:
- `samples_dir` (required): folder containing additional sample data.

Options:
- `--db-dir` (optional, default `./database`): existing database directory to update.

Examples:

```bash
MacroMapFF add-samples /data/new_samples
MacroMapFF add-samples /data/new_samples --db-dir /data/macro_db
```

### 3. `parameterize`

Generate a parameterized LAMMPS data file from one `.mol` structure.

```bash
MacroMapFF parameterize <molecule.mol> [--db-dir <db_dir>] [--out <output.lmp>]
```

Arguments:
- `molecule.mol` (required): input molecule in `.mol` format.

Options:
- `--db-dir` (optional, default `./database`): database directory to read.
- `--out` (optional, default `<molecule_stem>_param.lmp`): output LAMMPS data file.

Examples:

```bash
MacroMapFF parameterize /data/target/PS-oDMS7POSS.mol
MacroMapFF parameterize /data/target/PS-oDMS7POSS.mol --db-dir /data/macro_db
MacroMapFF parameterize /data/target/PS-oDMS7POSS.mol --db-dir /data/macro_db --out /data/out/PS-oDMS7POSS_param.lmp
```

## Typical Workflow

1. Build an initial database:

```bash
MacroMapFF build-db examples/ps_odms7poss_legacy --db-dir ./database
```

2. Add new samples over time:

```bash
MacroMapFF add-samples /path/to/new_samples --db-dir ./database
```

3. Parameterize a new molecule:

```bash
MacroMapFF parameterize /path/to/new_molecule.mol --db-dir ./database --out /path/to/new_molecule_param.lmp
```

## Expected Outputs

After `build-db` / `add-samples` in `--db-dir`:
- `Global_AtomMap.csv`
- `hop_env/hop2_KeyMap.csv`
- `hop_env/hop1_KeyMap.csv`
- `hop_env/hop0_KeyMap.csv`
- `Global_BondedTerms.csv`

After `parameterize`:
- output LAMMPS file (default `<molecule_stem>_param.lmp`)
- `parameterize.log`
- `atom_index_key_types.csv`

## Validation

```bash
pytest -q
```

## Developer Documentation

- Internal architecture and module details:
	- `docs/DEVELOPER_GUIDE.md`

## License

MIT (see LICENSE).
