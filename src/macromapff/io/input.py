#!/usr/bin/env python3
import csv
from collections import defaultdict
from pathlib import Path

from rdkit import Chem

from macromapff.domain import ENV_INDEX_COLUMNS


SECTION_NAMES = {
    "Masses",
    "Pair Coeffs",
    "Bond Coeffs",
    "Angle Coeffs",
    "Dihedral Coeffs",
    "Improper Coeffs",
    "Atoms",
    "Bonds",
    "Angles",
    "Dihedrals",
    "Impropers",
    "Velocities",
}


class DSU:
    """Simple disjoint-set union for grouping equivalent key ids."""

    def __init__(self):
        """Initialize parent mapping for DSU nodes."""
        self.parent = {}

    def find(self, x):
        """Find representative of a node with path compression."""
        if x not in self.parent:
            self.parent[x] = x
            return x
        if self.parent[x] != x:
            self.parent[x] = self.find(self.parent[x])
        return self.parent[x]

    def union(self, a, b):
        """Union two DSU sets by linking one root to another."""
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[rb] = ra


def load_input_structure(structure_path: Path) -> Chem.Mol:
    """Load a .mol structure for parameterization with robust sanitization."""
    suffix = structure_path.suffix.lower()
    if suffix != ".mol":
        raise ValueError(f"Only .mol input is supported; got: {structure_path}")

    mol = Chem.MolFromMolFile(
        str(structure_path), removeHs=False, sanitize=True, strictParsing=False
    )
    if mol is None:
        mol = Chem.MolFromMolFile(
            str(structure_path), removeHs=False, sanitize=False, strictParsing=False
        )
        if mol is None:
            raise ValueError(f"RDKit failed to read .mol file: {structure_path}")
        mol.UpdatePropertyCache(strict=False)
        try:
            Chem.SanitizeMol(mol)
        except Exception:
            Chem.SanitizeMol(
                mol,
                sanitizeOps=(
                    Chem.SanitizeFlags.SANITIZE_ALL
                    ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE
                ),
            )
    return mol


def load_structure_any(structure_path: Path):
    """Load supported structure formats (.mol/.mol2/.pdb) with fallbacks."""
    suffix = structure_path.suffix.lower()

    def _sanitize_or_raise(mol, source_path: Path, label: str):
        """Sanitize an RDKit molecule or raise a descriptive read error."""
        if mol is None:
            raise ValueError(f"RDKit failed to read {label} file: {source_path}")
        try:
            Chem.SanitizeMol(mol)
        except Exception as exc:
            raise ValueError(
                f"RDKit read {label} but sanitization failed: {source_path} ({exc})"
            )
        return mol

    def _load_pdb(pdb_path: Path):
        """Load pdb file with sanitization fallback and optional soft failure."""
        if not pdb_path.exists():
            return None
        mol = Chem.MolFromPDBFile(str(pdb_path), removeHs=False, sanitize=True)
        if mol is not None:
            return mol
        mol = Chem.MolFromPDBFile(str(pdb_path), removeHs=False, sanitize=False)
        if mol is None:
            return None
        return _sanitize_or_raise(mol, pdb_path, "pdb")

    if suffix == ".mol":
        mol = Chem.MolFromMolFile(str(structure_path), removeHs=False, sanitize=True)
        if mol is not None:
            return mol

        mol = Chem.MolFromMolFile(str(structure_path), removeHs=False, sanitize=False)
        if mol is not None:
            return _sanitize_or_raise(mol, structure_path, "mol")

        pdb_fallback = structure_path.with_suffix(".pdb")
        mol = _load_pdb(pdb_fallback)
        if mol is not None:
            print(
                f"[INFO] Failed to read mol; automatically fell back to pdb input: {pdb_fallback}",
                flush=True,
            )
            return mol

        raise ValueError(
            f"RDKit failed to read mol file: {structure_path}; fallback to same-name pdb also failed: {pdb_fallback}"
        )

    if suffix == ".pdb":
        mol = _load_pdb(structure_path)
        if mol is None:
            raise ValueError(f"RDKit failed to read pdb file: {structure_path}")
        return mol

    if suffix == ".mol2":
        mol = Chem.MolFromMol2File(str(structure_path), removeHs=False, sanitize=True)
        if mol is not None:
            return mol

        mol = Chem.MolFromMol2File(str(structure_path), removeHs=False, sanitize=False)
        if mol is not None:
            return _sanitize_or_raise(mol, structure_path, "mol2")

        return mol

    raise ValueError(
        f"Unsupported structure format: {structure_path}. Supported formats: .mol/.mol2/.pdb"
    )


def _parse_atom_line(toks):
    """Parse one LAMMPS Atoms row under supported atom-style variants."""
    if len(toks) < 6:
        return None

    atom_id = int(toks[0])
    candidates = []
    if len(toks) >= 7:
        candidates.append((2, 3, 4, 5, 6))
    if len(toks) >= 6:
        candidates.append((1, 2, 3, 4, 5))

    for t_col, q_col, x_col, y_col, z_col in candidates:
        try:
            atom_type = int(toks[t_col])
            charge = float(toks[q_col])
            x = float(toks[x_col])
            y = float(toks[y_col])
            z = float(toks[z_col])
            return atom_id, atom_type, charge, x, y, z
        except Exception:
            continue

    return None


def parse_lammps_sections(lmp_path: Path):
    """Parse Masses, Pair Coeffs, and Atoms sections from LAMMPS data."""
    lines = lmp_path.read_text(encoding="utf-8", errors="ignore").splitlines()

    current = None
    masses = {}
    pair_coeffs = {}
    atoms = {}

    for line in lines:
        stripped = line.strip()
        if not stripped:
            continue

        if stripped in SECTION_NAMES:
            current = stripped
            continue

        if stripped.startswith("#"):
            continue

        if current == "Masses":
            toks = stripped.split()
            if len(toks) < 2:
                continue
            if not toks[0].lstrip("+-").isdigit():
                continue
            type_id = int(toks[0])
            masses[type_id] = float(toks[1])

        elif current == "Pair Coeffs":
            toks = stripped.split()
            if len(toks) < 3:
                continue
            if not toks[0].lstrip("+-").isdigit():
                continue
            type_id = int(toks[0])
            pair_coeffs[type_id] = {
                "epsilon": float(toks[1]),
                "sigma": float(toks[2]),
            }

        elif current == "Atoms":
            toks = stripped.split()
            if not toks[0].lstrip("+-").isdigit():
                continue
            parsed = _parse_atom_line(toks)
            if parsed is None:
                continue
            atom_id, atom_type, charge, x, y, z = parsed
            atoms[atom_id] = {
                "lmp_type": atom_type,
                "charge": charge,
                "x": x,
                "y": y,
                "z": z,
            }

    return masses, pair_coeffs, atoms


def parse_lammps_masses(lmp_path: Path):
    """Parse only Masses section from a LAMMPS data file."""
    masses, _, _ = parse_lammps_sections(lmp_path)
    if not masses:
        raise ValueError(f"No Masses section parsed from: {lmp_path}")
    return masses


class LammpsDataParser:
    """Small object wrapper around LAMMPS section parsing helpers."""

    def __init__(self, lmp_path: Path) -> None:
        """Bind parser instance to one LAMMPS data file path."""
        self.lmp_path = lmp_path

    def parse_sections(self):
        """Parse Masses, Pair Coeffs, and Atoms sections."""
        return parse_lammps_sections(self.lmp_path)

    def parse_masses(self):
        """Parse only the Masses section for type-to-mass mapping."""
        return parse_lammps_masses(self.lmp_path)


def infer_atomic_num_from_mass(mass: float, max_z: int = 36, tol: float = 0.2):
    """Infer atomic number by nearest periodic-table mass within tolerance."""
    ptable = Chem.GetPeriodicTable()
    best = None
    for z in range(1, max_z + 1):
        w = ptable.GetAtomicWeight(z)
        d = abs(w - mass)
        if best is None or d < best[0]:
            best = (d, z)
    if best is None or best[0] > tol:
        raise ValueError(
            f"Failed to infer element from mass {mass} (min delta {best[0] if best else 'N/A'})"
        )
    return best[1]


def parse_lammps_data(lmp_path: Path):
    """Parse LAMMPS data and enrich atoms with inferred element/LJ info."""
    masses, pair_coeffs, atoms = parse_lammps_sections(lmp_path)

    if not masses:
        raise ValueError(f"No Masses section parsed from: {lmp_path}")
    if not pair_coeffs:
        raise ValueError(f"No Pair Coeffs section parsed from: {lmp_path}")
    if not atoms:
        raise ValueError(f"No Atoms section parsed from: {lmp_path}")

    type_info = {}
    ptable = Chem.GetPeriodicTable()
    for type_id, mass in masses.items():
        atomic_num = infer_atomic_num_from_mass(mass)
        type_info[type_id] = {
            "mass": mass,
            "atomic_num": atomic_num,
            "type_name": f"LMP_{type_id}",
            "element": ptable.GetElementSymbol(atomic_num),
        }

    for atom_id, info in atoms.items():
        lmp_type = info["lmp_type"]
        if lmp_type not in pair_coeffs:
            raise ValueError(f"LAMMPS atom type {lmp_type} is missing Pair Coeffs")
        if lmp_type not in type_info:
            raise ValueError(f"LAMMPS atom type {lmp_type} is missing Masses")
        info["sigma"] = pair_coeffs[lmp_type]["sigma"]
        info["epsilon"] = pair_coeffs[lmp_type]["epsilon"]
        info["atomic_num"] = type_info[lmp_type]["atomic_num"]
        info["type_name"] = type_info[lmp_type]["type_name"]
        info["mass"] = type_info[lmp_type]["mass"]

    return atoms, type_info


def _parse_int_tokens(line: str):
    """Split a topology/coeff line into tokens after stripping comments."""
    body = line.split("#", 1)[0].strip()
    if not body:
        return []
    toks = body.split()
    if not toks:
        return []
    return toks


def parse_lammps_topology_and_coeffs(lmp_path: Path):
    """Parse topology terms and coefficient sections from LAMMPS data."""
    lines = lmp_path.read_text(encoding="utf-8", errors="ignore").splitlines()

    coeff_sections = {
        "Bond Coeffs": "bond",
        "Angle Coeffs": "angle",
        "Dihedral Coeffs": "dihedral",
        "Improper Coeffs": "improper",
    }
    topo_sections = {
        "Bonds": ("bond", 2),
        "Angles": ("angle", 3),
        "Dihedrals": ("dihedral", 4),
        "Impropers": ("improper", 4),
    }
    known_sections = (
        set(coeff_sections)
        | set(topo_sections)
        | {
            "Masses",
            "Pair Coeffs",
            "Atoms",
            "Velocities",
        }
    )

    coeffs = {"bond": {}, "angle": {}, "dihedral": {}, "improper": {}}
    terms = {"bond": [], "angle": [], "dihedral": [], "improper": []}

    current = None
    for raw in lines:
        stripped = raw.strip()
        if not stripped:
            continue
        if stripped in known_sections:
            current = stripped
            continue
        if stripped.startswith("#"):
            continue

        if current in coeff_sections:
            toks = _parse_int_tokens(stripped)
            if len(toks) < 2:
                continue
            if not toks[0].lstrip("+-").isdigit():
                continue
            kind = coeff_sections[current]
            type_id = int(toks[0])
            params = toks[1:]
            coeffs[kind][type_id] = params
            continue

        if current in topo_sections:
            kind, n_atoms = topo_sections[current]
            toks = _parse_int_tokens(stripped)
            if len(toks) < 2 + n_atoms:
                continue
            if not toks[0].lstrip("+-").isdigit():
                continue
            term_id = int(toks[0])
            term_type = int(toks[1])
            atom_ids = [int(x) for x in toks[2 : 2 + n_atoms]]
            terms[kind].append(
                {
                    "term_id": term_id,
                    "term_type": term_type,
                    "atom_ids": atom_ids,
                }
            )

    return coeffs, terms


def load_atom_env(atom_env_csv: Path):
    """Load atom_env CSV into module name and atom-index mapping."""
    atom_map = {}
    module_name = ""
    with atom_env_csv.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            atom_index = int(row["atom_index"])
            atom_map[atom_index] = {
                "env_key": row["env_key"],
                "lmp_type": int(row["opls_type_id"]),
            }
            if not module_name:
                module_name = row.get("module", "")
    if not atom_map:
        raise ValueError(f"atom_env is empty or unreadable: {atom_env_csv}")
    return module_name, atom_map


def load_hop_param_db(hop_csv: Path):
    """Load hop-level atom parameter database with env and structured indexes."""
    if not hop_csv.exists():
        raise FileNotFoundError(f"Hop parameter database not found: {hop_csv}")

    structured_index = {}
    env_index = {}
    with hop_csv.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            source_ids = [x for x in str(row.get("source_key_ids", "")).split(";") if x]
            key_ids = sorted({int(x) for x in source_ids}) if source_ids else [-1]
            payload = {
                "key_ids": key_ids,
                "charge": float(row["charge_mean"]),
                "sigma": float(row["sigma_mean"]),
                "epsilon": float(row["epsilon_mean"]),
                "mass": float(row["mass_mean"]),
            }
            idx_key = tuple(str(row.get(c, "") or "") for c in ENV_INDEX_COLUMNS)
            structured_index[idx_key] = payload

            env_key = str(row.get("env_key", "") or "").strip()
            if env_key:
                env_index[env_key] = payload

    if not structured_index:
        raise ValueError(f"Hop parameter database is empty: {hop_csv}")
    return {
        "structured": structured_index,
        "env": env_index,
    }


def _parse_json(raw: str):
    """Parse JSON text and return None when decoding fails."""
    import json

    try:
        return json.loads(raw)
    except Exception:
        return None


def load_multiatom_db(multiatom_csv: Path):
    """Load multi-atom DB and build forward/inverted matching indexes."""
    if not multiatom_csv.exists():
        raise FileNotFoundError(f"Multi-atom parameter database not found: {multiatom_csv}")

    reversible_kinds = {"bond", "angle", "dihedral"}
    idx_rev_patterns = defaultdict(list)
    idx_imp_patterns = []
    idx_rev_inverted = defaultdict(lambda: defaultdict(lambda: defaultdict(set)))
    idx_imp_center_inverted = defaultdict(set)

    def _to_allowed_tuple(raw_list):
        """Normalize raw slot values to tuple of frozenset key-id slots."""
        allowed = []
        for item in raw_list:
            if isinstance(item, list):
                allowed.append(frozenset(int(x) for x in item))
            else:
                allowed.append(frozenset([int(item)]))
        return tuple(allowed)

    with multiatom_csv.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            kind = row["interaction_kind"].strip().lower()
            key_tuple_raw = _parse_json(row["key_type_tuple"])
            coeff_sets_raw = _parse_json(row["coeff_param_sets"])
            if not isinstance(key_tuple_raw, list) or not isinstance(
                coeff_sets_raw, list
            ):
                continue
            allowed = _to_allowed_tuple(key_tuple_raw)
            coeff_sets = {
                tuple(str(v) for v in one)
                for one in coeff_sets_raw
                if isinstance(one, list)
            }
            if not coeff_sets:
                continue

            if kind in reversible_kinds:
                pid = len(idx_rev_patterns[kind])
                idx_rev_patterns[kind].append(
                    {
                        "allowed": allowed,
                        "coeff_candidates": set(coeff_sets),
                    }
                )
                arity = len(allowed)
                inv = idx_rev_inverted[kind][arity]
                for pos, allowed_slot in enumerate(allowed):
                    for key_type in allowed_slot:
                        inv[(pos, key_type)].add(pid)
            elif kind == "improper":
                if len(allowed) != 4:
                    continue
                pid = len(idx_imp_patterns)
                idx_imp_patterns.append(
                    {
                        "center_allowed": allowed[0],
                        "other_allowed": (allowed[1], allowed[2], allowed[3]),
                        "coeff_candidates": set(coeff_sets),
                    }
                )
                for center_key in allowed[0]:
                    idx_imp_center_inverted[center_key].add(pid)

    return (
        idx_rev_patterns,
        idx_imp_patterns,
        idx_rev_inverted,
        idx_imp_center_inverted,
    )


def load_type_to_keyid(final_env_csv: Path):
    """Build mapping from global type token to canonical key_id."""
    if not final_env_csv.exists():
        raise FileNotFoundError(f"Final env CSV not found: {final_env_csv}")

    mapping = {}
    with final_env_csv.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            key_id = int(row["key_id"])
            raw = str(row.get("global_type_ids", "") or "")
            for global_type in [x.strip() for x in raw.split(";") if x.strip()]:
                if global_type not in mapping or key_id < mapping[global_type]:
                    mapping[global_type] = key_id

    if not mapping:
        raise ValueError(f"Final env CSV has no global_type_ids: {final_env_csv}")
    return mapping


def load_hop0_key_classes(hop0_csv: Path):
    """Load hop0 connectivity classes and expand each key to its class list."""
    if not hop0_csv.exists():
        raise FileNotFoundError(f"hop0 CSV not found: {hop0_csv}")

    dsu = DSU()
    seen_ids = set()
    with hop0_csv.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            raw = str(row.get("source_key_ids", "")).strip()
            if not raw:
                continue
            ids = [int(x) for x in raw.split(";") if x]
            if not ids:
                continue
            for key_id in ids:
                seen_ids.add(key_id)
                dsu.find(key_id)
            head = ids[0]
            for key_id in ids[1:]:
                dsu.union(head, key_id)

    groups = defaultdict(set)
    for key_id in seen_ids:
        groups[dsu.find(key_id)].add(key_id)

    key_to_class = {}
    for _, members in groups.items():
        cls = sorted(members)
        for key_id in members:
            key_to_class[key_id] = cls
    return key_to_class