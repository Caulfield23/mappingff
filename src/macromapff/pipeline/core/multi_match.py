#!/usr/bin/env python3
import csv
from collections import defaultdict
from itertools import combinations
from itertools import product
from pathlib import Path


def _parse_json(raw: str):
    import json

    try:
        return json.loads(raw)
    except Exception:
        return None


def _norm_rev_tuple(values):
    vals = tuple(values)
    rev = tuple(reversed(vals))
    return vals if vals <= rev else rev


def _append_build_log(log_path: Path, lines):
    if log_path is None:
        return
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("a", encoding="utf-8") as f:
        for line in lines:
            f.write(str(line) + "\n")


def load_multiatom_db(multiatom_csv: Path):
    if not multiatom_csv.exists():
        raise FileNotFoundError(f"Multi-atom parameter database not found: {multiatom_csv}")

    reversible_kinds = {"bond", "angle", "dihedral"}
    idx_rev_patterns = defaultdict(list)
    idx_imp_patterns = []
    idx_rev_inverted = defaultdict(lambda: defaultdict(lambda: defaultdict(set)))
    idx_imp_center_inverted = defaultdict(set)

    def _to_allowed_tuple(raw_list):
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


def _tuple_matches_allowed(key_tuple, allowed):
    if len(key_tuple) != len(allowed):
        return False
    return all(key_tuple[i] in allowed[i] for i in range(len(key_tuple)))


def _improper_matches_pattern(key_tuple, pattern):
    if len(key_tuple) != 4:
        return False
    if key_tuple[0] not in pattern["center_allowed"]:
        return False

    vals = [key_tuple[1], key_tuple[2], key_tuple[3]]
    slots = pattern["other_allowed"]

    for perm in ((0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 2, 0), (2, 0, 1), (2, 1, 0)):
        if all(vals[perm[i]] in slots[i] for i in range(3)):
            return True
    return False


def enumerate_terms(mol):
    bonds = []
    for bond in mol.GetBonds():
        i = bond.GetBeginAtomIdx() + 1
        j = bond.GetEndAtomIdx() + 1
        bonds.append((i, j))

    angles = []
    for center in range(mol.GetNumAtoms()):
        nbs = [n.GetIdx() + 1 for n in mol.GetAtomWithIdx(center).GetNeighbors()]
        for i, k in combinations(sorted(nbs), 2):
            angles.append((i, center + 1, k))

    dihedral_set = set()
    for bond in mol.GetBonds():
        j = bond.GetBeginAtomIdx() + 1
        k = bond.GetEndAtomIdx() + 1
        j_neighbors = [
            x.GetIdx() + 1
            for x in mol.GetAtomWithIdx(j - 1).GetNeighbors()
            if x.GetIdx() + 1 != k
        ]
        k_neighbors = [
            x.GetIdx() + 1
            for x in mol.GetAtomWithIdx(k - 1).GetNeighbors()
            if x.GetIdx() + 1 != j
        ]
        for i in j_neighbors:
            for l in k_neighbors:
                t = (i, j, k, l)
                dihedral_set.add(_norm_rev_tuple(t))
    dihedrals = sorted(dihedral_set)

    impropers = []
    for center in range(mol.GetNumAtoms()):
        nbs = [n.GetIdx() + 1 for n in mol.GetAtomWithIdx(center).GetNeighbors()]
        if len(nbs) < 3:
            continue
        for a, b, c in combinations(sorted(nbs), 3):
            impropers.append((center + 1, a, b, c))

    return bonds, angles, dihedrals, impropers


def assign_multiatom_params(
    kind: str,
    terms,
    atom_key_sets,
    idx_rev_patterns,
    idx_imp_patterns,
    idx_rev_inverted,
    idx_imp_center_inverted,
    strict_missing=True,
    build_log_path: Path = None,
):
    records = []
    missing = []
    ambiguous = []
    tuple_match_cache = {}
    cache_hit = 0
    cache_miss = 0

    reversible = kind in {"bond", "angle", "dihedral"}
    patterns = idx_rev_patterns[kind] if reversible else None

    def _candidate_ids(tuple_for_lookup):
        arity_inv = idx_rev_inverted.get(kind, {}).get(len(tuple_for_lookup))
        if not arity_inv:
            return None

        ids = None
        for pos, key_type in enumerate(tuple_for_lookup):
            hit = arity_inv.get((pos, key_type))
            if not hit:
                return set()
            if ids is None:
                ids = set(hit)
            else:
                ids &= hit
            if not ids:
                return set()
        return ids if ids is not None else set()

    def _normalize_cache_key(key_tuple):
        if reversible:
            rev = tuple(reversed(key_tuple))
            return key_tuple if key_tuple <= rev else rev
        return (key_tuple[0],) + tuple(sorted(key_tuple[1:]))

    def _resolve_coeff_candidates(key_tuple):
        if reversible:
            rev = tuple(reversed(key_tuple))
            coeff_candidates = set()
            cands_fwd = _candidate_ids(key_tuple)
            cands_rev = _candidate_ids(rev)
            if cands_fwd is None or cands_rev is None:
                pattern_ids = range(len(patterns))
            else:
                pattern_ids = cands_fwd | cands_rev

            for pid in pattern_ids:
                pattern = patterns[pid]
                if _tuple_matches_allowed(
                    key_tuple, pattern["allowed"]
                ) or _tuple_matches_allowed(rev, pattern["allowed"]):
                    coeff_candidates.update(pattern["coeff_candidates"])
            return tuple(sorted(coeff_candidates))

        coeff_candidates = set()
        if idx_imp_center_inverted is None:
            pattern_ids = range(len(idx_imp_patterns))
        else:
            pattern_ids = idx_imp_center_inverted.get(key_tuple[0], set())
        for pid in pattern_ids:
            pattern = idx_imp_patterns[pid]
            if _improper_matches_pattern(key_tuple, pattern):
                coeff_candidates.update(pattern["coeff_candidates"])
        return tuple(sorted(coeff_candidates))

    for atom_ids in terms:
        key_options = [sorted(atom_key_sets[x]) for x in atom_ids]
        if any(len(options) == 0 for options in key_options):
            missing.append((atom_ids, tuple()))
            continue

        matched_tuples = {}
        for key_tuple in product(*key_options):
            cache_key = _normalize_cache_key(key_tuple)
            coeff_candidates = tuple_match_cache.get(cache_key)
            if coeff_candidates is None:
                coeff_candidates = _resolve_coeff_candidates(key_tuple)
                tuple_match_cache[cache_key] = coeff_candidates
                cache_miss += 1
            else:
                cache_hit += 1

            if coeff_candidates:
                matched_tuples[key_tuple] = list(coeff_candidates)

        if not matched_tuples:
            combo_sample = []
            for i, combo in enumerate(product(*key_options)):
                combo_sample.append(combo)
                if i >= 19:
                    break
            missing.append((atom_ids, tuple(combo_sample)))
            records.append(
                {
                    "atom_ids": atom_ids,
                    "key_type_tuple": (
                        tuple(combo_sample[0]) if combo_sample else tuple()
                    ),
                    "coeff": None,
                }
            )
            continue

        sorted_matches = sorted(matched_tuples.items(), key=lambda x: x[0])
        chosen_key_tuple, chosen_coeff_list = sorted_matches[0]
        coeff = chosen_coeff_list[0]

        unique_coeffs = sorted(
            {tuple(c) for values in matched_tuples.values() for c in values}
        )
        if len(unique_coeffs) > 1:
            ambiguous.append(
                {
                    "atom_ids": atom_ids,
                    "matched_key_tuples": [k for k, _ in sorted_matches],
                    "matched_coeffs": unique_coeffs,
                    "chosen_key_tuple": chosen_key_tuple,
                    "chosen_coeff": coeff,
                }
            )

        records.append(
            {"atom_ids": atom_ids, "key_type_tuple": chosen_key_tuple, "coeff": coeff}
        )

    if missing:
        sample = ", ".join(f"atom_ids={x[0]} candidates={x[1]}" for x in missing[:10])
        _append_build_log(
            build_log_path,
            [
                (
                    f"[{kind.upper()}][ERROR] missing_count={len(missing)}"
                    if strict_missing
                    else f"[{kind.upper()}][WARN] missing_count={len(missing)}"
                ),
                (
                    f"[{kind.upper()}][ERROR] sample_missing={sample}"
                    if strict_missing
                    else f"[{kind.upper()}][WARN] sample_missing={sample}"
                ),
            ],
        )
        if strict_missing:
            raise ValueError(
                f"{kind} has {len(missing)} terms not found in parameter table. Sample key_type_tuple: {sample}"
            )

    if ambiguous:
        lines = [f"[{kind.upper()}][AMBIGUOUS] count={len(ambiguous)}"]
        for item in ambiguous[:200]:
            lines.append(
                f"[{kind.upper()}][AMBIGUOUS] atom_ids={item['atom_ids']} "
                f"chosen_key_tuple={item['chosen_key_tuple']} chosen_coeff={item['chosen_coeff']} "
                f"all_key_tuples={item['matched_key_tuples']} all_coeffs={item['matched_coeffs']}"
            )
        _append_build_log(build_log_path, lines)

    total_cache_req = cache_hit + cache_miss
    hit_rate = (cache_hit / total_cache_req) if total_cache_req else 0.0
    _append_build_log(
        build_log_path,
        [
            f"[{kind.upper()}][CACHE] hit={cache_hit} miss={cache_miss} size={len(tuple_match_cache)} hit_rate={hit_rate:.4f}"
        ],
    )

    return records, missing, ambiguous


def build_type_map(records):
    coeff_to_type = {}
    for rec in records:
        coeff = rec["coeff"]
        if coeff is None:
            continue
        if coeff not in coeff_to_type:
            coeff_to_type[coeff] = len(coeff_to_type) + 1

    for rec in records:
        if rec["coeff"] is None:
            rec["type_id"] = 0
        else:
            rec["type_id"] = coeff_to_type[rec["coeff"]]

    type_rows = sorted(
        [(tid, coeff) for coeff, tid in coeff_to_type.items()], key=lambda x: x[0]
    )
    return type_rows
