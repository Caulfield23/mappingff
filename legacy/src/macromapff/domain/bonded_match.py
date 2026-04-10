#!/usr/bin/env python3
"""Slot-based matching of bonded terms against global bonded patterns."""

from itertools import product


def _tuple_matches_allowed(key_tuple, allowed):
    """Check whether each slot value belongs to allowed slot sets."""
    if len(key_tuple) != len(allowed):
        return False
    return all(key_tuple[i] in allowed[i] for i in range(len(key_tuple)))


def _improper_matches_pattern(key_tuple, pattern):
    """Match improper tuple with fixed center and permutable outer atoms."""
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


def assign_bonded_params(
    kind: str,
    terms,
    atom_key_sets,
    idx_rev_patterns,
    idx_imp_patterns,
    idx_rev_inverted,
    idx_imp_center_inverted,
):
    """Assign best-matched coeffs for each multi-atom term candidate set."""
    records = []
    missing = []
    ambiguous = []
    tuple_match_cache = {}
    cache_hit = 0
    cache_miss = 0

    reversible = kind in {"bond", "angle", "dihedral"}
    patterns = idx_rev_patterns[kind] if reversible else None

    def _candidate_ids(tuple_for_lookup):
        """Find candidate pattern ids by intersecting inverted indexes."""
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
        """Normalize tuple orientation for stable cache reuse."""
        if reversible:
            rev = tuple(reversed(key_tuple))
            return key_tuple if key_tuple <= rev else rev
        return (key_tuple[0],) + tuple(sorted(key_tuple[1:]))

    def _resolve_coeff_candidates(key_tuple):
        """Resolve all coefficient candidates for one key tuple."""
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

    return records, missing, ambiguous, cache_hit, cache_miss, len(tuple_match_cache)


def build_type_map(records):
    """Assign compact local type ids for unique coefficient tuples."""
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
