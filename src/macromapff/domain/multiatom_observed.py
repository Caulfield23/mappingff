#!/usr/bin/env python3

INTERACTION_ORDER = {"bond": 0, "angle": 1, "dihedral": 2, "improper": 3}


def canonicalize_tuple(kind: str, env_tuple: tuple, type_tuple: tuple):
    """Canonicalize term tuples so equivalent permutations share one key."""
    if kind in {"bond", "angle", "dihedral"}:
        fwd = (env_tuple, type_tuple)
        rev = (tuple(reversed(env_tuple)), tuple(reversed(type_tuple)))
        return fwd if fwd <= rev else rev
    if kind == "improper" and len(env_tuple) == 4 and len(type_tuple) == 4:
        center_env = env_tuple[0]
        center_type = type_tuple[0]
        others = sorted(zip(env_tuple[1:], type_tuple[1:]))
        env_new = (center_env,) + tuple(x[0] for x in others)
        type_new = (center_type,) + tuple(x[1] for x in others)
        return env_new, type_new
    return env_tuple, type_tuple


def build_observed_mapping(coeffs, terms, atom_map):
    """Build observed multi-atom mapping rows and grouped summaries."""
    observed_records = {}
    env_summary = {}

    for kind in ["bond", "angle", "dihedral", "improper"]:
        for term in terms[kind]:
            atom_ids = term["atom_ids"]
            term_type = term["term_type"]
            coeff_params = coeffs[kind].get(term_type, [])

            try:
                env_tuple_raw = tuple(atom_map[a]["env_key"] for a in atom_ids)
                obs_types_raw = tuple(atom_map[a]["lmp_type"] for a in atom_ids)
            except KeyError as exc:
                raise ValueError(
                    f"atom_env is missing atom index {exc.args[0]}; cannot build multi-atom mapping."
                )

            can_env_obs, can_types_obs = canonicalize_tuple(
                kind, env_tuple_raw, obs_types_raw
            )

            summary_key = (kind, can_env_obs)
            if summary_key not in env_summary:
                env_summary[summary_key] = {
                    "interaction_kind": kind,
                    "n_atoms": len(can_env_obs),
                    "env_key_tuple": list(can_env_obs),
                    "n_terms": 0,
                    "observed_type_tuples": set(),
                    "observed_term_types": set(),
                }
            env_summary[summary_key]["n_terms"] += 1
            env_summary[summary_key]["observed_type_tuples"].add(can_types_obs)
            env_summary[summary_key]["observed_term_types"].add(term_type)

            rec_key = (kind, can_env_obs, can_types_obs)
            if rec_key not in observed_records:
                observed_records[rec_key] = {
                    "interaction_kind": kind,
                    "n_atoms": len(can_env_obs),
                    "env_key_tuple": list(can_env_obs),
                    "lmp_type_tuple": list(can_types_obs),
                    "term_count": 0,
                    "source_term_types": set(),
                    "coeff_param_sets": set(),
                }

            rec = observed_records[rec_key]
            rec["term_count"] += 1
            rec["source_term_types"].add(term_type)
            rec["coeff_param_sets"].add(tuple(coeff_params))

    observed_list = []
    for _, rec in sorted(
        observed_records.items(),
        key=lambda x: (
            INTERACTION_ORDER.get(x[1]["interaction_kind"], 99),
            tuple(x[1]["env_key_tuple"]),
            tuple(x[1]["lmp_type_tuple"]),
        ),
    ):
        observed_list.append(
            {
                "interaction_kind": rec["interaction_kind"],
                "n_atoms": rec["n_atoms"],
                "env_key_tuple": rec["env_key_tuple"],
                "lmp_type_tuple": rec["lmp_type_tuple"],
                "term_count": int(rec["term_count"]),
                "source_term_types": sorted(rec["source_term_types"]),
                "coeff_param_sets": [list(x) for x in sorted(rec["coeff_param_sets"])],
            }
        )

    summary_list = []
    for _, s in sorted(
        env_summary.items(),
        key=lambda x: (
            INTERACTION_ORDER.get(x[1]["interaction_kind"], 99),
            tuple(x[1]["env_key_tuple"]),
        ),
    ):
        summary_list.append(
            {
                "interaction_kind": s["interaction_kind"],
                "n_atoms": s["n_atoms"],
                "env_key_tuple": s["env_key_tuple"],
                "n_terms": int(s["n_terms"]),
                "n_observed_type_tuples": len(s["observed_type_tuples"]),
                "observed_type_tuples": [list(x) for x in sorted(s["observed_type_tuples"])],
                "observed_term_types": sorted(s["observed_term_types"]),
            }
        )

    return observed_list, summary_list
