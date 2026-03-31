#!/usr/bin/env python3
import json
from collections import defaultdict

from rdkit import Chem

from macromapff.domain.env import json_cell
from macromapff.domain.env_features import EnvFeatureBuilder


ENV_INDEX_COLUMNS = [
    "z",
    "formal_charge",
    "aromatic",
    "hybridization",
    "degree",
    "total_hs",
    "in_ring",
    "ring_count",
    "hop1_shell",
    "hop2_shell",
    "neighbor_sig",
    "bond_kinds",
]


def _env_index_key_from_obj(obj: dict):
    """Convert env feature dict into the structured index tuple key."""
    key = []
    for col in ENV_INDEX_COLUMNS:
        key.append(json_cell(obj.get(col, "")))
    return tuple(key)


def _env_key_from_obj(obj: dict):
    """Serialize env feature dict into canonical compact JSON key."""
    norm = {}
    for col in ENV_INDEX_COLUMNS:
        val = obj.get(col, "")
        if val == "" or val is None:
            continue
        norm[col] = val
    return json.dumps(norm, ensure_ascii=False, separators=(",", ":"), sort_keys=True)


def _env_features_at_hop(features: dict, hop: int):
    """Project features to a reduced representation for a target hop depth."""
    reduced = dict(features)
    for key in list(reduced.keys()):
        if key.startswith("hop") and key.endswith("_shell"):
            try:
                num = int(key[3:-6])
            except Exception:
                continue
            if num > hop:
                reduced[key] = ""
    return reduced


def parse_fallback_hops(raw: str):
    """Parse comma-separated fallback hop list into unique ordered tuple."""
    if raw is None:
        return tuple()
    items = []
    for token in str(raw).split(","):
        token = token.strip()
        if not token:
            continue
        hop = int(token)
        if hop < 0:
            raise ValueError(f"Fallback hop must be >= 0; got: {hop}")
        if hop not in items:
            items.append(hop)
    return tuple(items)


def build_atom_types_core(
    mol: Chem.Mol,
    hop2_env_to_atom_param: dict,
    hop1_env_to_atom_param: dict,
    hop0_env_to_atom_param: dict,
    hop_depth: int,
    fallback_hops=(),
):
    """Assign atom parameter payloads by exact and fallback env matching."""
    env_builder = EnvFeatureBuilder()
    atom_context = env_builder.precompute_atom_context(mol)

    atom_records = []
    missing = []
    fallback_hit_counter = defaultdict(int)
    struct_key_cache = {}

    def _lookup_params(db, features: dict, hop: int):
        """Lookup one atom's params from env-index and env-key maps."""
        cache_key = (id(features), hop)
        reduced = struct_key_cache.get(cache_key)
        if reduced is None:
            reduced = _env_features_at_hop(features, hop)
            struct_key_cache[cache_key] = reduced

        env_key = _env_key_from_obj(reduced)
        env_hit = db["env"].get(env_key)
        if env_hit is not None:
            return env_hit

        idx_key = _env_index_key_from_obj(reduced)
        return db["structured"].get(idx_key)

    for atom_idx in range(mol.GetNumAtoms()):
        atom = mol.GetAtomWithIdx(atom_idx)
        env_key_raw, env_features = env_builder.make_env_key(
            mol,
            atom,
            hop_depth=hop_depth,
            atom_ctx=atom_context.get(atom_idx),
        )
        params = _lookup_params(hop2_env_to_atom_param, env_features, 2)
        fallback_candidates = []

        matched_hop = 2
        if params is None:
            for fb_hop in fallback_hops:
                if fb_hop == 1:
                    params = _lookup_params(hop1_env_to_atom_param, env_features, 1)
                    fallback_candidates.append((1, "[structured-index]"))
                elif fb_hop == 0:
                    params = _lookup_params(hop0_env_to_atom_param, env_features, 0)
                    fallback_candidates.append((0, "[structured-index]"))
                else:
                    continue

                if params is not None:
                    matched_hop = fb_hop
                    fallback_hit_counter[fb_hop] += 1
                    break

        if params is None:
            missing.append(
                {
                    "atom_index": atom_idx + 1,
                    "symbol": atom.GetSymbol(),
                    "env_key_exact": env_key_raw,
                    "fallback_candidates": fallback_candidates,
                }
            )
            atom_records.append(
                {
                    "atom_id": atom_idx + 1,
                    "global_key_id": 0,
                    "global_key_ids": [0],
                    "charge": 0.0,
                    "sigma": 0.0,
                    "epsilon": 0.0,
                    "mass": 0.0,
                    "matched_hop_depth": -1,
                }
            )
            continue

        atom_records.append(
            {
                "atom_id": atom_idx + 1,
                "global_key_id": params["key_ids"][0],
                "global_key_ids": sorted(set(params["key_ids"])),
                "charge": params["charge"],
                "sigma": params["sigma"],
                "epsilon": params["epsilon"],
                "mass": params["mass"],
                "matched_hop_depth": matched_hop,
            }
        )

    used_global = sorted({a["global_key_id"] for a in atom_records})
    g2l = {gid: i + 1 for i, gid in enumerate(used_global)}

    for a in atom_records:
        a["atom_type"] = g2l[a["global_key_id"]]

    local_type_params = []
    for gid in used_global:
        any_atom = next(a for a in atom_records if a["global_key_id"] == gid)
        local_type_params.append(
            {
                "local_type": g2l[gid],
                "global_key_id": gid,
                "mass": any_atom["mass"],
                "sigma": any_atom["sigma"],
                "epsilon": any_atom["epsilon"],
            }
        )
    local_type_params.sort(key=lambda x: x["local_type"])

    return atom_records, local_type_params, dict(sorted(fallback_hit_counter.items())), missing
