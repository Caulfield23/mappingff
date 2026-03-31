#!/usr/bin/env python3
import json

ENV_KEY_PRIORITY = [
    "z",
    "formal_charge",
    "aromatic",
    "hybridization",
    "degree",
    "total_hs",
    "in_ring",
    "ring_count",
    "neighbor_sig",
    "bond_kinds",
]

SHELL_ITEM_KEY_PRIORITY = ["z", "fc", "ar", "deg", "h", "ring"]

ENV_SPLIT_COLUMNS = [
    "z",
    "formal_charge",
    "aromatic",
    "hybridization",
    "degree",
    "total_hs",
    "in_ring",
    "ring_count",
    "neighbor_sig",
    "bond_kinds",
    "hop1_shell",
    "hop2_shell",
]


def ordered_env_key_obj(obj: dict):
    """Reorder env-key dict fields into a stable, preferred key order."""
    if not isinstance(obj, dict):
        return obj

    out = {}
    for key in ENV_KEY_PRIORITY:
        if key in obj:
            out[key] = obj[key]

    for key in sorted(k for k in obj.keys() if k not in out):
        out[key] = obj[key]

    return out


def _stable_sort_list(values):
    """Sort JSON-like values deterministically by serialized representation."""
    return sorted(
        values,
        key=lambda x: json.dumps(
            x, ensure_ascii=False, sort_keys=True, separators=(",", ":")
        ),
    )


def _normalize_env_obj(obj, parent_key: str = ""):
    """Recursively normalize env objects for deterministic serialization."""
    if isinstance(obj, dict):
        if parent_key == "":
            ordered = ordered_env_key_obj(obj)
            out = {}
            for key, val in ordered.items():
                out[key] = _normalize_env_obj(val, parent_key=key)
            return out

        if all(k in SHELL_ITEM_KEY_PRIORITY for k in obj.keys()):
            out = {}
            for key in SHELL_ITEM_KEY_PRIORITY:
                if key in obj:
                    out[key] = _normalize_env_obj(obj[key], parent_key=key)
            for key in sorted(k for k in obj.keys() if k not in out):
                out[key] = _normalize_env_obj(obj[key], parent_key=key)
            return out

        out = {}
        for key in sorted(obj.keys()):
            out[key] = _normalize_env_obj(obj[key], parent_key=key)
        return out

    if isinstance(obj, list):
        norm = [_normalize_env_obj(v, parent_key=parent_key) for v in obj]
        if parent_key in {"neighbor_sig", "bond_kinds"} or (
            parent_key.startswith("hop") and parent_key.endswith("_shell")
        ):
            return _stable_sort_list(norm)
        return norm

    return obj


def canonicalize_env_key(env_key_raw: str, deep_normalize: bool = True) -> str:
    """Normalize raw env-key JSON into a stable canonical string."""
    try:
        obj = json.loads(env_key_raw)
        if deep_normalize:
            obj = _normalize_env_obj(obj)
        else:
            obj = ordered_env_key_obj(obj)
        return json.dumps(
            obj,
            ensure_ascii=False,
            sort_keys=False,
            separators=(",", ":"),
        )
    except Exception:
        return (env_key_raw or "").strip()


def json_cell(v):
    """Convert scalar/list/dict values into CSV-friendly string cells."""
    if isinstance(v, (list, dict)):
        return json.dumps(v, ensure_ascii=False, separators=(",", ":"), sort_keys=True)
    if v is None:
        return ""
    return str(v)


def split_env_key_columns(canonical_env_key: str):
    """Split canonical env-key JSON into fixed structured columns."""
    cols = {k: "" for k in ENV_SPLIT_COLUMNS}
    try:
        obj = json.loads(canonical_env_key)
    except Exception:
        return cols

    if not isinstance(obj, dict):
        return cols

    for k in ENV_SPLIT_COLUMNS:
        if k in obj:
            cols[k] = json_cell(obj[k])
    return cols
