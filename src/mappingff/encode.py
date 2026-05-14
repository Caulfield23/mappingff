"""Environment key encoding for mappingff.

Public API preserved by this module:

    get_hop0_subgraph(mol, atom) -> dict
    get_hop3_subgraph(mol, atom) -> dict
    compute_graph_hop_keys(mol, atom) -> tuple[str, str, str, str]

The encoder generates hop1/hop2/hop3 keys from rooted induced molecular
subgraphs.  Original RDKit atom indices are used only as temporary internal
handles and are never serialized or hashed.

Performance note
----------------
compute_graph_hop_keys() builds the radius-3 local environment once, then reuses
that local graph to derive hop1, hop2 and hop3 keys.  This avoids running three
separate BFS passes and avoids scanning all molecule bonds for every center atom.
"""

from __future__ import annotations

from collections import Counter, deque
from dataclasses import dataclass
import hashlib
import json
from typing import Any

from rdkit import Chem
from rdkit.Chem import rdchem

FINGERPRINT_VERSION = "mappingff-rooted-ego-v2"

_BOND_TYPE_CODE = {
    Chem.rdchem.BondType.SINGLE: "S",
    Chem.rdchem.BondType.DOUBLE: "D",
    Chem.rdchem.BondType.TRIPLE: "T",
    Chem.rdchem.BondType.AROMATIC: "A",
}


@dataclass(slots=True)
class _EgoData:
    """Reusable local graph data for one center atom."""

    center_idx: int
    distances: dict[int, int]
    labels: dict[int, dict[str, Any]]
    frozen_labels: dict[int, Any]
    edges: list[tuple[int, int, dict[str, Any], Any]]


# ── serialization helpers ────────────────────────────────────────────────────


def _canonical_json(obj: Any) -> str:
    return json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=True)


def _sha256_json(obj: Any) -> str:
    return hashlib.sha256(_canonical_json(obj).encode("utf-8")).hexdigest()


def _freeze(value: Any) -> Any:
    """Convert nested dict/list structures to deterministic tuples."""

    if isinstance(value, dict):
        return tuple((key, _freeze(value[key])) for key in sorted(value))
    if isinstance(value, (list, tuple)):
        return tuple(_freeze(item) for item in value)
    return value


def _digest_frozen(value: Any) -> str:
    return hashlib.sha256(repr(value).encode("utf-8")).hexdigest()


# ── atom and bond invariants ──────────────────────────────────────────────────


def _bond_type_code(bond: rdchem.Bond) -> str:
    return _BOND_TYPE_CODE.get(bond.GetBondType(), "U")


def _bond_label(bond: rdchem.Bond) -> dict[str, int | str]:
    return {
        "bt": _bond_type_code(bond),
        "ar": int(bond.GetIsAromatic()),
        "conj": int(bond.GetIsConjugated()),
        "ring": int(bond.IsInRing()),
    }


def _ring_count(mol: rdchem.Mol, atom_idx: int) -> int:
    try:
        return int(mol.GetRingInfo().NumAtomRings(atom_idx))
    except RuntimeError:
        return 0


def _atom_base_props(atom: rdchem.Atom, mol: rdchem.Mol) -> dict[str, int]:
    """Atom properties used by all rooted subgraph labels.

    The original atom index is intentionally excluded.
    """

    return {
        "z": int(atom.GetAtomicNum()),
        "fc": int(atom.GetFormalCharge()),
        "ar": int(atom.GetIsAromatic()),
        "deg": int(atom.GetTotalDegree()),
        "h": int(atom.GetTotalNumHs(includeNeighbors=True)),
        "ring": int(atom.IsInRing()),
        "ring_count": _ring_count(mol, atom.GetIdx()),
    }


def _atom_label(
    atom: rdchem.Atom,
    mol: rdchem.Mol,
    *,
    distance: int,
    is_root: bool,
) -> dict[str, int]:
    label = _atom_base_props(atom, mol)
    label["dist"] = int(distance)
    label["root"] = int(is_root)
    return label


# ── public API: hop0 descriptor ───────────────────────────────────────────────


def get_hop0_subgraph(mol: rdchem.Mol, atom: rdchem.Atom) -> dict[str, Any]:
    """Extract the hop0 descriptor centered on atom.

    hop0 is intentionally coarse: center atom properties plus sorted first-neighbor
    signatures and sorted bond kinds.  It remains useful as the weakest fallback
    level and as the bonded-term atom key in the existing database workflow.
    """

    center_idx = atom.GetIdx()
    ring_info = mol.GetRingInfo()

    center_props: dict[str, Any] = {
        "z": int(atom.GetAtomicNum()),
        "formal_charge": int(atom.GetFormalCharge()),
        "aromatic": int(atom.GetIsAromatic()),
        "hybridization": str(atom.GetHybridization()),
        "degree": int(atom.GetTotalDegree()),
        "total_hs": int(atom.GetTotalNumHs(includeNeighbors=True)),
        "in_ring": int(atom.IsInRing()),
        "ring_count": int(ring_info.NumAtomRings(center_idx)),
    }

    neighbor_sig: list[str] = []
    bond_kinds: list[str] = []

    for neighbor in atom.GetNeighbors():
        bond = mol.GetBondBetweenAtoms(center_idx, neighbor.GetIdx())
        bt_code = _bond_type_code(bond)
        neighbor_sig.append(
            f"{neighbor.GetAtomicNum()}:{bt_code}:{neighbor.GetFormalCharge()}"
        )
        bond_kinds.append(bt_code)

    center_props["neighbor_sig"] = sorted(neighbor_sig)
    center_props["bond_kinds"] = sorted(bond_kinds)

    return {
        "version": FINGERPRINT_VERSION,
        "kind": "hop0",
        "radius": 0,
        "center": center_props,
    }


# ── local rooted graph construction ───────────────────────────────────────────


def _shortest_distances_within_radius(
    mol: rdchem.Mol,
    center_idx: int,
    radius: int,
) -> dict[int, int]:
    distances: dict[int, int] = {center_idx: 0}
    queue: deque[int] = deque([center_idx])

    while queue:
        current_idx = queue.popleft()
        current_distance = distances[current_idx]
        if current_distance >= radius:
            continue

        current_atom = mol.GetAtomWithIdx(current_idx)
        for neighbor in current_atom.GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            if neighbor_idx in distances:
                continue
            distances[neighbor_idx] = current_distance + 1
            queue.append(neighbor_idx)

    return distances


def _prepare_ego_data(
    mol: rdchem.Mol, atom: rdchem.Atom, max_radius: int = 3
) -> _EgoData:
    """Build reusable rooted local graph data up to max_radius."""

    center_idx = atom.GetIdx()
    distances = _shortest_distances_within_radius(mol, center_idx, max_radius)
    included = set(distances)

    labels: dict[int, dict[str, Any]] = {}
    frozen_labels: dict[int, Any] = {}

    for idx, distance in distances.items():
        label = _atom_label(
            mol.GetAtomWithIdx(idx),
            mol,
            distance=distance,
            is_root=(idx == center_idx),
        )
        labels[idx] = label
        frozen_labels[idx] = _freeze(label)

    edges: list[tuple[int, int, dict[str, Any], Any]] = []
    seen_pairs: set[tuple[int, int]] = set()

    # Only inspect neighbors of atoms in the local environment.  This avoids a
    # full mol.GetBonds() scan for every atom in large molecules.
    for idx in included:
        rd_atom = mol.GetAtomWithIdx(idx)
        for neighbor in rd_atom.GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            if neighbor_idx not in included:
                continue

            u, v = sorted((idx, neighbor_idx))
            pair = (u, v)
            if pair in seen_pairs:
                continue
            seen_pairs.add(pair)

            bond = mol.GetBondBetweenAtoms(u, v)
            label = _bond_label(bond)
            edges.append((u, v, label, _freeze(label)))

    return _EgoData(
        center_idx=center_idx,
        distances=distances,
        labels=labels,
        frozen_labels=frozen_labels,
        edges=edges,
    )


def _filter_radius(
    data: _EgoData,
    radius: int,
) -> tuple[set[int], list[tuple[int, int, dict[str, Any], Any]]]:
    included = {idx for idx, distance in data.distances.items() if distance <= radius}
    edges = [edge for edge in data.edges if edge[0] in included and edge[1] in included]
    return included, edges


# ── canonicalization ──────────────────────────────────────────────────────────


def _rank_signatures(signatures: dict[int, Any]) -> dict[int, int]:
    rank_by_signature = {
        signature: rank
        for rank, signature in enumerate(sorted(set(signatures.values())))
    }
    return {idx: rank_by_signature[signature] for idx, signature in signatures.items()}


def _adjacency_from_edges(
    included: set[int],
    edges: list[tuple[int, int, dict[str, Any], Any]],
) -> dict[int, list[tuple[int, Any]]]:
    adjacency: dict[int, list[tuple[int, Any]]] = {idx: [] for idx in included}
    for u, v, _label, frozen_label in edges:
        adjacency[u].append((v, frozen_label))
        adjacency[v].append((u, frozen_label))
    return adjacency


def _refine_colors(
    data: _EgoData,
    included: set[int],
    edges: list[tuple[int, int, dict[str, Any], Any]],
) -> tuple[dict[int, int], dict[int, list[tuple[int, Any]]], int]:
    """Deterministic WL-style color refinement for stable serialization."""

    adjacency = _adjacency_from_edges(included, edges)
    colors = _rank_signatures({idx: data.frozen_labels[idx] for idx in included})
    rounds = 0

    for round_idx in range(max(1, len(included))):
        signatures: dict[int, Any] = {}
        for idx in included:
            neighbor_terms = tuple(
                sorted(
                    (edge_label, colors[neighbor_idx])
                    for neighbor_idx, edge_label in adjacency[idx]
                )
            )
            signatures[idx] = (data.frozen_labels[idx], neighbor_terms)

        new_colors = _rank_signatures(signatures)
        rounds = round_idx + 1
        if new_colors == colors:
            break
        colors = new_colors

    return colors, adjacency, rounds


def _final_node_signatures(
    data: _EgoData,
    included: set[int],
    adjacency: dict[int, list[tuple[int, Any]]],
    colors: dict[int, int],
) -> dict[int, Any]:
    signatures: dict[int, Any] = {}
    for idx in included:
        neighbor_terms = tuple(
            sorted(
                (edge_label, colors[neighbor_idx])
                for neighbor_idx, edge_label in adjacency[idx]
            )
        )
        signatures[idx] = (data.frozen_labels[idx], neighbor_terms)
    return signatures


def _canonicalize_radius(data: _EgoData, radius: int) -> dict[str, Any]:
    """Canonicalize the rooted induced subgraph at a given radius."""

    included, edges = _filter_radius(data, radius)
    colors, adjacency, rounds = _refine_colors(data, included, edges)
    final_signatures = _final_node_signatures(data, included, adjacency, colors)

    class_key_by_idx: dict[int, tuple[str, Any]] = {}
    class_entries: dict[tuple[str, Any], dict[str, Any]] = {}

    for idx in included:
        color_digest = _digest_frozen(final_signatures[idx])
        class_key = (color_digest, data.frozen_labels[idx])
        class_key_by_idx[idx] = class_key

        if class_key not in class_entries:
            class_entries[class_key] = {
                "id": "",
                "color": color_digest,
                "label": data.labels[idx],
                "count": 0,
            }
        class_entries[class_key]["count"] += 1

    ordered_class_items = sorted(
        class_entries.items(),
        key=lambda item: (
            0 if item[1]["label"].get("root", 0) else 1,
            item[1]["label"].get("dist", 0),
            _freeze(item[1]["label"]),
            item[1]["color"],
        ),
    )

    class_id_by_key: dict[tuple[str, Any], str] = {}
    class_order_by_id: dict[str, int] = {}
    ordered_classes: list[dict[str, Any]] = []

    for order, (class_key, entry) in enumerate(ordered_class_items):
        class_id = f"n{order}"
        entry["id"] = class_id
        class_id_by_key[class_key] = class_id
        class_order_by_id[class_id] = order
        ordered_classes.append(entry)

    edge_counter: Counter[tuple[str, str, Any]] = Counter()
    edge_label_by_frozen: dict[Any, dict[str, Any]] = {}

    for u, v, label, frozen_label in edges:
        u_class = class_id_by_key[class_key_by_idx[u]]
        v_class = class_id_by_key[class_key_by_idx[v]]

        if class_order_by_id[u_class] <= class_order_by_id[v_class]:
            left, right = u_class, v_class
        else:
            left, right = v_class, u_class

        edge_counter[(left, right, frozen_label)] += 1
        edge_label_by_frozen[frozen_label] = label

    canonical_edges = [
        {
            "u": left,
            "v": right,
            "label": edge_label_by_frozen[frozen_label],
            "count": count,
        }
        for (left, right, frozen_label), count in sorted(
            edge_counter.items(),
            key=lambda item: (
                class_order_by_id[item[0][0]],
                class_order_by_id[item[0][1]],
                item[0][2],
            ),
        )
    ]

    return {
        "version": FINGERPRINT_VERSION,
        "kind": "rooted_induced_subgraph",
        "radius": int(radius),
        "root": "n0",
        "node_classes": ordered_classes,
        "edges": canonical_edges,
        "stats": {
            "nodes": len(included),
            "bonds": len(edges),
            "color_rounds": rounds,
        },
    }


def _hop_key(data: _EgoData, radius: int) -> str:
    return _sha256_json(_canonicalize_radius(data, radius))


# ── public API: hop3 descriptor and all keys ──────────────────────────────────


def get_hop3_subgraph(mol: rdchem.Mol, atom: rdchem.Atom) -> dict[str, Any]:
    """Return the canonical radius-3 rooted induced subgraph for atom."""

    data = _prepare_ego_data(mol, atom, max_radius=3)
    return _canonicalize_radius(data, 3)


def compute_graph_hop_keys(
    mol: rdchem.Mol,
    atom: rdchem.Atom,
) -> tuple[str, str, str, str]:
    """Compute and return (hop3_key, hop2_key, hop1_key, hop0_key)."""

    data = _prepare_ego_data(mol, atom, max_radius=3)

    hop3_key = _hop_key(data, 3)
    hop2_key = _hop_key(data, 2)
    hop1_key = _hop_key(data, 1)
    hop0_key = _sha256_json(get_hop0_subgraph(mol, atom))

    return hop3_key, hop2_key, hop1_key, hop0_key
