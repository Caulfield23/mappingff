#!/usr/bin/env python3
import json

from rdkit import Chem

from macromapff.domain.env import ordered_env_key_obj


class EnvFeatureBuilder:
    """Builds atom-centered environment features and canonical env keys."""

    def bond_type_code(self, bond: Chem.Bond) -> str:
        """Map RDKit bond types to compact one-letter codes."""
        bt = bond.GetBondType()
        if bt == Chem.rdchem.BondType.SINGLE:
            return "S"
        if bt == Chem.rdchem.BondType.DOUBLE:
            return "D"
        if bt == Chem.rdchem.BondType.TRIPLE:
            return "T"
        if bt == Chem.rdchem.BondType.AROMATIC:
            return "A"
        return str(bt)

    def precompute_atom_context(self, mol: Chem.Mol):
        """Precompute reusable atom context such as distance matrix."""
        dist_matrix = Chem.GetDistanceMatrix(mol)

        context = {}
        for atom in mol.GetAtoms():
            idx = atom.GetIdx()
            context[idx] = {
                "dist_matrix": dist_matrix,
            }

        return context

    def hop_shell_signatures(
        self, mol: Chem.Mol, atom_idx: int, max_hop: int = 2, dist_matrix=None
    ):
        """Collect sorted neighbor shell signatures up to a hop depth."""
        max_hop = min(int(max_hop), 2)
        dist = dist_matrix if dist_matrix is not None else Chem.GetDistanceMatrix(mol)
        shells = {hop: [] for hop in range(1, max_hop + 1)}
        for idx in range(mol.GetNumAtoms()):
            hop = int(dist[atom_idx][idx])
            if hop < 1 or hop > max_hop:
                continue
            atom = mol.GetAtomWithIdx(idx)
            token = {
                "z": atom.GetAtomicNum(),
                "fc": atom.GetFormalCharge(),
                "ar": int(atom.GetIsAromatic()),
                "deg": atom.GetDegree(),
                "h": atom.GetTotalNumHs(includeNeighbors=True),
                "ring": int(atom.IsInRing()),
            }
            shells[hop].append(token)

        out = {}
        for hop in range(1, max_hop + 1):
            out[f"hop{hop}_shell"] = sorted(
                shells[hop],
                key=lambda x: (
                    x["z"],
                    x["fc"],
                    x["ar"],
                    x["deg"],
                    x["h"],
                    x["ring"],
                ),
            )
        return out

    def make_env_key(
        self,
        mol: Chem.Mol,
        atom: Chem.Atom,
        hop_depth: int = 2,
        atom_ctx: dict | None = None,
    ):
        """Create canonical env-key string and raw feature dict for an atom."""
        atom_idx = atom.GetIdx()
        nbr_sigs = []
        bond_kinds = []

        for bond in atom.GetBonds():
            nbr = bond.GetOtherAtom(atom)
            bcode = self.bond_type_code(bond)
            bond_kinds.append(bcode)
            nbr_sigs.append(f"{nbr.GetAtomicNum()}:{bcode}:{int(nbr.GetIsAromatic())}")

        ring_info = mol.GetRingInfo()
        ring_count = ring_info.NumAtomRings(atom_idx)
        ctx = atom_ctx or {}
        dist_matrix = ctx.get("dist_matrix")

        features = {
            "z": atom.GetAtomicNum(),
            "formal_charge": atom.GetFormalCharge(),
            "aromatic": int(atom.GetIsAromatic()),
            "hybridization": str(atom.GetHybridization()),
            "in_ring": int(atom.IsInRing()),
            "ring_count": ring_count,
            "degree": atom.GetDegree(),
            "total_hs": atom.GetTotalNumHs(includeNeighbors=True),
            "neighbor_sig": sorted(nbr_sigs),
            "bond_kinds": sorted(bond_kinds),
        }
        features.update(
            self.hop_shell_signatures(
                mol, atom_idx, max_hop=hop_depth, dist_matrix=dist_matrix
            )
        )

        key = json.dumps(
            ordered_env_key_obj(features),
            ensure_ascii=False,
            sort_keys=False,
            separators=(",", ":"),
        )
        return key, features
