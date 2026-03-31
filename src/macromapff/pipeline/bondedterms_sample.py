#!/usr/bin/env python3
"""Extract one sample's observed bonded-term mappings from a LAMMPS data file."""

from pathlib import Path

from macromapff.domain import build_observed_mapping
from macromapff.io import load_atom_env
from macromapff.io import parse_lammps_topology_and_coeffs
from macromapff.io import write_observed_csv


class BondedTermsSampleExtractor:
    """Extracts one sample's observed bonded-term mapping CSV data."""

    def __init__(self, lmp_path: Path, atom_env_csv: Path) -> None:
        """Bind extractor to one LAMMPS data file and one atom_env CSV."""
        self.lmp_path = lmp_path
        self.atom_env_csv = atom_env_csv

    def extract(self, out_prefix: Path):
        """Build observed mappings and write them as CSV output."""
        coeffs, terms = parse_lammps_topology_and_coeffs(self.lmp_path)
        _, atom_map = load_atom_env(self.atom_env_csv)
        observed_list, summary_list = build_observed_mapping(
            coeffs=coeffs,
            terms=terms,
            atom_map=atom_map,
        )

        out_csv = out_prefix.with_suffix(".csv")
        write_observed_csv(out_csv, observed_list)
        return out_csv, observed_list, summary_list
