#!/usr/bin/env python3
"""Extract one sample's observed bonded-term mappings from a LAMMPS data file."""

from pathlib import Path

from macromapff.domain.bonded_observed import build_observed_mapping
from macromapff.io.input import load_atom_env
from macromapff.io.input import parse_lammps_topology_and_coeffs
from macromapff.io.output import write_observed_csv


def extract_bondedterms_sample(lmp_path: Path, atom_env_csv: Path, out_prefix: Path):
    """Extract one sample's observed bonded-term mapping CSV data."""
    coeffs, terms = parse_lammps_topology_and_coeffs(lmp_path)
    _, atom_map = load_atom_env(atom_env_csv)
    observed_list, summary_list = build_observed_mapping(
        coeffs=coeffs,
        terms=terms,
        atom_map=atom_map,
    )

    out_csv = out_prefix.with_suffix(".csv")
    write_observed_csv(out_csv, observed_list)
    return out_csv, observed_list, summary_list
