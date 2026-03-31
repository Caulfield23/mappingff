#!/usr/bin/env python3
"""Unified domain facade for pure business logic modules."""

from macromapff.domain.atom_typing_core import ENV_INDEX_COLUMNS
from macromapff.domain.atom_typing_core import build_atom_types_core
from macromapff.domain.atom_typing_core import parse_fallback_hops
from macromapff.domain.env import ENV_SPLIT_COLUMNS
from macromapff.domain.env import canonicalize_env_key
from macromapff.domain.env import json_cell
from macromapff.domain.env import split_env_key_columns
from macromapff.domain.env_features import EnvFeatureBuilder
from macromapff.domain.keymap_merge import build_final_map
from macromapff.domain.keymap_merge import finalize_records
from macromapff.domain.multiatom_master_merge import INTERACTION_ORDER
from macromapff.domain.multiatom_master_merge import build_master
from macromapff.domain.multiatom_match_core import assign_multiatom_params_core
from macromapff.domain.multiatom_match_core import build_type_map
from macromapff.domain.multiatom_observed import build_observed_mapping
from macromapff.domain.term_enumeration import enumerate_terms

__all__ = [
	"ENV_INDEX_COLUMNS",
	"ENV_SPLIT_COLUMNS",
	"EnvFeatureBuilder",
	"INTERACTION_ORDER",
	"assign_multiatom_params_core",
	"build_atom_types_core",
	"build_final_map",
	"build_master",
	"build_observed_mapping",
	"build_type_map",
	"canonicalize_env_key",
	"enumerate_terms",
	"finalize_records",
	"json_cell",
	"parse_fallback_hops",
	"split_env_key_columns",
]
