from __future__ import annotations

import argparse
import runpy
from typing import Dict

MODULE_MAP: Dict[str, str] = {
    "build-envkey": "macromapff.pipeline.build_envkey_mapping",
    "build-final-keymap": "macromapff.pipeline.build_final_keymap",
    "build-hop-keymap": "macromapff.pipeline.build_hop_keymap",
    "extract-multiatom": "macromapff.pipeline.extract_multiatom_terms",
    "build-multiatom-master": "macromapff.pipeline.build_multiatom_master",
    "generate-lammps": "macromapff.pipeline.generate_lammps_data_from_mol2",
}


def main() -> None:
    parser = argparse.ArgumentParser(
        prog="MacroMapFF",
        description="MacroMapFF molecular parameterization pipeline CLI",
    )
    parser.add_argument("command", choices=sorted(MODULE_MAP.keys()))
    args, unknown = parser.parse_known_args()

    module = MODULE_MAP[args.command]
    # Forward all remaining CLI args to the migrated script module.
    import sys

    sys.argv = [f"MacroMapFF {args.command}", *unknown]
    runpy.run_module(module, run_name="__main__")


if __name__ == "__main__":
    main()
