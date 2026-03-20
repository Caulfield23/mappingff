#!/usr/bin/env python3
from pathlib import Path
import runpy
import sys

THIS_FILE = Path(__file__).resolve()
PROJECT_ROOT = None
for p in THIS_FILE.parents:
    if (p / "src" / "macromapff" / "pipeline").exists():
        PROJECT_ROOT = p
        break

if PROJECT_ROOT is None:
    raise RuntimeError("Cannot locate project root containing src/macromapff/pipeline")

PIPELINE_DIR = PROJECT_ROOT / "src" / "macromapff" / "pipeline"
if str(PIPELINE_DIR) not in sys.path:
    sys.path.insert(0, str(PIPELINE_DIR))

if __name__ == "__main__":
    target = PIPELINE_DIR / (Path(__file__).stem + ".py")
    runpy.run_path(str(target), run_name="__main__")
