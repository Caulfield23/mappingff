#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

PYTHON_BIN="${PYTHON_BIN:-/home/zhangzhao/anaconda3/bin/python}"
STRUCTURE="${1:-../PS-oDMS7POSS.mol}"
OUT="${2:-../outputs/$(basename "${STRUCTURE%.mol}")_from_db_hop2.lammps.lmp}"
HOP_DEPTH="${HOP_DEPTH:-2}"
FALLBACK_HOPS="${FALLBACK_HOPS:-1,0}"
HOP2_DB="${HOP2_DB:-../outputs/hop2_env_keymap.csv}"
HOP1_DB="${HOP1_DB:-../outputs/hop1_env_keymap.csv}"
HOP0_DB="${HOP0_DB:-../outputs/hop0_env_keymap.csv}"

mkdir -p "$(dirname "$OUT")"

"$PYTHON_BIN" build_hop_keymap.py \
  --final-env-csv ../outputs/final_env_keymap.csv \
  --hop2-out "$HOP2_DB" \
  --hop1-out "$HOP1_DB" \
  --hop0-out "$HOP0_DB"

"$PYTHON_BIN" generate_lammps_data_from_mol2.py \
  --structure "$STRUCTURE" \
  --hop-depth "$HOP_DEPTH" \
  --fallback-hops "$FALLBACK_HOPS" \
  --final-env-csv ../outputs/final_env_keymap.csv \
  --hop2-env-csv "$HOP2_DB" \
  --hop1-env-csv "$HOP1_DB" \
  --hop0-env-csv "$HOP0_DB" \
  --multiatom-csv ../outputs/multiatom_master_keytype.csv \
  --out "$OUT"

echo "生成完成: $OUT"
