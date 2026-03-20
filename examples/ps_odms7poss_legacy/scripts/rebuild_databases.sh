#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PYTHON_BIN="${PYTHON_BIN:-/home/zhangzhao/anaconda3/bin/python}"
OUT_DIR="${OUT_DIR:-${ROOT_DIR}/outputs}"
HOP_DEPTH="${HOP_DEPTH:-2}"
FROM_EXISTING=1
REFRESH_MULTIATOM="${REFRESH_MULTIATOM:-1}"

MODULE_ARGS=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --from-existing)
      FROM_EXISTING=1
      shift
      ;;
    --reextract-segments)
      FROM_EXISTING=0
      shift
      ;;
    --no-refresh-multiatom)
      REFRESH_MULTIATOM=0
      shift
      ;;
    *)
      MODULE_ARGS+=("$1")
      shift
      ;;
  esac
done

if [[ ${#MODULE_ARGS[@]} -gt 0 ]]; then
  MODULES=("${MODULE_ARGS[@]}")
else
  MODULES=(segment1 segment2 segment3 segment4)
fi

mkdir -p "${OUT_DIR}"

FINAL_SPECS=()
MULTI_SPECS=()

for module in "${MODULES[@]}"; do
  structure_path="${ROOT_DIR}/${module}/${module}.mol"
  lmp_path="${ROOT_DIR}/${module}/${module}/${module}.lammps.lmp"
  env_out_dir="${OUT_DIR}/${module}_envdb"
  atom_env_csv="${env_out_dir}/${module}_atom_env.csv"
  multi_prefix="${OUT_DIR}/${module}_multiatom_observed"

  echo "[RUN] ${module}"

  if [[ ${FROM_EXISTING} -eq 1 ]]; then
    if [[ ! -f "${atom_env_csv}" ]]; then
      echo "[ERROR] --from-existing 模式下找不到: ${atom_env_csv}" >&2
      exit 1
    fi
    echo "[SKIP] build_envkey_mapping (${module})，使用已有 ${atom_env_csv}"
  else
    "${PYTHON_BIN}" "${ROOT_DIR}/scripts/build_envkey_mapping.py" \
      --mol "${structure_path}" \
      --lmp "${lmp_path}" \
      --module "${module}" \
      --outdir "${env_out_dir}" \
      --hop-depth "${HOP_DEPTH}"
  fi

  if [[ ${REFRESH_MULTIATOM} -eq 1 ]]; then
    "${PYTHON_BIN}" "${ROOT_DIR}/scripts/extract_multiatom_terms.py" \
      --lmp "${lmp_path}" \
      --atom-env-csv "${atom_env_csv}" \
      --out-prefix "${multi_prefix}"
  else
    if [[ ! -f "${multi_prefix}.csv" ]]; then
      echo "[ERROR] --no-refresh-multiatom 下找不到: ${multi_prefix}.csv" >&2
      exit 1
    fi
    echo "[SKIP] extract_multiatom_terms (${module})，使用已有 ${multi_prefix}.csv"
  fi

  FINAL_SPECS+=("${module}::${atom_env_csv}::${lmp_path}")
  MULTI_SPECS+=("${module}::${multi_prefix}.csv")
done

FINAL_CMD=(
  "${PYTHON_BIN}" "${ROOT_DIR}/scripts/build_final_keymap.py"
  --out-prefix "${OUT_DIR}/final_env_keymap"
)
for spec in "${FINAL_SPECS[@]}"; do
  FINAL_CMD+=(--module-spec "${spec}")
done
"${FINAL_CMD[@]}"

"${PYTHON_BIN}" "${ROOT_DIR}/scripts/build_hop_keymap.py" \
  --final-env-csv "${OUT_DIR}/final_env_keymap.csv" \
  --hop2-out "${OUT_DIR}/hop2_env_keymap.csv" \
  --hop1-out "${OUT_DIR}/hop1_env_keymap.csv" \
  --hop0-out "${OUT_DIR}/hop0_env_keymap.csv"

MULTI_CMD=(
  "${PYTHON_BIN}" "${ROOT_DIR}/scripts/build_multiatom_master.py"
  --type-stats-csv "${OUT_DIR}/final_env_keymap_type_stats.csv"
  --hop0-env-csv "${OUT_DIR}/hop0_env_keymap.csv"
  --out-prefix "${OUT_DIR}/multiatom_master_keytype"
)
for spec in "${MULTI_SPECS[@]}"; do
  MULTI_CMD+=(--multiatom-spec "${spec}")
done
"${MULTI_CMD[@]}"

echo "[DONE] ${OUT_DIR}/final_env_keymap.csv"
echo "[DONE] ${OUT_DIR}/multiatom_master_keytype.csv"
echo "[DONE] ${OUT_DIR}/hop2_env_keymap.csv"
echo "[DONE] ${OUT_DIR}/hop1_env_keymap.csv"
echo "[DONE] ${OUT_DIR}/hop0_env_keymap.csv"
