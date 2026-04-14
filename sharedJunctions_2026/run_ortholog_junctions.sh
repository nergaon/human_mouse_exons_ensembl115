#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_SCRIPT="/gpfs0/tals/projects/Analysis/human_mouse_exons/scripts/ensembl115/sharedJunctions_2026/filter_ortholog_junctions.py"
POINTS_FILE="/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl115/unique_points_HN6.txt"
OUT_DIR="/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl115/sharedJunctions_2026"

INPUT_GSE116177="/gpfs0/tals/projects/Analysis/human_mouse_exons/GSE116177/leafcutter_0.2.9/GSE116177_JSR_junction_counts_no_strand.tsv"
INPUT_GSE180020="/gpfs0/tals/projects/Analysis/human_mouse_exons/GSE180020/leafcutter_0.2.9/GSE180020_JSR_junction_counts_no_strand.tsv"
INPUT_GSE60424="/gpfs0/tals/projects/Analysis/human_mouse_exons/GSE60424/leafcutter_0.2.9/GSE60424_JSR_junction_counts_no_strand.tsv"
INPUT_GSE115736="/gpfs0/tals/projects/Analysis/human_mouse_exons/GSE115736/leafcutter_0.2.9/GSE115736_JSR_junction_counts_no_strand.tsv"

DATASETS="both"      # GSE116177 | GSE180020 | GSE60424 | GSE115736 | both
DIRECTION="m2h"      # m2h | h2m | human-filter | both
MIN_READS=10
MIN_SAMPLES=2
DRY_RUN=0

usage() {
  cat <<EOF
Usage:
  $(basename "$0") [options]

Options:
  --datasets VALUE      GSE116177 | GSE180020 | GSE60424 | GSE115736 | both (default: both)
  --direction VALUE     m2h | h2m | human-filter | both (default: m2h)
  --points PATH         Path to unique_points_HN6.txt
  --python-script PATH  Path to filter_ortholog_junctions.py
  --out-dir PATH        Output directory for generated TSVs
  --min-reads N         Minimum reads per sample for express filter (default: 10)
  --min-samples N       Minimum number of samples passing --min-reads (default: 2)
  --dry-run             Print commands without executing
  -h, --help            Show this help

Direction meaning:
  m2h = source position_m -> target position_h
  h2m = source position_h -> target position_m
  human-filter = keep only junctions whose two breakpoints are in position_h

Examples:
  $(basename "$0")
  $(basename "$0") --datasets GSE116177 --direction m2h
  $(basename "$0") --datasets both --direction both
  $(basename "$0") --datasets GSE60424 --direction human-filter
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --datasets)
      DATASETS="$2"
      shift 2
      ;;
    --direction)
      DIRECTION="$2"
      shift 2
      ;;
    --points)
      POINTS_FILE="$2"
      shift 2
      ;;
    --python-script)
      PYTHON_SCRIPT="$2"
      shift 2
      ;;
    --out-dir)
      OUT_DIR="$2"
      shift 2
      ;;
    --min-reads)
      MIN_READS="$2"
      shift 2
      ;;
    --min-samples)
      MIN_SAMPLES="$2"
      shift 2
      ;;
    --dry-run)
      DRY_RUN=1
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown option: $1" >&2
      usage
      exit 1
      ;;
  esac
done

case "$DATASETS" in
  GSE116177|GSE180020|GSE60424|GSE115736|both) ;;
  *)
    echo "Invalid --datasets: $DATASETS" >&2
    exit 1
    ;;
esac

case "$DIRECTION" in
  m2h|h2m|human-filter|both) ;;
  *)
    echo "Invalid --direction: $DIRECTION" >&2
    exit 1
    ;;
esac

if [[ ! -f "$PYTHON_SCRIPT" ]]; then
  echo "Python script not found: $PYTHON_SCRIPT" >&2
  exit 1
fi

if [[ ! -f "$POINTS_FILE" ]]; then
  echo "Points file not found: $POINTS_FILE" >&2
  exit 1
fi

mkdir -p "$OUT_DIR"

run_one() {
  local dataset="$1"
  local direction="$2"
  local input_file mode source_col target_col output_file express_output_file

  if [[ "$dataset" == "GSE116177" ]]; then
    input_file="$INPUT_GSE116177"
  elif [[ "$dataset" == "GSE180020" ]]; then
    input_file="$INPUT_GSE180020"
  elif [[ "$dataset" == "GSE60424" ]]; then
    input_file="$INPUT_GSE60424"
  elif [[ "$dataset" == "GSE115736" ]]; then
    input_file="$INPUT_GSE115736"
  else
    echo "Unknown dataset: $dataset" >&2
    exit 1
  fi

  if [[ "$direction" == "m2h" ]]; then
    mode="remap"
    source_col="position_m"
    target_col="position_h"
    output_file="${OUT_DIR}/${dataset}_orthologs_junctions.tsv"
    express_output_file="${OUT_DIR}/${dataset}_express_orthologs_junctions.tsv"
  elif [[ "$direction" == "h2m" ]]; then
    mode="remap"
    source_col="position_h"
    target_col="position_m"
    output_file="${OUT_DIR}/${dataset}_mouseMapped_junctions.tsv"
    express_output_file="${OUT_DIR}/${dataset}_express_mouseMapped_junctions.tsv"
  else
    mode="filter-only"
    source_col="position_h"
    target_col="position_h"
    output_file="${OUT_DIR}/${dataset}_orthologs_junctions.tsv"
    express_output_file="${OUT_DIR}/${dataset}_express_orthologs_junctions.tsv"
  fi

  cmd=(
    python3 "$PYTHON_SCRIPT"
    --points "$POINTS_FILE"
    --mode "$mode"
    --source-column "$source_col"
    --inputs "$input_file"
    --outputs "$output_file"
    --express-outputs "$express_output_file"
    --min-reads "$MIN_READS"
    --min-samples "$MIN_SAMPLES"
  )

  if [[ "$mode" == "remap" ]]; then
    cmd+=(--target-column "$target_col")
  fi

  echo "Running: dataset=${dataset}, direction=${direction}"
  if [[ $DRY_RUN -eq 1 ]]; then
    printf '  %q ' "${cmd[@]}"
    printf '\n'
  else
    "${cmd[@]}"
  fi
}

is_pair_compatible() {
  local dataset="$1"
  local direction="$2"

  case "$direction" in
    m2h)
      [[ "$dataset" == "GSE116177" || "$dataset" == "GSE180020" ]]
      ;;
    h2m|human-filter)
      [[ "$dataset" == "GSE60424" || "$dataset" == "GSE115736" ]]
      ;;
    *)
      return 1
      ;;
  esac
}

datasets_to_run=()
if [[ "$DATASETS" == "both" ]]; then
  if [[ "$DIRECTION" == "m2h" ]]; then
    datasets_to_run=("GSE116177" "GSE180020")
  elif [[ "$DIRECTION" == "human-filter" ]]; then
    datasets_to_run=("GSE60424" "GSE115736")
  else
    datasets_to_run=("GSE116177" "GSE180020" "GSE60424" "GSE115736")
  fi
else
  datasets_to_run=("$DATASETS")
fi

directions_to_run=()
if [[ "$DIRECTION" == "both" ]]; then
  directions_to_run=("m2h" "h2m" "human-filter")
else
  directions_to_run=("$DIRECTION")
fi

ran_jobs=0
for ds in "${datasets_to_run[@]}"; do
  for dir in "${directions_to_run[@]}"; do
    if ! is_pair_compatible "$ds" "$dir"; then
      echo "Skipping incompatible pair: dataset=${ds}, direction=${dir}" >&2
      continue
    fi
    run_one "$ds" "$dir"
    ((ran_jobs+=1))
  done
done

if [[ ${ran_jobs:-0} -eq 0 ]]; then
  echo "No compatible dataset/direction pairs were selected." >&2
  echo "For m2h use: GSE116177 or GSE180020." >&2
  echo "For h2m/human-filter use: GSE60424 or GSE115736." >&2
  exit 1
fi

echo "Done. Outputs are in: $OUT_DIR"
