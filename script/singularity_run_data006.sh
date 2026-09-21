#!/bin/bash
# Build p-p, p-d and d-d SE/ME/QA histograms from one schema-3 reduced-tree ROOT file.
# Usage: ./script/singularity_run_data006.sh MAINCONF TREE.root OUT.root [maxEvents] [data006Config] [firstEvent]

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd -P)"

MAINCONF="${1:?Usage: $0 MAINCONF TREE.root OUT.root [maxEvents] [data006Config] [firstEvent]}"
TREE_FILE="${2:?}"
OUT_FILE="${3:?}"
MAX_EVENTS="${4:--1}"
DATA006_CONFIG="${5:-config/maker/data006_auau3p85fxt_pp_pd_dd.yaml}"
FIRST_EVENT="${6:-0}"

exec "$PROJECT_ROOT/script/singularity_run_anaFemtoPhiTreeDownstreamV3.sh" \
  "$MAINCONF" "$TREE_FILE" "$OUT_FILE" \
  -1 bufferAll "" -1 -1 0 0 "$MAX_EVENTS" 150 0.8 3.8 0 "$DATA006_CONFIG" "$FIRST_EVENT"
