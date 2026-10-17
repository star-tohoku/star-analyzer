#!/bin/bash
# Run StCandidateFemtoMaker (reconstruction parts + femto pairing).
# Host entry (needs an SL7 host, e.g. interactive sl7; on AL9 login nodes use script/singularity_run_anaCandidateFemto.sh).
# Usage: ./script/run_anaCandidateFemto.sh <inputFile|list> <outputFile> <jobid> <nEvents|-1> <mainconf>
# Every argument is required (no default mainconf; docs/ai/AGENT_RULES.md).

set -euo pipefail

if [[ $# -ne 5 ]]; then
  sed -n 2,5p "$0" >&2
  exit 1
fi
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$PROJECT_ROOT"

INPUT_FILE="$1"
OUTPUT_FILE="$2"
JOBID="$3"
NEVENTS="$4"
CONFIG_PATH="$5"

source ./script/setup.sh "$CONFIG_PATH"
export LD_LIBRARY_PATH="$PROJECT_ROOT/lib:${LD_LIBRARY_PATH:-}"
mkdir -p "$(dirname "$OUTPUT_FILE")"

echo "=== anaCandidateFemto ==="
echo "Input:   $INPUT_FILE"
echo "Output:  $OUTPUT_FILE"
echo "JobID:   $JOBID"
echo "nEvents: $NEVENTS"
echo "Config:  $CONFIG_PATH"

root4star -b -q \
  "analysis/run_anaCandidateFemto.C(\"$INPUT_FILE\",\"$OUTPUT_FILE\",\"$JOBID\",$NEVENTS,\"$CONFIG_PATH\")"
