#!/bin/bash
# Data-driven BDT (common/macro/datadriven_bdt, Python + xgboost) for one YAML recipe.
# Usage (from anywhere; paths in the recipe are relative to the repository root):
#   ./script/run_datadriven_bdt.sh <recipe.yaml> train
#   ./script/run_datadriven_bdt.sh <recipe.yaml> apply
#   ./script/run_datadriven_bdt.sh <recipe.yaml> study                         (training-size study)
#   ./script/run_datadriven_bdt.sh <recipe.yaml> report [label=apply.root ...]
#   ./script/run_datadriven_bdt.sh <recipe.yaml> model <dir>                   (model file for StCandidateFemtoMaker)
#   ./script/run_datadriven_bdt.sh <recipe.yaml> check <trees.root>            (wildcard in the file name allowed; quote it)
#   ./script/run_datadriven_bdt.sh <recipe.yaml> inspect                       (input size and memory needed; reads nothing)
# Runs on the host Python (xgboost, shap, uproot, scipy, matplotlib); STAR's PYTHONPATH is removed because it breaks
# those imports. Rules: docs/ai/AGENT_RULES.md "Candidate BDT"; procedure: docs/ai/skills/candidate-bdt.md.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT_REAL="$(cd "$SCRIPT_DIR/.." && pwd -P)"
[[ $# -ge 2 ]] || { sed -n '2,11p' "$0" >&2; exit 1; }
resolve() { if [[ "$1" = /* ]]; then echo "$1"; else echo "$PWD/$1"; fi; }
RECIPE="$(resolve "$1")"; shift
[[ -f "$RECIPE" ]] || { echo "ERROR: recipe not found: $RECIPE" >&2; exit 1; }
PYTHON="$(command -v python3)"
ARGS=("$@")
# a path argument of model/check is resolved against the caller's directory too
if [[ "${ARGS[0]}" == model || "${ARGS[0]}" == check ]] && [[ ${#ARGS[@]} -ge 2 ]]; then ARGS[1]="$(resolve "${ARGS[1]}")"; fi
echo "datadriven_bdt: recipe $RECIPE, action ${ARGS[*]}, python $PYTHON, repository $PROJECT_ROOT_REAL"
cd "$PROJECT_ROOT_REAL"
exec env -u PYTHONPATH "$PYTHON" common/macro/datadriven_bdt/datadriven_bdt.py "$RECIPE" "${ARGS[@]}"
