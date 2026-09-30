#!/bin/bash
set -euo pipefail
if [[ $# -ne 4 ]]; then
  echo "Usage: $0 direct|singularity ROOTFILE MAINCONF OUTPUT_PDF" >&2; exit 2
fi
LQ_MODE="$1"; LQ_INPUT="$2"; LQ_CONFIG="$3"; LQ_OUTPUT="$4"
LQ_SCRIPT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
LQ_PROJECT="$(cd "$LQ_SCRIPT/.." && pwd -P)"
cd "$LQ_PROJECT"
absolute_path() { case "$1" in /*) printf '%s' "$1";; *) printf '%s/%s' "$LQ_PROJECT" "$1";; esac; }
LQ_INPUT="$(absolute_path "$LQ_INPUT")"
LQ_CONFIG="$(absolute_path "$LQ_CONFIG")"
LQ_OUTPUT="$(absolute_path "$LQ_OUTPUT")"
[[ -f "$LQ_INPUT" && -f "$LQ_CONFIG" && "$LQ_OUTPUT" == *.pdf ]] || { echo "ERROR: ROOT/config required and PDF suffix required" >&2; exit 2; }
[[ ! -e "$LQ_OUTPUT" ]] || { echo "ERROR: refusing existing PDF" >&2; exit 2; }
cpp_quote() {
  local value="$1"
  value="${value//\\/\\\\}"; value="${value//\"/\\\"}"
  value="${value//$'\n'/\\n}"; value="${value//$'\r'/\\r}"
  printf '"%s"' "$value"
}
LQ_CALL="analysis/run_checkHistAnaFemtoLambda.C($(cpp_quote "$LQ_INPUT"),$(cpp_quote "$LQ_CONFIG"),$(cpp_quote "$LQ_OUTPUT"))"
source "$LQ_SCRIPT/femtoLambda_environment.sh"
femto_lambda_exec_root "$LQ_MODE" "$LQ_PROJECT" "$LQ_CONFIG" "$LQ_CALL" "$LQ_INPUT" "$LQ_CONFIG" "$LQ_OUTPUT"
