#!/bin/bash
# Internal environment wrapper. Public entry scripts supply species and execution mode.
set -euo pipefail
if [[ $# -ne 7 ]]; then
  echo "Usage: run_femtoLambda_common.sh direct|singularity SPECIES MAINCONF INPUT OUTPUT JOBID NEVENTS" >&2
  exit 2
fi
FL_MODE="$1"; FL_SPECIES="$2"; FL_MAINCONF="$3"; FL_INPUT="$4"
FL_OUTPUT="$5"; FL_JOBID="$6"; FL_NEVENTS="$7"
case "$FL_SPECIES" in d|t|3He|4He) ;; *) echo "ERROR: unsupported Lambda entry" >&2; exit 2;; esac
case "$FL_MODE" in direct|singularity) ;; *) echo "ERROR: unsupported execution mode" >&2; exit 2;; esac
if [[ -z "$FL_MAINCONF" || -z "$FL_INPUT" || -z "$FL_OUTPUT" || -z "$FL_JOBID" ||
      ! "$FL_NEVENTS" =~ ^(-1|[1-9][0-9]*)$ ]]; then
  echo "ERROR: explicit nonempty arguments required; NEVENTS is positive or -1" >&2; exit 2
fi
FL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
FL_PROJECT="$(cd "$FL_SCRIPT_DIR/.." && pwd -P)"
cd "$FL_PROJECT"
absolute_path() {
  case "$1" in /*|*://*) printf '%s' "$1";; *) printf '%s/%s' "$FL_PROJECT" "$1";; esac
}
FL_MAINCONF="$(absolute_path "$FL_MAINCONF")"
FL_INPUT="$(absolute_path "$FL_INPUT")"
FL_OUTPUT="$(absolute_path "$FL_OUTPUT")"
[[ -f "$FL_MAINCONF" ]] || { echo "ERROR: missing mainconf $FL_MAINCONF" >&2; exit 2; }
[[ -e "$FL_OUTPUT" ]] && { echo "ERROR: refusing existing output $FL_OUTPUT" >&2; exit 2; }
case "$FL_INPUT" in *://*) ;; *) [[ -f "$FL_INPUT" ]] || { echo "ERROR: missing input $FL_INPUT" >&2; exit 2; };; esac
cpp_quote() {
  local value="$1"
  value="${value//\\/\\\\}"
  value="${value//\"/\\\"}"
  value="${value//$'\n'/\\n}"
  value="${value//$'\r'/\\r}"
  printf '"%s"' "$value"
}
FL_ROOT_CALL="analysis/run_anaFemtoLambda_${FL_SPECIES}.C($(cpp_quote "$FL_INPUT"),$(cpp_quote "$FL_OUTPUT"),$(cpp_quote "$FL_JOBID"),$FL_NEVENTS,$(cpp_quote "$FL_MAINCONF"))"
echo "[mainconf] (argument) $FL_MAINCONF"
echo "[input] (argument) $FL_INPUT"
echo "[output] (argument) $FL_OUTPUT"
echo "[run] species=$FL_SPECIES jobid=$FL_JOBID requested=$FL_NEVENTS mode=$FL_MODE"
# Do not source setup.sh: this local run must not activate/change .current_mainconf.
source "$FL_SCRIPT_DIR/femtoLambda_environment.sh"
FL_INPUT_BINDS=("$FL_INPUT" "$FL_OUTPUT" "$FL_MAINCONF")
if [[ "$FL_INPUT" != *.root && -f "$FL_INPUT" ]]; then
  while IFS= read -r line || [[ -n "$line" ]]; do
    line="${line%%#*}"
    read -r path unused <<< "$line" || true
    [[ -z "${path:-}" ]] || FL_INPUT_BINDS+=("$path")
  done < "$FL_INPUT"
fi
femto_lambda_exec_root "$FL_MODE" "$FL_PROJECT" "$FL_MAINCONF" "$FL_ROOT_CALL" "${FL_INPUT_BINDS[@]}"
