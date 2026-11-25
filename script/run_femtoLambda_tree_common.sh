#!/bin/bash
# Environment and explicit arguments for the reduced-tree analysis entry.
set -euo pipefail
if [[ $# -ne 6 ]]; then
  echo "Usage: $0 direct|singularity MAINCONF INPUT OUTPUT JOBID NEVENTS" >&2
  exit 2
fi
FLT_MODE="$1"; FLT_MAINCONF="$2"; FLT_INPUT="$3"
FLT_OUTPUT="$4"; FLT_JOBID="$5"; FLT_NEVENTS="$6"
case "$FLT_MODE" in direct|singularity) ;; *) echo "ERROR: unsupported execution mode" >&2; exit 2;; esac
if [[ -z "$FLT_MAINCONF" || -z "$FLT_INPUT" || -z "$FLT_OUTPUT" || -z "$FLT_JOBID" ||
      ! "$FLT_NEVENTS" =~ ^(-1|[1-9][0-9]*)$ ]]; then
  echo "ERROR: explicit nonempty arguments required; NEVENTS is positive or -1" >&2
  exit 2
fi
FLT_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
FLT_PROJECT="$(cd "$FLT_SCRIPT_DIR/.." && pwd -P)"
cd "$FLT_PROJECT"
absolute_path() {
  case "$1" in /*|*://*) printf '%s' "$1";; *) printf '%s/%s' "$FLT_PROJECT" "$1";; esac
}
FLT_MAINCONF="$(absolute_path "$FLT_MAINCONF")"
FLT_INPUT="$(absolute_path "$FLT_INPUT")"
FLT_OUTPUT="$(absolute_path "$FLT_OUTPUT")"
[[ -f "$FLT_MAINCONF" ]] || { echo "ERROR: missing mainconf $FLT_MAINCONF" >&2; exit 2; }
[[ -e "$FLT_OUTPUT" ]] && { echo "ERROR: refusing existing output $FLT_OUTPUT" >&2; exit 2; }
case "$FLT_INPUT" in *://*) ;; *) [[ -f "$FLT_INPUT" ]] || { echo "ERROR: missing input $FLT_INPUT" >&2; exit 2; };; esac
cpp_quote() {
  local value="$1"
  value="${value//\\/\\\\}"
  value="${value//\"/\\\"}"
  value="${value//$'\n'/\\n}"
  value="${value//$'\r'/\\r}"
  printf '"%s"' "$value"
}
FLT_ROOT_CALL="analysis/run_anaFemtoLambda_4He_tree.C($(cpp_quote "$FLT_INPUT"),$(cpp_quote "$FLT_OUTPUT"),$(cpp_quote "$FLT_JOBID"),$FLT_NEVENTS,$(cpp_quote "$FLT_MAINCONF"))"
echo "[mainconf] (argument) $FLT_MAINCONF"
echo "[input] (argument) $FLT_INPUT"
echo "[output] (argument) $FLT_OUTPUT"
echo "[run] entry=anaFemtoLambda_4He_tree jobid=$FLT_JOBID requested=$FLT_NEVENTS mode=$FLT_MODE"
source "$FLT_SCRIPT_DIR/femtoLambda_environment.sh"
FLT_INPUT_BINDS=("$FLT_INPUT" "$FLT_OUTPUT" "$FLT_MAINCONF")
if [[ "$FLT_INPUT" != *.root && -f "$FLT_INPUT" ]]; then
  while IFS= read -r line || [[ -n "$line" ]]; do
    line="${line%%#*}"
    read -r path unused <<< "$line" || true
    [[ -z "${path:-}" ]] || FLT_INPUT_BINDS+=("$path")
  done < "$FLT_INPUT"
fi
femto_lambda_exec_root "$FLT_MODE" "$FLT_PROJECT" "$FLT_MAINCONF" "$FLT_ROOT_CALL" "${FLT_INPUT_BINDS[@]}"
