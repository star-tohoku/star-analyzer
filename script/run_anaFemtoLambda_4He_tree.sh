#!/bin/bash
set -euo pipefail
if [[ $# -ne 5 ]]; then
  echo "Usage: $0 MAINCONF INPUT OUTPUT JOBID NEVENTS" >&2
  exit 2
fi
FLT_ENTRY_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
exec bash "$FLT_ENTRY_DIR/run_femtoLambda_tree_common.sh" direct "$@"
