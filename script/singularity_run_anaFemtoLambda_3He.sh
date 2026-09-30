#!/bin/bash
# Explicit CLI: MAINCONF INPUT OUTPUT JOBID NEVENTS. Old implicit-mainconf mode is removed.
set -euo pipefail
if [[ $# -ne 5 ]]; then
  echo "Usage: $0 MAINCONF INPUT OUTPUT JOBID NEVENTS" >&2
  exit 2
fi
FL_ENTRY_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
exec bash "$FL_ENTRY_DIR/run_femtoLambda_common.sh" singularity 3He "$@"
