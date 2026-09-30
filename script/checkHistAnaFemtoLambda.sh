#!/bin/bash
set -euo pipefail
if [[ $# -ne 3 ]]; then echo "Usage: $0 ROOTFILE MAINCONF OUTPUT_PDF" >&2; exit 2; fi
LQ_ENTRY="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
exec bash "$LQ_ENTRY/checkHistAnaFemtoLambda_common.sh" direct "$@"
