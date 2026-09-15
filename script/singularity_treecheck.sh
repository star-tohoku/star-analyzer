#!/bin/bash
# Run one of tools/treecheck/*.C inside the STAR SL7 container, with the project's headers and
# libraries available. These are the checks the Step 7 runbook calls for.
# Usage: ./script/singularity_treecheck.sh MACRO.C 'arg1,arg2'
#   e.g. ./script/singularity_treecheck.sh MergeCheck.C '"rootfile/x/merged.root"'
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd -P)"
MACRO="${1:?Usage: $0 MACRO.C 'args'}"
ARGS="${2:-}"
MAINCONF="${STAR_ANA_MAINCONF:-config/mainconf/main_auau3p85fxt_anaFemtoPhiTree_prod.yaml}"
IMG=/cvmfs/singularity.opensciencegrid.org/star-bnl/star-sw:latest

cd "$PROJECT_ROOT"
LIBRARY_TAG=$(python3 script/analysis_info_helper.py --library-tag --mainconf "$MAINCONF" | xargs)
source ./script/setup.sh "$MAINCONF" > /dev/null
export LD_LIBRARY_PATH="$PROJECT_ROOT/lib:${LD_LIBRARY_PATH:-}"

CMD="
export STAR=/star/nfs4/AFS/star/packages/$LIBRARY_TAG
export STAR_HOST_SYS=$STAR_HOST_SYS
export STAR_LIB=\${STAR}/.\${STAR_HOST_SYS}/lib
export STAR_BIN=\${STAR}/.\${STAR_HOST_SYS}/bin
export PATH=\${STAR_BIN}:\${STAR}/mgr:\$PATH
export LD_LIBRARY_PATH=\${STAR_LIB}:\$LD_LIBRARY_PATH
cd $PROJECT_ROOT
root4star -b -q tools/treecheck/prelude.C 'tools/treecheck/$MACRO+($ARGS)'
"
exec singularity exec -B /gpfs:/gpfs -B /star/u:/star/u -B /star/nfs4/AFS:/star/nfs4/AFS \
  -B /home/starlib:/home/starlib --env LD_LIBRARY_PATH="$LD_LIBRARY_PATH" "$IMG" /bin/bash -lc "$CMD"
