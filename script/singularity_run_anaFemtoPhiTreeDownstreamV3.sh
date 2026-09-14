#!/bin/bash
# Run tree downstream SE/ME from a reduced ROOT file inside singularity.
# Usage: ./script/singularity_run_anaFemtoPhiTreeDownstream.sh MAINCONF TREE.root OUT.root [bufferSize] [mixingMode] [nSigmaD] [dcaD] [sigMin] [sigMax]

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
PROJECT_ROOT_REAL="$(cd "$PROJECT_ROOT" && pwd -P)"
PYTHON="$(command -v python3 2>/dev/null || command -v python 2>/dev/null)"
SINGULARITY_IMAGE="/cvmfs/singularity.opensciencegrid.org/star-bnl/star-sw:latest"

MAINCONF="${1:?Usage: $0 MAINCONF TREE OUT [bufferSize] [mixingMode] [nSigmaD] [dcaD] [sigMin] [sigMax]}"
TREE_FILE="${2:?}"
OUT_FILE="${3:?}"
BUFFER="${4:--1}"
MODE="${5:-bufferAll}"
NSIGD="${6:--1}"
DCAD="${7:--1}"
SIGMIN="${8:--1}"
SIGMAX="${9:--1}"
NDEDX="${10:--1}"
RECOPID="${11:-0}"
RECOMIX="${12:-0}"

mkdir -p "$(dirname "$OUT_FILE")"
LIBRARY_TAG=$(cd "$PROJECT_ROOT_REAL" && "$PYTHON" script/analysis_info_helper.py --library-tag --mainconf "$MAINCONF" | xargs)
source "$SCRIPT_DIR/setup.sh" "$MAINCONF"
export LD_LIBRARY_PATH="$PROJECT_ROOT_REAL/lib:${LD_LIBRARY_PATH:-}"

CONTAINER_CMD=$(cat <<EOF
export STAR=/star/nfs4/AFS/star/packages/$LIBRARY_TAG
export STAR_HOST_SYS=$STAR_HOST_SYS
export STAR_LIB=\${STAR}/.\${STAR_HOST_SYS}/lib
export STAR_BIN=\${STAR}/.\${STAR_HOST_SYS}/bin
export PATH=\${STAR_BIN}:\${STAR}/mgr:\$PATH
export LD_LIBRARY_PATH=\${STAR_LIB}:\$LD_LIBRARY_PATH
cd "$PROJECT_ROOT_REAL"
root4star -b -q 'tools/run_anaFemtoPhiTreeDownstreamV3.C("$TREE_FILE","$OUT_FILE","$MAINCONF",$BUFFER,"$MODE",$NSIGD,$DCAD,$SIGMIN,$SIGMAX,$NDEDX,$RECOPID,$RECOMIX)'
EOF
)

exec singularity exec \
  -B /gpfs:/gpfs \
  -B /star/u:/star/u \
  -B /star/nfs4/AFS:/star/nfs4/AFS \
  -B /home/starlib:/home/starlib \
  --env LD_LIBRARY_PATH="$LD_LIBRARY_PATH" \
  "$SINGULARITY_IMAGE" \
  /bin/bash -lc "$CONTAINER_CMD"
