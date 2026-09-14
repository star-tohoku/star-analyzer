#!/bin/bash
# Dump tree schema/branch sizes inside singularity.
# Usage: ./script/singularity_dumpFemtoPhiTree.sh MAINCONF TREE.root OUT_DIR
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
PROJECT_ROOT_REAL="$(cd "$PROJECT_ROOT" && pwd -P)"
PYTHON="$(command -v python3 2>/dev/null || command -v python 2>/dev/null)"
SINGULARITY_IMAGE="/cvmfs/singularity.opensciencegrid.org/star-bnl/star-sw:latest"
MAINCONF="${1:?}"
TREE_FILE="${2:?}"
OUT_DIR="${3:?}"
mkdir -p "$OUT_DIR"
LIBRARY_TAG=$(cd "$PROJECT_ROOT_REAL" && "$PYTHON" script/analysis_info_helper.py --library-tag --mainconf "$MAINCONF" | xargs)
source "$SCRIPT_DIR/setup.sh" "$MAINCONF"
CONTAINER_CMD=$(cat <<EOF
export STAR=/star/nfs4/AFS/star/packages/$LIBRARY_TAG
export STAR_HOST_SYS=$STAR_HOST_SYS
export STAR_LIB=\${STAR}/.\${STAR_HOST_SYS}/lib
export STAR_BIN=\${STAR}/.\${STAR_HOST_SYS}/bin
export PATH=\${STAR_BIN}:\${STAR}/mgr:\$PATH
export LD_LIBRARY_PATH=\${STAR_LIB}:\$LD_LIBRARY_PATH
cd "$PROJECT_ROOT_REAL"
root4star -b -q 'tools/dumpFemtoPhiTreeInfo.C("$TREE_FILE","$OUT_DIR")'
EOF
)
exec singularity exec \
  -B /gpfs:/gpfs -B /star/u:/star/u -B /star/nfs4/AFS:/star/nfs4/AFS -B /home/starlib:/home/starlib \
  --env LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}" \
  "$SINGULARITY_IMAGE" /bin/bash -lc "$CONTAINER_CMD"
