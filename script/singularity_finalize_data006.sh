#!/bin/bash
# Convert a hadd-merged DATA-006 downstream ROOT file into ROOT/CSV/QA products.
# Usage: ./script/singularity_finalize_data006.sh MAINCONF MERGED.root OUTDIR DATA006_CONFIG PRODUCTION_TAG GIT_COMMIT

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd -P)"
PYTHON="$(command -v python3 2>/dev/null || command -v python 2>/dev/null)"
SINGULARITY_IMAGE="/cvmfs/singularity.opensciencegrid.org/star-bnl/star-sw:latest"

MAINCONF="${1:?Usage: $0 MAINCONF MERGED.root OUTDIR DATA006_CONFIG PRODUCTION_TAG GIT_COMMIT}"
MERGED_ROOT="${2:?}"
OUT_DIR="${3:?}"
DATA006_CONFIG="${4:?}"
PRODUCTION_TAG="${5:?}"
GIT_COMMIT="${6:?}"

if [[ -z "$PYTHON" ]]; then
  echo "ERROR: python3 or python is required." >&2
  exit 1
fi
if ! command -v singularity >/dev/null 2>&1; then
  echo "ERROR: singularity command not found." >&2
  exit 1
fi

resolve_path() {
  "$PYTHON" - "$PROJECT_ROOT" "$1" <<'PY'
import os, sys
p = sys.argv[2]
print(os.path.realpath(p if os.path.isabs(p) else os.path.join(sys.argv[1], p)))
PY
}

MAINCONF_REAL="$(resolve_path "$MAINCONF")"
MERGED_REAL="$(resolve_path "$MERGED_ROOT")"
OUT_REAL="$(resolve_path "$OUT_DIR")"
DATA006_REAL="$(resolve_path "$DATA006_CONFIG")"
for required in "$MAINCONF_REAL" "$MERGED_REAL" "$DATA006_REAL"; do
  if [[ ! -f "$required" ]]; then
    echo "ERROR: required input not found: $required" >&2
    exit 1
  fi
done

LIBRARY_TAG=$(cd "$PROJECT_ROOT" && "$PYTHON" script/analysis_info_helper.py --library-tag --mainconf "$MAINCONF_REAL" | xargs)
STAR_ROOT="/star/nfs4/AFS/star/packages/$LIBRARY_TAG"
RESOLVED_STAR_HOST_SYS=""
for candidate in sl73_x8664_gcc485 sl74_x8664_gcc485 sl73_gcc485 sl74_gcc485; do
  if [[ -d "$STAR_ROOT/.$candidate" ]]; then RESOLVED_STAR_HOST_SYS="$candidate"; break; fi
done
if [[ -z "$RESOLVED_STAR_HOST_SYS" ]]; then
  echo "ERROR: no SL7 STAR_HOST_SYS found under $STAR_ROOT" >&2
  exit 1
fi

CONTAINER_CMD=$(cat <<EOF
module load root-5.34.38 mysql-5.6.43 libiconv-1.16 libxml2-2.9.13 log4cxx-0.10.0 gsl-2.7.1 fastjet-3.3.4 rave-2020-08-11 genfit-b496504a-root-5.34.38 kitrack-01-10-star-1-root-5.34.38 vc-0.7.4
export STAR=$STAR_ROOT
export STAR_HOST_SYS=$RESOLVED_STAR_HOST_SYS
export STAR_LIB=\${STAR}/.\${STAR_HOST_SYS}/lib
export STAR_BIN=\${STAR}/.\${STAR_HOST_SYS}/bin
export PATH=\${STAR_BIN}:\${STAR}/mgr:\$PATH
export LD_LIBRARY_PATH=\${STAR_LIB}:$PROJECT_ROOT/lib:\$LD_LIBRARY_PATH
cd "$PROJECT_ROOT"
root4star -b -q 'tools/run_finalizeData006.C("$MERGED_REAL","$DATA006_REAL","$OUT_REAL","$PRODUCTION_TAG","$GIT_COMMIT")'
EOF
)

exec singularity exec -B /gpfs:/gpfs -B /star/u:/star/u -B /star/nfs4/AFS:/star/nfs4/AFS \
  -B /home/starlib:/home/starlib "$SINGULARITY_IMAGE" /bin/bash -lc "$CONTAINER_CMD"
