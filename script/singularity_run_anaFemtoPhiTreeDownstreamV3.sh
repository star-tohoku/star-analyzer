#!/bin/bash
# Run tree downstream SE/ME from a reduced ROOT file inside singularity.
# Usage: ./script/singularity_run_anaFemtoPhiTreeDownstreamV3.sh MAINCONF TREE.root OUT.root [bufferSize] [mixingMode] [variation] [sigMin] [sigMax] [recomputePid] [recomputeMixBin] [maxEvents] [pairMtBins] [pairMtMin] [pairMtMax] [data005DiagnosticsOnly] [data006Config] [firstEvent] [threeEventEnabled] [threeEventMaxRawSamplesPerEvent] [threeEventSamplingSeed]

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
PROJECT_ROOT_REAL="$(cd "$PROJECT_ROOT" && pwd -P)"
PYTHON="$(command -v python3 2>/dev/null || command -v python 2>/dev/null)"
SINGULARITY_IMAGE="/cvmfs/singularity.opensciencegrid.org/star-bnl/star-sw:latest"

resolve_from_project() {
  "$PYTHON" - "$PROJECT_ROOT_REAL" "$1" <<'PY'
import os
import sys

project_root = sys.argv[1]
path = sys.argv[2]
full = path if os.path.isabs(path) else os.path.join(project_root, path)
print(os.path.realpath(full))
PY
}

resolve_star_host_sys() {
  local star_root="$1"
  local candidate
  local candidates=(sl73_x8664_gcc485 sl74_x8664_gcc485 sl73_gcc485 sl74_gcc485)
  for candidate in "${candidates[@]}"; do
    if [[ -d "$star_root/.$candidate" ]]; then
      printf '%s\n' "$candidate"
      return 0
    fi
  done
  return 1
}

MAINCONF="${1:?Usage: $0 MAINCONF TREE OUT [bufferSize] [mixingMode] [variation] [sigMin] [sigMax] [recomputePid] [recomputeMixBin] [maxEvents] [pairMtBins] [pairMtMin] [pairMtMax] [data005DiagnosticsOnly] [data006Config] [firstEvent] [threeEventEnabled] [threeEventMaxRawSamplesPerEvent] [threeEventSamplingSeed]}"
TREE_FILE="${2:?}"
OUT_FILE="${3:?}"
BUFFER="${4:--1}"
MODE="${5:-bufferAll}"
VARIATION="${6:-}"
SIGMIN="${7:--1}"
SIGMAX="${8:--1}"
RECOPID="${9:-0}"
RECOMIX="${10:-0}"
MAXEVENTS="${11:--1}"
PAIR_MT_BINS="${12:-240}"
PAIR_MT_MIN="${13:-0.8}"
PAIR_MT_MAX="${14:-3.2}"
DATA005_DIAGNOSTICS_ONLY="${15:-0}"
DATA006_CONFIG="${16:-}"
FIRST_EVENT="${17:-0}"
THREE_EVENT_ENABLED="${18:-0}"
THREE_EVENT_MAX_RAW_SAMPLES="${19:-0}"
THREE_EVENT_SAMPLING_SEED="${20:-0}"

if [[ -z "$PYTHON" ]]; then
  echo "ERROR: python3 or python is required." >&2
  exit 1
fi
if ! command -v singularity >/dev/null 2>&1; then
  echo "ERROR: singularity command not found." >&2
  exit 1
fi

MAINCONF_REAL="$(resolve_from_project "$MAINCONF")"
if [[ ! -f "$MAINCONF_REAL" ]]; then
  echo "ERROR: mainconf not found: $MAINCONF_REAL" >&2
  exit 1
fi
if [[ ! -f "$TREE_FILE" ]]; then
  echo "ERROR: input tree not found: $TREE_FILE" >&2
  exit 1
fi
if [[ -n "$DATA006_CONFIG" ]]; then
  DATA006_CONFIG_REAL="$(resolve_from_project "$DATA006_CONFIG")"
  if [[ ! -f "$DATA006_CONFIG_REAL" ]]; then
    echo "ERROR: DATA-006 config not found: $DATA006_CONFIG_REAL" >&2
    exit 1
  fi
else
  DATA006_CONFIG_REAL=""
fi

mkdir -p "$(dirname "$OUT_FILE")"
LIBRARY_TAG=$(cd "$PROJECT_ROOT_REAL" && "$PYTHON" script/analysis_info_helper.py --library-tag --mainconf "$MAINCONF_REAL" | xargs)
STAR_ROOT="/star/nfs4/AFS/star/packages/$LIBRARY_TAG"
if [[ ! -d "$STAR_ROOT" ]]; then
  echo "ERROR: STAR package directory not found: $STAR_ROOT" >&2
  exit 1
fi
if ! RESOLVED_STAR_HOST_SYS="$(resolve_star_host_sys "$STAR_ROOT")"; then
  echo "ERROR: could not resolve an SL7 STAR_HOST_SYS under $STAR_ROOT." >&2
  exit 1
fi
CONTAINER_CMD=$(cat <<EOF
module load root-5.34.38 mysql-5.6.43 libiconv-1.16 libxml2-2.9.13 log4cxx-0.10.0 gsl-2.7.1 fastjet-3.3.4 rave-2020-08-11 genfit-b496504a-root-5.34.38 kitrack-01-10-star-1-root-5.34.38 vc-0.7.4
export STAR=$STAR_ROOT
export STAR_HOST_SYS=$RESOLVED_STAR_HOST_SYS
export STAR_LIB=\${STAR}/.\${STAR_HOST_SYS}/lib
export STAR_BIN=\${STAR}/.\${STAR_HOST_SYS}/bin
export PATH=\${STAR_BIN}:\${STAR}/mgr:\$PATH
export LD_LIBRARY_PATH=\${STAR_LIB}:$PROJECT_ROOT_REAL/lib:\$LD_LIBRARY_PATH
cd "$PROJECT_ROOT_REAL"
root4star -b -q 'tools/run_anaFemtoPhiTreeDownstreamV3.C("$TREE_FILE","$OUT_FILE","$MAINCONF_REAL",$BUFFER,"$MODE","$VARIATION",$SIGMIN,$SIGMAX,$RECOPID,$RECOMIX,$MAXEVENTS,$PAIR_MT_BINS,$PAIR_MT_MIN,$PAIR_MT_MAX,$DATA005_DIAGNOSTICS_ONLY,"$DATA006_CONFIG_REAL",$FIRST_EVENT,$THREE_EVENT_ENABLED,$THREE_EVENT_MAX_RAW_SAMPLES,$THREE_EVENT_SAMPLING_SEED)'
EOF
)

# /cvmfs carries the STAR spack modules loaded above. It is automounted on the login nodes but
# has to be bound explicitly for the same wrapper to work inside a Condor job.
CVMFS_BIND=()
if [[ -d /cvmfs ]]; then CVMFS_BIND=(-B /cvmfs:/cvmfs); fi

exec singularity exec \
  -B /gpfs:/gpfs \
  ${CVMFS_BIND[@]+"${CVMFS_BIND[@]}"} \
  -B /star/u:/star/u \
  -B /star/nfs4/AFS:/star/nfs4/AFS \
  -B /home/starlib:/home/starlib \
  "$SINGULARITY_IMAGE" \
  /bin/bash -lc "$CONTAINER_CMD"
