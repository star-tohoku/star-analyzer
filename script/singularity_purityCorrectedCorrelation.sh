#!/bin/bash
# k*-differential purity-corrected CF of one StCandidateFemtoMaker channel, in the batch singularity context.
# Usage: ./script/singularity_purityCorrectedCorrelation.sh <root_file> <mainconf> <channel> <out_prefix> <kstarRebin> [cent9List=""] [options=""]
#   channel     e.g. k0s_xi_minus, kaon_plus_xi_minus (channel names of the candidateFemto YAML)
#   kstarRebin  fine k* bins (kstarAxis of the candidateFemto YAML) per CF bin, e.g. 2 (required)
#   cent9List   e.g. "4,5,6,7,8"; "" = all centralities
#   options     ';' separated: purityA=config|fit|unit; windowA=lo,hi; scoreMinA=x; scoreMinB=x; slideHalfWidth=N
# Writes <out_prefix>.pdf and <out_prefix>.root.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT_REAL="$(cd "$SCRIPT_DIR/.." && pwd -P)"
PYTHON="$(command -v python3 2>/dev/null || command -v python 2>/dev/null)"
if [[ $# -lt 5 ]]; then
  sed -n 2,8p "$0" >&2
  exit 1
fi
ROOT_FILE_REAL="$(realpath "$1")"
MAINCONF_REAL="$(realpath "$2")"
CHANNEL="$3"
mkdir -p "$(dirname "$4")"
OUT_PREFIX="$(realpath -m "$4")"
KSTAR_REBIN="$5"
[[ "$KSTAR_REBIN" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: kstarRebin must be a positive integer: $5" >&2; exit 1; }
CENT_LIST="${6:-}"
OPTIONS="${7:-}"
[[ -f "$ROOT_FILE_REAL" ]] || { echo "ERROR: ROOT file not found: $1" >&2; exit 1; }
[[ -f "$MAINCONF_REAL" ]] || { echo "ERROR: mainconf not found: $2" >&2; exit 1; }
LIBRARY_TAG=$(cd "$PROJECT_ROOT_REAL" && "$PYTHON" script/analysis_info_helper.py --library-tag --mainconf "$MAINCONF_REAL" | xargs)
source "$SCRIPT_DIR/setup.sh" "$MAINCONF_REAL"
export LD_LIBRARY_PATH="$PROJECT_ROOT_REAL/lib:$LD_LIBRARY_PATH"
CONTAINER_CMD=$(cat <<EOC
export STAR=/star/nfs4/AFS/star/packages/$LIBRARY_TAG
export STAR_HOST_SYS=$STAR_HOST_SYS
export STAR_LIB=\${STAR}/.\${STAR_HOST_SYS}/lib
export STAR_BIN=\${STAR}/.\${STAR_HOST_SYS}/bin
export PATH=\${STAR_BIN}:\${STAR}/mgr:\$PATH
export LD_LIBRARY_PATH=\${STAR_LIB}:\$LD_LIBRARY_PATH
cd "$PROJECT_ROOT_REAL"
root4star -b -q 'analysis/run_purityCorrectedCorrelation.C("$ROOT_FILE_REAL","$MAINCONF_REAL","$CHANNEL","$OUT_PREFIX",$KSTAR_REBIN,"$CENT_LIST","$OPTIONS")'
EOC
)
singularity exec -B /gpfs:/gpfs -B /star/nfs4/AFS:/star/nfs4/AFS -B /direct:/direct -B /star:/star \
  --env LD_LIBRARY_PATH="$LD_LIBRARY_PATH" --pwd "$PROJECT_ROOT_REAL" \
  /cvmfs/singularity.opensciencegrid.org/star-bnl/star-sw:latest /bin/bash -lc "$CONTAINER_CMD"
