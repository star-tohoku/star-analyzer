#!/bin/bash
# Run common/macro/candidate_femto/checkHistCandidateFemto.C (StCandidateFemtoMaker QA / CF PDF)
# inside the same singularity context as batch jobs.
# Usage: ./script/singularity_checkHistCandidateFemto.sh <root_file> <mainconf> <out_pdf> <centGroups>
#   root_file   ROOT file written by StCandidateFemtoMaker (hPairSE_<channel> / hPairME_<channel>)
#   mainconf    main config with the candidateFemto key (species, channels, windows, cfNormMin/Max)
#   out_pdf     output PDF; a ROOT file with the CFs is written next to it (<out_pdf without .pdf>.root)
#   centGroups  cent9 groups of the centrality pages, ';' between groups, ',' inside (e.g. "8;6,7;4,5"; required)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
PROJECT_ROOT_REAL="$(cd "$PROJECT_ROOT" && pwd -P)"
PYTHON="$(command -v python3 2>/dev/null || command -v python 2>/dev/null)"

ROOT_FILE="${1:-}"
MAINCONF="${2:-}"
OUT_PDF="${3:-}"
CENT_GROUPS="${4:-}"

if [[ -z "$ROOT_FILE" || -z "$MAINCONF" || -z "$OUT_PDF" || -z "$CENT_GROUPS" ]]; then
  echo "Usage: $0 <root_file> <mainconf> <out_pdf> <centGroups>" >&2
  echo "  root_file   ROOT file written by StCandidateFemtoMaker" >&2
  echo "  mainconf    main config (e.g. config/mainconf/main_<anaName>.yaml)" >&2
  echo "  out_pdf     output PDF path (a .root with the CFs is written next to it)" >&2
  echo "  centGroups  cent9 groups, e.g. \"8;6,7;4,5\"" >&2
  exit 1
fi
if [[ "$OUT_PDF" != *.pdf ]]; then
  echo "ERROR: out_pdf must end in .pdf: $OUT_PDF" >&2
  exit 1
fi

# Absolute path (relative paths are relative to the current directory, as usual for a CLI).
abs_path() {
  "$PYTHON" -c 'import os, sys; print(os.path.realpath(sys.argv[1]))' "$1"
}

ROOT_FILE_REAL="$(abs_path "$ROOT_FILE")"
MAINCONF_REAL="$(abs_path "$MAINCONF")"
OUT_PDF_REAL="$(abs_path "$OUT_PDF")"

if [[ ! -f "$ROOT_FILE_REAL" ]]; then
  echo "ERROR: ROOT file not found: $ROOT_FILE_REAL" >&2
  exit 1
fi
if [[ ! -f "$MAINCONF_REAL" ]]; then
  echo "ERROR: mainconf not found: $MAINCONF_REAL" >&2
  exit 1
fi
mkdir -p "$(dirname "$OUT_PDF_REAL")"

LIBRARY_TAG=$(cd "$PROJECT_ROOT_REAL" && "$PYTHON" script/analysis_info_helper.py --library-tag --mainconf "$MAINCONF_REAL") || exit 1
LIBRARY_TAG="$(echo "$LIBRARY_TAG" | xargs)"

source "$SCRIPT_DIR/setup.sh" "$MAINCONF_REAL"
export LD_LIBRARY_PATH="$PROJECT_ROOT_REAL/lib:$LD_LIBRARY_PATH"

# Stale ACLiC products next to the macro (the runner builds in tmp/kf-aclic/<host>-<pid>).
rm -f "$PROJECT_ROOT_REAL"/common/macro/candidate_femto/checkHistCandidateFemto_C.* \
      "$PROJECT_ROOT_REAL"/common/macro/candidate_femto/checkHistCandidateFemto_C_ACLiC_dict* 2>/dev/null || true

echo "[singularity_checkHistCandidateFemto] root file: $ROOT_FILE_REAL"
echo "[singularity_checkHistCandidateFemto] mainconf:  $MAINCONF_REAL"
echo "[singularity_checkHistCandidateFemto] out pdf:   $OUT_PDF_REAL"
echo "[singularity_checkHistCandidateFemto] centGroups: $CENT_GROUPS"

CONTAINER_CMD=$(cat <<EOF
export STAR=/star/nfs4/AFS/star/packages/$LIBRARY_TAG
export STAR_HOST_SYS=$STAR_HOST_SYS
export STAR_LIB=\${STAR}/.\${STAR_HOST_SYS}/lib
export STAR_BIN=\${STAR}/.\${STAR_HOST_SYS}/bin
export PATH=\${STAR_BIN}:\${STAR}/mgr:\$PATH
export LD_LIBRARY_PATH=\${STAR_LIB}:\$LD_LIBRARY_PATH
cd "$PROJECT_ROOT_REAL"
root4star -b -q 'analysis/run_checkHistCandidateFemto.C("$ROOT_FILE_REAL","$MAINCONF_REAL","$OUT_PDF_REAL","$CENT_GROUPS")'
EOF
)

singularity exec \
  -B /gpfs:/gpfs \
  -B /star/nfs4/AFS:/star/nfs4/AFS \
  -B /home/starlib:/home/starlib \
  -B /direct:/direct \
  -B /star:/star \
  --env LD_LIBRARY_PATH="$LD_LIBRARY_PATH" \
  /cvmfs/singularity.opensciencegrid.org/star-bnl/star-sw:latest \
  /bin/bash -lc "$CONTAINER_CMD"
