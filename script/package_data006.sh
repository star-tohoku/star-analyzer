#!/bin/bash
# Finalize and provenance-package a merged DATA-006 result. Refuses to overwrite OUTDIR.
# Usage: ./script/package_data006.sh MAINCONF MERGED.root OUTDIR [productionTag] [data006Config]

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd -P)"

MAINCONF="${1:?Usage: $0 MAINCONF MERGED.root OUTDIR [productionTag] [data006Config]}"
MERGED_ROOT="${2:?}"
OUT_DIR="${3:?}"
PRODUCTION_TAG="${4:-P24iy:production_3p85GeV_fixedTarget_2021}"
DATA006_CONFIG="${5:-config/maker/data006_auau3p85fxt_pp_pd_dd.yaml}"

if [[ -e "$OUT_DIR" ]]; then
  echo "ERROR: output already exists; refusing to overwrite: $OUT_DIR" >&2
  exit 1
fi
GIT_COMMIT=$(cd "$PROJECT_ROOT" && git rev-parse HEAD)

"$PROJECT_ROOT/script/singularity_finalize_data006.sh" \
  "$MAINCONF" "$MERGED_ROOT" "$OUT_DIR" "$DATA006_CONFIG" "$PRODUCTION_TAG" "$GIT_COMMIT"

mkdir -p "$OUT_DIR/provenance/config" "$OUT_DIR/provenance/source"
cp "$PROJECT_ROOT/$DATA006_CONFIG" "$OUT_DIR/provenance/config/"
cp "$PROJECT_ROOT/$MAINCONF" "$OUT_DIR/provenance/config/"

# Snapshot every direct config referenced by the mainconf. Missing files are fatal.
while IFS= read -r rel; do
  [[ -z "$rel" ]] && continue
  src="$PROJECT_ROOT/config/$rel"
  if [[ ! -f "$src" ]]; then
    echo "ERROR: mainconf-referenced config not found: $src" >&2
    exit 1
  fi
  cp "$src" "$OUT_DIR/provenance/config/"
done < <(awk -F: '/^[[:space:]]*(event|track|pid|mixing|centrality|nuclearid|maker|hist|analysis)[[:space:]]*:/ {v=$2; sub(/#.*/,"",v); gsub(/^[[:space:]]+|[[:space:]]+$/,"",v); print v}' "$PROJECT_ROOT/$MAINCONF")

cp "$PROJECT_ROOT/include/Data006Config.h" "$OUT_DIR/provenance/source/"
cp "$PROJECT_ROOT/tools/anaFemtoPhiTreeDownstreamV3.C" "$OUT_DIR/provenance/source/"
cp "$PROJECT_ROOT/tools/finalizeData006.C" "$OUT_DIR/provenance/source/"
cp "$PROJECT_ROOT/tools/validateData006Package.C" "$OUT_DIR/provenance/source/"
cp "$PROJECT_ROOT/tools/plotData006Correlations.C" "$OUT_DIR/provenance/source/"
cp "$PROJECT_ROOT/tools/run_plotData006Correlations.C" "$OUT_DIR/provenance/source/"
cp "$PROJECT_ROOT/script/singularity_plot_data006.sh" "$OUT_DIR/provenance/source/"
cp "$PROJECT_ROOT/tools/run_anaFemtoPhiTreeDownstreamV3.C" "$OUT_DIR/provenance/source/"
cp "$PROJECT_ROOT/tools/run_finalizeData006.C" "$OUT_DIR/provenance/source/"
cp "$PROJECT_ROOT/docs/others/DATA-006-pp-pd-dd-source-calibration.md" "$OUT_DIR/REQUEST.md"
cp "$PROJECT_ROOT/docs/others/DATA-006-implementation.md" "$OUT_DIR/IMPLEMENTATION.md"

(cd "$PROJECT_ROOT" && git status --short) > "$OUT_DIR/provenance/git_status.txt"
(cd "$PROJECT_ROOT" && git diff --binary) > "$OUT_DIR/provenance/git_diff.patch"
printf '%s\n' "$GIT_COMMIT" > "$OUT_DIR/provenance/git_commit.txt"
printf '%s\n' "$PRODUCTION_TAG" > "$OUT_DIR/provenance/production_tag.txt"
MERGED_REAL=$(readlink -f "$MERGED_ROOT")
MERGED_SIZE=$(stat -c %s "$MERGED_REAL")
MERGED_SHA256=$(sha256sum "$MERGED_REAL" | awk '{print $1}')
printf 'role,path,size_bytes,sha256\nmerged_data006_raw,%s,%s,%s\n' \
  "$MERGED_REAL" "$MERGED_SIZE" "$MERGED_SHA256" > "$OUT_DIR/INPUTS.csv"
printf '%q ' "$PROJECT_ROOT/script/package_data006.sh" "$MAINCONF" "$MERGED_ROOT" "$OUT_DIR" "$PRODUCTION_TAG" "$DATA006_CONFIG" > "$OUT_DIR/REGENERATE.command"
printf '\n' >> "$OUT_DIR/REGENERATE.command"

cat > "$OUT_DIR/README.md" <<'EOF'
# DATA-006 p-p / p-d / d-d source-calibration package

`data006_correlations.root` and `correlations.csv` contain matched SE, ME and normalized CF
objects. Merged centrality classes are formed by summing native SE and ME yields before a new
normalization; finished CFs are never averaged. ROOT k* is GeV/c and CSV k* is MeV/c.

The baseline inherits the frozen phi reduced-tree event/centrality/PID decisions. Matter p and d
only are used. Identical SE pairs are unordered (`i<j`); identical ME uses current x buffer once.
The p-d channel uses both mixing directions. The baseline close-pair veto is disabled to match the
current phi production, while SE/ME delta-eta/delta-phi-star QA is retained. SE and ME constituent
kinematics plus normalization-region pair-rapidity overlays are included for acceptance closure.
PID purity, primary fraction, feed-down, deuteron material background, corrected Nch and their
systematic uncertainties are `NOT_AVAILABLE`; they are not silently estimated.

See `CLOSURE_STATUS.txt`, `invalid_bins.csv`, `object_manifest.csv`, and `provenance/` before fit
ingestion.
EOF

(cd "$OUT_DIR" && find . -type f ! -name SHA256SUMS -print0 | sort -z | xargs -0 sha256sum) > "$OUT_DIR/SHA256SUMS"
echo "DATA-006 package: $OUT_DIR"
cat "$OUT_DIR/CLOSURE_STATUS.txt"
