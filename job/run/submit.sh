#!/bin/bash
# Submit from job/run/. Replaces __PROJECT_ROOT__ in the XML with the project root.
# Usage: cd job/run && ./submit.sh [joblist.xml]
# Example: ./submit.sh ../joblist/joblist_run_anaLambda.xml

set -e
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

"$SCRIPT_DIR/setup_symlinks.sh"

PROJECT_ROOT="$(cd .. && cd .. && pwd)"
TEMPLATE="${1:-../joblist/joblist_run_anaLambda.xml}"
OUTPUT="$(basename "$TEMPLATE")"

if [[ ! -f "$TEMPLATE" ]]; then
  echo "ERROR: Template not found: $TEMPLATE" >&2
  exit 1
fi

sed "s|__PROJECT_ROOT__|$PROJECT_ROOT|g" "$TEMPLATE" > "$OUTPUT"
echo "Submitting with PROJECT_ROOT=$PROJECT_ROOT"

# Capture star-submit output to get jobid for configlog
submit_out=$(mktemp)
trap "rm -f '$submit_out'" EXIT
star-submit "$OUTPUT" 2>&1 | tee "$submit_out"
submit_rc=${PIPESTATUS[0]}

# anaName from joblist filename (e.g. joblist_auau3p85fxt_anaPhi.xml -> auau3p85fxt_anaPhi)
anaName="${OUTPUT#joblist_}"
anaName="${anaName%.xml}"

# Extract jobid: from "Wrote scheduling report to : .../anaNameJOBID.report"
jobid=""
if report_line=$(grep "Wrote scheduling report to" "$submit_out" 2>/dev/null); then
  report_path=$(echo "$report_line" | sed 's/.*: *//' | tr -d '\r' | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
  if [[ -n "$report_path" ]]; then
    report_basename=$(basename "$report_path" .report)
    jobid="${report_basename#$anaName}"
  fi
fi
# Fallback: 32-char hex from "Writting files for process HEX_XX..done."
if [[ -z "$jobid" ]]; then
  jobid=$(grep -oE '[0-9A-Fa-f]{32}' "$submit_out" 2>/dev/null | head -1)
fi

if [[ -n "$jobid" ]]; then
  "$SCRIPT_DIR/snapshot_config.sh" "$anaName" "$jobid" || true
  "$SCRIPT_DIR/snapshot_joblist.sh" "$anaName" "$jobid" "$OUTPUT" || true
else
  echo "WARNING: Could not extract jobid; skipping configlog and joblistlog." >&2
fi

exit "$submit_rc"
