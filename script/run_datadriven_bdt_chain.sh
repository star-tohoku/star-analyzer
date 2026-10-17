#!/bin/bash
# Data-driven BDT for one recipe, from the candidate trees to the report: inspect, train, apply, study, report.
# Stops at the first failing step. The model file for StCandidateFemtoMaker is written only after you have read the report
# and set "validation" in the recipe (the command is printed at the end).
# Usage: ./script/run_datadriven_bdt_chain.sh <recipe.yaml>
# Procedure and rules: docs/ai/skills/candidate-bdt.md, docs/ai/candidate-bdt-guide.md.
set -euo pipefail
[[ $# -eq 1 ]] || { sed -n '2,6p' "$0" >&2; exit 1; }
RUN="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/run_datadriven_bdt.sh"
for step in inspect train apply study report; do
  echo "=== $step"
  "$RUN" "$1" "$step"
done
echo "=== done. Read the report (output directory of the recipe), set validation: <report path> in the recipe, then:"
echo "    $RUN $1 model <model dir, e.g. config/bdt/<anaName>>"
