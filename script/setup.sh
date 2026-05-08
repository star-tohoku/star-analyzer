#!/bin/bash
# Run from project root: ./script/setup.sh MAINCONF_PATH
# MAINCONF_PATH is required (e.g. config/mainconf/main_auau19_anaLambda.yaml).
# Runs starver with libraryTag from analysis info (mainconf -> analysis YAML).

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
MAINCONF="${1:?Usage: ./script/setup.sh MAINCONF_PATH (e.g. config/mainconf/main_auau19_anaLambda.yaml)}"

LIBRARY_TAG=$(cd "$PROJECT_ROOT" && python script/analysis_info_helper.py --library-tag --mainconf "$MAINCONF") || exit 1
LIBRARY_TAG=$(echo "$LIBRARY_TAG" | xargs)
if ! command -v starver >/dev/null 2>&1; then
  echo "ERROR: starver command not found. Load STAR environment first (for this project, run 'sl7' before setup/make)." >&2
  exit 1
fi
starver "$LIBRARY_TAG"
echo "LIBRARY_TAG: $LIBRARY_TAG"
