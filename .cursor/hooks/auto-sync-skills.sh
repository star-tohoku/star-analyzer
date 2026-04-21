#!/usr/bin/env bash
set -euo pipefail

INPUT_JSON="$(cat)"
export HOOK_INPUT_JSON="$INPUT_JSON"

# Exit quickly unless docs/ai/skills markdown was edited.
SHOULD_SYNC="$(
python3 - <<'PY'
import json
import os
import re

pattern = re.compile(r"docs/ai/skills/[^\"'\s]+\.md$")
raw = os.environ.get("HOOK_INPUT_JSON", "")

try:
    data = json.loads(raw)
except Exception:
    print("0")
    raise SystemExit(0)

def walk(obj):
    if isinstance(obj, dict):
        for v in obj.values():
            yield from walk(v)
    elif isinstance(obj, list):
        for v in obj:
            yield from walk(v)
    elif isinstance(obj, str):
        yield obj

matched = any(pattern.search(s.replace("\\\\", "/")) for s in walk(data))
print("1" if matched else "0")
PY
)"

if [ "$SHOULD_SYNC" != "1" ]; then
  exit 0
fi

script/sync_and_check_skills.sh

echo '{ "additional_context": "Auto-synced .cursor/skills wrappers after docs skill update." }'
