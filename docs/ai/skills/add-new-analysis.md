# Skill: add-new-analysis

Adds a new analysis with Maker, macros, run script, config, and joblist.

## Required outputs

1. `StMaker/StXXXMaker/` implementation (`.h` / `.cxx`)
2. Makefile target to produce `lib/libStXXXMaker.so`
3. `analysis/run_anaXxx.C` and `analysis/anaXxx.C`
4. `script/run_anaXxx.sh`
5. `config/mainconf/main_<anaName>.yaml` plus referenced YAML files
6. `config/analysis/analysis_info_*.yaml` with synchronized fields
7. `job/joblist/joblist_<anaName>.xml` from `script/generate_joblist.sh` (`analysis.anaName`)

## Rules

- Keep event loop in chain pattern only.
- Keep physics/cut logic in Maker, not runner macro.
- Prefer existing naming conventions and YAML aliases to keep `anaName` consistent.
- If scripts/workflow changed, follow `update-readme-scripts.md`.

## Reference

- `../../../docs/REFERENCE.md`
- `../../../PHILOSOPHY.md`
