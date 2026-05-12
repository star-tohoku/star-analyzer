# Agent rules

This document is the source of truth for project-specific agent workflow rules.
Framework-level design principles remain in `../../PHILOSOPHY.md`.

## Editing and implementation rules

- **Macros**: Do not implement custom event loops. Build `StChain`, run `Init()`, `Make(i)`, and `Finish()`. The Maker-using macro must be ACLiC-compiled and linked through `AddLinkedLibs`.
- **Makers**: Implement `Init()`, `Make()`, `Clear()`, and `Finish()`. Keep histogram lifecycle in `DeclareHistograms()` and `WriteHistograms()` / `Finish()`.
- **Config**: Prefer YAML updates over recompilation when changing cuts and histogram values. Main entry is `config/mainconf/main_<anaName>.yaml`.
- **No hardcoded analysis parameters**: Put thresholds and cuts in YAML configs.
- **analysis_info sync**: Keep `anaName`, macro base names, output naming, and `starTag` consistent.
- **Build environment for debug/repro**: For debug work and any build used for farm submission, run `make` inside an SL7 environment first (use `sl7`, then `source ./script/setup.sh <mainconf>` and `make`). Do not rely on host-only builds for heap/exit-crash diagnosis.

## Adding analyses and scripts

- **New analysis** requires: Maker code, Makefile target, `run_anaXxx.C`, `anaXxx.C`, `script/run_anaXxx.sh`, mainconf, referenced YAMLs, and analysis_info.
- **Joblist generation**: Use `script/generate_joblist.sh <mainconf>` and produce `job/joblist/joblist_<anaName>.xml` (`analysis.anaName`); the batch runtime mainconf is the one passed to `generate_joblist.sh`.
- **Skill sync after docs edits**: If a task adds/updates `docs/ai/skills/*.md`, run `script/sync_and_check_skills.sh` before finishing.
- **Scripts and workflow docs**: When script or run flow changes, update:
  - `../../docs/REFERENCE.md`
  - `../../INSTALL.md` if onboarding changed
  - `../../README.md` when index/table changes
  - `../../job/run/README.md` for submit/cleanup/archive details

## Naming and terminology

- **anaName**: `{system}_{anaId}[_condition]` and use it consistently across configs, outputs, and job artifacts.
- **jobid**: 32-char hex from star-submit, used in configlog/output naming.
- **mainconf**: Canonical configuration entry point under `config/mainconf/`.
- **baseRunMacro/baseAnaMacro**: Names without `.C`; joblist and ACLiC output names derive from them.

## Directory and path conventions

- `analysis/`: paired `run_anaXxx.C` and `anaXxx.C`
- `StMaker/`: compiled Maker source to `lib/libStXXXMaker.so`
- `config/`: concern-separated YAMLs; mainconf references others relative to `config/`
- `script/`: setup, local run, joblist generation helpers
- `job/run/`: submit and run-artifact management tools
- `lib/`: generated binaries, not committed

## Output conventions

- ROOT output: `rootfile/<anaName>/`
- Submitted config snapshot: `job/run/configlog/config_<anaName>_<jobid>.txt`
- QA PDF: `share/figure/<anaName>/<anaName>_checkHistAnaPhi[_<jobid>].pdf`
