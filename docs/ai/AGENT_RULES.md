# Agent rules

This document is the source of truth for project-specific agent workflow rules.
Framework-level design principles remain in `../../PHILOSOPHY.md`.

## Editing and implementation rules

- **File Modification Approval**: The agent MUST NOT modify, create, or delete any source, configuration, or documentation files in the repository without explicit prior approval from the user.
- **Macros**: Do not implement custom event loops. Build `StChain`, run `Init()`, `Make(i)`, and `Finish()`. The Maker-using macro must be ACLiC-compiled and linked through `AddLinkedLibs`.
- **Makers**: Implement `Init()`, `Make()`, `Clear()`, and `Finish()`. Keep histogram lifecycle in `DeclareHistograms()` and `WriteHistograms()` / `Finish()`.
- **Config**: Prefer YAML updates over recompilation when changing cuts and histogram values. Main entry is `config/mainconf/main_<anaName>.yaml`.
- **No hardcoded analysis parameters**: Put thresholds and cuts in YAML configs.
- **analysis_info sync**: Keep `anaName`, macro base names, output naming, and `starTag` consistent.
- **Build environment for debug/repro**: For debug work and any build used for farm submission, use a **batch-matched STAR toolchain** before submit. **Either** (a) interactive SL7: `sl7`, then `source ./script/setup.sh <mainconf>` and `make`, **or** (b) from any host with Singularity: **`./script/singularity_make.sh <mainconf>`** — option (b) runs `make` inside `star-bnl/star-sw:latest` with the same `sl73_*` / `sl74_*` `STAR_HOST_SYS` resolution as batch and **counts as the SL7-class build** for agents and reproducibility. Do not rely on host-only `make` for farm-bound libraries or heap/exit-crash diagnosis. For `root4star` / QA on fragile login nodes, use the other **`script/singularity_*`** wrappers in `../../docs/REFERENCE.md`.

## Batch safety: how values reach a job, and what happens when they do not

These rules exist because three defects of the same shape reached a 241-job farm pilot on
2026-09-15, each invisible to every local test
(`analysisnote/auau3p85fxt_anaFemtoPhi/root-tree-study/results/pilot-farm-20260915.md`).

- **A value the analysis needs arrives as an argument, not through the environment.** The
  environment may *override* for local convenience; it may not *supply*. Locally the wrapper
  scripts export `STAR_ANA_*`; the batch joblist calls `root4star` directly and exports nothing, so
  anything read only from the environment is empty on the farm — and empty was silently accepted
  for the job id, for the mainconf path, and for every tree schema key in turn.
- **Never continue with a default when a configuration cannot be read.** `if (missing) return
  kTRUE;` and "fall back to `.current_mainconf`" both turn a broken job into a job that produces
  plausible, wrong output. A job that cannot find its configuration must fail.
- **No built-in default mainconf.** A macro that is given no configuration stops.
- **Resolve relative paths against `gSystem->WorkingDirectory()`, never `gSystem->Getenv("PWD")`.**
  A SUMS job is a `csh` script that `cd`s into its runtime bundle, and csh does not update `PWD`.
- **Print the resolved configuration path and where it came from, unconditionally**, so that a job
  running with the wrong configuration is visible in the first lines of its log rather than in the
  output months later.
- **Required keys are checked, not defaulted.** A YAML read of the form "use the key if present,
  else keep the default" makes a typo silent. Keys that decide what an output *contains* (schema
  version, which species are stored, the storage envelope) must be present or the maker fails.
- **Provenance that was not supplied is recorded as `unrecorded`, never as an empty string,** and a
  warning is printed. An empty `gitRevision` reads like "no revision", which is what the pilot
  trees carried.

## Adding analyses and scripts

- **New analysis** requires: Maker code under `StMaker/StXXXMaker/` (auto-built; no Makefile edit), `run_anaXxx.C` that Loads/links `libStCommon` before the Maker `.so`, `anaXxx.C`, `script/run_anaXxx.sh`, mainconf, referenced YAMLs, and analysis_info.
- **Joblist generation**: Use `script/generate_joblist.sh <mainconf>` and produce `job/joblist/joblist_<anaName>.xml` (`analysis.anaName`); the batch runtime mainconf is the one passed to `generate_joblist.sh`.
- **Skill sync after docs edits**: If a task adds/updates `docs/ai/skills/*.md`, run `script/sync_and_check_skills.sh` before finishing.
- **Scripts and workflow docs**: When script or run flow changes, update:
  - `../../docs/REFERENCE.md`
  - `../../INSTALL.md` if onboarding changed
  - `../../README.md` when index/table changes
  - `../../job/run/README.md` for submit/cleanup/archive details

## StRoot vendoring

- **When to vendor:** only when `$STAR/lib` is insufficient and the code cannot live in `StMaker/`; see [`docs/ai/skills/reuse-star-stroot.md`](skills/reuse-star-stroot.md).
- **Canonical example:** [`StRoot/StRefMultCorr/`](../../StRoot/StRefMultCorr/) with `PROVENANCE.md`, `Makefile` target, and repo-local `.so` load in `run_ana*.C` (not `$STAR/lib` duplicate).
- **Agent procedure:** follow `reuse-star-stroot.md` when searching `$STAR/StRoot` or past trees and copying into `StRoot/<Package>/`.

## Centrality work

- **Bin conventions:** [`StRoot/StRefMultCorr/README.md`](../../StRoot/StRefMultCorr/README.md) is canonical (cent9: larger index = more central; 0–60% → cent9 2–8).
- **Agent procedure:** follow [`docs/ai/skills/centrality-strefmultcorr.md`](skills/centrality-strefmultcorr.md) when changing centrality YAML, `CentralityHelper`, cent axes, or `cfCent9Min`/`Max`.

## Femto φ–p checkHist / CF

- **CF and sideband logic live in** `common/macro/checkHistAnaFemtoPhiProton.C` (not in Maker). Skill: [`docs/ai/skills/femto-species-naming.md`](skills/femto-species-naming.md).
- **Dual PDF:** `singularity_checkHistAnaFemtoPhiProton.sh` writes QA PDF + `{anaName}_checkHistAnaFemtoPhiProton_CF_{jobid}.pdf` under `share/figure/<anaName>/`.
- **YAML (maker):** `cfCentSlices` (defaults in `FemtoConfig::SetDefaults`), `cfCentSlicesQaPdfInclude`, `cfPdfExcludeQaSlices`, `sidebandSubtractAlpha`, `sidebandAlphaMode`, `negativeBinPolicy`.
- **YAML (mixing):** `mixingMode`, `maxMixedPairsPerEvent`, `mixBothDirections`, `bufferSize` for ME statistics (`StFemtoMaker::FillMixedEventPairs`). Mixing bin = same Vz × **same cent9** × EP; ME is cross-event only.
- **SB-LR:** no new maker channel; checkHist `combineSidebandLR` on projected leftSB + rightSB histograms.
- **Sub CF keys (checkHist cache):** `CF_sig_raw`, `CF_sig_sub_SBL`, `CF_sig_sub_SBR`, `CF_sig_sub_SBLR` per slice id.
- **Plan / ops record:** `analysisnote/YYYYMMDD/femto_cent_sb_cf_plan.md` (when present).

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
- Submitted joblist snapshot: `job/run/joblistlog/joblist_<anaName>_<jobid>.xml`
- Submit-time runmeta manifest: `job/run/runmeta/runmeta_<anaName>_<jobid>.json`
- Submit-time runmeta sidecars: `job/run/runmeta/gitstatus_<anaName>_<jobid>.txt`, `gitdiff_<anaName>_<jobid>.patch`, `gitsubmodules_<anaName>_<jobid>.txt`, `runtime_bundle_<anaName>_<jobid>.tar.gz`, `sums_artifacts_<anaName>_<jobid>.tar.gz`, `submit_stdout_<anaName>_<jobid>.txt`
- Watch-merge log: `job/run/watchmerge/watchmerge_<anaName>_<jobid>.pid` / `.log`; completion metadata in runmeta `postProcess.watchMerge`
- QA PDF: `share/figure/<anaName>/<anaName>_checkHistAnaPhi[_<jobid>].pdf`

## Submit reproducibility artifacts

- Treat `jobid` as the key used to join config, submitted XML, runmeta, output ROOT, and QA.
- `submit.sh` should save per-`jobid` reproducibility artifacts immediately after jobid extraction, not as a later manual step.
- `runmeta` is the machine-readable index for submit-time provenance; use it before scanning loose `job/run/<anaName><jobid>*` files.
