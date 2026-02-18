---
name: add-new-analysis
description: Adds a new analysis (StMaker, macros, run script, config, joblist) following the project pattern. Use when the user wants to add a new analysis, create a new StMaker, or add an anaXxx (e.g. new physics channel or dataset).
---

# Add new analysis

Each analysis has: a Maker, a shared library, two macros (run_anaXxx.C + anaXxx.C), a run script, mainconf and referenced configs, and analysis_info. Use this checklist and refer to README.md "Adding a new analysis" and "Creating your own config" for detail.

## Before you start

Read existing StMaker and macros as reference so the new analysis follows the same patterns and signatures. For example:
- **StMaker**: `StMaker/StPhiMaker/`, `StMaker/StLambdaMaker/` (header and source: Init, Make, Clear, Finish, DeclareHistograms, WriteHistograms, config usage).
- **Macros**: `analysis/run_anaPhi.C`, `analysis/anaPhi.C` (or run_anaLambda.C, anaLambda.C). Check how the runner loads libs, compiles the analysis macro, and calls it; how the analysis macro loads config and builds the chain.
- **Run script**: `script/run_anaPhi.sh` or `script/run_anaLambda.sh`.

Then proceed with the checklist below.

## Checklist

1. **Maker**
   - Create `StMaker/StXXXMaker/StXXXMaker.h` and `StXXXMaker.cxx`. Subclass `StMaker`; implement `Init()`, `Make()`, `Clear()`, `Finish()`.
   - Histograms: create in `DeclareHistograms()`, write in `WriteHistograms()` (e.g. in `Finish()`).
   - Use `ConfigManager::GetInstance()` for cuts and `GetHistConfigPath()` for hist YAML if the Maker uses config. Use `StPicoDstMaker` to access PicoDst.

2. **Build**
   - Add a Makefile target for `lib/libStXXXMaker.so` (same pattern as `libStPhiMaker.so` / `libStLambdaMaker.so`). Link `libStarAnaConfig` if the Maker uses config.
   - Run `make` so `lib/libStXXXMaker.so` exists.

3. **Runner macro**
   - Copy `analysis/run_anaPhi.C` (or `run_anaLambda.C`) to `analysis/run_anaXxx.C`.
   - Replace: function name, `libStXXXMaker.so`, macro `anaXxx.C+`, call `anaXxx(...)`. Keep `libStarAnaConfig.so` and pass config path as 5th argument if used.

4. **Analysis macro**
   - Copy `analysis/anaPhi.C` or `analysis/anaLambda.C` to `analysis/anaXxx.C`.
   - Load main config: `ConfigManager::GetInstance().LoadConfig(mainConfigPath)` (from 5th argument or default mainconf path).
   - Replace Maker type and variable names; keep the same signature (including optional 5th `configPath`).

5. **Run script**
   - Copy `script/run_anaPhi.sh` or `script/run_anaLambda.sh` to `script/run_anaXxx.sh`.
   - Point to `run_anaXxx.C` and set default input, output, and config paths.

6. **Config**
   - **mainconf**: Copy an existing `config/mainconf/main_*.yaml` to `config/mainconf/main_<anaName>.yaml`. Set keys: event, track, pid, v0, mixing, maker key (e.g. phi/lambda), hist, analysis (path to analysis_info).
   - **analysis_info**: Create or copy `config/analysis/analysis_info_*.yaml`; set anaName, workDir, baseRunMacro, baseAnaMacro, mainConf, jobName, scratchSubdir, outputFileStem, nFiles, starTag. Use YAML alias so anaName stays in sync.
   - **maker/hist/cuts**: Add or copy YAMLs under `config/maker/`, `config/hist/`, `config/cuts/` as referenced by mainconf.

7. **Joblist**
   - Run `./script/generate_joblist.sh config/mainconf/main_<anaName>.yaml` to write `job/joblist/joblist_<baseRunMacro>.xml`.

8. **README**
   - Follow the **update-readme-scripts** skill if new scripts or run examples were added (e.g. run_anaXxx.sh in Directory layout / How to run).

## Naming

- **anaName**: `{system}_{anaId}[_condition]` (e.g. `auau19_anaLambda`, `auau3p85fxt_anaPhi`). Use consistently for mainconf, analysis_info, output paths, joblist.
- **baseRunMacro** / **baseAnaMacro**: e.g. `run_anaPhi`, `anaPhi`. Joblist file is `joblist_<baseRunMacro>.xml`; ACLiC builds `<baseAnaMacro>_C.so`.

## Reference

- Full steps and code patterns: README.md "Adding a new analysis (new StMaker)" and "Creating your own config".
