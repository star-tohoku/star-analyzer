# Reference — STAR analysis framework

Operational detail: prerequisites, setup, analysis_info, build, running locally and in batch, QA, adding analyses, config, and joblists.

- **Overview and doc map:** [README.md](../README.md)
- **First-time setup:** [INSTALL.md](../INSTALL.md)
- **Design principles and two-macro rationale:** [PHILOSOPHY.md](../PHILOSOPHY.md)
- **Submit / cleanup on the farm:** [job/run/README.md](../job/run/README.md)

For **why** `run_anaXxx.C` and `anaXxx.C` are separate, see [PHILOSOPHY.md](../PHILOSOPHY.md) (two macros section).

## Directory layout

Only **template/sample** content under `config/` and `job/joblist/` is tracked; built artifacts under `lib/` are git-ignored.

| Directory | Description |
|-----------|-------------|
| **analysis/** | ROOT macros: `run_anaXxx.C` (runner: loads libs, compiles `anaXxx.C+`, calls analysis) and `anaXxx.C` (StChain + event loop). One pair per analysis (e.g. Lambda, Phi). |
| **config/** | YAML configs. **Templates/samples only** tracked. Subdirs: `mainconf/` (main YAML that includes the rest), `maker/`, `hist/`, `cuts/` (event, track, pid, v0reco, mixing), `analysis/` (e.g. **analysis_info_temp.yaml** — used by setup.sh and joblist generator), `picoDstList/` (input file lists; user lists are typically untracked). |
| **include/** | Framework headers: `ConfigManager.h`, `HistManager.h`, cut configs (`cuts/*.h`). Used by StMaker and `src/`. |
| **job/** | Job submission: `job/joblist/` = **template** job XMLs (tracked); `job/run/` = submit directory (`submit.sh`, `cleanup_job_run.sh`, `archive_job_run.sh`, generated/copied files). On successful submit, joblist and config are saved to `job/run/joblistlog/joblist_<anaName>_<jobid>.xml` and `job/run/configlog/config_<anaName>_<jobid>.txt`. Files under `job/run/*.xml` and SUMS outputs are git-ignored. |
| **lib/** | Built shared libraries (`libStarAnaConfig.so`, `libStXXXMaker.so`). **Contents git-ignored**; produced by `make`. |
| **StMaker/** | One subdir per Maker (e.g. `StLambdaMaker/`, `StPhiMaker/`). Each has `.h` and `.cxx`; built into `lib/libStXXXMaker.so`. |
| **script/** | Environment and run scripts: `setup.sh` (starver from analysis info), `generate_joblist.sh` (joblist XML from mainconf), `run_anaLambda.sh`, `run_anaPhi.sh`, `checkHistAnaPhi.sh` (QA PDF from run_anaPhi output ROOT), `analysis_info_helper.py` (libraryTag + joblist generation), `sync_cursor_skills.py` / `check_cursor_skill_sync.py` (sync and validate `docs/ai/skills/*.md` ↔ `.cursor/skills/*/SKILL.md` parity), `sync_and_check_skills.sh` (single command to run sync + parity check), `time_NYT_to_JST.py` (NY time → JST), `time_now_NY_to_JST.py` (current NY server time → JST), and helpers (e.g. `get_file_list_*.sh`). |

## Prerequisites and setup

- **STAR environment**: use `starver` (e.g. `starver SL24y`). The version is taken from **analysis info** (see below), not hardcoded.
- **Build tools**: ROOT, gcc, and CMake (for yaml-cpp). The framework is intended to be buildable with any starver environment (subject to future development).

**Setup** reads the main config (mainconf) and the analysis info YAML it points to, then runs `starver` with the `libraryTag` from that file. From the project root, pass the mainconf path as the first argument:

```bash
./script/setup.sh config/mainconf/main_auau19_anaLambda.yaml
```

For another analysis (e.g. Phi):

```bash
./script/setup.sh config/mainconf/main_auau19_anaPhi.yaml
```

If the helper script or the analysis info is missing or invalid, `setup.sh` will exit with an error.

## Python 3 environment (for joblist generation)

Generating a joblist with `./script/generate_joblist.sh` uses **Python 3** and **PyYAML**. Prepare once:

1. **Python 3**  
   Use the system or module `python3` (e.g. `module load python` if your site provides it).

2. **PyYAML**  
   Install for Python 3 (user install if you lack admin rights):
   ```bash
   python3 -m pip install --user pyyaml
   ```

3. **Run the script**  
   From the project root:
   ```bash
   ./script/generate_joblist.sh config/mainconf/main_auau19_anaLambda.yaml
   ```
   The script prefers `python3` when available and writes e.g. `job/joblist/joblist_run_anaLambda.xml`.

## Analysis info (analysis_info_temp.yaml)

The file **config/analysis/analysis_info_temp.yaml** (or the one referenced by your mainconf’s `analysis:` key) holds metadata used by `setup.sh` and by the **joblist generator**. After cloning, copy or edit this file so that paths and tags match your environment.

**Main keys:**

| Section | Keys | Purpose |
|--------|------|--------|
| **starTag** | `libraryTag`, `triggerSets`, `productionTag`, `filetype`, `filenameFilter`, `storageExclude` | `setup.sh` uses `libraryTag` for `starver`. The rest are used to build the SUMS catalog URL when generating a joblist. |
| **dataset** | `allPicoDstList`, `runRange`, etc. | Dataset description; can be used by run scripts or docs. |
| **analysis** | `anaName`, `name`, `workDir`, `baseRunMacro`, `baseAnaMacro`, `mainConf`, `jobName`, `scratchSubdir`, `outputFileStem`, `nFiles` | **anaName** is the canonical name for this analysis (see Naming conventions below). Use YAML alias: define `anaName: &anaName "auau19_anaLambda_temp"` and set `name`, `jobName`, `scratchSubdir`, `outputFileStem` to `*anaName` so they stay in sync. **baseRunMacro** / **baseAnaMacro** are the macro base names without `.C` (e.g. `run_anaLambda`, `anaLambda`); the joblist generator builds the run macro as `baseRunMacro + ".C"` and writes `joblist_<baseRunMacro>.xml`. **workDir** is the base path for log/err/output. |
| **analyst** | `name`, `institute`, `email` | For documentation. |

**Script that uses it:** `script/analysis_info_helper.py`

- **`--library-tag`**: Reads mainconf → analysis info and prints `starTag.libraryTag`. Used by `setup.sh`. Works without PyYAML (minimal grep fallback).
- **`--generate-joblist`**: Reads mainconf → analysis info and fills **job/joblist/job_template_from_conf.xml**, then writes e.g. **job/joblist/joblist_run_anaLambda.xml**. Requires PyYAML; use the **Python 3 environment** above.

## AI skill wrapper synchronization

To keep Antigravity (source docs) and Cursor wrappers aligned, treat `docs/ai/skills/*.md` as source-of-truth.

In Cursor, project hook config at `.cursor/hooks.json` runs `.cursor/hooks/auto-sync-skills.sh` after file edits and automatically syncs when a matching path under `docs/ai/skills/*.md` changes.

Single manual command:

```bash
script/sync_and_check_skills.sh
```

- `sync_and_check_skills.sh` runs sync then parity check.
- `sync_cursor_skills.py` writes/updates `.cursor/skills/<skill>/SKILL.md` wrappers from docs sources.
- `check_cursor_skill_sync.py` exits non-zero when wrappers and docs sources are out of sync (useful for CI/pre-commit).

### Naming conventions

- **anaName** = `{system}_{anaId}[_condition]` (e.g. `auau19_anaLambda`, `auau19_anaLambda_mid`). Define it once per analysis and tie output files, mainconf, and analysis_info to it in a 1:1 way. Use a YAML anchor in analysis_info: `anaName: &anaName "auau19_anaLambda_temp"` and reference it with `*anaName` for `name`, `jobName`, `scratchSubdir`, `outputFileStem`.
- **baseRunMacro** / **baseAnaMacro**: Macro base names without extension (e.g. `run_anaLambda`, `anaLambda`). The runner macro file is `baseRunMacro + ".C"`; the joblist file is `joblist_<baseRunMacro>.xml`. ACLiC builds `anaLambda_C.so` from the analysis macro.

## Submodules

This repository uses a git submodule for YAML parsing:

- **src/third_party/yaml-cpp** — [jbeder/yaml-cpp](https://github.com/jbeder/yaml-cpp); required to build `libStarAnaConfig.so`.

**After cloning**, initialize and update the submodule:

```bash
git submodule update --init --recursive
```

If you clone with recursive submodules, you can skip the above:

```bash
git clone --recurse-submodules <repository-url>
```

Without the submodule populated, `make` will fail when building the config library.

## Build

From the project root:

```bash
./script/setup.sh config/mainconf/main_auau19_anaLambda.yaml
make
```

This builds `lib/libStarAnaConfig.so`, `lib/libStPhiMaker.so`, and `lib/libStLambdaMaker.so`. The Makefile uses `$STAR` and `root-config`; other Makers need their own targets (see "Adding a new analysis" below).

## How to run

### Lambda analysis example (local with root4star)

Run from the **project root**:

```bash
./script/run_anaLambda.sh [inputList] [outputRoot] [jobid] [nEvents] [configPath]
```

Defaults:

- `inputList`   = `config/picoDstList/auau19GeV_lambda.list`
- `outputRoot`  = `rootfile/auau19_anaLambda_temp/auau19_anaLambda_temp.root`
- `jobid`       = `0`
- `nEvents`     = `-1` (all events)
- `configPath`  = (default main config; omit to use `config/mainconf/main_auau19_anaLambda.yaml`)

Example (first 100 events):

```bash
./script/run_anaLambda.sh config/picoDstList/auau19GeV_lambda.list rootfile/auau19_anaLambda_temp/out.root 0 100
```

With a custom main config (5th argument):

```bash
./script/run_anaLambda.sh config/picoDstList/auau19GeV_lambda.list rootfile/auau19_anaLambda_temp/out.root 0 -1 config/mainconf/main_auau19_anaLambda.yaml
```

The script sets `LD_LIBRARY_PATH` and runs:

```bash
root4star -b -q "analysis/run_anaLambda.C(\"$INPUT\",\"$OUTPUT\",\"$JOBID\",$NEVENTS[, \"$CONFIG\"])"
```

`run_anaLambda.C` loads STAR libs, `libStarAnaConfig.so`, and `libStLambdaMaker.so`, compiles `anaLambda.C+`, and calls `anaLambda(...)`.

### Result QA (Phi): checkHistAnaPhi.sh

After running the Phi analysis (locally or after merging batch output), you can produce a histogram QA PDF from the output ROOT file:

```bash
./script/checkHistAnaPhi.sh <root_file> <mainconf_path>
```

Example (single or merged ROOT file):

```bash
./script/checkHistAnaPhi.sh rootfile/auau3p85fxt_anaPhi/auau3p85fxt_anaPhi_CCBCC32EA67793F5A24B5F6BA44EE413_merge.root config/mainconf/main_auau3p85fxt_anaPhi.yaml
```

If the input filename has the form `anaName_jobid_merge.root` (32‑char hex jobid), the PDF is written as `share/figure/<anaName>/<anaName>_checkHistAnaPhi_<jobid>.pdf`; otherwise as `share/figure/<anaName>/<anaName>_checkHistAnaPhi.pdf`. When config is loaded, cut regions (event, track, phi) are overlaid on pre-cut histograms as red dashed lines.

### Batch (star-submit)

First-time flow (after git clone): customize analysis info → setup → build → generate joblist → submit.

1. **Set your paths in analysis info**  
   Edit **config/analysis/analysis_info_temp.yaml** (or the file your mainconf’s `analysis:` points to). At least set **analysis.workDir** to your project directory (e.g. `/star/u/$USER/star-analyzer`). This will be used as the base for log, err, and output in the generated joblist.

2. **Setup and build** (from project root):
   ```bash
   ./script/setup.sh config/mainconf/main_auau19_anaLambda.yaml
   make
   ```
   `setup.sh` reads the mainconf you pass and the analysis info, then runs `starver` with the `libraryTag` from that file.

3. **Generate the joblist** (from project root). Pass the mainconf path:
   ```bash
   ./script/generate_joblist.sh config/mainconf/main_auau19_anaLambda.yaml
   ```
   For another analysis (e.g. Phi):
   ```bash
   ./script/generate_joblist.sh config/mainconf/main_auau19_anaPhi.yaml
   ```
   This writes e.g. **job/joblist/joblist_run_anaLambda.xml** from **job/joblist/job_template_from_conf.xml** and the analysis info. Requires the **Python 3 environment** (see above).

4. **Submit** from `job/run/`:
   ```bash
   cd job/run
   ./submit.sh ../joblist/joblist_run_anaLambda.xml
   ```
   If your generated XML uses `__PROJECT_ROOT__`, `submit.sh` will replace it with the actual project path. Log and output URLs come from **analysis.workDir** in the analysis info. On successful submit, the mainconf and all referenced config files are written to **job/run/configlog/config_<anaName>_<jobid>.txt** (one text file per job) for reproducibility. See `job/run/README.md` for more.

5. **Cleaning up job/run**  
   After submission, `job/run/` is filled with many files named `anaName+jobid+*` (`.csh`, `.list`, etc.). To remove them (avoids "Argument list too long"): `cd job/run && ./cleanup_job_run.sh <anaName+jobid>`. To move them into an archive instead: `cd job/run && ./archive_job_run.sh <anaName+jobid>` (files go to `job/run/joblog/<anaName>/`; the directory is created if needed).

## Adding a new analysis (new StMaker)

Each analysis has:

1. A Maker in `StMaker/StXXXMaker/` (e.g. `StLambdaMaker`).
2. A shared library `lib/libStXXXMaker.so` (e.g. `libStLambdaMaker.so`).
3. Two macros: `analysis/run_anaXxx.C` and `analysis/anaXxx.C`.
4. A run script: `script/run_anaXxx.sh`.

### 1. Maker code

- Create `StMaker/StXXXMaker/StXXXMaker.h` and `StMaker/StXXXMaker/StXXXMaker.cxx`.
- Subclass `StMaker`; implement `Init()`, `Make()`, `Clear()`, `Finish()`.
- Create histograms in `DeclareHistograms()` and write them in `WriteHistograms()` (e.g. in `Finish()`).
- Use `ConfigManager::GetInstance()` for cuts and `GetHistConfigPath()` for the hist YAML path if your Maker uses config.
- Use `StPicoDstMaker` to access PicoDst; get the chain from it as needed.

### 2. Build the shared library

- In the **Makefile**, add a target for `lib/libStXXXMaker.so` (same pattern as `libStPhiMaker.so` / `libStLambdaMaker.so`). If the Maker uses config, link against `libStarAnaConfig`.
- Run `make` so that `lib/libStXXXMaker.so` exists.

### 3. Runner macro (run_anaXxx.C)

- Copy `analysis/run_anaPhi.C` (or `run_anaLambda.C`) to `analysis/run_anaXxx.C`.
- Replace the analysis name in: function name, library `libStXXXMaker.so`, macro `anaXxx.C+`, and call `anaXxx(...)`. If using config, keep loading `libStarAnaConfig.so` and pass the config path as the 5th argument where applicable.

### 4. Analysis macro (anaXxx.C)

- Copy `analysis/anaPhi.C` or `analysis/anaLambda.C` to `analysis/anaXxx.C`.
- Before building the chain, load the main config: `ConfigManager::GetInstance().LoadConfig(mainConfigPath)` (resolve `mainConfigPath` from the 5th argument or default to e.g. `config/mainconf/main_XXX.yaml`).
- Replace the Maker type and variable names; keep the same signature so the runner can call it (including optional 5th parameter `configPath` if used).

### 5. Run script (run_anaXxx.sh)

- Copy `script/run_anaPhi.sh` or `script/run_anaLambda.sh` to `script/run_anaXxx.sh`.
- Point to `run_anaXxx.C` and set default input, output, and config paths.

## Creating your own config

- **Main config**: Copy `config/mainconf/main_auau19_anaLambda.yaml` (or `main_auau19_anaPhi.yaml`) to e.g. `config/mainconf/main_myanalysis.yaml`. It references (paths are relative to `config/`):
  - **Cuts**: `event`, `track`, `pid`, `v0`, `mixing`, and optionally analysis-specific keys (e.g. `phi`, `lambda`).
  - **Maker**: e.g. `lambda: maker/maker_lambda.yaml`.
  - **Hist**: `hist: hist/hist_lambda.yaml`.
  - **Analysis info**: `analysis: analysis/analysis_info_temp.yaml` (or your own file). This file is used by `setup.sh` and by `script/analysis_info_helper.py --generate-joblist`.
- **Maker config**: Add e.g. `config/maker/maker_my.yaml` and reference it in the main config under the key your Maker expects. Makers read cuts via `ConfigManager::GetInstance().GetXXXCuts()` and the hist path via `GetHistConfigPath()`.
- **Hist config**: Add e.g. `config/hist/hist_my.yaml` with the same structure as existing hist YAMLs (`axes`, `histograms`). Set the `hist` key in the main config to this file.
- **New cut type**: If you need a new cut category, add a new key in the main YAML, a new `XxxCutConfig` in `include/cuts/` and `src/cuts/`, and register it in `ConfigManager`. For a new analysis that only uses existing event/track/pid/v0/mixing and maker keys, copying and editing the existing YAMLs under `config/cuts/`, `config/maker/`, and `config/hist/` is enough.

## Creating a joblist

**Recommended: generate from analysis info**

1. Ensure **config/analysis/analysis_info_temp.yaml** (or the file referenced by your mainconf’s `analysis:` key) has the keys described in **Analysis info** above, including `anaName`, `workDir`, `baseRunMacro`, `baseAnaMacro`, `mainConf`, `jobName`, `scratchSubdir`, `outputFileStem`, `nFiles`, and the `starTag` fields for the catalog URL.
2. From the project root, run (requires the **Python 3 environment**):
   ```bash
   ./script/generate_joblist.sh config/mainconf/main_myanalysis.yaml
   ```
   This fills **job/joblist/job_template_from_conf.xml** from the analysis info and writes e.g. **job/joblist/joblist_run_anaLambda.xml**.
3. Submit from `job/run`: `./submit.sh ../joblist/joblist_run_anaLambda.xml`.

**Manual edit**

- **Templates** in `job/joblist/`: **job_template_from_conf.xml** (for script-based generation) and concrete joblists (e.g. `joblist_run_anaLambda.xml`).
- To edit by hand: copy a joblist to a new name, then set **command** (macro, `$FILELIST`, `$SCRATCH` path, config path), **stdout/stderr/output toURL** (e.g. from your `workDir`), **input** (catalog URL, `nFiles`), and **SandBox** `File` entries as needed.
- Submit from `job/run`: `./submit.sh ../joblist/YourJoblist.xml`. Ensure `make` has been run so `lib/` contains the required `.so` files.
