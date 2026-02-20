# STAR analysis framework

STAR PicoDst-based analysis using the StChain/StMaker pattern: analysis macros drive the event loop, and custom Makers (compiled as shared libraries) perform the physics logic.

## Philosophy

- **No custom event loop in macros.** The macro builds a `StChain`, calls `chain->Init()`, then runs `chain->Make(i)` in a loop and `chain->Finish()`. The framework controls the loop.
- **Physics in StMaker subclasses.** Each analysis has a Maker (e.g. `StPhiMaker`) that implements `Init()`, `Make()`, `Clear()`, `Finish()`. Histograms are created in `DeclareHistograms()` and written in `WriteHistograms()`.
- **Makers are compiled.** Maker code lives in `StMaker/StXXXMaker/` and is built into `lib/libStXXXMaker.so`. This keeps heavy logic out of the interpreter and allows reuse.
- **Macros run under root4star.** The entry point is a ROOT macro invoked with `root4star -b -q "..."`. A small shell script (`script/run_anaXxx.sh`) sets the environment and calls the macro.
- **YAML-driven config.** Cuts and histogram definitions are read from YAML via `ConfigManager`, so you can change them without recompiling Makers when using existing cut/hist keys.
- **No hardcoding.** Do not hardcode concrete values in code; put all cut values, thresholds, and analysis parameters in YAML config so they can be changed without recompiling.
- **One config per concern.** Use one maker config per StMaker (referenced from mainconf). Use one cut config per cut type/source (e.g. one config for event cuts, one for track, one for PID, one for v0, one for mixing); mainconf references each.
- **Mainconf as the single entry point.** The main config (`config/mainconf/main_<anaName>.yaml`) references all other configs (paths relative to `config/`). Analysis should be fully reproducible using only the files referenced there. Scripts and workflows take mainconf as the primary argument so that setup, execution, I/O, and outputs are tied to one mainconf.

## After cloning

Do the following once after you clone the repository.

1. **Submodules** (required for build):
   ```bash
   git submodule update --init --recursive
   ```
   Without this, `make` will fail when building the config library.

2. **Directories not in git** — The framework expects these at runtime; they are git-ignored. Create them if you use batch submission or local ROOT output:
   - **log/** — Batch job stdout (e.g. `log/stdout.$JOBID.out`).
   - **err/** — Batch job stderr (e.g. `err/stderr.$JOBID.err`).
   - **rootfile/** — Analysis output ROOT files (local runs and batch copy-back). Use subdirs per analysis, e.g. `rootfile/auau19_anaPhi/`.
   - **lib/** — Filled by `make` (shared libraries). No need to create by hand.
   - **share/figure/** — QA PDFs from `checkHistAnaPhi.sh` (e.g. `share/figure/<anaName>/`). Created by the macro when needed.

   Example (create log, err, and rootfile so batch jobs can write outputs):
   ```bash
   mkdir -p log err rootfile
   ```

3. **Your analysis info and config** — Create analysis-specific config from the template:
   - Copy `config/analysis/analysis_info_temp.yaml` to `config/analysis/analysis_info_<anaName>.yaml`. For **anaName** (e.g. `{system}_{anaId}[_condition]`) see [Analysis info](#analysis-info-analysis_info_tempyaml) and [Naming conventions](#naming-conventions). Edit at least **analysis.workDir** and other keys as needed.
   - Run `script/setup_config_from_analysisinfo.py` with your analysis_info so that configs and mainconf for that analysis are created:
     ```bash
     python script/setup_config_from_analysisinfo.py config/analysis/analysis_info_<anaName>.yaml
     ```
     This creates the mainconf and cut/maker configs referenced by it (e.g. `config/mainconf/main_<anaName>.yaml`).

4. **PicoDst file list (optional)** — Once mainconf exists, you can build an input file list using the FileCatalog (`get_file_list.pl`):
   ```bash
   python script/generate_picodst_list.py config/mainconf/main_<anaName>.yaml
   ```
   Use `--dry-run` to print the command and output path without running. Requires STAR environment. Output path is taken from analysis_info (e.g. `config/picoDstList/<anaName>.list`). See the script docstring for details.

5. **Setup and build** — From the project root (use the mainconf created in step 3):
   ```bash
   ./script/setup.sh config/mainconf/main_<anaName>.yaml
   make
   ```

6. **Run analysis** — With joblist and mainconf you can run locally or submit batch jobs. See [How to run](#how-to-run) (local) and [Batch (star-submit)](#batch-star-submit) (batch).

7. **Job submission** — When submitting jobs, see [Batch (star-submit)](#batch-star-submit) and [Creating a joblist](#creating-a-joblist).

8. **Cursor (optional)** — To develop with the agent, open this directory in Cursor as the project folder (see [Development with Cursor](#development-with-cursor-agent-ai)).


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
| **script/** | Environment and run scripts: `setup.sh` (starver from analysis info), `generate_joblist.sh` (joblist XML from mainconf), `run_anaLambda.sh`, `run_anaPhi.sh`, `checkHistAnaPhi.sh` (QA PDF from run_anaPhi output ROOT), `analysis_info_helper.py` (libraryTag + joblist generation), `time_NYT_to_JST.py` (NY time → JST), `time_now_NY_to_JST.py` (current NY server time → JST), and helpers (e.g. `get_file_list_*.sh`). |

## Development with Cursor (agent AI)

When developing or modifying this framework, open the repository in [Cursor](https://cursor.com/) as the **project directory** (File → Open Folder → select the cloned repo root). The project includes rules and skills so that the agent follows this README.

- **After cloning**, open the cloned directory in Cursor as the project folder. The agent will then have full context of the codebase and the conventions below.
- **`.cursor/rules/`** — Project rules (e.g. README philosophy) are applied automatically. The agent is instructed to follow the Philosophy above (no custom event loop, physics in StMaker, YAML config, no hardcoding, one config per concern, mainconf as entry point). Do not remove or ignore `.cursor` if you want the agent to respect these conventions.
- **`.cursor/skills/`** — For common tasks (adding a new analysis, adding histograms to a StMaker, daily analysis log, updating README for scripts), the project provides skills; asking the agent for those tasks will trigger the corresponding workflow.
- If you change the **Philosophy** section in this README, update `.cursor/rules/README-philosophy.mdc` so the rule stays in sync.

## Prerequisites and setup

- **STAR environment**: use `starver` (e.g. `starver SL24y`). The version is taken from **analysis info** (see below), not hardcoded.
- **Build tools**: ROOT, gcc, and CMake (for yaml-cpp). The framework is intended to be buildable with any starver environment (subject to future development).

**Setup** reads the main config (mainconf) and the analysis info YAML it points to, then runs `starver` with the `libraryTag` from that file. From the project root, pass the mainconf path as the first argument:

```bash
source script/setup.sh config/mainconf/main_auau19_anaLambda.yaml
```

For another analysis (e.g. Phi):

```bash
source script/setup.sh config/mainconf/main_auau19_anaPhi.yaml
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
- **`--generate-joblist`**: Reads mainconf → analysis info and fills **job/joblist/job_template_from_conf.xml**, then writes e.g. `job/joblist/joblist_run_anaLambda.xml`. Requires PyYAML; use the **Python 3 environment** above.

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
source script/setup.sh config/mainconf/main_auau19_anaLambda.yaml
make
```

This builds `lib/libStarAnaConfig.so`, `lib/libStPhiMaker.so`, and `lib/libStLambdaMaker.so`. The Makefile uses `$STAR` and `root-config`; other Makers need their own targets (see "Adding a new analysis" below).

## analysis code (run_anaXxx.C and anaXxx.C)

ROOT's CINT interpreter does **not** resolve symbols from dynamically loaded shared libraries. So a single macro that does `gSystem->Load("libStPhiMaker.so")` and then `new StPhiMaker(...)` fails with "declared but not defined."

The fix is to **compile** the macro that uses the Maker with ACLiC (`.L anaPhi.C+`) and link it against the Maker library. That requires:

1. **Before** compilation: load STAR libs, load `libStXXXMaker.so`, and tell ACLiC to link against it (`gSystem->AddLinkedLibs(...)`).
2. **Then** compile and load the analysis macro (`.L anaXxx.C+`).
3. **Then** call the analysis function (e.g. `anaPhi(...)`).

Because `root4star -q` accepts only **one** macro invocation, that sequence is implemented in a **runner** macro, and the actual chain/event loop lives in a **second** macro that gets compiled:

- **run_anaXxx.C** — Loads STAR and project libs, sets `AddLinkedLibs`, compiles `anaXxx.C+`, and calls `anaXxx(...)`.
- **anaXxx.C** — Builds `StChain`, adds `StPicoDstMaker` and your Maker, runs the event loop. This file is compiled by ACLiC and linked to `libStXXXMaker.so`.

So each new analysis uses this two-file pattern on purpose; the runner is analysis-specific and not shared.


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
   source script/setup.sh config/mainconf/main_auau19_anaLambda.yaml
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
