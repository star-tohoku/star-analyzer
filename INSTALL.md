# Install and first run

Step-by-step: **clone → configure → build → run locally** (optional: batch). For design background see [PHILOSOPHY.md](PHILOSOPHY.md). For command-line details, tables, and extending the framework see [docs/REFERENCE.md](docs/REFERENCE.md). A short overview is in [README.md](README.md).

---

## Step 0: Prerequisites

Confirm the following **before** cloning (each line: what / why).

| Requirement | What | Why |
|-------------|------|-----|
| **STAR stack** | You can run `starver` and use `root4star`. | Makers and macros assume STAR ROOT and libraries. |
| **Build** | `gcc`, `make`, CMake; `root-config` on `PATH` after `starver`. | Builds `libStarAnaConfig.so` and Maker `.so` files. |
| **Git** | Clone and submodules. | `yaml-cpp` lives in a submodule. |
| **Python 3** | `python3` available. | Used by `setup_config_from_analysisinfo.py` and optional list/joblist scripts. |
| **PyYAML** | `python3 -m pip install --user pyyaml` (if you lack it). | **Required** for `./script/generate_joblist.sh` (batch). Optional for local-only runs if you do not generate joblists. |

---

## Step 1: Clone

**What:** Get a copy of the repository. **Why:** Local tree for build and config.

```bash
git clone <repository-url>
cd <repository-directory>
```

Use your real URL and directory name instead of the placeholders above.

---

## Step 2: Submodules

**What:** Populate `src/third_party/yaml-cpp`. **Why:** Without it, `make` fails building `libStarAnaConfig.so`.

```bash
git submodule update --init --recursive
```

Or clone with: `git clone --recurse-submodules <repository-url>` and skip the command.

**Gotcha:** Forgetting this step is the most common first build failure.

---

## Step 3: Runtime directories

**What:** Create writable dirs that are git-ignored. **Why:** Batch and local jobs write logs and ROOT output here.

```bash
mkdir -p log err rootfile
```

- **log/** — job stdout (e.g. `log/stdout.$JOBID.out`).
- **err/** — job stderr.
- **rootfile/** — output ROOT files (use subdirs per analysis, e.g. `rootfile/auau19_anaLambda/`).
- **lib/** — produced by `make` (no need to create).
- **share/figure/** — QA PDFs from `checkHistAnaPhi.sh`; often created on demand.

---

## Step 4: analysis_info and mainconf

**What:** Point the framework at your STAR version, paths, and analysis name. **Why:** `setup.sh` and joblists read this metadata.

**Option A — use a mainconf already in the repo (e.g. Lambda / Phi)**

1. Edit `config/analysis/analysis_info_temp.yaml` (or whatever your mainconf’s `analysis:` key references).
2. Set at least **analysis.workDir** to your project root (e.g. `/star/u/$USER/work/star-analysis`).
3. Adjust **starTag.libraryTag** if your site needs a different `starver` tag.
4. Use that mainconf in Step 6 (e.g. `config/mainconf/main_auau19_anaLambda.yaml`).

**Option B — new analysis from the template**

1. Copy `config/analysis/analysis_info_temp.yaml` to `config/analysis/analysis_info_<anaName>.yaml`.
2. Edit `anaName`, `workDir`, macro names, and related keys (see [docs/REFERENCE.md](docs/REFERENCE.md) — Analysis info).
3. Ensure `config/mainconf/mainconf.yaml` exists and uses the `__ANANAME__` placeholder where analysis-specific names are required.
4. Generate configs:

```bash
python script/setup_config_from_analysisinfo.py config/analysis/analysis_info_<anaName>.yaml
```

This creates `config/mainconf/main_<anaName>.yaml` and referenced cut/maker YAMLs.

Replace `<anaName>` with your token (e.g. `auau19_anaLambda`).

---

## Step 5 (optional): PicoDst file list

**What:** Build an input `.list` from the STAR FileCatalog. **Why:** Local or batch runs need an input list.

```bash
python script/generate_picodst_list.py config/mainconf/main_<anaName>.yaml
```

Use `--dry-run` to print the command without running. Requires a working STAR environment. Output path comes from analysis_info (often `config/picoDstList/<anaName>.list`).

---

## Step 6: Setup and build

**What:** Run `starver` from analysis_info and compile libraries. **Why:** `make` needs `$STAR` and ROOT from the correct release.

From the **project root**, run:

```bash
./script/setup.sh config/mainconf/main_<anaName>.yaml
make
```

**Note:** This project standardizes on `./script/setup.sh <mainconf>` for setup commands.

---

## Step 7: Run locally

**What:** Execute one analysis interactively or in batch mode on the login node. **Why:** Validate build and config before farm submission.

Example (Lambda, first 100 events):

```bash
./script/run_anaLambda.sh config/picoDstList/auau19GeV.list rootfile/auau19_anaLambda_temp/out.root 0 100
```

Pass a fifth argument to override mainconf (see [docs/REFERENCE.md](docs/REFERENCE.md) — How to run). For Phi, use `./script/run_anaPhi.sh` and the matching list/mainconf.

The run scripts set **`LD_LIBRARY_PATH`** (and related env) before `root4star`. If you invoke `root4star` by hand, you must reproduce that environment or linking may fail.

---

## Step 8 (optional): Batch (star-submit)

**What:** Generate a joblist XML and submit. **Why:** Process many files on the farm.

Requires **Python 3 + PyYAML** (Step 0).

```bash
./script/generate_joblist.sh config/mainconf/main_<anaName>.yaml
cd job/run
./submit.sh ../joblist/joblist_run_anaLambda.xml
```

Use the XML name your generator produced (`joblist_<baseRunMacro>.xml`). After submit, see [job/run/README.md](job/run/README.md) for `configlog`, `cleanup_job_run.sh`, and `archive_job_run.sh`.

---

## Troubleshooting (short)

| Symptom | Things to check |
|---------|-------------------|
| `make` fails in yaml-cpp / config lib | Submodule: Step 2. |
| Wrong STAR / missing `root-config` | Re-run **`./script/setup.sh ...`** (Step 6) before `make`. |
| Library load errors at runtime | Run via **`script/run_anaXxx.sh`** or match its `LD_LIBRARY_PATH` setup. |
| Joblist script errors | Install PyYAML for the same `python3` you use. |
| Batch paths wrong | **analysis.workDir** and paths in analysis_info (Step 4). |

---

## Quick reference (existing Lambda mainconf)

After editing `config/analysis/analysis_info_temp.yaml` (**analysis.workDir**):

```bash
git submodule update --init --recursive
mkdir -p log err rootfile
./script/setup.sh config/mainconf/main_auau19_anaLambda.yaml
make
./script/run_anaLambda.sh config/picoDstList/auau19GeV_lambda.list rootfile/auau19_anaLambda_temp/out.root 0 -1
# Optional batch:
./script/generate_joblist.sh config/mainconf/main_auau19_anaLambda.yaml
cd job/run && ./submit.sh ../joblist/joblist_run_anaLambda.xml
```

---

## More information

- [README.md](README.md) — Overview and documentation index.
- [PHILOSOPHY.md](PHILOSOPHY.md) — Design principles and two-macro rationale.
- [docs/REFERENCE.md](docs/REFERENCE.md) — Full reference (analysis_info, batch detail, adding analyses, config, joblists, QA).
- [job/run/README.md](job/run/README.md) — Submit directory and cleanup.
