# star-submit run directory

Run `star-submit` from **this directory**. SUMS will write generated files (`.csh`, `.list`, `.condor`, `.report`, etc.) here.

## Steps

1. **Build at the project root** (batch-matched STAR toolchain — **either** path is valid)
   ```bash
   cd /path/to/star-analysis
   ./script/singularity_make.sh config/mainconf/main_auau19_anaLambda.yaml
   ```
   **Or** interactive SL7:
   ```bash
   cd /path/to/star-analysis
   sl7
   source ./script/setup.sh config/mainconf/main_auau19_anaLambda.yaml
   make
   ```
   `singularity_make.sh` runs `make` inside `star-bnl/star-sw:latest` with the same `sl73_*` / `sl74_*` toolchain as farm jobs; use it when you do not want to enter `sl7` manually (for example on AL9 login nodes). For `csh` / `tcsh`, use `source ./script/setup.csh ...` instead when building inside `sl7`.
   For debug/repro jobs (especially heap/exit issues), always use one of these batch-like builds before submit — not host-only `make` on a mismatched OS.

2. **Move to this directory and submit**
   ```bash
   cd job/run
   ./submit.sh ../joblist/joblist_auau19_anaLambda_test.xml
   ```
   Pass the joblist XML you want (for example the file produced by `./script/generate_joblist.sh`, named **`joblist_<anaName>.xml`** where **`anaName`** is `analysis.anaName` in your analysis_info). By default, `submit.sh` uses **`../joblist/joblist_auau19_anaLambda_temp.xml`** if you omit the argument (adjust the default in `submit.sh` if your primary analysis differs). Example for Phi: `./submit.sh ../joblist/joblist_auau19_anaPhi.xml`.

  `submit.sh` replaces `__PROJECT_ROOT__` in the template with the actual project path, so the same template works for any user.
  Before submit it runs a preflight check: it infers **`anaName`** from the joblist basename (`joblist_<anaName>.xml` → `<anaName>`), then extracts the embedded **`config/mainconf/...yaml`** path from the joblist and uses that as the single source of truth for `libraryTag` resolution, rebuilds, ELF class consistency of `lib/*.so` vs `root4star`, and runtime linker sanity in the singularity context. This fails fast on mismatch and avoids secondary errors like missing `fromScratch` output after an early crash.

  After preflight succeeds, **`submit.sh` creates the stdout, stderr, and ROOT output directories parsed from `<stdout>`, `<stderr>`, and `<output toURL>`** if they are missing, and verifies each directory is writable. Generated joblists use `analysis.workDir` for ROOT output and optional `analysis.logDir` / `analysis.errDir` for stdout/stderr. This prevents SUMS from scheduling a job that later fails when copying output.

   If preflight fails and you want the script to try recovery once, use:
   ```bash
   ./submit.sh --rebuild-if-needed ../joblist/joblist_auau19_anaLambda_test.xml
   ```
   This runs `source ./script/setup.sh <embedded-mainconf> && make`, then retries preflight with the same embedded mainconf. If host `make` is not batch-matched (for example AL9), run **`./script/singularity_make.sh <embedded-mainconf>`** at the project root first, then re-run `submit.sh` without relying on this recovery path.

  On successful submit, reproducibility artifacts are saved per `jobid`:
  - **joblistlog/joblist_<anaName>_<jobid>.xml** — submitted XML after `__PROJECT_ROOT__` replacement
  - **configlog/config_<anaName>_<jobid>.txt** — embedded mainconf plus all referenced YAMLs
  - **runmeta/runmeta_<anaName>_<jobid>.json** — machine-readable manifest tying together the saved artifacts below
  - **runmeta/gitstatus_<anaName>_<jobid>.txt** — `git status --porcelain=v2 --branch` at submit time
  - **runmeta/gitdiff_<anaName>_<jobid>.patch** — `git diff --binary HEAD` at submit time
  - **runmeta/gitsubmodules_<anaName>_<jobid>.txt** — `git submodule status --recursive` at submit time
  - **runmeta/runtime_bundle_<anaName>_<jobid>.tar.gz** — submit-time code/runtime bundle for replay (`analysis/`, `config/`, `include/`, `StMaker/`, `StRoot/`, `lib/`, and build/source metadata when present)
  - **runmeta/sums_artifacts_<anaName>_<jobid>.tar.gz** — stable snapshot of SUMS-generated `anaName+jobid+*` files such as `.list`, `.csh`, `.condor`, `.report`, and `.session.xml`
  - **runmeta/submit_stdout_<anaName>_<jobid>.txt** — captured `star-submit` console output

## Lambda–He4 reduced-tree jobs

`anaFemtoLambda_4He_tree` writes target-event/candidate/track information and
selected non-pair QA into the same producer ROOT. Its `reducedTree.mode: write`
stage does not make SE/ME or a mixing pool. The `replay` stage reads completed
producer ROOTs and performs SE/ME without PicoDst or KF reconstruction. It is a
separate format from the `candidateStore:` cache described below.

Use `config/mainconf/main_auau3p85_anaFemtoLambda_4He_tree.yaml` for production,
`main_auau3p85_anaFemtoLambda_4He_tree_replay.yaml` for original-mixing closure,
and `main_auau3p85_anaFemtoLambda_4He_tree_replay_ep1.yaml` for the one-EP-bin
control on the same stores. Inputs are explicit: only the producer uses the
covered 2021 P24iy PicoDst list. Replay must receive an ordered list of producer
ROOT files, never that PicoDst list.

`main_auau3p85_anaFemtoLambda_4He_tree_replay_pairY5.yaml` is a further replay
control with five rather than ten pair-rapidity bins and otherwise identical
EP8 mixing. Replay outputs retain the per-plane k* distributions; producer QA
does not acquire these pair histograms. These controls use the same stored
input and require no new reconstruction-stage job.

```bash
python3 script/generate_joblist_from_list.py \
  --mainconf config/mainconf/main_auau3p85_anaFemtoLambda_4He_tree.yaml \
  --input <covered-2021-P24iy-PicoDst-list> --files-per-subjob <n> \
  --label treeWrite --output-dir <new-producer-directory>
python3 script/generate_joblist_from_list.py \
  --mainconf config/mainconf/main_auau3p85_anaFemtoLambda_4He_tree_replay.yaml \
  --input <ordered-producer-ROOT-list> --files-per-subjob <n> \
  --label treeReplay --output-dir <new-pair-directory>
```

The producer has one ROOT output, so `--tree-output-dir` is not needed.
The standard writer/replay/EP1/pairY5 profiles use compact storage with
`storePrePidControl: false`. They retain the producer's 108 non-pair QA
histograms and selected-pair PID/rigidity QA; only nonselected pre-PID pair
controls are outside the saved envelope. Use
`main_auau3p85_anaFemtoLambda_4He_tree_fullControl.yaml` and
`main_auau3p85_anaFemtoLambda_4He_tree_replay_fullControl.yaml` with the explicit
true capability for full-control studies or existing full-control pilots.
Do not mix their size estimates with compact production estimates, and do not
ask a compact replay to recover controls that were not stored.
Rebuild with the batch-matched toolchain and verify write/read, one-pass
closure, mixing boundaries, merging geometry and measured quota before a full
production submission. The uncalibrated 81 runs are not enabled by these
entries. Producer QA remains associated with its full original input and is
not silently added once per replay condition or per split/repacked copy.
Use an ordered list of completed producer files for replay. Ordered hadd or
repack is supported only when complete producer segments and their metadata
remain contiguous. Arbitrary splitting within a segment is not supported.
Retain the original per-producer source/input manifests and source/config
checksums; a hadded file-level metadata object does not preserve all original
source-list mappings. An additional hadded tree copy must be included in the
quota estimate.

## Lambda–He4 candidate cache jobs

The `anaFemtoLambda_4He` entry supports candidate persistence and replay through
mainconf `candidateStore:`. The writer reads PicoDst and retains the original
reconstruction/PID/pairing; its ROOT output also contains `femtoCandidateEvents`
and layout, selection, original-mixing and file-UUID metadata. The replay stage
reads those exact candidates without PicoDst/KF reconstruction. This is the
`StFemtoMaker` cache format, separate from the reconstruction-part store used by
`StCandidateFemtoMaker`. The explicit cache species include d/t along with
Lambda/He3/He4: legacy d/t candidates affect original pool admission even when
only He3/He4 pair channels are enabled. Do not filter these collections out of
the cache or change their PID. Selection provenance also guards effective
dE/dx calibration content.

Build and run a local closure test before farm submission. The write/replay
mainconfs are
`main_auau3p85_anaFemtoLambda_4He_KFParticle_candidateStoreWrite.yaml` and
`main_auau3p85_anaFemtoLambda_4He_KFParticle_candidateStoreReplay.yaml` under
`config/mainconf/`. Both keep the production selection references. The standard
He4 dataset is 2021 P24iy; do not substitute archived 2019 validation input.

The [no-Sparse writer mainconf](../../config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle_candidateStoreWriteNoJointSparse.yaml)
keeps candidate writing and physical cuts unchanged. Its dedicated nuclear pair
QA YAML uses explicit `jointSparse: false`, with the full 8D definition commented
out. This removes only `hNuclearPidPairJoint`; dense QA and pre-PID candidate
caching remain enabled. A complete `jointSparse` mapping enables the 8D output;
missing, null or `true` values are errors.

From the project root, create explicit file-list jobs with separate destinations:

```bash
python3 script/generate_joblist_from_list.py \
  --mainconf <candidateStoreWrite-mainconf> --input <2021-P24iy-PicoDst-list> \
  --files-per-subjob <n> --label candidateStore --output-dir <new-cache-directory>
python3 script/generate_joblist_from_list.py \
  --mainconf <candidateStoreReplay-mainconf> --input <ordered-original-store-list> \
  --files-per-subjob <n> --label candidateReplay --output-dir <new-pair-directory>
```

Submit the reviewed XML from `job/run` in the normal way. The command retains
explicit input, output, job ID, event limit and mainconf. In replay, `nEvents`
counts **stored accepted events**, not producer PicoDst input events; `-1` reads
all store rows. Save the submitted runtime/config snapshots and the exact ordered
cache list for each variant.

Keep candidate caches as their **original completed per-job files**. Do not use
`--watch-merge`, hadd, repacking, or generic histogram cleanup on the cache stage:
ROOT UUID and producer pool boundaries form part of the replay contract.
Byte-for-byte copies retain the UUID, but duplicate UUIDs/files in one replay list
are rejected. Inventory completion and expected coverage before replay; scheduler
exit 0 alone does not establish a valid completed producer. Submit cache jobs
without `--watch-merge`, including when the 8D histogram is disabled. Use
`script/watch_job_notify.py --watch` for completion email without a merge. Start
it on the submit host after runmeta is available, in tmux or another persistent
session. tmux survives a client disconnect, not a submit-host reboot.

Use explicit `resetPoolAtFileBoundary: true` for original per-job closure. It
starts an empty pool at each original store file; grouping several files into one
replay job therefore does not accidentally mix across source-job boundaries.
Setting it to `false` deliberately continues across files, changes the ME
reference and must be recorded as a variation with a fixed input order.
`requireOriginalMixing: true` guards the producer's mixing configuration;
`false` permits a deliberate mixing override while selection/PID provenance stays
fixed. The supplied `candidateStoreRemixVz` mainconf varies only Vz bins 10→20.

The cache retains original `psi2` separately from auxiliary unflattened TPC Q2 and
validity. Storing Q2 does not change the writer's pairing event plane or establish
a calibrated physics event plane. Pair-only replay does not regenerate
PicoDst-dependent reconstruction/helix QA, so compare core pair histograms rather
than expecting all producer QA to be refilled. Replay histogram outputs can use
the normal guarded merge after completion; their input caches remain untouched.

With the full-QA configuration (`jointSparse` mapping enabled), the ROOT5
limited-statistics check stores/replays 3,858 accepted rows from 4,027
PicoDst events and reproduces 1,675 dense pair histograms plus one Sparse exactly,
including pre-PID QA. The Vz 10→20 variation leaves SE unchanged and changes the
signal He4 ME population from 251 to 235. This confirms the cache/replay contract
on the sample, not full-statistics physics behavior. One auxiliary TPC plane is
invalid, so EP-binned replay stops explicitly under the current validity rule.

The validation runner can reuse an explicit successful reconstruction record after
strict recorded input/config/library checks, without new PicoDst/KF commands;
original ROOT fingerprints and the comparison-macro change remain in the record.
The exact local pair comparator and stage commands are in the [reference](../../docs/REFERENCE.md#cached-lambdahe4-pairing),
with the full [Maker storage contract](../../StMaker/StFemtoMaker/README.md#lambdahe4-candidate-event-store-and-replay).

## Auto merge after batch completion (`--watch-merge`)

To poll subjob ROOT output and run `merge_root_files.csh` automatically when all subjobs finish:

```bash
cd job/run
./submit.sh --watch-merge ../joblist/joblist_auau3p85fxt_anaFemtoPhiProton.xml
```

Unified phi femto (all bachelors in one pass):

```bash
./submit.sh --watch-merge ../joblist/joblist_auau3p85fxt_anaFemtoPhi.xml
```

After merge, QA PDF:

```bash
./script/singularity_checkHistAnaFemtoPhi.sh \
  rootfile/auau3p85fxt_anaFemtoPhi/auau3p85fxt_anaFemtoPhi_<jobid>_merge.root \
  config/mainconf/main_auau3p85fxt_anaFemtoPhi.yaml
```

Unified K⁻ femto (all bachelors in one pass):

```bash
./submit.sh --watch-merge ../joblist/joblist_auau3p85fxt_anaFemtoKaon.xml
```

After merge, QA + CF PDF:

```bash
./script/singularity_checkHistAnaFemtoKaon.sh \
  rootfile/auau3p85fxt_anaFemtoKaon/auau3p85fxt_anaFemtoKaon_<jobid>_merge.root \
  config/mainconf/main_auau3p85fxt_anaFemtoKaon.yaml
```

- **`--watch-merge`**: after a successful submit, starts `script/watch_job_and_merge.sh` in the background (`nohup`).
- **`--watch-merge-foreground`**: same watcher, but blocks until merge finishes (debug).
- **Log**: `job/run/watchmerge/watchmerge_<anaName>_<jobid>.log`
- **Completion rule**: count subjob ROOT files under `rootfile/<anaName>/` (excluding `*_merge.root`) until it matches the number of SUMS `{anaName}{jobid}_*.list` files in `job/run/`; requires two consecutive polls with the same count (GPFS settle). Default poll interval: 300 s (`WATCH_MERGE_POLL_SEC`). Default timeout: 72 h (`WATCH_MERGE_TIMEOUT_SEC` or `--timeout-hours` on the watch script).
- **runmeta update**: on finish, `runmeta_<anaName>_<jobid>.json` gets a `postProcess.watchMerge` block (`status`, paths, counts).
- **Skip**: if `*_merge.root` already exists, the watcher exits without re-merging (use `watch_job_and_merge.sh --force-merge` to override).
- **Bad-root exclusion (manual list)**: `watch_job_and_merge.sh --exclude-bad-roots <list.txt>` forwards to `merge_root_files.csh --exclude-list=<list.txt>`.
- **QA PDF**: not run automatically; after merge, run the usual `singularity_checkHistAna*` script on `*_merge.root` manually.

Manual re-run (e.g. after a timeout):

```bash
./script/watch_job_and_merge.sh --runmeta job/run/runmeta/runmeta_<anaName>_<jobid>.json
```

Manual re-run with exclusion list:

```bash
./script/watch_job_and_merge.sh \
  --runmeta job/run/runmeta/runmeta_<anaName>_<jobid>.json \
  --force-merge \
  --exclude-bad-roots /tmp/bad_subjob_roots_<anaName>_<jobid>.txt
```

## Manual bad-root scan and re-merge

Use this when `*_merge.root` looks suspiciously small or `hadd` failed due to bad subjob ROOT files.

```bash
# 1) Build exclusion list from a sample subjob file
./script/scan_bad_subjob_roots.sh \
  --sample rootfile/<anaName>/<anaName>_<jobid>_0.root \
  --output /tmp/bad_subjob_roots_<anaName>_<jobid>.txt

# 2) Re-merge with exclusions (auto-overwrites *_merge.root)
./script/merge_root_files.csh \
  --exclude-list=/tmp/bad_subjob_roots_<anaName>_<jobid>.txt \
  rootfile/<anaName>/<anaName>_<jobid>_0.root
```

Notes:
- `merge_root_files.csh` now runs `scan_bad_subjob_roots.sh` automatically by default.
- Use `--skip-bad-scan` to disable automatic scan.
- Use `--exclude-list=<file>` to force a specific exclusion list.

## Quick Check

- Generate from the intended mainconf: `./script/generate_joblist.sh config/mainconf/main_auau19_anaPhi_test.yaml`
- Confirm the joblist embeds the same path: `python script/analysis_info_helper.py --mainconf-from-joblist job/joblist/joblist_auau19_anaPhi_test.xml`
- Confirm setup metadata resolves from that same path: `python script/analysis_info_helper.py --library-tag --mainconf config/mainconf/main_auau19_anaPhi_test.yaml`
- For short farm smoke tests, set `analysis.maxEvents` (for example `100`) in the referenced `analysis_info`; the generated batch command will pass that as the 4th `root4star` argument
- Submit and check the preflight log prints the same mainconf path before `star-submit`
- After submit, verify `runmeta/runmeta_<anaName>_<jobid>.json` exists and points to the expected sidecars

## Cleaning up this directory

After submission, SUMS leaves many files named `anaName+jobid+*` (e.g. `auau3p85fxt_anaPhiCCBCC32EA67793F5A24B5F6BA44EE413_0.csh`). To avoid "Argument list too long" when removing them:

- **Delete one job's loose files:** `./script/cleanup_job_run.sh <anaName+jobid>` (e.g. `./script/cleanup_job_run.sh auau3p85fxt_anaPhiCCBCC32EA67793F5A24B5F6BA44EE413`).
- **Delete all loose SUMS files:** `./script/cleanup_job_run.sh --all`. This removes only regular files directly under
  `job/run/` whose names contain a 32-hex jobid. It preserves `runmeta/`, `configlog/`, `joblog/`, and all other
  subdirectories.
- **Archive** them under `joblog/<anaName>/`: `./archive_job_run.sh <anaName+jobid>` (creates `joblog/<anaName>/` if needed).

`runmeta/sums_artifacts_<anaName>_<jobid>.tar.gz` already preserves a stable submit-time copy of those SUMS-generated files, so cleanup/archive of the loose originals is now an operational convenience rather than the only reproducibility path.

## Notes

- Always `cd` into **job/run/** before running `./submit.sh`.
- To use a different joblist template, run e.g. `./submit.sh ../joblist/YourJoblist.xml`.
- Submit-time reproducibility artifacts now live in `job/run/runmeta/`; `script/archive_all_job_logs.sh` archives that directory alongside `configlog` and `joblistlog`.
- Batch command now clears `analysis/<baseAnaMacro>_C.*` before `root4star` so stale ACLiC outputs (`.so/.d/.pcm` etc.) do not mix across environments.
- Batch runtime now copies `analysis/`, `config/`, `lib/`, `include/`, `StMaker/`, `StRoot/`, and optional build artifacts into a scratch-local runtime bundle before `singularity exec`, so moved or symlinked repositories do not break macro lookup inside the container.
- Current `auau19_anaLambda` joblists run `root4star` via `singularity exec ... star-bnl/star-sw:latest` with `-B /star/nfs4/AFS`, `-B /home/starlib:/home/starlib`, and inherited `LD_LIBRARY_PATH` to satisfy `libgfortran.so.3`.
- The same container strategy is available for local Phi QA via `script/singularity_checkHistAnaPhi.sh` when host `root4star` cannot start due to missing runtime libraries.
- For local LL/KP CF fitting on AL9, use `script/singularity_run_fitCorrelation.sh <root_file> <mainconf> [hist_name]` (compiled `fit_correlation` binary inside Singularity; input must be a TH1D CF histogram, not raw maker merge ROOT).
- Spack `root-config` / batch `root4star` builds here omit `Netx`/`RFIO` plugins, so `root://` and `rfio://` URLs in the SUMS `.list` cannot be opened. The test joblist rewrites the list to POSIX paths (`sed` strips `root://host:port//` and `rfio://`) before `root4star`, then reads `/home/starlib/...` as local files inside the container.

## root4star exit abort (`corrupted size vs. prev_size`)

If analysis finishes and ROOT files copy successfully but stderr shows glibc heap errors at process exit, see **`TROUBLESHOOTING_root4star_exit.md`** in this directory. It lists **debug joblists** (minimal exit, `maxEvt=1` / `100`, ACLiC vs interpreter, `gSystem->Exit(0)`), **Pr3** `Make()` stage bisect (`joblist_auau19_anaLambda_debug_pr3_max1_no_hist_s[0-4].xml`), **Pr4** Phi control and optional **`MALLOC_CHECK_=3`** joblist, a **reproduction matrix** template, **hypotheses**, and **operational success criteria** when a full fix is not available. For a **local symbolic backtrace** (no gdb in `star-sw:latest`), use **`./local_root4star_anaLambda_backtrace.sh`**. The matrix also includes **`joblist_auau19_anaLambda_debug_max1_no_hist.xml`** (Lambda without `HistManager`); see troubleshooting row **B2**.

## Aggregate completion checks and guarded merging

The job template accepts SUMS one-column paths or two-column `path count` input
lists. It builds a path-only ROOT input list and checks scratch paths, input
readability, process status and produced ROOT output before copying results.
Scheduler exit code 0 alone does not establish successful analysis. Completion
checks also require expected outputs and successful analysis completion records.

`script/watch_job_notify.py --runmeta RUNMETA --event-log CONDOR_LOG --email ADDRESS --check`
is read-only. `--watch` waits for all subjobs to terminate, an empty queue and
a settle interval, then sends one aggregate completion email. It checks each ROOT
file's header and the analysis log's matching requested/completed counts and
successful status. It does not read every ROOT object or validate the complete
candidate payload. Failed or unconfirmed analysis is reported as requiring
attention, not successful completion. State and one-time delivery records live
under `job/run/notify/<anaName>_<jobid>/`; restart on the same submit host and
inspect uncertain delivery records before retrying mail.

Add `--email-content minimal` to send only completion or attention text, without
the analysis name, SUMS ID, paths, counts or diagnostic details. Its Message-ID
is opaque; full diagnostics remain in the local state files. The selected
`emailContent` is recorded in settings and status, including the sent record.
The default `--email-content detailed` preserves the aggregate diagnostic email.
Use `--state-dir PATH` to place the lock, settings, status and one-time delivery
records on an explicit filesystem. The resolved location is recorded as
`stateDirectory`; the default remains `job/run/notify/<anaName>_<jobid>/`.
`--check` reports that location without creating it.

For histogram-only outputs, validate and merge after all subjobs terminate:

```bash
python3 script/watch_job_merge_notify.py \
  --runmeta RUNMETA --event-log CONDOR_LOG \
  --root-runner script/singularity_root_batch.sh --mainconf MAINCONF \
  --email ADDRESS --watch
```

Use `--check` for a read-only status inspection and `--no-email` when email is not
wanted. Full successful input coverage is required unless `--allow-partial` is
explicit. The tool records included/excluded jobs, validates ROOT inputs and
merged output, refuses an existing destination and does not retry an uncertain
merge. State is local under `job/run/merge_notify/`; no study artifacts are
needed in Git. Review failures before starting another merge.

For comparisons that require the same source population, add
`--allow-partial --include-job-ids COMMON_IDS.list` to each merge command.
The selection file contains one unique, known nonnegative subjob number per
line, with no blank lines or comments. IDs are sorted numerically. The tool
records the file path, SHA256, IDs, and count in its settings, result, and input
provenance, and rejects a changed selection on restart or before merging.
Otherwise valid jobs outside that selection are recorded as
`outside_requested_selection`; their ROOT files remain outside the merged
input archive. Missing or incomplete jobs retain their own exclusion reasons.
The default still requires full coverage, and a corrupt ROOT header stops the
merge even in partial mode. `--check` only inspects scheduler completion; the
watch path performs the ROOT and analysis-log input checks.

After a validated publication, `python3 script/archive_merged_subjobs.py --runmeta RUNMETA`
relocates exactly the validated inputs into `subjobs/<SUMS id>` with a resumable
plan, without deleting ROOT data or overwriting a destination. It does not merge
or email. ROOT basenames follow the submitted XML's `fromScratch` template;
the submission name may contain a separate label used for state identification.
Oversized THnSparse objects require a separate, explicitly partitioned
merge: preserve same-axis Sparse sidecars, then sum low-dimensional projections.
The generic hadd workflow does not automatically convert such objects.
