# star-submit run directory

Run `star-submit` from **this directory**. SUMS will write generated files (`.csh`, `.list`, `.condor`, `.report`, etc.) here.

## Steps

1. **Build at the project root**
   ```bash
   cd /path/to/star-analysis
   sl7
   source ./script/setup.sh config/mainconf/main_auau19_anaLambda.yaml
   make
   ```
   For `csh` / `tcsh`, use `source ./script/setup.csh ...` instead.
   For debug/repro jobs (especially heap/exit issues), always build in SL7 before submit.

2. **Move to this directory and submit**
   ```bash
   cd job/run
   ./submit.sh ../joblist/joblist_auau19_anaLambda_test.xml
   ```
   Pass the joblist XML you want (for example the file produced by `./script/generate_joblist.sh`, named **`joblist_<anaName>.xml`** where **`anaName`** is `analysis.anaName` in your analysis_info). By default, `submit.sh` uses **`../joblist/joblist_auau19_anaLambda_temp.xml`** if you omit the argument (adjust the default in `submit.sh` if your primary analysis differs). Example for Phi: `./submit.sh ../joblist/joblist_auau19_anaPhi.xml`.

  `submit.sh` replaces `__PROJECT_ROOT__` in the template with the actual project path, so the same template works for any user.
  Before submit it runs a preflight check: it infers **`anaName`** from the joblist basename (`joblist_<anaName>.xml` → `<anaName>`), then extracts the embedded **`config/mainconf/...yaml`** path from the joblist and uses that as the single source of truth for `libraryTag` resolution, rebuilds, ELF class consistency of `lib/*.so` vs `root4star`, and runtime linker sanity in the singularity context. This fails fast on mismatch and avoids secondary errors like missing `fromScratch` output after an early crash.

   If preflight fails and you want the script to try recovery once, use:
   ```bash
   ./submit.sh --rebuild-if-needed ../joblist/joblist_auau19_anaLambda_test.xml
   ```
   This runs `source ./script/setup.sh <embedded-mainconf> && make`, then retries preflight with the same embedded mainconf.

   On successful submit, the joblist and config used for that run are saved: **joblistlog/joblist_<anaName>_<jobid>.xml** (submitted XML) and **configlog/config_<anaName>_<jobid>.txt** (the embedded mainconf and all referenced YAMLs in one text file).

## Quick Check

- Generate from the intended mainconf: `./script/generate_joblist.sh config/mainconf/main_auau19_anaPhi_test.yaml`
- Confirm the joblist embeds the same path: `python script/analysis_info_helper.py --mainconf-from-joblist job/joblist/joblist_auau19_anaPhi_test.xml`
- Confirm setup metadata resolves from that same path: `python script/analysis_info_helper.py --library-tag --mainconf config/mainconf/main_auau19_anaPhi_test.yaml`
- For short farm smoke tests, set `analysis.maxEvents` (for example `100`) in the referenced `analysis_info`; the generated batch command will pass that as the 4th `root4star` argument
- Submit and check the preflight log prints the same mainconf path before `star-submit`

## Cleaning up this directory

After submission, SUMS leaves many files named `anaName+jobid+*` (e.g. `auau3p85fxt_anaPhiCCBCC32EA67793F5A24B5F6BA44EE413_0.csh`). To avoid "Argument list too long" when removing them:

- **Delete** all matching files: `./cleanup_job_run.sh <anaName+jobid>` (e.g. `./cleanup_job_run.sh auau3p85fxt_anaPhiCCBCC32EA67793F5A24B5F6BA44EE413`).
- **Archive** them under `joblog/<anaName>/`: `./archive_job_run.sh <anaName+jobid>` (creates `joblog/<anaName>/` if needed).

## Notes

- Always `cd` into **job/run/** before running `./submit.sh`.
- To use a different joblist template, run e.g. `./submit.sh ../joblist/YourJoblist.xml`.
- Batch command now clears `analysis/<baseAnaMacro>_C.*` before `root4star` so stale ACLiC outputs (`.so/.d/.pcm` etc.) do not mix across environments.
- Current `auau19_anaLambda` joblists run `root4star` via `singularity exec ... star-bnl/star-sw:latest` with `-B /star/nfs4/AFS`, `-B /home/starlib:/home/starlib`, and inherited `LD_LIBRARY_PATH` to satisfy `libgfortran.so.3`.
- Spack `root-config` / batch `root4star` builds here omit `Netx`/`RFIO` plugins, so `root://` and `rfio://` URLs in the SUMS `.list` cannot be opened. The test joblist rewrites the list to POSIX paths (`sed` strips `root://host:port//` and `rfio://`) before `root4star`, then reads `/home/starlib/...` as local files inside the container.

## root4star exit abort (`corrupted size vs. prev_size`)

If analysis finishes and ROOT files copy successfully but stderr shows glibc heap errors at process exit, see **`TROUBLESHOOTING_root4star_exit.md`** in this directory. It lists **debug joblists** (minimal exit, `maxEvt=1` / `100`, ACLiC vs interpreter, `gSystem->Exit(0)`), **Pr3** `Make()` stage bisect (`joblist_auau19_anaLambda_debug_pr3_max1_no_hist_s[0-4].xml`), **Pr4** Phi control and optional **`MALLOC_CHECK_=3`** joblist, a **reproduction matrix** template, **hypotheses**, and **operational success criteria** when a full fix is not available. For a **local symbolic backtrace** (no gdb in `star-sw:latest`), use **`./local_root4star_anaLambda_backtrace.sh`**. The matrix also includes **`joblist_auau19_anaLambda_debug_max1_no_hist.xml`** (Lambda without `HistManager`); see troubleshooting row **B2**.
