# star-submit run directory

Run `star-submit` from **this directory**. SUMS will write generated files (`.csh`, `.list`, `.condor`, `.report`, etc.) here.

## Steps

1. **Build at the project root**
   ```bash
   cd /path/to/star-analysis
   source script/setup.sh config/mainconf/main_auau19_anaLambda.yaml && make
   ```

2. **Move to this directory and submit**
   ```bash
   cd job/run
   ./submit.sh
   ```
   By default, `submit.sh` uses the template `../joblist/joblist_run_anaLambda.xml`. You can pass another template, e.g. `./submit.sh ../joblist/joblist_run_anaPhi.xml`.

   `submit.sh` replaces `__PROJECT_ROOT__` in the template with the actual project path, so the same template works for any user.

   On successful submit, the joblist and config used for that run are saved: **joblistlog/joblist_<anaName>_<jobid>.xml** (submitted XML) and **configlog/config_<anaName>_<jobid>.txt** (mainconf and all referenced YAMLs in one text file).

## Cleaning up this directory

After submission, SUMS leaves many files named `anaName+jobid+*` (e.g. `auau3p85fxt_anaPhiCCBCC32EA67793F5A24B5F6BA44EE413_0.csh`). To avoid "Argument list too long" when removing them:

- **Delete** all matching files: `./cleanup_job_run.sh <anaName+jobid>` (e.g. `./cleanup_job_run.sh auau3p85fxt_anaPhiCCBCC32EA67793F5A24B5F6BA44EE413`).
- **Archive** them under `joblog/<anaName>/`: `./archive_job_run.sh <anaName+jobid>` (creates `joblog/<anaName>/` if needed).

## Notes

- Always `cd` into **job/run/** before running `./submit.sh`.
- To use a different joblist template, run e.g. `./submit.sh ../joblist/YourJoblist.xml`.
