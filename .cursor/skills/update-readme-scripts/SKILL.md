---
name: update-readme-scripts
description: Updates README.md and job/run/README.md when adding or changing scripts so that someone who clones the repo can run the analysis from the docs. Use when adding scripts under script/ or job/run/, changing submit/run/QA workflow, or when the user asks to document scripts or update README for onboarding.
---

# README update for scripts and run workflow

Goal: README should let someone **clone → setup → build → run/submit → check results → clean up** without guessing.

## When to apply

- New or changed scripts in **script/** or **job/run/**.
- Submit/run/QA flow changes (e.g. configlog, new helpers).
- User asks to document scripts or improve "getting started" / onboarding in README.

## Files to touch

| File | Update when |
|------|-------------|
| **README.md** | Any script or workflow change that affects "How to run" or directory layout. |
| **job/run/README.md** | Changes under job/run (submit, configlog, cleanup, archive). |

## Checklist: README.md

1. **Directory layout table**
   - **script/**  
     List each user-facing script and one-line role (e.g. `checkHistAnaPhi.sh` = QA PDF from run_anaPhi output). Keep helpers as "e.g. get_file_list_*.sh".
   - **job/**  
     Mention `job/run/` scripts (`submit.sh`, `cleanup_job_run.sh`, `archive_job_run.sh`) and that config is written to `job/run/configlog/config_<anaName>_<jobid>.txt` on successful submit.

2. **How to run**
   - **Local run**: Keep or add run script usage (e.g. `run_anaLambda.sh`, `run_anaPhi.sh`) with optional args and one example.
   - **Result QA**: If there is a script that turns output ROOT into a report/PDF (e.g. `checkHistAnaPhi.sh`), add a short subsection: usage, one example, and where the output is written (path pattern).
   - **Batch (star-submit)**:  
     - Step 4 (submit): Add one sentence if submit now does something extra (e.g. "On successful submit, config is written to job/run/configlog/config_<anaName>_<jobid>.txt").  
     - If there are cleanup/archive scripts: add a short step or paragraph (e.g. step 5) with `cleanup_job_run.sh` and `archive_job_run.sh`: when to use (many anaName+jobid+* files), one example each (delete vs move to joblog/anaName).

3. **Tone**
   - One sentence per script: what it does and when to use it (e.g. "after submit" / "for QA" / "to avoid argument list too long").

## Checklist: job/run/README.md

1. **Submit**
   - After describing `submit.sh`, add one sentence if submit writes config (e.g. configlog path and that it's one text file per job).

2. **Cleaning up**
   - Add a short section (e.g. "Cleaning up this directory") if job/run gets many anaName+jobid+* files:
     - **cleanup_job_run.sh** `<anaName+jobid>`: delete all matching files (avoids "Argument list too long").
     - **archive_job_run.sh** `<anaName+jobid>`: move them to joblog/<anaName>/ (directory created if needed).
   - One example each is enough.

3. **Notes**
   - Keep: always `cd job/run` before submit; how to pass a different joblist.

## Principles

- **Onboarding first**: Prefer "clone → run analysis" over listing every option.
- **One sentence per script**: Name + what it does + when (e.g. "after submit", "for QA").
- **No duplication**: If job/run/README.md has the detail, main README can point there ("See job/run/README.md for more").
- **Consistent names**: Use the same anaName/jobid/configlog/cleanup/archive wording in both READMEs.
