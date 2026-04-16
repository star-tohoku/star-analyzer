---
name: update-readme-scripts
description: Updates docs/REFERENCE.md, INSTALL.md, README.md (index), and job/run/README.md when adding or changing scripts so that someone who clones the repo can run the analysis from the docs. Use when adding scripts under script/ or job/run/, changing submit/run/QA workflow, or when the user asks to document scripts or update onboarding.
---

# Documentation update for scripts and run workflow

Goal: **[INSTALL.md](../../../INSTALL.md)** and **[docs/REFERENCE.md](../../../docs/REFERENCE.md)** let someone **clone → setup → build → run/submit → check results → clean up** without guessing. **[README.md](../../../README.md)** stays a short index; only update it when the doc map or top-level layout summary changes.

## When to apply

- New or changed scripts in **script/** or **job/run/**.
- Submit/run/QA flow changes (e.g. configlog, new helpers).
- User asks to document scripts or improve getting started / onboarding.

## Files to touch

| File | Update when |
|------|-------------|
| **docs/REFERENCE.md** | Directory layout (script list), How to run, batch steps, QA — the detailed reference. |
| **INSTALL.md** | First-time or prerequisite steps change (e.g. new mandatory setup command). |
| **README.md** | Documentation table or repository layout summary needs a new link or one-line note. |
| **job/run/README.md** | Changes under job/run (submit, configlog, cleanup, archive). |

## Checklist: docs/REFERENCE.md

1. **Directory layout (long table)**  
   - **script/** — List each user-facing script and one-line role (e.g. `checkHistAnaPhi.sh` = QA PDF from run_anaPhi output). Keep helpers as "e.g. get_file_list_*.sh".  
   - **job/** — Mention `job/run/` scripts (`submit.sh`, `cleanup_job_run.sh`, `archive_job_run.sh`) and that config is written to `job/run/configlog/config_<anaName>_<jobid>.txt` on successful submit.

2. **How to run**  
   - **Local run**: Run script usage (e.g. `run_anaLambda.sh`, `run_anaPhi.sh`) with optional args and one example.  
   - **Result QA**: Script that turns output ROOT into a report/PDF — usage, one example, output path pattern.  
   - **Batch (star-submit)**: Submit sentence if behavior changes; cleanup/archive (`cleanup_job_run.sh`, `archive_job_run.sh`) when many `anaName+jobid+*` files appear.

3. **Tone**  
   - One sentence per script: what it does and when to use it.

## Checklist: INSTALL.md

- If a new script is part of the **minimal** clone-to-first-run path, add it to the numbered steps or Quick reference.  
- Keep **`source script/setup.sh`** wording consistent with REFERENCE.

## Checklist: README.md

- Usually no change unless a **new top-level doc** or a major layout change needs a row in the documentation table or repository layout table.

## Checklist: job/run/README.md

1. **Submit** — One sentence if submit writes config (configlog path; one text file per job).  
2. **Cleaning up** — `cleanup_job_run.sh` / `archive_job_run.sh` with one example each.  
3. **Notes** — `cd job/run` before submit; how to pass a different joblist.

## Principles

- **Onboarding first**: INSTALL + REFERENCE favor "clone → run" over listing every option.  
- **No duplication**: job/run detail lives in job/run/README.md; REFERENCE points there for farm cleanup.  
- **Consistent names**: Same anaName / jobid / configlog / cleanup wording across INSTALL, REFERENCE, and job/run/README.
