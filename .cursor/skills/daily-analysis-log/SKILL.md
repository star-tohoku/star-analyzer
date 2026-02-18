---
name: daily-analysis-log
description: Generates end-of-day analysis logs for reproducibility. Creates analysisnote/YYYYMMDD/summaryYYYYMMDD.md and auto-fills development (git) and analyses run (configlog). Use when the user asks for a daily summary, end-of-day log, analysis note, or "今日の作業をまとめ".
---

# Daily analysis log

Goal: One markdown file per day under **analysisnote/YYYYMMDD/summaryYYYYMMDD.md** so that development, studies, and runs are recorded for **reproducibility**. The directory **analysisnote/** is in `.gitignore`; logs are local only.

## When to apply

- User asks for a **daily summary**, **end-of-day log**, **analysis note**, or similar.
- User says "今日の作業をまとめ" or "今日の開発を記録".

## Output

| Item | Value |
|------|--------|
| Directory | `analysisnote/YYYYMMDD/` (create if missing) |
| File | `analysisnote/YYYYMMDD/summaryYYYYMMDD.md` |
| Date | Today (default) or user-specified date as YYYYMMDD |

## Auto-generation workflow

1. **Resolve date**
   - If the user gave a date, use it as YYYYMMDD.
   - Otherwise use today in project timezone (e.g. `date +%Y%m%d` from repo root).

2. **Create directory**
   - `mkdir -p analysisnote/YYYYMMDD`

3. **Gather data (from repo root)**
   - **Development**:  
     `git log --since="YYYY-MM-DD 00:00" --until="YYYY-MM-DD 23:59" --oneline --name-status`  
     Optionally for each commit: `git show <hash> --stat` or a one-line summary.
   - **Configlog**:  
     If `job/run/configlog/` exists, list files: `config_<anaName>_<jobid>.txt`.  
     Optionally filter by date: `find job/run/configlog -name "*.txt" -newermt "YYYY-MM-DD" ! -newermt "YYYY-MM-DD 23:59:59"` (or list all and note "possibly from this day" if no mtime filter).

4. **Write summary file**
   - Use the template below.
   - **Development (code/config)**: Fill from git log (commits and changed files).
   - **Analyses run**: Fill from configlog listing: for each `config_<anaName>_<jobid>.txt` relevant to the day, add a line with mainconf path (`config/mainconf/main_<anaName>.yaml`), jobid, and configlog path (`job/run/configlog/config_<anaName>_<jobid>.txt`). If no runs found, leave a single line like "*(none recorded)*".
   - **Studies** and **Notes**: Leave placeholder lines so the user can add content.

5. **Reproducibility**
   - In "Analyses run", always include: **mainconf** path, **jobid**, and **configlog** path. Optionally the exact submit command or script if known (e.g. `submit.sh` with anaName).

## Template for summaryYYYYMMDD.md

```markdown
# Daily summary — YYYY-MM-DD

## Date
YYYY-MM-DD

## Development (code / config)
<!-- Auto-filled from git log for this day -->
- (commits and changed files)

## Studies
<!-- User: papers, notes, discussions -->

## Analyses run
<!-- Auto-filled from job/run/configlog when available -->
- mainconf: `config/mainconf/main_<anaName>.yaml` | jobid: `<jobid>` | configlog: `job/run/configlog/config_<anaName>_<jobid>.txt`

## Notes
<!-- User: follow-ups, TODOs -->
```

## Rules

- **Do not commit** `analysisnote/`; it is already in `.gitignore`.
- Prefer **overwriting** the summary file when re-running for the same day so the latest git/configlog state is reflected; the user can copy to backup before re-run if needed.
- If the user only wants to **open or edit** an existing summary, do that without re-gathering git/configlog unless they ask to refresh.
