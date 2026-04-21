# Skill: daily-analysis-log

Creates an end-of-day reproducibility note under `analysisnote/YYYYMMDD/summaryYYYYMMDD.md`.

## Workflow

1. Resolve date (`YYYYMMDD`) from user input or current date.
2. Ensure directory exists: `analysisnote/YYYYMMDD/`.
3. Gather development summary from git history for the day.
4. Gather run evidence from `job/run/configlog/` when present.
5. Write or refresh summary file with date, development, studies, analyses run, and notes.

## Reproducibility requirements

- In "Analyses run", include:
  - mainconf path
  - jobid
  - configlog path
- Keep `analysisnote/` uncommitted (`.gitignore` managed).

## Notes

- Prefer overwrite-on-refresh for the same day unless the user asks to append.
