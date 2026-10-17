# AI docs index

This directory contains AI guidance documents managed as the primary source for agent behavior.

## Source-of-truth documents

- `../../PHILOSOPHY.md`: Framework philosophy and non-negotiable architecture rules.
- `AGENT_RULES.md`: Project workflow, naming, layout, and editing conventions for agents (including **batch-matched builds**: interactive `sl7` + `make`, or **`./script/singularity_make.sh <mainconf>`** as the SL7-equivalent Singularity build).
- `skills/`: Task-oriented procedural playbooks used by agent tools.
- `candidate-bdt-guide.md`: Design principles of the candidate BDTs and how to add one for a new particle (for people and agents).
- `MIGRATION_PLAN.md`: Migration design and rollout policy for this structure.

## Tool mapping

- Antigravity / Claude Code reads `../../CLAUDE.md`, which points to these docs.
- Cursor rules and skills under `.cursor/` are wrappers that reference this directory.
