---
name: singularity-local-build-run
description: Build and debug locally through the batch-like Singularity image (`star-bnl/star-sw:latest`). **`./script/singularity_make.sh <mainconf>` is the standard SL7-equivalent `make`** (same STAR `sl73_*` toolchain as batch) when not inside interactive `sl7`; also use when host `make` or `root4star` is unreliable (AL9, missing `libgfortran.so.3`, etc.). Use when the user asks for this skill or Singularity local build/run.
---

# Skill wrapper: singularity-local-build-run

Canonical procedure:

- [docs/ai/skills/singularity-local-build-run.md](../../../docs/ai/skills/singularity-local-build-run.md)

Related source-of-truth:

- [docs/ai/AGENT_RULES.md](../../../docs/ai/AGENT_RULES.md)
- [PHILOSOPHY.md](../../../PHILOSOPHY.md)
