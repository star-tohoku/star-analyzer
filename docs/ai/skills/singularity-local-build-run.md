# Skill: singularity-local-build-run

Build and debug locally through the batch-like Singularity runtime when host `make` or `root4star` is unreliable (for example AL9 login nodes, missing `libgfortran.so.3`, or mixed 32-bit ROOT plugins for `root://`).

## When to use

- Host `make` or `root4star` fails before the analysis macro runs.
- Farm preflight or batch jobs already use `singularity exec ... star-bnl/star-sw:latest`, and local behavior should match that runtime.
- After StMaker or Config C++ changes, rebuild `lib/*.so` before local or smoke validation.

## Build

From the project root:

```bash
./script/singularity_make.sh config/mainconf/main_<anaName>.yaml
```

- Default: `make clean && make` with `BUILD_BITS=64`.
- Use `--no-clean` when `src/third_party/yaml-cpp/build` already exists and you only need to rebuild project libraries.
- After `make clean`, CMake is required to rebuild `yaml-cpp`; the wrapper prepends cvmfs `cmake` when available.
- Prefer SL7 + `source ./script/setup.sh <mainconf> && make` when that host environment is healthy; use this skill as the fallback.

## Run and QA

- **Lambda:** `./script/singularity_run_anaLambda.sh` — same arguments as `run_anaLambda.sh`.
- **Phi:** `./script/singularity_run_anaPhi.sh MAINCONF [inputFile] [outputFile] [jobid] [nEvents]` — same arguments as `run_anaPhi.sh`.
- **Phi QA:** `./script/singularity_checkHistAnaPhi.sh <root_file> <mainconf_path>` when `checkHistAnaPhi.sh` fails for the same runtime reason.

## Debug workflow

1. Rebuild with `singularity_make.sh` after Maker or Config C++ changes.
2. Run a short local job with the matching `singularity_run_ana*.sh` wrapper and the same `mainconf` used for the build.
3. For Phi, generate QA with `singularity_checkHistAnaPhi.sh` when host QA fails.
4. For farm validation, continue with `generate_joblist.sh` and `job/run/submit.sh`; do not treat host-only `root4star` success as proof of batch readiness.

## Rules

- Keep physics and cut logic in Makers and YAML; these scripts only change the execution environment.
- Use the same `mainconf` for build, local run, and batch joblist generation unless the task explicitly compares variants.
- If script names or workflow change, follow `update-readme-scripts.md`.

## Reference

- `../../../INSTALL.md`
- `../../../docs/REFERENCE.md`
- `../../../docs/ai/AGENT_RULES.md`
