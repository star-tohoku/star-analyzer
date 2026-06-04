# Femto species and channel naming

Use when adding femto particles, channels, or histograms in `StFemtoMaker`.

## Canonical rules location

**`StMaker/StFemtoMaker/README.md`** — source of truth for:

- species keys (`lower_snake_case`)
- `particleKey` / `builderType` conventions
- channel naming (`partA_partB`)
- flat femto YAML key layout

Headers `include/FemtoCandidate.h` and `include/cuts/FemtoConfig.h` point to this file; do not duplicate rules elsewhere.

## Required steps when extending femto

1. Read `StMaker/StFemtoMaker/README.md` and follow naming exactly.
2. Update `config/cuts/femto/femto_<anaName>.yaml` with new species/channel keys.
3. Extend builder dispatch in `StFemtoMaker` for new `particleKey` values (generic builders, not analysis-specific classes).
4. Add histogram names using channel suffix: `hKstarSE_<channel>`, `hKstarME_<channel>`, `hCF_<channel>`.
5. Update `config/hist/hist_anaFemto*.yaml` and checkHist macro pages if QA changes.
6. Document new keys in `StMaker/StFemtoMaker/README.md`.

## Guardrails

- Do not modify `StPhiMaker` / `StLambdaMaker` production behavior for femto work.
- Keep cuts in YAML; no hardcoded thresholds in maker code.
- k* uses PRF relative momentum (`|q*|/2` after pair boost), not lab 3-momentum subtraction.

## Related docs

- `docs/femto_track_sharing_concerns.md` — overlap / sharing risks
- `.cursor/skills/add-new-analysis/SKILL.md` — general analysis entry pattern
