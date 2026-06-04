# StMaker/common

Shared compiled helpers linked into multiple Maker libraries (same object file pattern as `CentralityHelper.o`).

## Phi KK reconstruction

**Source of truth:** [`StPhiKKReconstruction.h`](StPhiKKReconstruction.h) / [`StPhiKKReconstruction.cxx`](StPhiKKReconstruction.cxx)

Helix DCA, φ→K⁺K⁻ invariant mass, opening angle, pair rapidity, strict TOF pair cuts, and TOF fill logic used by `StPhiMaker` and `StFemtoMaker` live here. Cut thresholds come from YAML via `ConfigManager` (`PhiCutConfig`, `PIDCutConfig`).

Makers keep event loops, histogram filling, and analysis-specific orchestration (e.g. femto `BuildResonanceCandidates`, SE/ME pairing).
