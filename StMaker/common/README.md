# StMaker/common

Shared compiled helpers linked into multiple Maker libraries (same object file pattern as `CentralityHelper.o`).

## Phi KK reconstruction

**Source of truth:** [`StPhiKKReconstruction.h`](StPhiKKReconstruction.h) / [`StPhiKKReconstruction.cxx`](StPhiKKReconstruction.cxx)

Helix DCA, φ→K⁺K⁻ invariant mass, opening angle, pair rapidity, strict TOF pair cuts, and TOF fill logic used by `StPhiMaker` and `StFemtoMaker` live here. Cut thresholds come from YAML via `ConfigManager` (`PhiCutConfig`, `PIDCutConfig`).

Makers keep event loops, histogram filling, and analysis-specific orchestration (e.g. femto `BuildResonanceCandidates`, SE/ME pairing).

## Nuclear ID helper

**Source of truth:** [`NuclearIdDeDxVsMom.h`](NuclearIdDeDxVsMom.h), [`StNuclearIdHelper.h`](StNuclearIdHelper.h) / [`StNuclearIdHelper.cxx`](StNuclearIdHelper.cxx)

LogPoly dE/dx vs p calibration tables (header-only constants) and light-nucleus identification (`IsHe4`, TPC nσ, TOF m², best-species). Operational cuts (`maxNSigmaNuclear`, `m2_selection`, etc.) come from `NuclearIdCutConfig` YAML. Used by `StFemtoMaker` for ⁴He bachelor candidates without loading `StNuclearIdMaker` on the chain.
