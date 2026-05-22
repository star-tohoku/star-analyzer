# StRefMultCorr (vendored)

Copied from STAR papers analysis tree:

- **Source**: `/gpfs/mnt/gpfs01/star/pwg/oura/papers/psn0836/9GeV/StRoot/StRefMultCorr/`
- **Date**: 2026-05-19
- **Files**: `StRefMultCorr.{h,cxx}`, `CentralityMaker.{h,cxx}`, `Param.{h,cxx}`, `BadRun.h`, `StRefMultCorrLinkDef.h`, `checkStRefMultCorr.C` (reference macro)

## Notes

- PicoDst I/O remains from `$STAR`; only centrality correction code is vendored.
- `Param.h` targets BES-II / low-energy runs; verify `init(runId)` for 19 GeV Au+Au run ranges. If failures appear, merge parameters from `papers/psn0836/19GeV/StRoot/StRefMultCorr/`.
- Load `lib/libStRefMultCorr.so` from this repo, not `$STAR/lib/libStRefMultCorr.so`, to avoid duplicate symbols.
- **Location:** `StRoot/StRefMultCorr/` (not under `StMaker/`; this code is not a `StMaker` subclass).
- `ClassDef` / `ClassImp` are stripped so the library loads without a ROOT dictionary (used only from compiled `CentralityHelper` / `StPhiMaker`).
- **FXT run coverage:** `fxtmult` index 0 run stop was extended to `22129999` so P24iy `production_3p85GeV_fixedTarget_2021` (e.g. run `22121033`) can initialize until dedicated 2021:3.85 parameters are vendored from STAR/papers (calibration still Run18 3.85 GeV tables).
