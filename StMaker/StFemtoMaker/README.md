# StFemtoMaker naming rules

Canonical location for femto **species keys**, **particle keys**, and **channel names**.
When adding particles or channels, follow these rules and update `FemtoConfig` YAML accordingly.

## Species keys (`species` map in femto YAML)

- Format: `lower_snake_case`
- One key per candidate source registered in the event store
- Examples: `proton`, `anti_proton`, `kaon_plus`, `phi`
- Each species entry requires:
  - `builderType`: `track` or `resonance`
  - `particleKey`: builder dispatch key (see below)

## Particle keys (`particleKey` field)

- Selects builder logic inside generic `TrackPidBuilder` / `ResonanceBuilder`
- Phase 1 supported values: `proton` (track), `phi` (resonance from KK)
- Background QA: `phi_rotation` (resonance builder dispatch → `BuildRotatedPhiCandidates`, stored as species `phi_rot`)
- Use PDG-inspired names; charge sign goes in species key when needed (`anti_proton`)

## Channel names

- Format: `{partA}_{partB}` using species keys
- Example: `phi_proton` with `partA: phi`, `partB: proton`
- Background channels (same `partB`, different `partA` / mass window): `phi_proton_signal`, `phi_proton_leftSB`, `phi_proton_rightSB`, `phi_rot_proton`
- Histogram suffixes use channel name: `hKstarSE_phi_proton`, `hCF_phi_proton`, `hKstarSEVsCent_phi_proton_signal`, etc.
- Mass window per channel: `signalMin` / `signalMax` on resonance `partA` invMass at pairing time (sidebands need no extra maker logic)

## Correlation function (CF)

- `StFemtoMaker` fills **SE/ME k* histograms only** (`hKstarSE_*`, `hKstarME_*`). It does not compute CF at `Finish()`.
- `hCF_*` may exist in hist YAML as empty shells; do not rely on CF in subjob or merged ROOT from Maker.
- **CF is computed in `checkHistAnaFemtoPhiProton.C`** from merged SE/ME after all events are summed:
  - SE/ME are rebinned first when `cfRebinFactor` > 1 in maker YAML (e.g. 100 k* bins / 5 → 20 CF bins).
  - **Cent slice CF**: project `hKstarSEVsCent_*` / `hKstarMEVsCent_*` over `cfCent9Min`–`cfCent9Max` (default 2–8 ≈ 0–60%), then same CF formula. Drawn as TGraphErrors with Poisson errors.
  - Bin-wise: `C(k*) = SE / (ME * seNorm/meNorm)` with `seNorm`, `meNorm` integrals over `normQMin`–`normQMax` from maker YAML.
  - Default norm region for `anaFemtoPhiProton`: **0.5–1.0 GeV/c** (high k*, where C→1 when SE≈ME).
- Run QA on `*_merge.root`: `./script/singularity_checkHistAnaFemtoPhiProton.sh <merge.root> <mainconf>`.

## Centrality cuts (event level)

- Centrality YAML (`centrality:` in mainconf) may set optional **per-bin refMult_corr cap**:
  - `cent9MaxRefMultCorrBin: 8`, `cent9MaxRefMultCorr: 220` → reject cent8 events with refMult_corr > 220 (auau3p85fxt tail trim).
  - Applied in `CentralityHelper::AcceptCentBin` after `StRefMultCorr` bin assignment; requires analysis re-run to affect ROOT histograms.

## Mainconf and maker YAML

- StFemtoMaker analyses use mainconf key **`maker:`** (one file per StMaker; see `PHILOSOPHY.md`).
- Path: `config/maker/maker_<anaName>.yaml` (e.g. `maker_auau3p85fxt_anaFemtoPhiProton.yaml`).
- The same flat YAML holds both:
  - **φ resonance builder** keys (`PhiCutConfig`: `nSigmaKaon`, `minPairRapidity`, …)
  - **femto species/channels** keys (`FemtoConfig`: `speciesKeys`, `channel_*`, `proton*`, `rotation*`)
- `ConfigManager` loads `maker:` into `GetPhiCuts()` and `GetFemtoConfig()` from one file.

## YAML flat-key layout

`YamlParser` is flat; the femto section of maker YAML uses prefixed keys:

```yaml
speciesKeys: proton,phi
species_proton_builderType: track
species_proton_particleKey: proton
species_phi_builderType: resonance
species_phi_particleKey: phi

channelName: phi_proton
partA: phi
partB: proton
signalMin: 1.01
signalMax: 1.03
normQMin: 0.5
normQMax: 1.0
```

## Proton femto bachelor cuts (Zhangwei-like)

YAML keys on femto config (not `PIDCutConfig`): `protonChargeMode`, `protonMaxDca`, `protonMinPtPre`, `protonMinPtPair`, `protonMaxPtPair`, `protonMaxAbsEta`, `protonMaxAbsNSigma`, `protonMinNHitsFit`, `protonMinNHitsRatio`, `protonTofMomentumThreshold`, `protonMinMass2`, `protonMaxMass2`, `protonMinRapidityCm`, `protonMaxRapidityCm`. Applied in `PassFemtoProtonCuts` after generic `PassProtonCuts` + `IsProton`.

## Rotation background

When `rotationEnabled: true`, species `phi_rot` with `particleKey: phi_rotation` fills rotated φ candidates (`rotationN` azimuth rotations per real KK pair; DCA/PID/TOF evaluated on real tracks). Pair with proton via channel `phi_rot_proton`.

## Agent / developer checklist

1. Add species block with unique species key
2. Implement or extend builder dispatch for new `particleKey` in `StFemtoMaker`
3. Add channel entry referencing existing species keys
4. Add histogram entries to `config/hist/hist_anaFemto*.yaml` using channel name suffix
5. Document new keys in this file

See also `docs/femto_track_sharing_concerns.md` for overlap risks when extending species.
