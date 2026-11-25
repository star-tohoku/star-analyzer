# StFemtoMaker naming rules

Canonical location for femto **species keys**, **particle keys**, and **channel names**.
When adding particles or channels, follow these rules and update `FemtoConfig` YAML accordingly.

## Species keys (`species` map in femto YAML)

- Format: `lower_snake_case`
- One key per candidate source registered in the event store
- Examples: `proton`, `anti_proton`, `kaon_plus`, `phi`
- Each species entry requires:
  - `builderType`: `track`, `resonance`, or `derived` (`derived` is populated by a parent
    background builder rather than dispatched independently)
  - `particleKey`: builder dispatch key (see below)

## Particle keys (`particleKey` field)

- Selects builder logic inside generic `TrackPidBuilder` / `ResonanceBuilder`
- Phase 1 supported values: `proton` (track), `he4` (track, via `StNuclearIdHelper`), `deuteron` (track, via `StNuclearIdHelper`), `kaon_minus` (track, K⁻ bachelor), `phi` (resonance from KK)
- Background QA: `phi_rotation` (resonance builder dispatch → `BuildRotatedPhiCandidates`, stored as species `phi_rot`)
- Use PDG-inspired names; charge sign goes in species key when needed (`anti_proton`)

## Channel names

- Format: `{partA}_{partB}` using species keys
- Example: `phi_proton` with `partA: phi`, `partB: proton`; or `phi_he4` with `partB: he4`; or `phi_deuteron` with `partB: deuteron`
- Background channels (same `partB`, different `partA` / mass window): `phi_proton_signal`, `phi_proton_leftSB`, `phi_proton_rightSB`, `phi_rot_proton` (or `phi_he4_*` / `phi_rot_he4` for ⁴He; `phi_deuteron_*` / `phi_rot_deuteron` for d)
- Histogram suffixes use channel name: `hKstarSE_phi_proton`, `hCF_phi_proton`, `hKstarSEVsCent_phi_proton_signal`, `hPhiMKK_vs_KstarSE_phi_proton_signal` (TH3: M_KK × k* × cent9), etc.
- Mass window per channel: `signalMin` / `signalMax` on resonance `partA` invMass at pairing time (sidebands need no extra maker logic)

## Correlation function (CF)

- `StFemtoMaker` fills **SE/ME k* histograms only** (`hKstarSE_*`, `hKstarME_*`). It does not compute CF at `Finish()`.
- `hCF_*` may exist in hist YAML as empty shells; do not rely on CF in subjob or merged ROOT from Maker.
- **CF is computed in `checkHistAnaFemtoPhiProton.C`** from merged SE/ME after all events are summed:
  - SE/ME are rebinned first when `cfRebinFactor` > 1 in maker YAML (e.g. 100 k* bins / 5 → 20 CF bins).
  - **Legacy cent page (Page 20)**: project over `cfCent9Min`–`cfCent9Max` (default 2–8 ≈ 0–60%); columns signal / leftSB / rightSB / rot.
  - **Multi-slice CF (`cfCentSlices`)**: default 15 slices (cent9\_0…8 + pct\_0\_10…60). Per slice: signal, SB-L, SB-R, **SB-LR** (left+right histogram sum in checkHist, no new maker channel).
  - **Sideband-subtracted CF**: count-level `N_corr = N_sig - α N_SB` on SE and ME separately, then same CF formula (`sidebandSubtractAlpha`, default 1.0; `negativeBinPolicy: zero|skip`).
  - Bin-wise: `C(k*) = SE / (ME * seNorm/meNorm)` with `seNorm`, `meNorm` integrals over `normQMin`–`normQMax` from maker YAML.
  - Default norm region for `anaFemtoPhiProton`: **0.5–1.0 GeV/c** (high k*, where C→1 when SE≈ME).
- **PDF output** (one checkHist run, `share/figure/<anaName>/`):
  - QA: `{anaName}_checkHistAnaFemtoPhiProton_{jobid}.pdf` — Pages 1–19 + Page 20 + representative slices (`cfCentSlicesQaPdfInclude`, default pct\_0\_10/20/30).
  - CF: `{anaName}_checkHistAnaFemtoPhiProton_CF_{jobid}.pdf` — remaining slices (`cfPdfExcludeQaSlices: true` avoids duplicate).
  - **Per-slice layout (representative + CF PDF):** 4 columns (signal / SB-L / SB-R / SB-LR) × 4 rows (SE k*, ME k*, raw `CF_sig_raw`, sub `CF_sig_sub_SBL|SBR|SBLR` on row 4).
- Run QA on `*_merge.root`: `./script/singularity_checkHistAnaFemtoPhiProton.sh <merge.root> <mainconf>`.

### Topic 3: k*-binned purity and CF_genuine (`anaFemtoPhi` unified)

- Maker fills **TH3** `hPhiMKK_vs_KstarSE/ME_<channel>_signal` (x = M_KK, y = k*, z = cent9) for signal channels only.
- checkHist (`checkHistAnaFemtoPhi.C`) post-processes per `cfCentSlices` slice:
  - scaled SE-ME subtraction on M_KK (per k* bin), Gaussian + constant background fit -> `lambda_sig(k*)`, `lambda_bkg = 1 - lambda_sig`.
  - `C_bkg` from ME mass shape: scale `hPhiMKK_vs_KstarME` in k* norm band, project signal mass window -> combinatorial CF (`cfBkgMode: me_mass`).
  - `C_genuine = 1 + [(C_meas-1) - lambda_bkg(C_bkg-1)] / lambda_sig` with `C_meas` = raw signal-channel CF.
  - `C_genuine` uncertainties use analytic propagation in `C_meas`, `C_bkg`, and `lambda_sig`; the `lambda_sig` error currently comes from the Gaussian-fit partial derivatives.
- YAML keys: `purityFitUseConstantBkg`, `purityFitGaussSigmaMin/Max`, `purityMinKstar/Max`, `purityMinEntriesPerBin`, `purityClampMin/Max`, `cfBkgMode`.
- QA PDF pages: `lambda_sig(k*)` and `C_meas` vs `C_genuine` for `cfCentSlicesQaPdfInclude` slices × all bachelor bases.
- Requires `legacyCfPagesEnabled: true`. Batch re-run after hist/Maker deploy so TH3 keys exist in merge ROOT.

### kstarMassFitCF (`anaFemtoPhi` unified)

Primary φ–p / φ–d correlation function. Per \(k^*\) bin, fit the full \(M_{KK}\) spectrum after subtracting a scaled ROT or MIX template:

\[
S = F - \alpha B,\qquad
\alpha = \frac{\int_{M_\alpha} F\,dM}{\int_{M_\alpha} B\,dM},\qquad
C_{\mathrm{raw}}(k^*) = \frac{Y_{\mathrm{SE}}(k^*)}{Y_{\mathrm{ME}}(k^*)}.
\]

\(Y\) is the gaus-only integral of \(S\) in the channel `signalMin`–`signalMax` window. \(C_{\mathrm{norm}}\) scales \(C_{\mathrm{raw}}\) to mean 1 in `normQMin`–`normQMax`.

- **Default background template is ROT** (`kstarMassFitCfTemplate: rot`). MIX is a cross-check (`kstarMassFitCfCrossCheck: true`) because older MIX ROOT files still use a superseded mixing/cap scheme.
- Maker fills **wide** TH3 `hPhiMKK_vs_KstarSE/ME_<base>_wide` (full \(M_{KK}\), no signal-window cut at fill). ROT/MIX templates use `hPhiMKK_vs_Kstar{SE,ME}_phi_{rot|mix}_{proton,deuteron}_wide`. Histogram names are unchanged.
- YAML (`kstarMassFitCf*`): `Enabled`, `Template` (`rot|mix`), `CrossCheck`, `FitMassMin/Max`, `KstarBinWidth` (0.050 GeV/\(c\) in `auau3p85fxt_anaFemtoPhi`), `LowKstarMergeBins` (2 → first point covers \(0<k^*<0.1\) GeV/\(c\)), `AlphaMassMin/Max` (1.04–1.06; must not overlap `phi_proton_signal`), `WriteSidecar`.
- \(\alpha\) window must not overlap the signal window (Validate error). Overlap with the fit range is a warning.
- Invalid \(k^*\) bins are recorded (`kmf_fitstatus_*`): 0=ok, 1=fit fail, 2=negative yield, 3=low stat, 4=norm fail, 5=non-finite. Failed bins are omitted from \(C_{\mathrm{raw}}\) (no fake zero yield). Negative \(S\) bins are **not** clipped.
- QA PDF: `share/figure/<anaName>/<anaName>_checkHistAnaFemtoPhi_kstarMassFitCf[_jobid].pdf`. Sidecar: `..._CFkmf[_jobid].root`.
- Graph names: `CF_kmf_{rot|mix}_{base}_{slice}_{raw|norm}`, `kmf_Y_{SE|ME}_{rot|mix}_{base}_{slice}`, `kmf_fitstatus_{rot|mix}_{base}_{slice}`.
- Simple SE/ME count ratio in the main QA PDF is diagnostic only.
- `legacyCfPagesEnabled: false` (default) skips Topic 3, old direct mass-fit, and Method 5 pages.
- Old names (one-release YAML aliases only): Method3 / Method3AllKstar / `cfDirectPurityMode: method3` / `method3BkgSub*` / `purityDirectFitMass*` / `CF_method3_bkgsub_*` / PDF `..._Method3AllKstar` / sidecar `..._CFmethod3`.
- Closure toy: `tools/test_cf_kstarmassfit.C`.
- Hist keys generated by `script/generate_hist_anaFemtoPhi.py` (`WIDE_FEMTO_CHANNELS`).

### Topic / method 5: sideband CF-subtraction (`anaFemtoPhi` unified)

- **Distinct from kstarMassFitCF and Topic 3 `C_genuine`.** Off unless `legacyCfPagesEnabled`. Method 5 uses sideband correlation functions, not ME-mass `C_bkg`.
- Formula:
  - `CF_sig = norm(SE_sig / ME_sig)`, `CF_SB = norm(SE_SB / ME_SB)` (SBLR = left+right SE/ME sum by default; no width α).
  - `CF_CFsub = [CF_sig - (1-P) CF_SB] / P` with constant per-slice purity `P` from φ M_KK gaus+const fit (`fit_slice`) or YAML `cfSubPurityFixed`.
- Cache keys (non-destructive): `CF_CFsub_SBLR_<base>`, `CF_CFsub_SBL_*`, `CF_CFsub_SBR_*`, `CF_sig_method5_*`, `CF_SB_*`.
- Existing products kept: raw CF, count-level `CF_sig_sub_SBL|SBR|SBLR`, Topic 3 `CF_genuine`.
- YAML: `cfSubtractionMode` (`none|method5`), `cfSubPurityMode` (`fixed|fit_slice`), `cfSubPurityFixed`, `cfSubSidebandCombine` (`sumLR|avgCF_LR`), `cfSubWriteSidecarRoot`, `cfSubLowStatsRebinExtra` (extra rebin for t/³He/⁴He).
- Sidecar ROOT: `share/figure/<anaName>/<anaName>_checkHistAnaFemtoPhi_CFsub_<jobid>.root` (graphs + purity parameters + meta).
- Unit / toy tests: `tools/test_cfsub_method5.C` (`P=1`, `CF_SB=1`, hand numeric).
- Reference job for development: `7BB546C8496FB964E74D00B21D3B1505`.

### h–K± two-body CFs and Kubo-rule background (`anaFemtoPhi` unified)

- **Additive; all existing paths unchanged.** Everything below is behind off-by-default flags.
- **New track species** `phikaon_plus` / `phikaon_minus` (`particleKey: phi_kaon_plus|phi_kaon_minus`). They use **production** φ-daughter PID (`PassPhiDaughterTofPid`, same as real/ROT/MIX production K), NOT the loose TPC collection and NOT the anaFemtoKaon `kaon_minus` bachelor cuts. Registered via `speciesKeys` in the maker YAML. Pool copies store `tofMatch`, `mass2`, and `deltaOneOverBeta` so MIX can apply the same daughter PID using each track's origin-event TOF.
- **Two-body h–K± channels** `phikaon_plus_<h>` / `phikaon_minus_<h>` for `h ∈ {proton, deuteron, triton, he3, he4}` (partA = kaon species, partB = hadron). Filled by the existing `FillSameEventPairs`/`FillMixedEventPairs` (a same-event same-track guard was added to `TracksOverlap`). Hist keys `hKstarSE/ME[VsCent]_<channel>`, empty `hCF_<channel>` shell; `*Kstar` for p,d and coarser `*KstarCoarse` (150 bins) for t/³He/⁴He. Gated by `enableHKaonTwoBody`; macro pages gated by the same flag (`isHKaonTwoBodyEnabled()`).
- **Kubo-rule 3-body background** (`FillKuboTripletBackground`, gated by `enableKuboTriplet`, bases `phi_proton`, `phi_deuteron`):
  - Legacy 1D (signal-window gated): `hKstarTripSEKp/SEKm/Mix_<base>`, `hMKKtriplet_<base>`.
  - Full-mass TH3 when `kuboStoreFullMass: true`: `hKuboMKK_vs_KstarSEK{p,m}/Mix/KK_<base>` (x=M_KK, y=k*, z=cent9). KK topology = ME hadron + SE K⁺K⁻.
  - Reject counter: `hKuboNRejectShared_<base>`.
- **ROT full-mass TH3**: `hPhiMKK_vs_Kstar{SE,ME}_phi_rot_{proton,deuteron}_wide` (fill before mass window; channels `phi_rot_proton`, `phi_rot_deuteron`).
  - MIX full-mass TH3 (standard event-mixed KK: current K×buffer opposite K, psn0585-style): `hPhiMKK_vs_Kstar{SE,ME}_phi_mix_{proton,deuteron}_wide` (channels `phi_mix_proton`, `phi_mix_deuteron`; YAML `fullyMixedEnabled`). kstarMassFitCF uses these as the MIX background template (cross-check; default is ROT).
  - Old vs Kubo background / genuine CF (checkHist): from 1D triplets; legacy direct-mass-fit × Kubo closure pages when `legacyCfPagesEnabled` and Kubo are both on.
- **YAML flags**: `enableHKaonTwoBody`, `enableKuboTriplet`, `enableKuboGenuine`, `kuboStoreFullMass` (default false).

### Phi–bachelor pair momentum angle QA (`anaFemtoPhi` unified)

- **Not used in CF** — QA histograms only; `FillSameEventPairs` / mixing unchanged.
- Angle θ: acos(p\_φ · p\_bachelor / |p\_φ||p\_bachelor|) in lab frame (distinct from KK opening angle).
- **Four φ QA stores / selections** for `phi_<bach>_signal` channels:
  - `hPhiPairMomAngle_phi_<bach>_signal` / `hPhiPairMomAngle_vs_MKK_*`: from `m_phiQaLoose` (passStage kinematics plus signal mass window, **no** `PassPairTofCut`); not stored in mixing pool.
  - `hPhiPairMomAngle_*_tofStrict` / `hPhiPairMomAngle_vs_MKK_*_tofStrict`: from `m_eventCandidates["phi"]` with `reso.betaGamma > 0` (both K daughters TOF-matched) plus signal mass window.
  - `hPhiPairMomAngle_vs_MKK_*_preMass`: all M_KK after DCA, track pT/eta, and nSigmaK cuts; no signal mass window, no KK opening-angle cut, no pair-rapidity cut.
  - `hPhiPairMomAngle_vs_MKK_*_preMass_tofStrict`: same preMass selection plus both K daughters pass the strict pair TOF requirement.
- Bachelor: existing `PassFemto*Cuts` candidates in `m_eventCandidates`.
- checkHist: `drawPhiBachelorPairAngleQaPage` — one PDF page per bachelor (p, d, t, ³He, ⁴He).
- Hist keys generated by `script/generate_hist_anaFemtoPhi.py` (30 histograms for 5 signal channels).

### Phi-near-track PID QA (`anaFemtoPhi` unified)

- **Not used in CF** — inclusive PID QA only (`FillPhiNearTrackPidQa`).
- Sample: all `PassTrackCuts` tracks (not bachelor-PID selected), excluding φ daughter indices.
- Selection: reconstructed φ in `m_eventCandidates["phi"]` with M_KK in `phiNearTrackSignalChannel` signal window; PRF `k*` vs track (mass hyp `phiNearTrackMassHyp`, default pion) below YAML thresholds.
- Histograms:
  - `hDedxVsP_PhiSignalNear_k1` / `hMass2ChargeVsP_PhiSignalNear_k1` — `k* < phiNearTrackMaxKstarLoose` (default 1.0 GeV/c)
  - `hDedxVsP_PhiSignalNear_k03` / `hMass2ChargeVsP_PhiSignalNear_k03` — `k* < phiNearTrackMaxKstarTight` (default 0.3 GeV/c)
- `m²×charge` y-axis covers ⁴He (`Mass2ChargeY`, ±30). TOF-matched tracks only for m² plots.
- Maker YAML: `phiNearTrackQaEnabled`, `phiNearTrackSignalChannel`, `phiNearTrackMaxKstarLoose/Tight`, `phiNearTrackMassHyp`.
- checkHist: one page (2×2) immediately before legacy CF pages (Page 15+).

### Phi-daughter production PID (`anaFemtoPhi` unified)

- **Loose K** (`kaonsPlus` / `kaonsMinus`): TPC quality + nσ (`PassKaonCuts`) plus collection filter `IsKaon` → `PassTofKaonPid` (`tofFallbackMode: tpcOnly` allows unmatched tracks). Used only for loose QA (`m_phiQaLoose`, `m_phiQaPreMassLoose`). `hNKaonPlus` / `hNKaonMinus` count this sample.
- **Production K**: `PassPhiDaughterTofPid` (full momentum \(p\), YAML `pMomKaonPID`). Same predicate for real φ, ROT, MIX, and `phikaon_*`. **Charge-asymmetric by default** (`phiDaughterKaonMinusRequireTof: true`):
  - K⁺: \(p \le pMomKaonPID\), no TOF allowed (TPC already applied); high \(p\) requires TOF
  - K⁻: **TOF match always required** (including \(p \le pMomKaonPID\)). Low-\(p\) TPC-only is K⁺ only
  - \(p > pMomKaonPID\), no TOF: rejected for both charges (`tpcOnly` must not pass high-\(p\) production K)
  - TOF-matched: `tofUseMass2Cut` / `tofUseDeltaInvBetaCut` as configured
- Do **not** make K⁺/K⁻ TOF policy charge-symmetric at FXT: low-\(p\) unmatched K⁻ is dominated by \(\pi^-\)/\(e^-\) and spoils \(\phi\) S/N and ROT \(M_{KK}\) templates.
- **Where production PID is required:** real φ candidates after `PassPairTofCut`, ROT, MIX (current and buffer daughters using stored TOF), and `phikaon_plus` / `phikaon_minus`.
- MIX does not apply a common DCA (cross-event). Pair topology and mixing sampler are unchanged.
- QA: `hPhiDauPid_*` (loose / prod / reject, K⁺/K⁻, category vs \(p\)/TOF) and `hPhiDauPidUsed_*_{real,rot,mix}_*` for daughters that entered production pairs. Generated by `script/generate_hist_anaFemtoPhi.py`.

## Event mixing (ME statistics)

- Config: `config/cuts/mixing/mixing_<ana>.yaml` via mainconf `mixing:` key (`MixingConfig`).
- Both modes use one eligible candidate-pair population: current A × buffered B plus,
  when `mixBothDirections` is true, buffered A × current B. The two directions are
  independent, so a missing current B does not suppress forward mixing and a missing
  current A does not suppress reverse mixing. Empty buffer-side candidate collections
  never enter the eligible population.
- `mixingMode`:
  - `randomSample` (default): uniformly sample without replacement from the eligible
    population, up to `maxMixedPairsPerEvent`.
  - `bufferAll`: exhaust the same eligible population (Zhangwei-like) and ignore
    `maxMixedPairsPerEvent`.
- P0-2 intentionally corrects the old `randomSample` behavior. The old implementation
  required both current A and B, derived its attempt count from the same-event
  multiplicities, sampled buffer events uniformly, and discarded draws with empty
  buffered B; it is no longer available as a physics mode.
- `maxMixedPairsPerEvent`: maximum attempted pairs per current event and channel for
  `randomSample` (default 500; 0 = unlimited).
- Weighting: every candidate-pair instance has equal inclusion probability in
  `randomSample`; `bufferAll` includes each once. A buffer event therefore contributes
  in proportion to its candidate-pair multiplicity, rather than receiving equal event
  weight. No additional multiplicity or event weight is applied. With a finite cap, a
  current event contributes `min(N_eligible, maxMixedPairsPerEvent)` attempts. Later CF
  normalization can absorb total ME scale but not a multiplicity-correlated shape change.
- QA: `hMixSamplerQA` is a channel-indexed counter matrix. Its status bins verify
  `eligible = attempted + skipped_by_cap`, `attempted = filled + skipped_by_pair_cut`,
  and `skipped = eligible - filled`. `selected_empty_buffer_direction` must remain zero;
  forward/reverse eligible and filled counts diagnose the two directions separately.
  These counters are raw counts; centrality/event weights are not applied.
- **Mixing bin:** `GetMixingBin(vz, cent9, psi2)` uses Vz bin × **cent9 as bin index** × EP bin. ME pairs only **different events** in the same bin (never same-event pairs; those are SE).
- **Shared-track provenance:** each resonance daughter is stored as `(eventIndex, trackIndex)`. Real and rotated φ
  daughters share the current event index; fully mixed φ daughters retain their distinct current/pool event indices.
  Overlap rejection requires both fields to match, so a reused local track index from another event is not rejected.
- **Centrality slices in checkHist:** `pct_0_10` etc. project cent9 ranges **after** fill; Maker mixing stays at cent9 resolution.
- Changing mixing requires new pair processing; a projection of existing `hKstar*VsCent_*` cannot rebuild ME. The optional Lambda–He4 candidate store below supports replay without rereading PicoDst. Analyses without a candidate cache require reconstruction again.
- **Rollout:** compare `randomSample` and `bufferAll` on a short joblist before full catalog; check high-k* `C→1`, ME shape, job time, merge ROOT size.
- Production mainconfs keep the `randomSample` mode name and receive its corrected semantics.

## Lambda–He4 reduced tree and producer QA

The dedicated `anaFemtoLambda_4He_tree` entry uses the same `StFemtoMaker`
candidate/pair implementation with a required `reducedTree:` feature YAML:

| Mode | Processing | Output |
|---|---|---|
| `write` | PicoDst event selection, centrality/EP, Lambda reconstruction and nuclear PID; no pairing or pool registration | One ROOT with the reduced tree and selected non-pair QA |
| `replay` | Stored candidates and raw relevant-track data; no PicoDst reader or Lambda provider | SE/ME and pair diagnostics with the requested downstream conditions |

The new format retains canonical Float pairing values and richer event,
Lambda-fit, daughter and raw-track information. Physical nuclear momentum is
separate from measured primary/global rigidity momentum and helix geometry.
The event TTree is `femtoLambdaReducedEvents`; `FemtoLambdaReducedMetadata`
records the persisted producer segments and their selection/layout/counts.
The reader uses those segment identities rather than requiring the current
transport file to retain its original ROOT UUID. It rejects incomplete coverage,
duplicate segments/events and noncontiguous segment order. Use completed
producer files in their recorded order; arbitrary partial splitting is not an
accepted input format.
Ordered hadd/repack of complete producer segments is supported. Keep each
original producer source-input manifest: file-level ROOT metadata alone does
not preserve every segment-to-source-list mapping after hadd.
Only target candidates and daughters are persisted in the standard compact
profile. The optional full-control profile also retains nonselected pre-PID
controls. Selected-pair PID/rigidity QA uses physical candidates and remains
available in both profiles;
all-track diagnostics can remain producer histograms without saving unrelated
particle rows. Events with no target pair, including events whose other legacy
candidates admitted them to the original pool, retain their processing and
pool-admission information.

Track `hasTof` means that the referenced BTOF traits object exists. `tofValid`
means only that beta is finite and positive; it does not certify all TOF fields
or apply a match-quality cut. Beta, time, match flag, local Y/Z and hit position
use native PicoDst accessors. The reserved `tofPathLength`, `tofDeltaY` and
`tofDeltaZ` fields are unavailable from the current native interface and remain
-999. Treat those fields as unavailable, including when `hasTof` is true.

The producer histogram references explicitly select non-pair QA. Pair histograms
and sparse outputs are not filled in that stage. Producer QA keeps its own
original population; it is not a replacement for a condition-dependent QA
histogram reconstructed during replay and is not automatically duplicated into
each replay output. Keep producer QA identity/coverage with any later merge.

The replay histogram profile includes optional k* × EPD-bin × cent9 outputs
for physical Lambda–He4 and Lambda–He3, SE/ME and all four mass classes.
The EPD coordinate labels the current processing event using the bin count of
the pair histogram's EPD axis. All such axes must have the same contiguous
integer bins. With the active mixing partition, both mixed events share that
class. A legacy EP1 profile may retain eight diagnostic classes, in which case
the label need not be the Lambda-side class. Summing the complete EPD axis
recovers the matching inclusive k* × cent9 distribution.
The five-bin pair-rapidity replay control keeps the same inclusive pair
physics and EP8 mixing as the ten-bin profile.

Feature settings are explicit: `mode`, `speciesKeys`,
`resetPoolAtSegmentBoundary`, `requireOriginalMixing`, `compressionLevel`,
`autoFlush`, `storePrePidControl`. The standard value of the last key is false:
the producer's 108 single-particle/event QA histograms and replay's selected
PID/rigidity diagnostics (540 histograms) remain, while only the 540
nonselected pre-PID pair controls are excluded. The full-control profiles set it true and retain that larger
bank. Capability metadata and selection fingerprints prevent requesting full
control QA from a compact input. The macro supplies producer identity and source input from its
mandatory jobid/input arguments. Original-mixing replay is the closure mode;
the separate EP1 control configuration relaxes original-mixing equality while
retaining the stored reconstruction selection. Changing downstream conditions
does not recover candidates outside that stored selection.
The embedded source revision is explicitly `unrecorded`; retain the source,
configuration and dependency checksums in the submit-time manifest alongside
the original ordered input list. A Git HEAD alone does not describe a dirty
working tree.

The entrypoints and complete commands are documented in
[`docs/REFERENCE.md`](../../docs/REFERENCE.md). The old `candidateStore:` format
below remains unchanged, and a mainconf cannot request both formats.

## Lambda–He4 candidate event store and replay

The `anaFemtoLambda_4He` entry supports an optional `candidateStore:` mainconf
reference. It uses the same `StFemtoMaker` pair functions for reconstruction and
replay; it does not introduce another pairing Maker. Without that key, the
PicoDst/KF path remains active.

| Feature mode | Input and processing | Output |
|---|---|---|
| absent | PicoDst, KF Lambda reconstruction and nuclear PID | Existing histograms and ledgers |
| `write` | The same reconstruction, PID and pairing | Existing output plus accepted events and exact candidate payloads |
| `replay` | Original candidate-store ROOT files; no PicoDst reader or Lambda provider | SE/ME and pair diagnostics from the stored candidates |

Storage lives in [`FemtoCandidateEventStore.h`](../../include/FemtoCandidateEventStore.h),
with the explicit settings and provenance policy in
[`FemtoCandidateEventStoreConfig.h`](../../include/FemtoCandidateEventStoreConfig.h).
The cache is part of the writer's analysis ROOT file. Its `femtoCandidateEvents`
TTree contains every accepted event in processing order, including events with
no selected candidates. Physical and pre-PID nuclear-QA collections retain their
presence flags, species, candidate order and QA selection flags. The cache species
list is explicit: `lambda,deuteron,triton,he3,he4`. Legacy nuclear collection
contains d/t even when the configured pair channels use only He3/He4. Those d/t
collections affect whether an event enters the original mixing pool, so dropping
them would change ME; caching them does not add d/t pair channels or alter PID.

Every `FemtoCandidate` field is stored: the original `Float_t` `px`, `py`, `pz`,
**energy**, derived kinematics, nuclear/track PID and resonance quantities, plus
integer charge, source, quality flags and track/daughter event references.
Replay restores these values directly and does not reconstruct energy with
`SetP4` or run PID again. Event records contain the original input index,
run/event IDs, cent9/cent16, raw/corrected multiplicity, track and TOF counts,
primary vertex, magnetic field, centrality weight and event-plane information.
Local event references are remapped together during replay so tracks from
different source events cannot acquire a false shared-track identity.

The accompanying ROOT metadata is:

- `FemtoCandidateEventStoreLayout`: schema and exact species/field layout;
- `FemtoCandidateEventStoreSelection`: selection settings, effective dE/dx calibration content and auxiliary event-plane definition;
- `FemtoCandidateEventStoreMixing`: the producer's mixing settings;
- `FemtoCandidateEventStoreFileUUID`: the original ROOT stream identity.

Use the original completed per-job files. Schema, species, branch types and
lengths, candidate references, selection provenance and stream identity are
checked before pairing. Duplicate files/UUIDs and changed entry counts are
errors. **Do not hadd or repack candidate caches:** changing their ROOT UUID or
losing producer boundaries prevents faithful replay. A byte-for-byte copy
preserves the UUID, but two copies of the same stream must not appear in one list.

### Configuration and execution

The standard He4 mainconfs are:

- `main_auau3p85_anaFemtoLambda_4He_KFParticle_candidateStoreWrite.yaml`;
- `main_auau3p85_anaFemtoLambda_4He_KFParticle_candidateStoreWriteNoJointSparse.yaml`;
- `main_auau3p85_anaFemtoLambda_4He_KFParticle_histOnlyNoJointSparse.yaml`;
- `main_auau3p85_anaFemtoLambda_4He_KFParticle_candidateStoreReplay.yaml`;
- `main_auau3p85_anaFemtoLambda_4He_KFParticle_candidateStoreRemixVz.yaml`.

They keep the production reconstruction/PID references. The write/replay pair
uses the original mixing; the Vz diagnostic changes only `nVzBins: 10` to `20`
over the same vertex range. The `candidateStoreWriteNoJointSparse` entry keeps
candidate writing and dense nuclear pair QA while disabling only the 8D PID
histogram through its dedicated QA YAML. Reconstruction, PID and mixing cuts
remain the same.

The [histogram-only mainconf](../../config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle_histOnlyNoJointSparse.yaml)
omits `candidateStore:` and uses the same dense-only nuclear pair QA. It writes
neither a candidate-event cache nor the joint 8D PID histogram; dense histograms,
PID, TOF selections and mixing match the no-Sparse candidate writer. Later mixing
changes require PicoDst reconstruction again because no replay cache is saved.

The candidate-store feature YAML requires explicit `mode`,
`speciesKeys`, `requireOriginalMixing`, `resetPoolAtFileBoundary`, `eventPlaneSource` and
`eventPlaneCalculation` settings.

For closure, use `requireOriginalMixing: true` and
`resetPoolAtFileBoundary: true`. The latter is the default boundary policy:
each producer file starts with a fresh pool, matching independent source jobs.
Setting it to `false` deliberately continues the pool across files and changes
the mixed-event reference. Preserve the ordered input list for reproducibility.
Setting `requireOriginalMixing: false` permits an explicit mixing-only variation;
it does not permit different cached selection/PID settings or candidate reselection.

The existing runner and wrapper select the stage from the mainconf:

```bash
bash script/singularity_run_anaFemtoLambda_4He.sh \
  config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle_candidateStoreWrite.yaml \
  <explicit-2021-P24iy-picoDst-or-list> <new-store.root> store_job -1
bash script/singularity_run_anaFemtoLambda_4He.sh \
  config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle_candidateStoreReplay.yaml \
  <original-store.root-or-list> <new-replay.root> replay_job -1
```

Replay lists have one ROOT path per line; relative paths are resolved against the
list's directory. Replay `nEvents` counts **stored accepted-event rows**, not
original PicoDst input events. `-1` processes every row; a positive request beyond
the available rows fails rather than truncating silently. Input, output, job ID
and mainconf remain explicit arguments, and existing outputs are refused.
Rebuild with the batch-matched toolchain after changing the Maker or store layout.

### Event-plane payload and validation scope

The original pairing `psi2` and its validity are separate from the auxiliary
TPC second-harmonic Q vector. The writer can calculate pT-weighted, **unflattened**
TPC Q2 from configured charged-track cuts, optionally excluding candidate tracks.
It stores Q2 components, angle, multiplicity and validity. This diagnostic angle
is not a calibrated or resolution-corrected physics event plane.
`eventPlaneSource: original` preserves the original pairing behavior, including
its `psi2=-1` state. `eventPlaneSource: stored` explicitly selects the auxiliary
angle. EP-binned replay requires a valid angle for every accepted event and stops
on an invalid value; it does not silently assign invalid angles to a valid bin.

Replay does not regenerate reconstruction, single-track or helix/merging QA
that requires PicoDst. Pair closure covers selected and cached pre-PID k*, mass,
centrality, mT, species/mass mappings and associated PID diagnostics.
[`compareFemtoCandidateReplay.C`](../../tests/compareFemtoCandidateReplay.C)
compares classes, axes, every bin including flows, Sumw2/errors, Entries,
statistics and sparse coordinates exactly. Its `sameEventOnly` mode verifies
unchanged SE for a mixing variation and reports whether the He4 ME probe changes.
### Verified replay behavior

With the full-QA configuration (`jointSparse` mapping enabled), the ROOT5
limited-statistics validation processes 4,027 PicoDst events, stores
3,858 accepted-event rows and replays all rows without a PicoDst reader or Lambda
provider. Baseline, writer and original-mixing replay agree exactly in 1,675 dense
pair histograms and one pre-PID QA Sparse: 39,963,446 dense cells including flows
and 89,815 occupied sparse cells, with identical errors, Entries and statistics.
The Vz-bin diagnostic preserves all 837 compared SE histograms exactly and changes
the He4 signal ME histogram from 251 to 235 entries, with 13 different cells.
This verifies persistence and controlled remix behavior; it does not establish a
physical CF baseline or normalization plateau.

The auxiliary TPC angle is invalid for one of the 3,858 stored events. EP-binned
replay therefore stops under the current strict validity rule. Comparing EP
classes requires an explicit policy for invalid-plane events; no silent dropping
or assignment to a valid bin is part of this workflow.

The [validation runner](../../tools/kfparticle_summary/validate_femto_candidate_store.py)
performs fresh baseline/write/replay checks or explicitly reuses successful
reconstruction records after checking the frozen input and recorded config/library
fingerprints. Reuse records the comparison-macro change and original ROOT
fingerprints; its execution log contains no new PicoDst reconstruction commands.

## QA PDF pre/post layout (`checkHistAnaFemtoPhiProton.C`)

- Pre-cut pages are immediately followed by post-cut pages (event, track, proton femto, phi candidate).
- **KK pair (A):** `hOpeningAngle_Raw`, `hPairRapidity_vs_Pt`, … filled after strict TOF, before opening/rapidity; **AfterCuts** on pass.
- **Femto φ (B):** `hPhi_*_PreCut` at candidate push; `hPhi_*` + `hPhi_PtVsY_PostCut` in `FillCandidateQA`.
- **Proton:** `hP_*_PreFemtoCut` in `FillProtonPreFemtoQa`; post-femto in `FillProtonFemtoQa` / `FillCandidateQA`.
- Re-run analysis after Maker/hist changes so new ROOT keys exist in merge output.

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

`YamlParser` is flat; the femto section of maker YAML uses prefixed keys.

**`nChannels`** must be at least one greater than the highest `channel_N_*` index present (e.g. `channel_0` … `channel_4` requires `nChannels: 5`). If `nChannels` is too small, sideband and rotation channels are silently skipped at load time.

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

YAML keys on femto config (not `PIDCutConfig`): `protonChargeMode`, `protonMaxDca`, `protonMinPtPre`, `protonMinPtPair`, `protonMaxPtPair`, `protonMaxAbsEta`, `protonMaxAbsNSigma`, `protonMinNHitsFit`, `protonMinNHitsRatio`, `protonTofMomentumThreshold`, `protonMinMass2`, `protonMaxMass2`, `protonMinRapidityCm`, `protonMaxRapidityCm`.

TOF/m² gate is **only** in `PassTofProtonPid` (also used by `IsProton` and `PassFemtoProtonCuts` — same predicate, no second rule):

- \(|p| <\) `protonTofMomentumThreshold`: TPC-only (no TOF match required; **no m² veto** even if TOF-matched).
- \(|p| \ge\) threshold: TOF match required and `protonMinMass2` ≤ m² ≤ `protonMaxMass2`.

Does **not** use global `PIDCutConfig.requireTOF` / `tofFallbackMode` (those remain for φ-daughter / collection paths).

## K⁻ bachelor (`anaFemtoKaon`)

- **Species key:** `kaon_minus` with `builderType: track`, `particleKey: kaon_minus`.
- **Femto cuts** (`kaonMinus*` keys in maker YAML): mirror proton femto layout — DCA, nσ, nHits, TOF momentum threshold + m² window, pair pT, CM rapidity. TPC pre-cuts use `kaonMinusMaxDca` / `kaonMinusMaxAbsNSigma` (not `PhiCutConfig`).
- **TOF:** `IsKaon` (`PIDCutConfig`) at ID time; femto TOF rule via `kaonMinusTofMomentumThreshold` + `kaonMinusMin/MaxMass2`.
- **Channels (5, no background):** `kaon_minus_proton`, `kaon_minus_deuteron`, `kaon_minus_triton`, `kaon_minus_he3`, `kaon_minus_he4`.
- **QA hist prefix:** `hKm_`; k* keys `hKstarSE/ME_<channel>`.
- **CF:** computed in `checkHistAnaFemtoKaon.C` from merged SE/ME (Maker fills histograms only). Common `normQMin`/`normQMax` and `cfCentSlices` across channels (default: cent9 8–8, 7–8, 6–8, 5–8).
- **checkHist:** `checkHistAnaFemtoKaon.C` — K⁻ QA + bachelor-loop + inclusive CF + cent slices (no SB/rotation). Per cent slice: one page with 5-channel SE + norm-scaled ME k* overlay, then one CF page (3×2 grid).

## Unified multi-bachelor (`auau3p85fxt_anaFemtoKaon`)

- **anaName:** `auau3p85fxt_anaFemtoKaon` — one picoDst pass for K⁻–(p, d, t, ³He, ⁴He).
- **Species keys:** `kaon_minus`, `proton`, `deuteron`, `triton`, `he3`, `he4`.
- **Channels (5):** `kaon_minus_<bachelor>` only (purity=1.0, no sidebands).
- **Nuclear ID:** same as unified `anaFemtoPhi` (`nuclearid_auau3p85fxt_anaFemtoPhi.yaml`).
- **Maker YAML:** minimal `PhiCutConfig` keys only (`rapidityFrame`, event-plane, `maxNTr`) — no φ reconstruction keys.

## Rotation background

When `rotationEnabled: true`, species `phi_rot` with `particleKey: phi_rotation` fills rotated φ candidates (`rotationN` azimuth rotations per real KK pair; DCA/PID/TOF evaluated on real tracks). Pair with proton via channel `phi_rot_proton`.

- `rotationChargeMode: legacyKp|legacyKm|bothSeparated`; production remains `legacyKp` until the
  charge-symmetric validation is accepted. `bothSeparated` additionally fills derived species
  `phi_rot_kp` and `phi_rot_km`; channels are named `phi_rot_kp_<bachelor>` and
  `phi_rot_km_<bachelor>`. The legacy `phi_rot` store stays identical to the K⁺ leg. Non-legacy
  modes are currently implemented only by `anaFemtoPhiTreeDownstreamV3`; direct `StFemtoMaker`
  rejects them until the validation gate is passed and the Maker/histogram rollout is completed.
- `rotationPairCutMode: preRotationLegacy|postRotation`. The former applies opening-angle and
  pair-rapidity cuts before rotation; the latter applies them independently to each rotated leg.

## Fully-mixed KK (`phi_mix`)

- **Species:** `phi_mix` / `particleKey: phi_fully_mixed` when `fullyMixedEnabled: true`.
- `mixChargeMode: combinedLegacy|bothSeparated`. The latter preserves the accepted directions as
  derived species `phi_mix_curkp` (current K⁺ × buffer K⁻) and `phi_mix_curkm` (buffer K⁺ ×
  current K⁻), while `phi_mix` remains their exact raw union under the same combined cap/sample.
- **Population:** production-PID daughters (`PassPhiDaughterTofPid`) in both directions combined — current K⁺ × buffered K⁻ and buffered K⁺ × current K⁻ — into one pair-index space. Eligible pairs also pass invMass>0, opening-angle, and pair-rapidity cuts.
- **Cap:** `fullyMixedMaxCandidates` is the stored-candidate count **per current event** on the combined set (not per direction). Production uses `2000`. `<=0` keeps every eligible pair (validation only; do not run uncapped on full production).
- **Sampling:** lazy Fisher–Yates permutation of 64-bit flat pair indices (`include/FemtoPhiMixSampler.h`). The first `cap` eligible pairs in that random order are stored, which is a uniform sample without replacement from the eligible set. No sampling weights / event reweighting.
- **RNG:** independent SplitMix64, **not** `gRandom`. YAML `fullyMixedSamplingSeed` (production `314159`). `>0` is fully reproducible for the same input; `0` is time-based and `Validate()` warns. Changing `rotationN` does not change the `phi_mix` sample.
- **QA:** `NPairPopulation` is the exact post-PID combination count. `NAttempted` is how many indices were evaluated. `EligibleLowerBound` is a lower bound when the scan stops at cap+1; `NEligibleExact` is filled only on a full scan. `hPhiMix_MKK` is stored candidates only (not a raw/unscanned MKK).
- **Tests:** `make test-phi-mix-sampler` (standalone). Does not change `FillMixedEventPairs`.

## Agent / developer checklist

1. Add species block with unique species key
2. Implement or extend builder dispatch for new `particleKey` in `StFemtoMaker`
3. Add channel entry referencing existing species keys
4. Add histogram entries to `config/hist/hist_anaFemto*.yaml` using channel name suffix
5. Document new keys in this file

## ⁴He bachelor (`anaFemtoPhi4He`)

- **Nuclear ID** (`nuclearid:` YAML → `StNuclearIdHelper`): TPC nσ + `requireBestSpecies`. Set `m2_selection: false` to skip TOF m² at ID time.
- **Femto TOF rule** (`he4TofMomentumThreshold` in maker YAML): tracks with `|p| >= threshold` require TOF m² in `[he4MinMass2, he4MaxMass2]`. Set threshold high (e.g. `99`) to disable.
- **Kinematics** (maker `he4*` keys): primary window `he4MinPMom` / `he4MaxPMom`; pT limits `he4MinPtPre`–`he4MaxPtPre` (pre QA) and `he4MinPtPair`–`he4MaxPtPair` (femto cut).
- **QA histograms:** `hHe4_Mass2VsP` (m² 0–16) plus `hHe4_Mass2VsP_*_wide` (m² 0–20) for TOF diagnostics when fills are enabled.

## Deuteron bachelor (`anaFemtoPhiDeuteron`)

- **Nuclear ID** (`nuclearid:` YAML → `StNuclearIdHelper::IsDeuteron`): TPC nσ + `requireBestSpecies` (best among d/t/³He/⁴He). Set `m2_selection: false` to skip TOF m² at ID time (initial default matches `anaFemtoPhi4He`).
- **Femto TOF rule** (`deuteronTofMomentumThreshold` in maker YAML): tracks with `|p| >= threshold` require TOF m² in `[deuteronMinMass2, deuteronMaxMass2]`. Set threshold high (e.g. `99`) to disable.
- **Kinematics** (maker `deuteron*` keys): primary window `deuteronMinPMom` / `deuteronMaxPMom`; pT limits `deuteronMinPtPre`–`deuteronMaxPtPre` (pre QA) and `deuteronMinPtPair`–`deuteronMaxPtPair` (femto cut).
- **QA histograms:** `hDeuteron_Mass2VsP` (m² 0–16) plus `hDeuteron_Mass2VsP_*_wide` (m² 0–20) for TOF diagnostics when fills are enabled.

See also `docs/femto_track_sharing_concerns.md` for overlap risks when extending species.

## Unified multi-bachelor (`anaFemtoPhi`)

- **anaName:** `auau3p85fxt_anaFemtoPhi` — one picoDst pass for φ–(p, d, t, ³He, ⁴He).
- **Variant `auau3p85fxt_anaFemtoPhi_dedxOnly`:** nuclei femto TOF m² off (`deuteron/triton/he3/he4TofMomentumThreshold: 99`, NuclearId `m2_selection: false`). **Proton:** `protonTofMomentumThreshold: 2.0` — low \(|p|\) TPC-only, high \(|p|\) TPC+TOF(m²). **Deuteron:** dE/dx-only. Mainconf: `main_auau3p85fxt_anaFemtoPhi_dedxOnly.yaml` (same macros/hist).
- **Species keys:** `proton`, `deuteron`, `triton`, `he3`, `he4`, `phi`, `phi_rot` (rotation for proton only).
- **Channels (21):** `phi_<bachelor>`, `_signal`, `_leftSB`, `_rightSB` per bachelor; plus `phi_rot_proton`.
- **Nuclear ID:** `StNuclearIdHelper::IsTriton()` / `IsHe3()` with `requireBestSpecies: true` (bestType 1=triton, 2=he3).
- **QA hist prefixes:** `hP_`, `hDeuteron_`, `hTriton_`, `hHe3_`, `hHe4_`; k* keys `hKstarSE/ME_<channel>`.
- **checkHist:** `checkHistAnaFemtoPhi.C` — bachelor-loop QA; single PDF (no CF-only PDF).

## DATA-006 p–p / p–d / d–d reduced-tree channels

The source-calibration channels use canonical track-species names
`proton_proton`, `proton_deuteron`, and `deuteron_deuteron`. They are built downstream from
the schema-3 `anaFemtoPhiTree` production, so event, centrality, track-quality, and PID decisions
are the stored decisions from the phi production. The producer mainconf and its 119-key
`FemtoFlagConfig` snapshot must match; the separate downstream-only configuration is
`config/maker/data006_auau3p85fxt_pp_pd_dd.yaml`.

- Identical SE pairs are unordered (`i < j`), excluding self-pairs and double counting.
- Identical ME pairs use current × buffer once; the reverse direction is not duplicated.
- p–d ME uses both current-p × buffer-d and buffer-p × current-d, following `mixBothDirections`.
- A shared `(eventUID, trackIndex)` is rejected when one physical track appears in both species.
- The baseline keeps the phi production's disabled close-pair veto. Uncut Δη–Δφ* SE/ME QA is
  still stored; enabling the DATA-006 switch applies the configured window to SE and ME alike.
- Native k* is 5 MeV/c. ROOT retains 0–3 GeV/c, fitter CSV covers 0–500 MeV/c, and normalization
  is channel-configured at 0.5–1.0 GeV/c.

Use `script/singularity_run_data006.sh` per reduced-tree block, merge its ROOT outputs with
`hadd`, then use `script/package_data006.sh`. Merged centralities add SE and ME before
normalization; finished CFs are never averaged.
For visual checks, run `script/singularity_plot_data006.sh` on the finalized package. Its display
rebinning recomputes the CF from rebinned SE and ME and never replaces the native 5 MeV/c data.

## KF Lambda–nucleus profile (`anaFemtoLambda_{d,t,3He,4He}`)

These entries use `StFemtoMaker` with an injected `FemtoLambdaProvider`. The neutral
interface carries no vendor KF/SIMD type. `createFemtoLambdaKfProvider()` lives in
`libStKfParticleCommon.so`; the Maker owns the provider. Phi/Kaon entries have no
provider and no mandatory KF linkage.

The basic profile has Lambda plus one nucleus and four mass-region channels
(inclusive, signal, leftSB, rightSB). `lambdaHe3ControlEnabled: true` opts into
exactly `lambda,he4,he3` and eight channels. Each accepted He3 pair also fills a
`he3_as_he4` comparison using He4 mass and unchanged physical He3 momentum. This
is a derived representation, not another species or mixing channel. It changes
k* and mean-mass mT; track provenance and Lambda mass-region membership remain
the same. The Maker applies no contamination coefficient or CF subtraction.

### PID and acceptance

The standard 2021 P24iy entry is
`config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle.yaml`.
Its explicit `nuclearid` selects a He4-only log10(R) polynomial calibration:
cubic mean, quadratic Gaussian width. d/t/He3 use `NuclearIdDeDxVsMom.h`
unchanged. The schema-2 YAML includes all coefficients, units, domain and
provenance; source paths are identification records, not runtime dependencies.
Missing or invalid explicit calibration fails Init. Outside its declared
calibration domain, the configured legacy policy applies.

TPC best-species selection precedes rigidity and TOF cuts. A rejected winner is
never reassigned. For all four nuclei, R<2 GeV/c requires no TOF; R>=2 requires
valid TOF plus the species' strict mass-squared window. He4 uses the deuteron
window because both have similar m/q. These gates use raw rigidity and (m/q)^2;
He physical momentum has its Z=2 correction applied exactly once. He4 requires
strictly 0.61<R<5 GeV/c with no He4 rigidity veto. d/t/He3 retain their legacy
minima and vetoes. All numerical cuts live in nuclearid YAML.

`lambdaEventPolicy: lambda_imp5_compat` uses finite PV/maxNTr and centrality;
it does not add generic Phi primary-track cuts. `qaVertexByMode` sets QA centers,
not event acceptance. KF Finder/Topo sees the full PicoDst and shares its final
selector with `StLambdaKFParticleMaker`. Daughter provenance uses Pico array
indices and event indices. `lambdaPairMassMode: fixed` uses the configured
Lambda mass for kinematics while `reso.invMass` retains unconstrained KF mass.
Signal is +/-3 configured sigmas, including both edges; sidebands exclude the
adjacent signal edge. Mixing uses the configured Vz/cent9/EP bins, pool depth,
Vz clamping and directions. Shared-track rejection remains mandatory.

Single-candidate acceptance maps fill once per selected candidate, independently
of pair multiplicity. Lambda signal and He4 have beta vs lab eta and lab Y;
beta is physical |p|/E, not measured TOF beta. He4 beta/eta, beta/Y and pT/Y maps
also have `lowR_noTof`, `lowR_tofMatch`, `lowR_tofM2`, `highR_tofM2` samples.
The three low-R samples are nested, so do not add them. For each variable,
lowR_noTof + highR_tofM2 equals inclusive, including flow bins.

### He4 lab-rapidity and track-merging diagnostics

The diagnostic mainconfs retain the physical Lambda/nuclear selections and add
pair-level QA:

| Mainconf suffix under `config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle_` | Purpose |
|---|---|
| `epd8_ylab_merging.yaml` | EPD selection and eight EP mixing bins, with the diagnostics below |
| `epd1_ylab_merging.yaml` | The EPD-selected one-bin mixing control, with the same diagnostics |
| `diagnostics_noEp.yaml` | Rapidity and merging diagnostics with no EPD event selection |

The EPD source and validity contract is separate from these diagnostic fills.
Rapidity splitting does not select whole events: one event can contain He4
candidates on both sides, and each Lambda-He4 pair follows its own He4 candidate.
The optional maker block is explicit:

```yaml
he4RapiditySplit:
  enabled: true
  yLabThreshold: -0.5
```

It uses the selected physical He4 candidate's stored four-vector in the lab
frame. SE and both ME directions use that same candidate-level definition.
Strictly smaller and strictly larger rapidities fill `_YLabBelow` and
`_YLabAbove`; exact equality and invalid rapidity fill neither diagnostic.
The existing inclusive pair fill is retained. The He3-as-He4 mass representation
does not enter this He4 diagnostic.

The 32 added keys follow
`hKstar{SE,ME}[VsCent]_lambda_he4{,_signal,_leftSB,_rightSB}_YLab{Below,Above}`.
Their k* axes match the corresponding inclusive spectra exactly, and `VsCent`
retains cent9. The ROOT output contains `TParameter<double>
he4RapiditySplitYLabThreshold` and `TNamed FemtoHe4RapiditySplitDefinition`.
An enabled block requires the complete histogram set and matching axes; putting
these histograms in a profile with the feature disabled is an Init error.

All Lambda profiles require the corrected merging geometry explicitly:

```yaml
mergingRadiusMeters: 1.4
mergingGeometry:
  method: starHelixFixedGlobalEta
  tpcHalfLengthCm: 210.0
  radiusToleranceCm: 0.001
  includePion: true
```

`includePion: false` permits the existing proton-only histogram inventory.
The method and numerical geometry fields are required; missing or unsupported
settings fail Init. Every production path uses `StPicoTrack::helix` with the
event's kG field and the global track momentum. A primary-track fit is not
required for a Lambda daughter. The helix already represents curvature with
the stored charge and rigidity, so no extra helium Z=2 scaling is applied.

The calculation takes the first finite nonnegative intersection with the common
cylinder, rejects STAR's no-solution value, checks the radial residual and the
global TPC z extent, then fills wrapped daughter-position phi minus
nucleus-position phi against their global-momentum eta difference. There is no
asin clamp or fixed-field approximation. The bare
`hDphiDeta_proton_<representation>` and `_EventField` names are compatibility
aliases of this same corrected definition. With pion output enabled,
`hDphiDeta_pion_<representation>` and its `_EventField` alias are also required.
The representations are `4He`, `3He` and, when configured, `3HeAs4He`.

These maps use signal SE pairs with the current PicoDst available. They do not
veto close pairs or alter physical SE/ME acceptance. Pair replay cannot rebuild
track geometry from candidate four-vectors. A cylinder intersection is a helix
prediction, not a measurement of shared hits or merging efficiency. The output
`TNamed FemtoLambdaMergingDefinition` records the geometry and alias contract.

The dedicated histogram profiles retain every baseline key, class and axis and
append the diagnostics. Existing merging-map bins change to the corrected
definition. To regenerate the dedicated profiles from their explicit maker
settings without replacing the base profiles, run
`tools/kfparticle_summary/generate_lambda_diagnostic_histograms.py` with required
`--base-mainconf` and `--mainconf` arguments. Existing generated output requires
an explicit `--overwrite`.

### Single-He4 rapidity histograms

The optional maker block retains a binned view of every accepted physical
Lambda-He4 pair without changing the pair or pool selection:

```yaml
he4RapidityBinning:
  enabled: true
  frame: lab
  nBins: 8
  min: -1.0
  max: 0.0
```

Every key is required when the block is present. The interval uses equal-width
ROOT bins `[low, high)`; finite rapidities below `min` or at/above `max` remain
in underflow/overflow. Invalid nucleus rapidity skips this diagnostic only.
The example gives eight intervals of width 0.125. Add adjacent pairs of raw
SE and ME intervals to obtain four intervals of width 0.25, then normalize
each selected interval. Do not average separately normalized CFs. Include
rapidity flow bins when checking closure against the inclusive spectrum.

Eight TH3D keys follow
`hKstar{SE,ME}VsHe4YLabVsCent_lambda_he4{,_signal,_leftSB,_rightSB}`.
Their axes are k*, physical-He4 lab rapidity and cent9. The k*/cent9 axes
match the inclusive pair spectra exactly. Rapidity comes from the nucleus
candidate's stored physical four-vector in SE and both ME directions, with
no extra charge correction or CM shift. Lambda-He3 and He3-as-He4 do not
enter this view. The existing strict two-region He4 split and the separate
pair-rapidity diagnostic retain their own definitions. In particular, an
exact threshold value is included in the upper native rapidity interval,
while the strict two-region split excludes exact equality from both sides.

An enabled block requires the complete eight-histogram profile and matching
axes; a disabled/absent block rejects profiles containing these keys.
`he4RapidityLabNBins`, `he4RapidityLabMin` and `he4RapidityLabMax` are typed
ROOT parameters in first-value merge mode.
`FemtoHe4RapidityBinningDefinition` records the physical definition and flow
policy. Existing reduced event trees retain the information needed for this
diagnostic; changing these histogram bins does not change reconstruction
provenance or require another PicoDst pass.

### Lambda-He4 pair-rapidity diagnostics

The pair-rapidity diagnostic uses the laboratory rapidity of
`lambda.P4() + he4.P4()`, from the selected physical candidates. It applies no
CM shift and does not average the two particle rapidities. The explicit maker
configuration is:

```yaml
pairRapiditySplit:
  enabled: true
  frame: lab
  nBins: 10
  min: -1.5
  max: 0.5
```

The ten equal-width intervals have width 0.2. The eight `TH3D` keys are
`hKstar{SE,ME}VsPairYVsCent_lambda_he4{,_signal,_leftSB,_rightSB}`, with axes
`(k*, pairYLab, cent9)`. Their k* and cent9 axes match the existing inclusive
spectra. SE and both ME directions compute rapidity from the same pair-level
four-vector definition. He3-as-He4 does not enter these histograms.

This diagnostic is independent of the single-He4 `_YLabBelow`/`_YLabAbove`
views: neither is an additional cut on the other. It leaves inclusive raw
SE/ME counts, pair acceptance and mixing unchanged. No pre-normalized CF is
stored in the pair-rapidity maps. Adjacent intervals can be added separately
for SE and ME before the chosen CF normalization procedure is applied.

ROOT bins use `[low, high)`. A rapidity equal to 0.5 enters overflow; finite
values outside the configured range remain in underflow or overflow. Nonfinite
sums, `E <= abs(Pz)`, and nonfinite computed rapidities skip only this
diagnostic. The output retains cent9 and k* flow bins, so a rapidity projection
can include or inspect the flow regions explicitly.

`FemtoHe4PairRapiditySplitDefinition` records the definition. The configuration
values `he4PairRapidityLabNBins`, `he4PairRapidityLabMin` and
`he4PairRapidityLabMax` are ROOT parameters with first-value merge mode, as is
the single-He4 `he4RapiditySplitYLabThreshold`. Worker merging preserves these
constants instead of adding them. Init requires the complete eight-map set,
`TH3D` storage and matching axes. A disabled or absent feature does not accept
pair-rapidity histogram definitions.

The QA macro draws raw k* spectra for all ten intervals and reports the
rapidity flow counts. Its `FemtoLambdaQa::ProjectPairRapidity` helper takes an
inclusive rapidity-bin range and inclusive cent9 values. Use it for an exact
underflow-only projection on ROOT 5, where the native projection range `(0,0)`
otherwise denotes the complete axis. Examples, including adjacent-bin merging,
are in [REFERENCE](../../docs/REFERENCE.md#lambda-he4-pair-rapidity-diagnostics).

The pair diagnostic regression passes 1,603,559 checks, including 1,140,649
pair-rapidity checks for definition, boundaries, flow bins, both ME directions,
He3 exclusion and metadata. A 10,000-event EPD-selected sample retains all
1,707 existing pair-physics histograms and 40,002,230 cells exactly, including errors,
entries and statistics. Its physical He4 signal has SE 13 and ME 403 entries;
the new pair-rapidity projections retain those counts, with no rapidity flow
entries in that sample. A two-worker merge retains the metadata constants.

### Optional nuclear pair PID QA

`nuclearPidPairQaHist` enables a separate dense-QA HistManager and pre-PID
QA-candidate storage. Absent key means no extra QA histograms or QA-candidate
storage. The referenced QA YAML must explicitly choose its 8D output:

- A complete `jointSparse` mapping creates `THnSparseD hNuclearPidPairJoint`.
- `jointSparse: false` disables only that histogram; dense QA and cached pre-PID
  candidates remain available.
- Missing, null or `true` values are errors; there is no implicit output default.

The [no-Sparse writer mainconf](../../config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle_candidateStoreWriteNoJointSparse.yaml)
references [the dense-only QA YAML](../../config/hist/hist_auau3p85_femtoLambda_he4_he3_pid_pair_qa_noJointSparse.yaml).
It sets `jointSparse: false` and retains the full 8D definition as comments.
Candidate writing stays enabled and physical cuts remain unchanged.

`selected` uses final nuclear PID. `prePID` removes only the species' own TPC
pull window and best-species competition, retaining track quality, permitted
rigidity domains and applicable TOF cuts. Selected is a subset of prePID; He3
and He3-as-He4 share the same tracks and pairs. Neither pair of representations
is statistically independent. Pull is nSigmaHe4 for He4 and nSigmaHe3 for both
He3 representations. Pair-QA rigidity comes from stored Float_t physical P4/2;
selection boundaries use the original double-precision track quantities.

QA payloads accompany the existing physical mixing pool. QA-only events do not
enter or evict that pool, or change CF eligibility, pair counters or randomness.
SE TH2D keys live in `true/`, ME in `mix/`. The 1,080 TH2D retain k* vs pull
and k* vs R by representation, mass region, centrality and PID stage.
When enabled, `hNuclearPidPairJoint` stores axes
`(Kstar,POverQ,NSigma,Cent9,Representation,MixingMode,MassClass,PidStage)`.
MassClass is exclusive: SBNeg, signalLow, signalHigh, SBPos. Full signal is the
sum of its two classes. The TH2 projections alone cannot recover the R/pull
correlation; use the Sparse for simultaneous cuts. Sparse size grows with
occupied bins and can exceed ROOT's single-object serialization limit. Preserve
it in separate input-group files and sum matching low-dimensional projections
when one merged object would exceed that limit.

The full-QA configuration defines 1,786 ordinary histograms plus one Sparse.
The no-Sparse configuration retains the ordinary histograms and omits the Sparse.
The [complete histogram catalog](../../mdfiles/femto/implementation/femto_lambda_histograms_20260930.md)
lists names, classes, axes, sample definitions, filling stages and metadata.
The calibration and selection policies are also saved as ROOT TNamed records.

Build with `script/singularity_make.sh MAINCONF`. Regression entry points under
`tests/` cover calibration validation, legacy coefficients, exact TOF boundaries,
He3 control k*, acceptance sample identities, Sparse projections and QA-pool
independence. Run commands are in [REFERENCE](../../docs/REFERENCE.md).


### Xialei EPD first-order event plane

The optional top-level `eventPlane` reference selects
`config/event_plane/lambda_he4_xialei_epd8.yaml`. This computes the corrected
East EPD AB (rows 1–8) first-order plane, with the source recentering and
20-term shift calibration. The mixing class reproduces the original Float_t
`psi1 / 2` conversion and eight-bin calculation. It is separate from the
existing auxiliary TPC second-order plane. Event-keyed tile sampling preserves
the caller's random stream; original event-by-event angles are not claimed to
be bitwise identical without the original random sequence.

The supplied calibration covers **388 runs, 22158032–22179022**, within the same
2021 P24iy fixed-target production. It does not cover the entire current
analysis input. Unknown runs, missing calibration objects, insufficient bin
counts and empty required calibration cells fail explicitly. Use covered-run
input for the EP profiles; additional calibration is required for other runs.
See [calibration provenance](../../config/event_plane/xialei_epd_2021P24iy/README.md).

The canonical stored/calibrated EPD class remains eight bins. Mixing uses the
physical first-order angle with `mixing.nEventPlaneBins`; one bin keeps all
EP-valid events in the same angular class. Partitions of 1, 8, 16 and 32 bins
use the same validity selection and preserve the original Float_t boundary
convention. They do not require new recenter/shift calibrations. Keep cent9,
Vz, both pool-depth limits and the admission policy identical when comparing
these partitions. The auxiliary TPC plane is never used as an EPD fallback.

Inclusive raw SE/ME histograms sum the classes; optional k* × EPD × cent9
histograms retain each diagnostic class. The active pool count is separate
from canonical eight-bin calibration QA. `hEpdPoolBinVsCent`, when configured,
must match the active partition. CF normalization and background processing
remain outside the Maker.

Outputs include `hEpdPsi1VsCent`, `hEpdMixingBinVsCent`, `hEpdStatus`, the
`FemtoEpdEventPlaneDefinition` record and seen/valid/rejected event counters.
`writeDetailedLedger: false` avoids per-event EP branches in the standard
histogram-only run. Candidate-store writing, when explicitly enabled with EPD,
uses schema 2 and separate EPD fields; ordinary no-EP schema 1 is unchanged.
Replay of EPD analysis requires those fields and matching calibration provenance.

### Lambda–He4 close-pair selection

The optional maker `mergingPairCut` block requires explicit `enabled`,
`nucleus: he4`, `daughter: proton`, `phiMax`, `etaMax`,
`unreachablePolicy: keep` and `writeQA` values. Widths are positive, finite
limits in radians and pseudorapidity. A physical Lambda–He4 pair is rejected
only when both strict inequalities hold. The same predicate precedes every
SE/ME mass-class and differential pair fill; it does not reject events,
candidates or pool registrations. Lambda–He3 and its alternative He4 mass
hypothesis remain independent of this selection.

The geometry uses the configured comparison radius, each track's origin-event
magnetic field, raw rigidity/charge and global momentum pseudorapidity.
Intersection positions define wrapped daughter-minus-He4 azimuth. Mixed pools
own the precomputed track points, including in reverse mixing. Unreachable
tracks are retained and counted; a missing track or mismatched event cache
fails explicitly. No fixed-field fallback is used for the selection.

With `writeQA: true`, proton and pion signal-pair maps are written before and
after selection for SE and ME. `hDphiDetaPairCutQA` records both populations,
rejections and unreachable tracks. Pion maps diagnose the effect of the proton
selection; pion proximity does not reject a pair. Existing SE merging aliases
retain their population before selection. Output metadata records the cut
parameters, geometry and reachability policy. Without the block, the original
pairing path and QA population are preserved.
