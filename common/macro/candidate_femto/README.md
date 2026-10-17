# candidate_femto: analysis macros for StCandidateFemtoMaker output

Macros that turn the pair histograms of `StMaker/StCandidateFemtoMaker` (THnSparse
`hPairSE_<channel>` / `hPairME_<channel>`) into correlation functions and QA. They read the same
`candidateFemto` YAML as the Maker (via `StMaker/StCandidateFemtoMaker/CandidateFemtoConfig.h`), so
windows, sidebands, normalisation range and purity-fit settings have one source. Correlation
functions are computed here, never in the Maker.

| File | Content |
|---|---|
| `CandidateFemtoProjection.h` | projections of the THnSparse with a cent9 list and BDT score thresholds: k* per mass-class combination, (k*, M_A, M_B) |
| `purityCorrectedCorrelation.C` | CF of one channel with the k*-differential purity correction, plus constant purity, sliding-window purity and merged k* bins where a fit fails (final result) |
| `checkHistCandidateFemto.C` | QA and CF pages for every channel (one combined PDF) |

Runners and wrappers: `analysis/run_purityCorrectedCorrelation.C` /
`script/singularity_purityCorrectedCorrelation.sh`, `analysis/run_checkHistCandidateFemto.C` /
`script/singularity_checkHistCandidateFemto.sh`. Every input (ROOT file, mainconf, output path) is an
argument; nothing is read from the environment. Analysis settings are keys of the `candidateFemto` YAML (no
numbers built into the macros): `cfNormMin/Max`, `qaCfKstarBinWidth`, `qaClosePairMaxDeltaEta`,
`purityFitMinEntries`, `slidingWindowHalfWidth`, `purityCorrectionMaxMergedBins`, and per species
`purityFitRange`, `purityFitSeed`, `qaPurityFitSigmaRange`, `qaPurityBackgroundOrder` (documented in
`CandidateFemtoConfig.h`). Particle masses come from `TDatabasePDG`. What remains in the code: fit start
values and limits, the cent9 → percentile table (a specification), and layout (text positions, plot
headroom, rebinning of QA displays).

## Purity correction

```bash
./script/singularity_purityCorrectedCorrelation.sh <file.root> <mainconf> <channel> <out/prefix> <kstarRebin> [cent9List] [options]
# e.g. 0-40 % (cent9 4-8), BDT thresholds from species_<s>_scoreMin of the YAML:
./script/singularity_purityCorrectedCorrelation.sh merged.root config/mainconf/main_<anaName>_bdtSelection_<thresholds>.yaml <channel> out/<channel> 2 "4,5,6,7,8"
```
Options (`;`-separated, scans only; a chosen value belongs in the YAML): `purityA=config|fit|unit`, `windowA=lo,hi`,
`slideHalfWidth=N`, `scoreMinA=x` / `scoreMinB=x` (BDT thresholds of species A / B; must be score-axis edges).

With BDTs, the job runs a `_bdt` mainconf (scores on the score axes, no thresholds). A `_bdtSelection_<thresholds>`
mainconf differs from it only by `species_<s>_scoreMin` in its candidateFemto YAML; give it to the macros together with the
output of the `_bdt` job.

## Choosing BDT thresholds

A threshold is an analysis decision, made here from the stored score axes, never in the Maker and never by rerunning a job.

1. Run the CF macro on the merged `_bdt` output for a grid of thresholds with `scoreMinA` / `scoreMinB` (every value must
   be a score-axis edge; a grid of 0.1 is enough), same centrality and k* binning as the final result.
2. For each point record: SE and ME purities of both species, the number of k* bins whose purity fit failed, the pairs in
   both signal windows, and the final CF (merged bins) averaged in a few k* ranges with its statistical error.
3. Choose the loosest thresholds that keep both SE purities high (≥ about 0.94), with few failed purity fits and close
   to the smallest statistical error, where the CF is stable against the neighbouring thresholds. Tighter cuts make SE and
   ME purities closer (less reliance on the correction) but cost pairs and make the per-k* purity fits fail; looser cuts
   make the result rely more on the correction.
4. Write the choice as a `_bdtSelection_<thresholds>` YAML/mainconf and keep the scan as a study (table, summary figure,
   driver scripts), so the spread over neighbouring thresholds can later enter the systematic uncertainty.

Normalisation: `purityCorrectedCorrelation.C` normalises each CF by the inverse-variance weighted mean of SE/ME in
`[cfNormMin, cfNormMax]`; `checkHistCandidateFemto.C` (QA) by ΣME/ΣSE in the same range. Purities are per channel
(each correlation function is computed independently; the Maker's shared-track rule acts per channel):
`purityCorrectedCorrelation.C` fits the mass projections of the channel's pair THnSparse (k*-differential and
integrated); `checkHistCandidateFemto.C` fits the mass of each species in the channel's SE pairs whose partner is in
its signal window, with the score thresholds of both species. The channel candidate histograms
`hChannelMassVsCent9[VsScore]_<channel>_<s>` are not used for purities: they hold every candidate of the events with
both species, without the partner's threshold, and their purity can be much lower than in the pairs. Without
`hPairSE`, or with a threshold that is not a pair score edge, checkHist falls back, with a warning and a "FALLBACK"
label, to the per-species spectra (`hCandidateMass*_<s>`), which the species pages show for information. One species and one
charge each (Ξ⁻ and Ξ⁺ are separate species); the reconstruction QA histograms of the KF part are split by charge
too (`hXiMinus_*`, `hXiPlus_*`).

- `species_<s>_purityModel` decides the fit of each species: `invariantMassGausPol2`, `invariantMassStudentTPol2`
  (Student-t peak, for peaks with non-Gaussian tails that a Gaussian counts as background),
  `invariantMassDoubleGausPol2` (purity macro only) or `tofMassSquaredKaonPion`. `purityByDefinition: true` sets the
  purity to 1 (a species whose selection gives a clean sample); the fits are then drawn as QA only. Species B needs an invariant-mass
  model. checkHist supports the Gaussian and Student-t models (its "±nσ" lines use the Student-t scale for the latter).
  The shape parameters (Student-t ν) come from the k*-integrated fit and are fixed in the k*-resolved fits.
- Fit ranges: `species_<s>_purityFitRange`.
- SE and ME purities are fitted separately: the ME events are not the SE events (see the Maker README, "ME purity
  differs from SE purity"), so one purity for both would bias the CF.
- TOF m² model (`purityA=fit` for a TOF-identified species): mean and width float in every k*-resolved fit (per bin, sliding window,
  merged bins), started from the k*-integrated fit of the same sample; a fit is valid only if it converged and
  its width is not at a limit (0.3 or 3 times the start value).
- Output names (`<prefix>.root`): `cfRaw`, `cfConstantPurity`, `cfPurityCorrected`,
  `cfSlidingWindowPurity`, `cfPurityCorrectedMergedBins` (final), purities `purityA_*` / `purityB_*`
  (SE, ME), `puritySlidingA_*` / `puritySlidingB_*`, pairs in both windows `pairsInWindows_*`, and every
  fitted projection (`fits/`, `fits_sliding/`, `fits_merged/`).
- Validation: on toy data with the real purities and statistics, the merged-bin correction recovers the input CF within
  its errors.

## QA pages (checkHist)

```bash
./script/singularity_checkHistCandidateFemto.sh <file.root> <mainconf> <out.pdf> "8;6,7;4,5"   # cent9 groups, required
```
