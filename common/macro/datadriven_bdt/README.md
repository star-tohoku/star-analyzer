# datadriven_bdt: data-driven BDT (mass window vs side bands), Python + xgboost

This toolkit trains, applies, evaluates and checks a BDT that selects candidates forming a mass peak (K0S, Λ, Ξ⁻, Ω⁻, ...) without Monte Carlo.

- **Training sample.** Candidates in the mass window (signal + background) are trained against candidates in the side bands (background only). This is the "CWoLa" approach.
- **Channel specifics.** Everything channel-specific lives in one YAML recipe; a new particle needs a new recipe, not new code.
- **Folds.** Two folds by run parity: fold f is trained on runs with `runId % 2 == f`, and every candidate is scored with the fold that did not see it. `StCandidateFemtoMaker` applies the model the same way.
- **Learner.** xgboost (binary logistic). The trained trees are exported to a plain text file that the Maker evaluates itself (`StMaker/StCandidateFemtoMaker/CandidateTreeEnsemble.h`), so the Maker does not need xgboost.

Scope:

- **Particles with a mass peak only** (decays reconstructed with KFParticle: K0S, Λ, Ξ⁻, ...). The labels come from the
  mass, so a species without a peak (an identified track such as K⁺) is not a case for this toolkit; select it with PID
  cuts (a BDT for it would need labels other than the mass).
- **Data only, no Monte Carlo.** The training needs nothing but the candidate trees of the analysis.
- **Every variable of the part as input** (`inputs.use: all`) unless a check shows a reason to drop one; the SHAP and
  sculpting pages of the report are where such a reason shows up.
- **Decide the training statistics with the `study` action** before a full-statistics training.

## Files

| File | Role |
|---|---|
| `datadriven_bdt.py` | Recipe reading and checking; actions `inspect`, `train`, `apply`, `study`, `report`, `model`, `check` |
| `candidate_expression.py` | The expression grammar of the Maker (`CandidateExpression` in `CandidateBdtScorer.h`), ported line by line, so the Maker and the toolkit compute the same numbers |
| `tree_ensemble.py` | Export of the xgboost trees to the tree-ensemble file, and the reference evaluation of that file (the score definition) |
| `datadriven_bdt_report.py` | One PDF with bookmarks: every evaluation item of the recipe |
| `examples/template_new_particle.yaml` | Starting point for a new particle: replace every `EDIT` (a recipe that still has one is refused) |
| `examples/{xi_minus,k0s}_candidate_trees.yaml` | Worked recipes (Ξ⁻, K0S of the 3 GeV analysis) |
| `requirements.txt` | The Python packages (tested versions); see "Python environment" |
| `../../../script/run_datadriven_bdt.sh` | Entry point (host Python; removes STAR's `PYTHONPATH`) |
| `../../../script/run_datadriven_bdt_chain.sh` | inspect → train → apply → study → report for one recipe, stopping at the first failure; prints the `model` command to run after reading the report |
| `../../../config/bdt/<anaName>/` | Production model files (`<name>.model.yaml`) and tree-ensemble files |

## Usage

```bash
./script/run_datadriven_bdt_chain.sh <recipe.yaml>   # the usual sequence in one go (below: the single actions)
R=common/macro/datadriven_bdt/examples/xi_minus_candidate_trees.yaml   # copy and edit; paths relative to the repository root
./script/run_datadriven_bdt.sh $R train            # both folds -> <output>/<name>_f{0,1}.{xgboost.json,trees.txt,training.json,shap.npz}
./script/run_datadriven_bdt.sh $R apply            # every candidate scored by the other fold -> <output>/<name>_apply.root
./script/run_datadriven_bdt.sh $R study            # training-size study (if evaluation.checks has trainingSize) -> <output>/<name>_training_size.json
./script/run_datadriven_bdt.sh $R report [label=other_apply.root ...]   # -> <output>/<name>_report.pdf
./script/run_datadriven_bdt.sh $R model <dir>      # model file for StCandidateFemtoMaker (after validation)
./script/run_datadriven_bdt.sh $R check '<dir>/*.tree.root'   # Maker scores of that model == this model's scores
./script/run_datadriven_bdt.sh $R inspect          # files, candidates and memory needed, before reading anything
```

The whole procedure (training trees → model → Maker) is `docs/ai/skills/candidate-bdt.md`; the design and how to add a BDT for a new particle: `docs/ai/candidate-bdt-guide.md`.

## Python environment

The toolkit runs on the host (login or interactive node), not in the STAR container and not as a batch job: batch jobs only
read the exported tree-ensemble files, which the Maker evaluates in C++. Set the environment up once:

```bash
python3 -m venv ~/.venv                      # Python 3.9 (tested); any location works
source ~/.venv/bin/activate
pip install -r common/macro/datadriven_bdt/requirements.txt   # numpy, xgboost, shap, uproot, awkward, scipy, matplotlib, PyYAML
```

Activate the environment before running `./script/run_datadriven_bdt.sh`: the script uses the `python3` found on `PATH` and
removes STAR's `PYTHONPATH`, which would otherwise shadow these packages.

**Run time and memory** (3 GeV Au+Au, full statistics, one merged tree file of 21 GB on GPFS, 8-core login node):

| | candidates read | train (2 folds + SHAP) | apply | study | peak memory |
|---|---|---|---|---|---|
| Ξ⁻, all events | 0.67 M | 1.3 min | 1 min | 0.5 min | 0.6 GB |
| K0S, `eventFraction: 0.01` | 2.9 M | 10 min | 4 min | 1.3 min | 2.0 GB |
| K0S, `eventFraction: 0.03` | 8.6 M | 18 min | 1 min | | 3.7 GB |

The K0S BDT reaches the same AUC at 0.01 and 0.03 (and as the production model). Reading the same Ξ⁻ candidates from 646
files took 2.4 min instead of 1.3 min, and from 34661 per-PicoDst-job files about 40 min: opening files costs, not data.

**Memory check.** Before reading each file the toolkit adds up the candidates it will keep and stops with a clear message if the
expected peak (1.5 GB + 4 × the candidates read, from the measurements above) exceeds the memory available at the start
(`MemAvailable`); `inspect` prints the same estimate. On a shared node keep a margin for other users.

So:

- Merge the training trees (`hadd`, as for histograms) or write them with the pair stage on candidate stores (one tree file
  per pair-stage job); do not read thousands of per-PicoDst-job files. One merged file per species is fine.
- Before a large training run `inspect` the recipe: it prints the number of files and candidates and the memory the read
  candidates take (training needs a few times that). If it is too much, or more than the training needs (check with `study`:
  the K0S BDT of the 3 GeV analysis reached its plateau at about 70k window + side-band candidates per fold), lower
  `data.eventFraction`.
- The candidates read by `train`, `apply` and `study` are cached in `<output>/<name>_candidates_cache.npz` and reused while
  the tree files (path, size, modification time), the tree and the needed branches are unchanged; anything else triggers a
  new read. The cache is about the size of the read arrays and can be deleted at any time.
- To use part of the input, lower `data.eventFraction` (not a subset of files: the sample then depends on the file layout).

## Recipe

Every key is required (no built-in defaults); an unknown key stops the program. Example: `examples/xi_minus_candidate_trees.yaml`.

| Key | Meaning |
|---|---|
| `name` | Name of the model and its files (no `%`, `/`, spaces) |
| `particle` | Decay, reconstruction and any training condition the trees cannot express (goes into the model file) |
| `validation` | Path of the report that validated the model, or `"none: <reason>"` (goes into the model file) |
| `data.trees` | Candidate-tree files; a wildcard only in the file name |
| `data.tree` | Tree name, e.g. `candidates_xi_minus` |
| `data.eventSelection`, `data.preselection` | Expressions; together they define the training domain. `"1"` = none. The Maker gives every candidate outside the domain the outside-training score |
| `data.trainingSelection` | Expression applied to the **training samples only** (train, training-size study), never to the domain the Maker scores; `""` = none. Use it to drop candidates the shared-track rule removed in the channel(s) the model serves, e.g. K0S: `sharedTrackRemoved_k0s_xi_minus==0`; Xi-: `sharedTrackRemoved_k0s_xi_minus==0 && sharedTrackRemoved_kaon_plus_xi_minus==0` (one Xi- model scores Xi- for every channel, so it is trained on Xi- that every channel keeps) |
| `data.folds` | `runParity` (the only scheme) |
| `data.eventFraction` | Fraction of the events of the input used by `train`, `apply` and `study` (`1` = all). Events are chosen by a hash of (`runId`, `eventId`): the same events whatever the file layout, all candidates of an event together, independent of the folds. Read in chunks, so only the sample is kept in memory. Not applied by `check` or by the Maker |
| `labels.method` | How signal and background are defined. `massSidebands` (the only method; a new method goes into `Labels` in `datadriven_bdt.py`) |
| `labels.mass`, `labels.window`, `labels.sidebands` | Mass expression, `[lo, hi]`, `[[lo, hi], ...]` (open intervals) |
| `inputs.use` | `all` = every source variable of the reconstruction part (all tree branches except `runId eventId cent9 mass px py pz score` and the flags `sharedTrackRemoved_*`), or a list of expressions. The flags are bookkeeping and can never be inputs |
| `inputs.exclude` | Inputs removed from `use` (e.g. `[massError]`) |
| `inputs.log10` | Inputs replaced by `log10(max(x,1e-6))` (heavy-tailed variables) |
| `inputs.massCorrelationWarning` | Warn when an input's correlation with the mass in the side bands exceeds this (sculpting risk) |
| `reweight.variables`, `reweight.bins` | Two expressions and `[[n, lo, hi], [n, lo, hi]]`: side-band candidates are reweighted to the window distribution in these, per fold. `none` / `none` = no reweighting |
| `learner.type` | `xgboost` |
| `learner.seed` | Seed of the train/test split, of the learner and of the samples drawn by `study` and the SHAP sample |
| `learner.testFraction` | Fraction of each fold kept out of the training (overtraining check) |
| `learner.parameters` | Passed to `xgboost.XGBClassifier` (objective fixed to `binary:logistic`; no `random_state`, use `seed`) |
| `evaluation.massRange`, `evaluation.binWidth` | Mass range and bin width of the apply output and the fits |
| `evaluation.purities` | Purities of the summary table, the spectra and the checks |
| `evaluation.sculptingPurity` | Purity of the ε(m) and equal-purity ratio pages |
| `evaluation.windowScan` | Fractions of the window half-width for the window-width scan |
| `evaluation.checks` | Any of `sculpting overtraining importance backgroundModels trainingSize` |
| `evaluation.trainingSize` | `{fractions: [...], purity: p}`: training on these fractions of each fold; S at purity p and AUC on the other fold |
| `evaluation.purityBins` | Bins of `pt`, `y` (laboratory rapidity), `cent9` for the purity-in-bins page (any subset; `{}` = none) |
| `evaluation.shapCandidates` | Number of test candidates per fold for the SHAP values |
| `output` | Output directory, outside the tracked directories (e.g. a study directory) |

`particle`, `validation` and `evaluation` may change after training; any other change requires training again (`model` refuses otherwise).

**Expressions** (selections and inputs) are those of the Maker: `+ - * /`, comparisons, `&& || !`, `log10 log exp sqrt abs max min`, numbers and tree variables; mixed `&&`/`||` need parentheses; comparisons do not chain. Edge values follow ROOT 5.34 TTreeFormula (`log` of x ≤ 0 = 0, `sqrt(|x|)`, `exp` clamped at ±700, `x/0 = 0`). Inputs may not use `mass`, `runId`, `eventId` or `score`; inputs using `px py pz` give a warning (keep pT out, see the lessons).

## Outputs and the model contract

- **`train`** writes per fold: the native xgboost model (`.xgboost.json`), the tree-ensemble file (`.trees.txt`, what the Maker reads), a record (`.training.json`: recipe, git revision, candidate counts, AUC on training and test candidates, score histograms, gain, mean |SHAP|, correlation of each input with the mass), and SHAP values of a test sample (`.shap.npz`). It stops if the exported file does not reproduce xgboost to 1e-5.
- **Cache.** `train`, `apply` and `study` keep the candidates they read in `<name>_candidates_cache.npz` (see Run time).
- **`apply`** scores with the saved xgboost model of each fold (fast) and checks a sample of 2000 candidates per fold against the
  tree-ensemble file, the score definition (stops above 1e-5; typical difference 1e-6). `check` uses the tree-ensemble file itself.
  It writes `<name>_apply.root`: `mass_vs_score` (TH2D) and the tree `candidates` (score, mass, runId, cent9, px, py, pz) of every candidate in the domain and the mass range.
- **Score.** margin = sum of the leaves + base margin; score = 2 / (1 + exp(−margin)) − 1 ∈ [−1, 1]. The score is defined by the tree-ensemble file (`tree_ensemble.py`); the Maker evaluates the same file (`CandidateTreeEnsemble.h`).
- **`model <dir>`** writes `<dir>/<name>.model.yaml` and copies `<dir>/<name>_f{0,1}.trees.txt`. Keys: `name particle foldBy weights eventSelection preselection trainingData trainingCandidates recipe gitRevision learner validation created` (`weights: <name>_f%d.trees.txt`).
- **`check <trees.root>`** (a file or a wildcard in the file name): candidate trees written by the Maker *with* this model. Every candidate in the domain must have this model's score (to 1e-6), every other candidate a score < −1. Prints `CHECK OK` or stops with exit code 1.

## What the report contains

1. **Summary:** recipe, candidates and AUC per fold, S, S/√(S+B) and threshold at each purity, compared with any other apply outputs given (`label=file.root`, e.g. another model of the same particle).
2. **Curves:** S/√(S+B) and S versus purity, ratio to this BDT.
3. **Spectra** at each purity (linear and log).
4. **Sculpting:** background efficiency ε(m) with a line fitted in the side bands and the pulls of the gap bins next to the window; mass-spectrum ratios at equal purity (between the BDTs of the report; with a single BDT the lower panel shows selected / before the score cut, scaled to 1 in the side bands).
5. **Background models:** S and purity with pol2 / pol4 / fit-free side-band count; window half-width scan.
6. **Purity in bins** of pT, y and cent9.
7. **Overtraining:** training versus test score distributions, AUC.
8. **Importance:** total gain and mean |SHAP| per input and fold; SHAP summary plots; correlation of each input with the mass.
9. **Training size:** S at fixed purity and AUC on the other fold versus the training fraction.

S and B come from a fit at each threshold: a Student-t peak plus pol2 (the K0S and Ξ⁻ peaks have tails; for K0S the deviance/ndf is 1.2 with a Student-t and 4.1 with two Gaussians); B is the integral of pol2 in the window.

## Lessons (Au+Au 3 GeV, Ξ⁻ and K0S)

- **Log-transform heavy-tailed inputs.**
  - With a learner that cuts on a uniform grid between the minimum and the maximum, log10 raised S/√(S+B) at equal purity by a factor 1.2 (Ξ⁻, purity 0.87).
  - xgboost bins by quantiles and is less sensitive, but log10 keeps the inputs readable and the SHAP plots meaningful.
- **Inputs matter more than the training method.**
  - Bigger trees, two-stage training and gradient boosting instead of AdaBoost changed nothing. More inputs did help (track nHitsFit, whether the Λ points back to the PV, ...).
  - xgboost with regularisation (600 trees, depth 3, rate 0.02, min_child_weight 50, subsample 0.7) gave 13% more Ξ⁻ signal than an AdaBoost BDT at equal purity with the same 23 inputs, 17% with every variable of the part; for K0S it was equal (10% sample). Without regularisation it overtrains strongly: the window is mostly background, so the training AUC rises while the test AUC does not.
- **Training statistics (10% sample).** Ξ⁻ (about 25,000 window + side-band candidates per fold) is close to the plateau; K0S reaches it at about 70,000 per fold (30% of 100 files). Run the `study` action before deciding how much to train on.
- **Reweight the side bands, and keep pT out of the inputs.**
  - The kinematics of the background change with the candidate mass. If pT is an input, or the side bands are not reweighted, the BDT learns the window itself.
  - Inputs built from daughter momenta (Λ mass, qT, a Λ-hypothesis mass of a K0S) need the same care. Watch the mass-correlation warnings and the sculpting pages.
- **Never cut only on the window.** Window and side bands must get identical selections; otherwise the BDT learns the cut.
- **Neither of these helped:**
  - Wide side bands made things worse; background far from the window differs from the background under it.
  - Negative-weight side-band subtraction made things worse.
- **Always check for sculpting:**
  - ratios at equal purity (no step at the window edges);
  - ε(m) linear in the side bands;
  - purity stable against the background model.
- **Non-Gaussian peaks.** A peak that one Gaussian cannot describe leaves positive residuals next to the window. Fit a peak with tails (the report uses a Student-t) before calling this sculpting. Example: the K0S gap pulls (+13..+15 at the 10% sample) were peak tails, not sculpting; for this reason the K0S side bands start ≥ 6 σ68 from the peak.
