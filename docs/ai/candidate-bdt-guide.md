# Candidate BDT: design and how to add a BDT for a new particle

This guide is for whoever wants a machine-learned candidate selection for a particle that does not have one yet
(e.g. Λ, φ, Ω). The K0S and Ξ⁻ BDTs of the 3 GeV K0S–Ξ⁻ / K⁺–Ξ⁻ analysis are the worked example. The binding rules are in
[`docs/ai/AGENT_RULES.md`](AGENT_RULES.md) ("Candidate BDT"), the step-by-step procedure in
[`docs/ai/skills/candidate-bdt.md`](skills/candidate-bdt.md), the tool reference in
[`common/macro/datadriven_bdt/README.md`](../../common/macro/datadriven_bdt/README.md) and the Maker side in
[`StMaker/StCandidateFemtoMaker/README.md`](../../StMaker/StCandidateFemtoMaker/README.md).

## 1. Design principles (and why)

| Principle | Why |
|---|---|
| **Labels from data, no Monte Carlo.** Signal-like = candidates in the mass window, background = candidates in the side bands. | The selection is defined by the analysis sample itself; no simulation tuning, the same recipe works for any particle with a mass peak. |
| **A score belongs to a candidate, a threshold to the analysis.** One model per species; the femto Maker stores the score on an axis of the pair THnSparse. | The same candidate has the same score in every channel; thresholds can be scanned and changed without rerunning jobs. |
| **Decisions outside the Maker.** Thresholds, windows and purity corrections are chosen in the CF macros. | Jobs only fill; every choice is reproducible from one production output. |
| **One toolkit, one YAML recipe per model.** `common/macro/datadriven_bdt` (Python + xgboost). | Channel specifics live in the recipe; a new particle needs a new recipe, not new training code. |
| **A model is a file with its provenance.** `<name>.model.yaml` + `<name>_f{0,1}.trees.txt`, written by the toolkit, never edited by hand. | Training data, recipe, git revision and validation travel with the model. |
| **Two folds by run parity.** A candidate is scored by the fold that did not see its run. | No candidate is scored by a model trained on it; events are never split. |
| **The Maker evaluates the trees itself** (C++98, no Python at job time) and must agree with the toolkit to 1e-6 (`check`). | Batch jobs run ROOT 5.34; one score definition on both sides. |
| **Validate before physics:** sculpting, background models, overtraining, SHAP, training size, Maker = offline. | A data-driven BDT can learn the mass window; these checks catch it. |

Environment: the toolkit runs in a Python virtual environment on the host (`pip install -r
common/macro/datadriven_bdt/requirements.txt`, see the toolkit README); batch jobs need no Python.

## 2. What is generic and what you decide

Generic (no code to write): training, applying, reports, the model export, scoring in `StCandidateFemtoMaker`, the score
axes, the threshold options of the CF macros, the purity correction.

You decide, per particle, in YAML:

| Decision | Where | Guidance |
|---|---|---|
| Mass variable, signal window, side bands | recipe `labels` | Window about ±3σ of the peak; side bands start well outside the peak tails (≥ 6σ68 for K0S) and stay close to the window (wide side bands made things worse). Fit the peak with a tailed shape before judging. |
| Preselection | the part's YAML (a loose preselection for the training trees) and recipe `data.preselection` | Loose, identical for window and side bands (never cut only in the window). |
| Scoring domain | recipe `data.eventSelection` | Candidates outside get the outside-training score (e.g. the centrality range of the analysis). |
| Training sample restrictions | recipe `data.trainingSelection` | E.g. drop candidates removed by the shared-track rule in the channels the model serves. Never in the domain keys. |
| Inputs | recipe `inputs` | `use: all` (every variable of the part), exclude what must not be learned (mass, pT, run, score), `log10` for heavy-tailed variables; watch the mass-correlation warnings. |
| Side-band reweighting | recipe `reweight` | Reweight in (pT, pz) so the BDT cannot learn the window from kinematics; keep pT out of the inputs. |
| Learner settings | recipe `learner` | Start from the templates (600 trees, depth 3, rate 0.02, min_child_weight 50, subsample 0.7); without regularisation it overtrains. |
| Input size and training statistics | merged tree file(s), recipe `data.eventFraction`, `evaluation.trainingSize` | Merge the trees (`hadd`); `inspect` shows the candidates and memory before reading; run the `study` action and use the fraction of events where the curve is flat (`eventFraction` < 1 samples whole events, independent of files and folds). |
| Score-axis edges | `candidateBdt` YAML | Fine steps where thresholds will be scanned (e.g. 0.05 above 0), coarse elsewhere; every studied threshold must be an edge. |
| Thresholds | `species_<s>_scoreMin` of a selection YAML | Scan with the CF macro options and choose by SE/ME purities, failed purity fits, statistical error and stability (`common/macro/candidate_femto/README.md`). |

## 3. Adding a BDT for a new particle

**A. The particle is produced by an existing part** (`KfParticleCandidateSource`: K0S, Ξ⁻, Ξ⁺). Only YAML and a recipe:
declare the species and channels in the `candidateFemto` YAML, write candidate trees (`writeCandidateTrees`; merge them with
`hadd`), copy `common/macro/datadriven_bdt/examples/template_new_particle.yaml` and replace every `EDIT`, run
`./script/run_datadriven_bdt_chain.sh <recipe>` (inspect → train → apply → study → report), read the report, then `model` →
`candidateBdt` YAML → `check` → threshold scan (every step is in the skill).

**B. A new KFParticle particle** (e.g. Λ, Ω): add it to `KfParticleCandidateSource` (`SpeciesKeys`, `DaughterNames`,
`VariableNames`, `ProcessEvent`; all KFParticle particles of a job come from one Finder pass) and to the kf YAML. The variables
you expose are the BDT inputs and tree branches: expose every topological and daughter quantity that could help. Then case A.

**C. A particle from another reconstruction** (e.g. φ → K⁺K⁻ from track pairs): write a new part implementing
`CandidateSource` (`include/CandidateSource.h`): `massValue` = invariant mass, `daughters` = the PicoDst track indices,
`variables` = the quantities the BDT may use (daughter PID, pair kinematics, …). The Maker, the trees, the store and the BDT
need no change. Then case A. If the side bands do not describe the background under the peak (e.g. a broad peak on a steep
combinatorial background), another way to define the labels is a new `labels.method` in the toolkit, not a new script.

**D. An analysis in another Maker** (e.g. `StFemtoMaker`): move `CandidateBdtScorer.h` and `CandidateTreeEnsemble.h` to a
shared place and call them from that Maker; do not write a second scorer or evaluator.

## 4. Before a model is used for physics

1. `datadriven_bdt check` on trees written with the model: `CHECK OK`.
2. The report: sculpting (ratio at equal purity, ε(m), pulls next to the window), background-model variation, overtraining,
   SHAP (no input that dominates for a non-physical reason), training size.
3. Record the report in the model's `validation`, add the model to the analysis's catalog (`bdt_models.md` of the version).
4. Purities for the CF from the pair THnSparse, SE and ME separately.

## 5. Extending the toolkit

New ways to label signal and background go into `Labels` (`datadriven_bdt.py`, `labels.method`); another learner type needs
an exporter in the toolkit, a reference evaluator and a C++98 evaluator behind `CandidateBdtScorer`. Keep the score definition
and the Maker = offline check for every extension.
