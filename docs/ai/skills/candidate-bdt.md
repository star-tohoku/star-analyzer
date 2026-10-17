# Candidate BDT: train, validate, register, apply

Use when training a new BDT for a femto species, changing a BDT, or wiring a BDT into
`StCandidateFemtoMaker`. Rules (the training/application contract): `docs/ai/AGENT_RULES.md`,
section "Candidate BDT". Tool details: `common/macro/datadriven_bdt/README.md`. Maker side:
`StMaker/StCandidateFemtoMaker/README.md`, `CandidateBdtScorer.h` and `CandidateTreeEnsemble.h`. Design and the decisions to
make for a new particle (new KFParticle particle, a new reconstruction part, another Maker): `docs/ai/candidate-bdt-guide.md`.

The toolkit is Python + xgboost; the Maker evaluates the exported tree-ensemble files itself (`CandidateTreeEnsemble.h`,
no Python or xgboost at job time).

Before you start: a BDT is for a species with a mass peak (decays); identified tracks (K⁺) use PID cuts. Training uses
data only (window vs side bands), every variable of the part as input unless a check says otherwise, and a
training-size study before any full-statistics training.

## 1. Training trees

1. Mainconf for the training run: the analysis mainconf with a candidateFemto YAML that has
   `writeCandidateTrees: true` and `candidateTreeSpecies` (the species that get a tree), and **no** `candidateBdt` key (example:
   `config/mainconf/main_auau3_anaCandidateFemto_k0sXiKaonXi_bdtTraining.yaml`).
2. Local test first (`script/singularity_run_anaCandidateFemto.sh`), then the farm. Generate the joblist with
   `script/generate_joblist_from_list.py ... --tree-output-dir <dir>` (it also copies `<stem>_$JOBID.tree.root`).
   If candidate stores exist (reconstruction stage, `writeCandidateStore: true`), run this as a pair stage on the stores:
   same mainconf, store files instead of PicoDsts, no PicoDst read (much faster), and one tree file per pair-stage job
   (50–100 stores) instead of one per PicoDst job: the toolkit's run time is dominated by the number of files it opens.
3. Trees: `candidates_<species>` in `<output>.tree.root` (one entry per usable candidate of an accepted
   event; branches `runId eventId cent9 mass px py pz score`, the part's variables and one flag per channel of the species,
   `sharedTrackRemoved_<channel>` = 1 if the shared-track rule removed the candidate from that channel).

## 2. Recipe

Copy `common/macro/datadriven_bdt/examples/template_new_particle.yaml` (replace every `EDIT`) or a worked recipe
(`examples/{xi_minus,k0s}_candidate_trees.yaml`); every key is required, the table is in the toolkit README. Set at least:

```yaml
name: <species>_<what>                       # readable; becomes the model and file names
particle: <decay and reconstruction, any condition the trees cannot express>
validation: "none: not validated yet"        # set after step 4
data: {trees: rootfile/<anaName>/.../*.tree.root, tree: candidates_<species>, eventSelection: ..., preselection: "1",
       trainingSelection: "sharedTrackRemoved_<channel>==0 [&& ... for every channel the model serves]", folds: runParity,
       eventFraction: 1}                    # < 1: use a deterministic fraction of the events (large merged inputs)
labels: {method: massSidebands, mass: mass, window: [lo, hi], sidebands: [[lo, hi], [lo, hi]]}
inputs: {use: all, exclude: [massError], log10: [...heavy-tailed variables...], massCorrelationWarning: 0.1}
reweight: {variables: ["sqrt(px*px+py*py)", pz], bins: [[30, 0, 3], [30, -3, 0.5]]}
learner: {type: xgboost, seed: ..., testFraction: 0.2, parameters: {...}}
evaluation: {...}                            # checks: sculpting overtraining importance backgroundModels trainingSize
output: share/figure/<data>_<anaName>/<major>/<version>/studies/<topic>/data/<name>
```

- `trainingSelection` acts on the training samples only; `eventSelection`/`preselection` define the scoring domain and must
  not use the flags. `inputs.use: all` never includes the flags.
- Expressions use only tree branches and the shared grammar (rules). Log10 heavy-tailed inputs, keep pT out of the
  inputs and reweight the side bands (lessons in the toolkit README).

## 3. Train

Once: a Python virtual environment with `common/macro/datadriven_bdt/requirements.txt` (toolkit README, "Python environment");
activate it before running the toolkit.

```bash
R=<recipe.yaml>
./script/run_datadriven_bdt.sh $R inspect    # size of the input and memory needed; lower data.eventFraction if too large
./script/run_datadriven_bdt.sh $R train      # both folds; prints the mass-correlation warnings and the AUC per fold
```

## 4. Validate

```bash
./script/run_datadriven_bdt.sh $R apply
./script/run_datadriven_bdt.sh $R study      # training-size study (learning curve), if evaluation.checks has trainingSize
./script/run_datadriven_bdt.sh $R report [label=<other>_apply.root ...]
```

Read the report: sculpting (ratio at equal purity, ε(m), background-model variation), overtraining, SHAP (no input that
dominates for a non-physical reason; act on mass-correlation warnings), training size (is the curve flat at the
statistics used?). Then set `validation: <report path>` (`particle`, `validation` and `evaluation` are the only keys that
may change between training and `model`).

## 5. Register the model

```bash
./script/run_datadriven_bdt.sh $R model <dir>
```

- Writes `<dir>/<name>.model.yaml` and `<dir>/<name>_f{0,1}.trees.txt`. Refuses folds trained with a different recipe.
- Production: `<dir>` = `config/bdt/<anaName>/`; study: the study's `models/`.
- A model used on the farm must be under `config/` (only `config/` and the code directories are shipped). A model elsewhere
  (e.g. in a study directory, referenced as `../share/...` relative to `config/`) works locally and fails on the farm.
- Add the model to the analysis catalog (`bdt_models.md`).

## 6. Apply in the Maker

candidateBdt YAML (referenced by the mainconf key `candidateBdt`):

```
speciesKeys: xi_minus
outsideTrainingScore: -2
xi_minus_model: bdt/<anaName>/<name>.model.yaml
xi_minus_scoreEdges: -3, -1, 0, 0.05, 0.1, 0.15, 0.2, 1
```

- The first score bin holds only the outside-training score; the other edges cover [-1, 1].
- Every threshold to be studied must be an edge: use a fine grid where thresholds will be scanned (e.g. 0.05 steps
  above 0) and coarse bins elsewhere, since each edge multiplies the THnSparse size.

## 7. Check the Maker scores

Run the Maker with the model **and** `writeCandidateTrees: true` on a few files, then:

```bash
./script/run_datadriven_bdt.sh $R check '<out dir>/*.tree.root'
```

It must print `CHECK OK`: every in-domain Maker score equals the offline score (≤ 1e-6), and the domain
decisions agree. Only then submit production.

## 8. Choose the thresholds

Thresholds are chosen in the CF macros from the score axes of one `_bdt` production output, never by rerunning a job.
Procedure and criteria: `common/macro/candidate_femto/README.md`, "Choosing BDT thresholds" (scan with the macro options
`scoreMinA` / `scoreMinB`; look at SE and ME purities, failed purity fits, statistical error and stability; keep the scan
as a study). The chosen values go into `species_<s>_scoreMin` of a `_bdtSelection_<thresholds>` candidateFemto YAML.
Purities for the CF come from the pair THnSparse (SE and ME separately), not from the per-species histograms.
