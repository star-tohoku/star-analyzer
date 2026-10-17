#!/usr/bin/env python3
# datadriven_bdt.py -- data-driven (no Monte Carlo) BDT for candidates that form a mass peak (K0S, Lambda, Xi, Omega, ...):
# candidates in the mass window (signal + background) are trained against candidates in the side bands (background),
# two folds by run parity, xgboost. Everything channel-specific is in one YAML recipe (README.md, examples/).
#
#   ./script/run_datadriven_bdt.sh <recipe.yaml> train                     both folds
#   ./script/run_datadriven_bdt.sh <recipe.yaml> apply                     score every candidate with the fold that did not see it
#   ./script/run_datadriven_bdt.sh <recipe.yaml> study                     training-size study (learning curve)
#   ./script/run_datadriven_bdt.sh <recipe.yaml> report [label=apply.root ...]   one PDF with every evaluation item
#   ./script/run_datadriven_bdt.sh <recipe.yaml> model <dir>               model file for StCandidateFemtoMaker
#   ./script/run_datadriven_bdt.sh <recipe.yaml> check <trees.root>        Maker scores == scores of this model
#   ./script/run_datadriven_bdt.sh <recipe.yaml> inspect                   size of the input and memory needed, reads nothing
#
# Run from the repository root (the wrapper does that); relative paths in the recipe are relative to it.
# Rules: docs/ai/AGENT_RULES.md "Candidate BDT"; procedure: docs/ai/skills/candidate-bdt.md.
import copy, datetime, glob, hashlib, json, os, subprocess, sys
import numpy as np
import uproot
import yaml

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from candidate_expression import CandidateExpression, ExpressionError  # noqa: E402
import tree_ensemble  # noqa: E402

# Branches every candidate tree of StCandidateFemtoMaker starts with (CandidateModelVariableNames + score).
EVENT_BRANCHES = ['runId', 'eventId', 'cent9', 'mass', 'px', 'py', 'pz']
MAKER_SCORE = 'score'
# Inputs that would let the BDT recognise the window or the run directly.
FORBIDDEN_INPUT_VARIABLES = ['mass', 'runId', 'eventId', MAKER_SCORE]
# Per-channel flag branches of the candidate trees (StCandidateFemtoMaker): 1 = removed by the shared-track rule in that
# channel. They are bookkeeping, never inputs; a recipe uses them only in data.trainingSelection.
FLAG_PREFIX = 'sharedTrackRemoved_'


def fatal(msg):
    print(f'DATADRIVEN_BDT FATAL: {msg}', flush=True)
    sys.exit(1)


# ----------------------------------------------------------------------------------------------- recipe
# Every key is required (no built-in defaults); unknown keys are errors. Lists of allowed keys per section:
SCHEMA = {
    None: ['name', 'particle', 'validation', 'data', 'labels', 'inputs', 'reweight', 'learner', 'evaluation', 'output'],
    'data': ['trees', 'tree', 'eventSelection', 'preselection', 'trainingSelection', 'folds', 'eventFraction'],
    'labels': ['method', 'mass', 'window', 'sidebands'],
    'inputs': ['use', 'exclude', 'log10', 'massCorrelationWarning'],
    'reweight': ['variables', 'bins'],
    'learner': ['type', 'seed', 'testFraction', 'parameters'],
    'evaluation': ['massRange', 'binWidth', 'purities', 'sculptingPurity', 'windowScan', 'checks', 'trainingSize', 'purityBins',
                   'shapCandidates'],
}
LABEL_METHODS = ['massSidebands']   # signal/background definitions; add new methods in Labels below
CHECKS = ['sculpting', 'overtraining', 'importance', 'backgroundModels', 'trainingSize']
# keys that may change between training and `model` (they do not change the trained model)
FREE_KEYS = ['particle', 'validation', 'evaluation']


def load_recipe(path):
    text = open(path).read()
    left = [f'{i}: {l.strip()}' for i, l in enumerate(text.splitlines(), 1) if 'EDIT' in l.split('#')[0]]
    if left:
        fatal(f'{path} still contains EDIT placeholders (examples/template_new_particle.yaml): ' + '; '.join(left))
    try:
        r = yaml.safe_load(text)
    except yaml.YAMLError as e:
        fatal(f'cannot parse {path}: {e}')
    if not isinstance(r, dict):
        fatal(f'{path}: not a YAML mapping')
    for sec, keys in SCHEMA.items():
        d = r if sec is None else r.get(sec)
        if not isinstance(d, dict):
            fatal(f'{path}: section "{sec}" missing or not a mapping')
        for k in keys:
            if k not in d:
                fatal(f'{path}: key "{k if sec is None else sec + "." + k}" missing (every key is required)')
        for k in d:
            if k not in keys:
                fatal(f'{path}: unknown key "{k if sec is None else sec + "." + k}"')
    f = r['data']['eventFraction']
    if isinstance(f, bool) or not isinstance(f, (int, float)) or not 0. < float(f) <= 1.:
        fatal(f'data.eventFraction must be a number in (0, 1], not {f!r}')
    if r['data']['folds'] != 'runParity':
        fatal('data.folds must be runParity (the only scheme StCandidateFemtoMaker applies)')
    if r['labels']['method'] not in LABEL_METHODS:
        fatal(f"labels.method must be one of {LABEL_METHODS}")
    if len(r['labels']['window']) != 2 or not all(len(s) == 2 for s in r['labels']['sidebands']):
        fatal('labels.window is [lo, hi], labels.sidebands a list of [lo, hi]')
    if r['learner']['type'] != 'xgboost':
        fatal('learner.type must be xgboost')
    if not 0. < float(r['learner']['testFraction']) < 1.:
        fatal('learner.testFraction must be in (0, 1)')
    p = r['learner']['parameters']
    if not isinstance(p, dict) or 'objective' in p and p['objective'] != 'binary:logistic':
        fatal('learner.parameters must be a mapping; the objective is binary:logistic')
    for k in ('random_state', 'seed'):
        if k in p:
            fatal(f'learner.parameters.{k}: use learner.seed')
    rw = r['reweight']
    if rw['variables'] != 'none' and (len(rw['variables']) != 2 or len(rw['bins']) != 2 or any(len(b) != 3 for b in rw['bins'])):
        fatal('reweight: two variables and bins [[n, lo, hi], [n, lo, hi]], or variables: none and bins: none')
    ev = r['evaluation']
    for c in ev['checks']:
        if c not in CHECKS:
            fatal(f'evaluation.checks: unknown check "{c}" (known: {CHECKS})')
    if 'trainingSize' in ev['checks'] and (not isinstance(ev['trainingSize'], dict) or set(ev['trainingSize']) != {'fractions', 'purity'}):
        fatal('evaluation.trainingSize: {fractions: [...], purity: p}')
    if set(ev['purityBins']) - {'pt', 'y', 'cent9'}:
        fatal('evaluation.purityBins: keys pt, y, cent9 (any subset)')
    if '%' in r['name'] or '/' in r['name'] or ' ' in r['name']:
        fatal('name: no %, / or spaces')
    root = os.getcwd()
    r['_path'] = os.path.abspath(path)
    r['_path_in_repo'] = os.path.relpath(r['_path'], root) if r['_path'].startswith(root + '/') else r['_path']
    r['_text'] = text
    out = r['output']
    print(f"DATADRIVEN_BDT recipe {r['_path']}\n  name {r['name']}, trees {r['data']['trees']} ({r['data']['tree']}), output {out}", flush=True)
    return r


def trained_part(r):
    """The recipe without the keys that may change after training (compared by `model`)."""
    d = {k: copy.deepcopy(v) for k, v in r.items() if not k.startswith('_') and k not in FREE_KEYS}
    return d


def git_revision():
    try:
        rev = subprocess.run(['git', 'rev-parse', '--short', 'HEAD'], capture_output=True, text=True).stdout.strip() or 'unrecorded'
        n = len([l for l in subprocess.run(['git', 'status', '--porcelain', '--untracked-files=no'], capture_output=True,
                                           text=True).stdout.splitlines() if l.strip()])
        return rev + (f' (+{n} modified tracked files)' if n else '')
    except OSError:
        return 'unrecorded'


# ----------------------------------------------------------------------------------------------- data
def input_files(pattern):
    d, b = os.path.split(pattern)
    if any(c in d for c in '*?['):
        fatal(f'wildcards only in the file name: {pattern}')
    files = sorted(glob.glob(pattern))
    if not files:
        fatal(f'no file matches {pattern}')
    return files


def event_fraction_mask(run_id, event_id, fraction):
    """Deterministic event sample: keeps an event if a hash of (runId, eventId) mapped to [0, 1) is below fraction.
    Independent of how the candidates are split into files, of the order of the entries and of the run-parity folds;
    all candidates of an event are kept or dropped together."""
    if fraction >= 1.:
        return np.ones(len(run_id), dtype=bool)
    x = (run_id.astype(np.uint64) << np.uint64(32)) ^ event_id.astype(np.uint64)
    with np.errstate(over='ignore'):                           # splitmix64 finaliser, arithmetic modulo 2^64
        x = x + np.uint64(0x9E3779B97F4A7C15)
        x = (x ^ (x >> np.uint64(30))) * np.uint64(0xBF58476D1CE4E5B9)
        x = (x ^ (x >> np.uint64(27))) * np.uint64(0x94D049BB133111EB)
        x = x ^ (x >> np.uint64(31))
    return (x >> np.uint64(11)).astype(np.float64) / float(1 << 53) < fraction


# Peak memory of train/apply/study: a fixed part plus a multiple of the candidates read. Measured with the 3 GeV trees
# (peak / read): K0S 2.0 GB / 0.28 GB and 3.7 GB / 0.83 GB, Xi- 0.6 GB / 0.08 GB; these constants stay above all three.
# Used to stop before the machine runs out of memory.
MEMORY_PER_READ_BYTE = 4.
MEMORY_FIXED_BYTES = 1.5e9


def available_memory_bytes():
    """MemAvailable of /proc/meminfo, or None where it does not exist."""
    try:
        for line in open('/proc/meminfo'):
            if line.startswith('MemAvailable:'):
                return int(line.split()[1]) * 1024
    except OSError:
        pass
    return None


def cache_key(files, tree, want, fraction):
    """Identifies what Candidates reads: the files (path, size, modification time), the tree, the branches and the
    event fraction."""
    h = hashlib.sha256()
    h.update(json.dumps([tree, want, float(fraction)]).encode())
    for f in files:
        s = os.stat(f)
        h.update(f'{os.path.abspath(f)}|{s.st_size}|{s.st_mtime_ns}\n'.encode())
    return h.hexdigest()


class Candidates:
    """Candidates of one tree (all files), the variables a recipe needs, as float64 / int arrays.
    cache: optional .npz path. Reading many tree files dominates the run time, so the arrays are saved there and reused by
    later actions as long as the files (path, size, modification time), the tree and the branches are unchanged; otherwise
    the files are read again and the cache is rewritten. The cache can be deleted at any time."""

    def __init__(self, files, tree, names_needed, with_maker_score=False, cache=None, fraction=1.):
        from concurrent.futures import ThreadPoolExecutor
        self.files = files
        branches = tree_branch_names(files, tree)
        if branches[:len(EVENT_BRANCHES)] != EVENT_BRANCHES:
            fatal(f'{files[0]}: {tree} does not start with {EVENT_BRANCHES} (not a candidate tree of StCandidateFemtoMaker)')
        missing = [n for n in names_needed if n not in branches]
        if missing:
            fatal(f'{tree} has no branch {missing}')
        want = [n for n in branches if n in names_needed or n in EVENT_BRANCHES or (with_maker_score and n == MAKER_SCORE)]
        self.branches = branches
        self.model_variables = [b for b in branches if b != MAKER_SCORE]   # CandidateModelVariableNames
        self.source_variables = [b for b in self.model_variables if b not in EVENT_BRANCHES]
        key = cache_key(files, tree, want, fraction) if cache else None
        if cache and os.path.exists(cache):
            try:
                with np.load(cache, allow_pickle=False) as z:
                    if str(z['__key__']) == key:
                        self.v = {k: z[k] for k in want}
                        self.n = len(self.v['mass'])
                        print(f'DATADRIVEN_BDT input: {len(files)} files, {self.n} candidates in {tree}'
                              + (f' (event fraction {fraction})' if fraction < 1. else '') + f' (from the cache {cache})', flush=True)
                        return
            except Exception as e:                            # damaged or foreign file: read the trees instead
                print(f'WARNING: cache {cache} not used ({type(e).__name__}: {e}); reading the tree files', flush=True)

        import threading
        lock = threading.Lock()
        budget = {'read': 0., 'available': available_memory_bytes(), 'error': None}

        class ReadError(Exception):
            pass

        def read(f):
            """One file; the first error stops the other reading threads and is reported once."""
            if budget['error']:
                return None
            try:
                return read_file(f)
            except ReadError as e:
                msg = str(e)
            except Exception as e:                            # a damaged or incomplete file: say which one
                msg = f'cannot read {tree} from {f}: {type(e).__name__}: {str(e).splitlines()[0]}'
            with lock:
                budget['error'] = budget['error'] or msg
            return None

        def reserve(t, f):
            """Stop before reading a file if the expected peak memory exceeds the memory available at the start."""
            per = sum(t[b].interpretation.numpy_dtype.itemsize for b in want)
            with lock:
                budget['read'] += t.num_entries * per * fraction
                need = MEMORY_FIXED_BYTES + MEMORY_PER_READ_BYTE * budget['read']
                if budget['available'] is not None and need > budget['available']:
                    raise ReadError(f'not enough memory: the candidates of {tree} up to {f} need about {need / 1e9:.1f} GB (to read: '
                          f'{budget["read"] / 1e9:.2f} GB, x {MEMORY_PER_READ_BYTE:g} for the training, + {MEMORY_FIXED_BYTES / 1e9:g} GB), '
                          f'{budget["available"] / 1e9:.1f} GB available; lower data.eventFraction (the inspect action shows '
                          'the size) or use a machine with more free memory')

        def read_file(f):
            with uproot.open(f) as fh:
                if tree not in fh:
                    raise ReadError(f'no tree {tree} in {f}')
                t = fh[tree]
                if [b.name for b in t.branches] != branches:
                    raise ReadError(f'{f}: branches of {tree} differ from the first file')
                reserve(t, f)
                if fraction >= 1.:
                    return t.arrays(want, library='np')
                # read in chunks and keep only the sampled events, so a large merged file needs memory for the sample only
                chunks = []
                for a in t.iterate(want, library='np', step_size='200 MB'):
                    m = event_fraction_mask(a['runId'], a['eventId'], fraction)
                    chunks.append({k: v[m] for k, v in a.items()})
                return {k: np.concatenate([c[k] for c in chunks]) for k in want}

        parts = []
        with ThreadPoolExecutor(max_workers=8) as pool:      # reading is I/O bound; order of the files is kept
            for i, a in enumerate(pool.map(read, files)):
                parts.append(a)
                if (i + 1) % 500 == 0 and not budget['error']:
                    print(f'  read {i + 1}/{len(files)} files', flush=True)
        if budget['error']:
            fatal(budget['error'])
        self.v = {k: np.concatenate([a[k] for a in parts]) for k in want}
        self.n = len(self.v['mass'])
        mb = sum(v.nbytes for v in self.v.values()) / 1e6
        print(f'DATADRIVEN_BDT input: {len(files)} files, {self.n} candidates in {tree}'
              + (f' (event fraction {fraction})' if fraction < 1. else '') + f', {mb:.0f} MB in memory', flush=True)
        if cache:
            try:
                os.makedirs(os.path.dirname(cache) or '.', exist_ok=True)
                tmp = cache + '.tmp.npz'
                np.savez(tmp, __key__=np.array(key), **self.v)
                os.replace(tmp, cache)
                print(f'DATADRIVEN_BDT wrote the cache {cache}', flush=True)
            except OSError as e:
                print(f'WARNING: cache {cache} not written ({e})', flush=True)

    def subset(self, mask):
        c = copy.copy(self)
        c.v = {k: v[mask] for k, v in self.v.items()}
        c.n = int(mask.sum())
        return c


def open_tree_file(f):
    """uproot.open with a one-line error naming the file (damaged or incomplete files)."""
    try:
        return uproot.open(f)
    except Exception as e:
        fatal(f'cannot open {f}: {type(e).__name__}: {str(e).splitlines()[0]}')


def tree_branch_names(files, tree):
    with open_tree_file(files[0]) as fh:
        if tree not in fh:
            fatal(f'no tree {tree} in {files[0]}')
        return [b.name for b in fh[tree].branches]


def input_expressions(r, branches):
    """The input expressions of the recipe, in order."""
    model_vars = [b for b in branches if b != MAKER_SCORE]
    source = [b for b in model_vars if b not in EVENT_BRANCHES and not b.startswith(FLAG_PREFIX)]
    spec = r['inputs']
    if spec['use'] == 'all':
        base = list(source)       # every source variable of the part (not the event/kinematic branches)
    elif isinstance(spec['use'], list):
        base = list(spec['use'])
    else:
        fatal('inputs.use: "all" or a list of expressions')
    for e in spec['exclude']:
        if e not in base:
            fatal(f'inputs.exclude: "{e}" is not an input')
    base = [e for e in base if e not in spec['exclude']]
    for e in spec['log10']:
        if e not in base:
            fatal(f'inputs.log10: "{e}" is not an input')
    exprs = [f'log10(max({e},1e-6))' if e in spec['log10'] else e for e in base]
    if len(set(exprs)) != len(exprs):
        fatal('duplicate inputs')
    for e in exprs:
        try:
            ce = CandidateExpression(e, model_vars)
        except ExpressionError as err:
            fatal(f'input {err}')
        bad = [v for v in ce.variables if v in FORBIDDEN_INPUT_VARIABLES or v.startswith(FLAG_PREFIX)]
        if bad:
            fatal(f'input "{e}" uses {bad}: the BDT would see the mass, the run or a bookkeeping flag directly')
        if set(ce.variables) & {'px', 'py', 'pz'}:
            print(f'DATADRIVEN_BDT WARNING: input "{e}" uses the candidate momentum; keep pT out of the inputs (README, lessons)')
    return exprs


def needed_variables(r, exprs, branches):
    model_vars = [b for b in branches if b != MAKER_SCORE]
    texts = list(exprs) + [r['data']['eventSelection'], r['data']['preselection'], r['labels']['mass']]
    if str(r['data']['trainingSelection']).strip():
        texts.append(r['data']['trainingSelection'])
    if r['reweight']['variables'] != 'none':
        texts += list(r['reweight']['variables'])
    need = set()
    for t in texts:
        try:
            need |= set(CandidateExpression(str(t), model_vars).variables)
        except ExpressionError as err:
            fatal(str(err))
    return sorted(need)


def evaluate(c, text):
    return CandidateExpression(str(text), c.model_variables).evaluate(c.v)


def domain(r, c):
    """Candidates inside the training domain (eventSelection && preselection), as in the Maker."""
    return (evaluate(c, r['data']['eventSelection']) != 0.) & (evaluate(c, r['data']['preselection']) != 0.)


def training_mask(r, c):
    """data.trainingSelection: applied to the training samples only (train, study), never to the domain the Maker scores.
    Empty = every domain candidate may be used for training."""
    t = str(r['data']['trainingSelection']).strip()
    return np.ones(c.n, bool) if not t else evaluate(c, t) != 0.


def input_matrix(c, exprs):
    """float32 (n, n_inputs): the value of each expression rounded to float32 (the Maker's mInputs)."""
    return np.column_stack([evaluate(c, e).astype(np.float32) for e in exprs]) if exprs else np.zeros((c.n, 0), np.float32)


def run_parity(c):
    return (c.v['runId'].astype(np.int64) % 2).astype(np.int64)


class Labels:
    """Signal-like (1) / background (0) / unused (-1) per candidate. Only massSidebands is implemented."""

    def __init__(self, r):
        self.method = r['labels']['method']
        self.mass_expr = r['labels']['mass']
        self.window = [float(x) for x in r['labels']['window']]
        self.sidebands = [[float(a), float(b)] for a, b in r['labels']['sidebands']]

    def mass(self, c):
        return evaluate(c, self.mass_expr)

    def in_window(self, m):
        return (m > self.window[0]) & (m < self.window[1])

    def in_sidebands(self, m):
        s = np.zeros(len(m), bool)
        for lo, hi in self.sidebands:
            s |= (m > lo) & (m < hi)
        return s

    def label(self, c):
        m = self.mass(c)
        lab = np.full(c.n, -1, np.int64)
        lab[self.in_sidebands(m)] = 0
        lab[self.in_window(m)] = 1
        return lab


def bin_index(x, n, lo, hi):
    """TH1-like bin index with under/overflow: 0 below lo, n+1 at or above hi."""
    w = (hi - lo) / n
    i = np.floor((x - lo) / w).astype(np.int64) + 1
    return np.clip(i, 0, n + 1)


def sideband_weights(r, c, lab):
    """Side-band candidates reweighted to the window distribution in the two reweight variables (TH2 with overflow),
    separately in each fold (each fold's BDT sees only its own runs)."""
    w = np.ones(c.n)
    if r['reweight']['variables'] == 'none':
        print('DATADRIVEN_BDT note: no side-band reweighting (reweight.variables: none)')
        return w
    (nx, xl, xh), (ny, yl, yh) = r['reweight']['bins']
    a = evaluate(c, r['reweight']['variables'][0]); b = evaluate(c, r['reweight']['variables'][1])
    k = bin_index(a, int(nx), float(xl), float(xh)) * (int(ny) + 2) + bin_index(b, int(ny), float(yl), float(yh))
    size = (int(nx) + 2) * (int(ny) + 2)
    par = run_parity(c)
    for fold in (0, 1):
        f = par == fold
        hw = np.bincount(k[f & (lab == 1)], minlength=size).astype(float); hb = np.bincount(k[f & (lab == 0)], minlength=size).astype(float)
        hw /= max(hw.sum(), 1.); hb /= max(hb.sum(), 1.)
        ratio = np.where(hb > 0, hw / np.where(hb > 0, hb, 1.), 0.)
        sb = f & (lab == 0)
        w[sb] = ratio[k[sb]]
    return w


def mass_correlations(r, c, lab, x, exprs):
    """Pearson correlation of every input with the mass, in the side bands (before training)."""
    m = Labels(r).mass(c)
    sb = lab == 0
    rows = []
    for i, e in enumerate(exprs):
        v = x[sb, i].astype(np.float64)
        rr = float(np.corrcoef(v, m[sb])[0, 1]) if sb.sum() > 2 and v.std() > 0 else 0.
        rows.append((e, rr))
    lim = float(r['inputs']['massCorrelationWarning'])
    for e, rr in sorted(rows, key=lambda t: -abs(t[1])):
        if abs(rr) > lim:
            print(f'DATADRIVEN_BDT WARNING: input "{e}" is correlated with the mass in the side bands (r = {rr:+.3f}, '
                  f'limit {lim}): it can sculpt the background; consider excluding it')
    return rows


# ----------------------------------------------------------------------------------------------- learner
def make_classifier(r, seed=None):
    import xgboost as xgb
    p = dict(r['learner']['parameters'])
    p['objective'] = 'binary:logistic'
    return xgb.XGBClassifier(random_state=int(r['learner']['seed'] if seed is None else seed), **p)


def train_fold(r, x, y, w, seed_offset=0):
    """Train one fold; returns (classifier, train mask over the given rows)."""
    rng = np.random.default_rng(int(r['learner']['seed']) + seed_offset)
    test = rng.random(len(y)) < float(r['learner']['testFraction'])
    ww = w.copy()
    # equal total weight for window and side bands, so the learner sees balanced classes
    for s in (0, 1):
        sel = ~test & (y == s)
        tot = ww[sel].sum()
        if tot <= 0:
            fatal(f'no training candidates in the {"window" if s == 1 else "side bands"} (training part of the fold: window '
                  f'{int((~test & (y == 1)).sum())}, side bands {int((~test & (y == 0)).sum())}); check labels.window / '
                  'labels.sidebands, the selections and data.eventFraction')
    norm = ww[~test & (y == 1)].sum() / ww[~test & (y == 0)].sum()
    ww[y == 0] *= norm
    clf = make_classifier(r, int(r['learner']['seed']) + seed_offset)
    clf.fit(x[~test], y[~test], sample_weight=ww[~test])
    return clf, ~test, ww


def weighted_auc(y, s, w):
    o = np.argsort(s)
    y, w = y[o], w[o]
    ws = w * (y == 1); wb = w * (y == 0)
    cb = np.cumsum(wb)
    return float((ws * (cb - 0.5 * wb)).sum() / max(ws.sum() * wb.sum(), 1e-300))


def paths(r):
    out, n = r['output'], r['name']
    return dict(model=lambda f: f'{out}/{n}_f{f}.xgboost.json', ensemble=lambda f: f'{out}/{n}_f{f}.trees.txt',
                training=lambda f: f'{out}/{n}_f{f}.training.json', shap=lambda f: f'{out}/{n}_f{f}.shap.npz',
                apply=f'{out}/{n}_apply.root',
                study=f'{out}/{n}_training_size.json', report=f'{out}/{n}_report.pdf',
                cache=f'{out}/{n}_candidates_cache.npz')


def prepare(r):
    """Read the candidates of the recipe; returns (candidates, input expressions, labels, domain mask)."""
    files = input_files(r['data']['trees'])
    branches = tree_branch_names(files, r['data']['tree'])
    exprs = input_expressions(r, branches)
    c = Candidates(files, r['data']['tree'], needed_variables(r, exprs, branches), cache=paths(r)['cache'],
                   fraction=float(r['data']['eventFraction']))
    return c, exprs


def action_train(r):
    import xgboost as xgb
    c, exprs = prepare(r)
    dom = domain(r, c)
    c = c.subset(dom)
    tsel = training_mask(r, c)
    print(f'DATADRIVEN_BDT {r["name"]}: training selection keeps {int(tsel.sum())} of {c.n} domain candidates', flush=True)
    c = c.subset(tsel)
    lab = Labels(r).label(c)
    use = lab >= 0
    c = c.subset(use); lab = lab[use]
    x = input_matrix(c, exprs)
    w = sideband_weights(r, c, lab)
    corr = mass_correlations(r, c, lab, x, exprs)
    par = run_parity(c)
    P = paths(r)
    os.makedirs(r['output'], exist_ok=True)
    rev = git_revision()
    for fold in (0, 1):
        sel = (par == fold) & (w > 0)
        xf, yf, wf = x[sel], lab[sel], w[sel]
        print(f"DATADRIVEN_BDT {r['name']} fold {fold}: window {int((yf == 1).sum())}, side bands {int((yf == 0).sum())}, "
              f'inputs {len(exprs)}', flush=True)
        clf, train, ww = train_fold(r, xf, yf, wf)
        booster = clf.get_booster()
        booster.save_model(P['model'](fold))
        ens = tree_ensemble.from_xgboost(booster, exprs, f"xgboost {xgb.__version__} binary:logistic")
        tree_ensemble.write(ens, P['ensemble'](fold))
        ens = tree_ensemble.read(P['ensemble'](fold))
        # the exported file reproduces xgboost (float32 margins inside xgboost)
        p_xgb = clf.predict_proba(xf)[:, 1].astype(np.float64)
        s_ens = tree_ensemble.score(ens, xf)
        diff = float(np.max(np.abs((2 * p_xgb - 1) - s_ens))) if len(xf) else 0.
        if diff > 1e-5:
            fatal(f'exported trees differ from xgboost by {diff:.2e} (fold {fold})')
        s = s_ens
        info = dict(name=r['name'], fold=fold, recipe=trained_part(r), recipePath=r['_path_in_repo'], recipeText=r['_text'],
                    gitRevision=rev, learner=f'xgboost {xgb.__version__}', created=datetime.datetime.now().isoformat(timespec='seconds'),
                    inputs=exprs, candidates=dict(window=int((yf == 1).sum()), sidebands=int((yf == 0).sum())),
                    exportedVsXgboostMaxDifference=diff,
                    auc=dict(train=weighted_auc(yf[train], s[train], ww[train]), test=weighted_auc(yf[~train], s[~train], ww[~train])),
                    massCorrelation=[[e, v] for e, v in corr],
                    scoreHistograms={k: np.histogram(s[m], bins=50, range=(-1, 1), weights=ww[m])[0].tolist() for k, m in (
                        ('trainWindow', train & (yf == 1)), ('trainSidebands', train & (yf == 0)),
                        ('testWindow', ~train & (yf == 1)), ('testSidebands', ~train & (yf == 0)))},
                    gain={exprs[int(k[1:])]: float(v) for k, v in booster.get_score(importance_type='total_gain').items()})
        if 'importance' in r['evaluation']['checks']:
            # SHAP values (log-odds units) of a random sample of this fold's test candidates, for the report
            import shap
            test_rows = np.flatnonzero(~train)
            rng = np.random.default_rng(int(r['learner']['seed']) + 77 + fold)
            pick = rng.choice(test_rows, size=min(int(r['evaluation']['shapCandidates']), len(test_rows)), replace=False)
            sv = shap.TreeExplainer(booster).shap_values(xf[pick])
            np.savez_compressed(P['shap'](fold), shap=np.asarray(sv, np.float32), x=xf[pick], label=yf[pick], inputs=np.array(exprs))
            info['meanAbsShap'] = {e: float(v) for e, v in zip(exprs, np.abs(sv).mean(axis=0))}
        json.dump(info, open(P['training'](fold), 'w'), indent=1)
        print(f"DATADRIVEN_BDT TRAINED {r['name']} fold {fold}: AUC train {info['auc']['train']:.4f} test {info['auc']['test']:.4f}; "
              f"exported trees == xgboost to {diff:.1e} -> {P['ensemble'](fold)}", flush=True)


def load_ensembles(r):
    P = paths(r)
    ens = []
    for f in (0, 1):
        if not os.path.exists(P['ensemble'](f)):
            fatal(f"fold {f} is not trained: {P['ensemble'](f)}")
        ens.append(tree_ensemble.read(P['ensemble'](f)))
    if ens[0]['inputs'] != ens[1]['inputs']:
        fatal('the two folds have different inputs')
    return ens


APPLY_REFERENCE_SAMPLE = 2000      # candidates per fold whose xgboost score is compared with the reference evaluation
APPLY_REFERENCE_TOLERANCE = 1e-5   # largest allowed |xgboost - reference| score difference


def score_candidates_xgboost(r, ens, c):
    """Fast scores for apply: the saved xgboost model of each fold (fold 1 - parity, as the Maker). The score definition is
    the exported tree-ensemble file; a sample of each fold is evaluated with that reference too and must agree to
    APPLY_REFERENCE_TOLERANCE (float32 inside xgboost), otherwise the action stops."""
    import xgboost as xgb
    P = paths(r)
    x = input_matrix(c, ens[0]['inputs'])
    par = run_parity(c)
    s = np.empty(c.n)
    rng = np.random.default_rng(int(r['learner']['seed']))
    for p in (0, 1):
        f = 1 - p
        if not os.path.exists(P['model'](f)):
            fatal(f"no xgboost model of fold {f} ({P['model'](f)}): train again")
        booster = xgb.Booster()
        booster.load_model(P['model'](f))
        m = np.flatnonzero(par == p)
        if not len(m):
            continue
        prob = booster.predict(xgb.DMatrix(x[m], missing=np.nan)).astype(np.float64)
        s[m] = 2. * prob - 1.
        pick = rng.choice(len(m), size=min(APPLY_REFERENCE_SAMPLE, len(m)), replace=False)
        diff = float(np.max(np.abs(s[m][pick] - tree_ensemble.score(ens[f], x[m][pick]))))
        if diff > APPLY_REFERENCE_TOLERANCE:
            fatal(f'fold {f}: xgboost and the tree-ensemble file differ by {diff:.1e} (> {APPLY_REFERENCE_TOLERANCE:g}); train again')
        print(f'DATADRIVEN_BDT apply fold {f}: {len(m)} candidates, xgboost == tree-ensemble file to {diff:.1e} '
              f'({len(pick)} checked)', flush=True)
    return s


def score_candidates(ens, c):
    """Reference scores (the tree-ensemble file, as the Maker): fold 1 - parity. Used where exactness matters (check)."""
    x = input_matrix(c, ens[0]['inputs'])
    par = run_parity(c)
    s = np.empty(c.n)
    for p in (0, 1):
        m = par == p
        s[m] = tree_ensemble.score(ens[1 - p], x[m])
    return s


def write_apply(path, r, mass, score, c):
    lo, hi = [float(v) for v in r['evaluation']['massRange']]
    nb = int(round((hi - lo) / float(r['evaluation']['binWidth'])))
    h, xe, ye = np.histogram2d(mass, score, bins=[nb, 4000], range=[[lo, hi], [-1., 1.]])
    with uproot.recreate(path) as o:
        o['mass_vs_score'] = (h, xe, ye)
        o.mktree('candidates', dict(score=np.float32, mass=np.float32, runId=np.int32, cent9=np.int32, px=np.float32, py=np.float32,
                                    pz=np.float32))   # a TTree (uproot would otherwise default to RNTuple)
        o['candidates'].extend(dict(score=score.astype(np.float32), mass=mass.astype(np.float32), runId=c.v['runId'].astype(np.int32),
                               cent9=c.v['cent9'].astype(np.int32), px=c.v['px'].astype(np.float32), py=c.v['py'].astype(np.float32),
                               pz=c.v['pz'].astype(np.float32)))


def action_apply(r):
    ens = load_ensembles(r)
    files = input_files(r['data']['trees'])
    branches = tree_branch_names(files, r['data']['tree'])
    c = Candidates(files, r['data']['tree'], needed_variables(r, ens[0]['inputs'], branches), cache=paths(r)['cache'],
                   fraction=float(r['data']['eventFraction']))
    c = c.subset(domain(r, c))
    lo, hi = [float(v) for v in r['evaluation']['massRange']]
    m = Labels(r).mass(c)
    c = c.subset((m > lo) & (m < hi))
    m = Labels(r).mass(c)
    s = score_candidates_xgboost(r, ens, c)
    P = paths(r)
    write_apply(P['apply'], r, m, s, c)
    print(f"DATADRIVEN_BDT APPLIED {r['name']}: {c.n} candidates -> {P['apply']}", flush=True)


def action_study(r):
    """Training-size study: train both folds on a fraction of their candidates, score the other fold, compare."""
    import datadriven_bdt_report as report
    ev = r['evaluation']
    if 'trainingSize' not in ev['checks']:
        fatal('evaluation.checks has no trainingSize')
    c, exprs = prepare(r)
    c = c.subset(domain(r, c))
    lo, hi = [float(v) for v in ev['massRange']]
    lab_all = Labels(r)
    m_all = lab_all.mass(c)
    c = c.subset((m_all > lo) & (m_all < hi))
    m = lab_all.mass(c)
    lab = lab_all.label(c)
    x = input_matrix(c, exprs)
    lab_train = np.where(training_mask(r, c) & (lab >= 0), lab, -1)   # training samples only (data.trainingSelection)
    w = sideband_weights(r, c, lab_train)
    par = run_parity(c)
    rows = []
    for frac in ev['trainingSize']['fractions']:
        frac = float(frac)
        s = np.empty(c.n)
        n_train = []
        auc = []
        for fold in (0, 1):
            sel = (par == fold) & (lab_train >= 0) & (w > 0)
            rng = np.random.default_rng(int(r['learner']['seed']) + 1000 + fold)
            sub = sel & (rng.random(c.n) < frac)
            clf, train, ww = train_fold(r, x[sub], lab_train[sub], w[sub], seed_offset=fold)
            other = par == 1 - fold
            s[other] = 2 * clf.predict_proba(x[other])[:, 1] - 1
            n_train.append(int(sub.sum()))
            held = other & (lab_train >= 0) & (w > 0)
            auc.append(weighted_auc(lab[held], s[held], w[held]))
        spec = report.Spectra.from_arrays(m, s, r)
        spec.scan()
        rows.append(dict(fraction=frac, trainingCandidates=n_train, aucOtherFold=auc,
                         signalAtPurity=spec.at(float(ev['trainingSize']['purity']), 1),
                         significanceAtPurity=spec.at(float(ev['trainingSize']['purity']), 4)))
        print(f"DATADRIVEN_BDT STUDY fraction {frac}: training candidates {n_train}, AUC (other fold) {auc[0]:.4f} {auc[1]:.4f}, "
              f"S at purity {ev['trainingSize']['purity']} = {rows[-1]['signalAtPurity']:.0f}", flush=True)
    P = paths(r)
    json.dump(dict(name=r['name'], purity=float(ev['trainingSize']['purity']), rows=rows), open(P['study'], 'w'), indent=1)
    print(f"DATADRIVEN_BDT STUDY -> {P['study']}")


def action_model(r, d):
    P = paths(r)
    for k in ('particle', 'validation'):
        if not str(r[k]).strip() or '#' in str(r[k]):
            fatal(f'model needs a non-empty "{k}" without "#"')
    for k in ('eventSelection', 'preselection'):
        if '#' in str(r['data'][k]):
            fatal(f"'#' cannot be stored in a model file: {r['data'][k]}")
    infos = []
    for f in (0, 1):
        if not os.path.exists(P['training'](f)):
            fatal(f"fold {f} is not trained: {P['training'](f)}")
        info = json.load(open(P['training'](f)))
        if info['recipe'] != json.loads(json.dumps(trained_part(r))):
            fatal(f"fold {f} was trained with a different recipe ({P['training'](f)}); train again "
                  f"(only {FREE_KEYS} may change after training)")
        infos.append(info)
    os.makedirs(d, exist_ok=True)
    for f in (0, 1):
        dest = f"{d}/{r['name']}_f{f}.trees.txt"
        open(dest, 'w').write(open(P['ensemble'](f)).read())
    fn = f"{d}/{r['name']}.model.yaml"
    counts = '; '.join(f"fold {i['fold']} (runParity): window {i['candidates']['window']}, side bands {i['candidates']['sidebands']}"
                       for i in infos)
    with open(fn, 'w') as o:
        o.write('# BDT model for StCandidateFemtoMaker, written by common/macro/datadriven_bdt (action model).\n'
                '# Do not edit by hand: train again and rewrite. The input expressions are in the tree-ensemble files.\n')
        o.write(f"name: {r['name']}\nparticle: {r['particle']}\nfoldBy: runParity\nweights: {r['name']}_f%d.trees.txt\n"
                f"eventSelection: {r['data']['eventSelection']}\npreselection: {r['data']['preselection']}\n"
                f"trainingData: {r['data']['trees']} (tree {r['data']['tree']})\ntrainingCandidates: {counts}\n"
                f"recipe: {r['_path_in_repo']}\ngitRevision: {infos[0]['gitRevision']}\nlearner: {infos[0]['learner']}, "
                + ', '.join(f'{k}={v}' for k, v in sorted(r['learner']['parameters'].items()))   # no braces: flat YAML
                + f", seed={r['learner']['seed']}\n"
                f"validation: {r['validation']}\ncreated: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    print(f"DATADRIVEN_BDT MODEL {r['name']} -> {fn}")


def action_check(r, pattern):
    """Candidate trees written by StCandidateFemtoMaker WITH this model: in-domain Maker scores == this model's scores
    (to 1e-6), out-of-domain Maker scores < -1."""
    ens = load_ensembles(r)
    files = input_files(pattern)
    branches = tree_branch_names(files, r['data']['tree'])
    if MAKER_SCORE not in branches:
        fatal(f'{files[0]}: no {MAKER_SCORE} branch')
    c = Candidates(files, r['data']['tree'], needed_variables(r, ens[0]['inputs'], branches), with_maker_score=True)
    maker = c.v[MAKER_SCORE].astype(np.float64)
    dom = domain(r, c)
    bad = ~np.isfinite(maker)
    mismatch = int(bad.sum()) + int(((~dom) & (maker >= -1.) & ~bad).sum()) + int((dom & (maker < -1.) & ~bad).sum())
    inside = c.subset(dom & ~bad & (maker >= -1.))
    s = score_candidates(ens, inside) if inside.n else np.zeros(0)
    md = float(np.max(np.abs(s - maker[dom & ~bad & (maker >= -1.)]))) if inside.n else 0.
    print(f"DATADRIVEN_BDT CHECK {r['name']} on {pattern} ({len(files)} files): {c.n} candidates, {int(dom.sum())} inside the training "
          f"domain, {int((~dom).sum())} outside; domain mismatches {mismatch}; max |maker - this model| = {md:.3g}")
    if c.n == 0 or mismatch or md > 1e-6:
        print('DATADRIVEN_BDT CHECK FAILED')
        sys.exit(1)
    print('DATADRIVEN_BDT CHECK OK')


def action_inspect(r):
    """Size of the input before reading it: files, candidates, needed branches and the memory they take."""
    from concurrent.futures import ThreadPoolExecutor
    files = input_files(r['data']['trees'])
    tree = r['data']['tree']
    branches = tree_branch_names(files, tree)
    want = [n for n in branches if n in needed_variables(r, input_expressions(r, branches), branches) or n in EVENT_BRANCHES]

    def entries(f):
        with open_tree_file(f) as fh:
            t = fh[tree]
            return t.num_entries, sum(t[b].interpretation.numpy_dtype.itemsize for b in want)

    with ThreadPoolExecutor(max_workers=8) as pool:
        res = list(pool.map(entries, files))
    n = sum(e for e, _ in res)
    per = res[0][1] if res else 0
    f = float(r['data']['eventFraction'])
    print(f'DATADRIVEN_BDT INSPECT {tree}: {len(files)} files, {n} candidates, {len(want)} branches read ({per} bytes per candidate)')
    size = lambda b: f'{b / 1e9:.2f} GB' if b >= 1e9 else f'{b / 1e6:.1f} MB'
    print(f'  memory for the candidates read: {size(n * per)} at eventFraction 1, {size(n * per * f)} at eventFraction {f}'
          ' (training needs a few times this for the inputs and the learner)')
    avail = available_memory_bytes()
    need = MEMORY_FIXED_BYTES + MEMORY_PER_READ_BYTE * n * per * f
    print(f'  expected peak memory of train/apply/study: about {size(need)}; available now: '
          + (size(avail) if avail is not None else 'unknown') + ('  -> lower data.eventFraction' if avail is not None and need > avail else ''))
    print('  reading time grows with the number of files opened (merge many small files with hadd first); the candidates are '
          'cached in the output directory after the first read')


def main(argv):
    if len(argv) < 3:
        fatal('usage: datadriven_bdt.py <recipe.yaml> <train|apply|study|report|model|check|inspect> [args]')
    r = load_recipe(argv[1])
    a = argv[2]
    if a == 'train' and len(argv) == 3:
        action_train(r)
    elif a == 'apply' and len(argv) == 3:
        action_apply(r)
    elif a == 'study' and len(argv) == 3:
        action_study(r)
    elif a == 'report':
        import datadriven_bdt_report as report
        report.make_report(r, argv[3:])
    elif a == 'model' and len(argv) == 4:
        action_model(r, argv[3])
    elif a == 'check' and len(argv) == 4:
        action_check(r, argv[3])
    elif a == 'inspect' and len(argv) == 3:
        action_inspect(r)
    else:
        fatal(f'bad action or arguments: {" ".join(argv[2:])}')


if __name__ == '__main__':
    main(sys.argv)
