# tree_ensemble.py -- the tree-ensemble file that StCandidateFemtoMaker reads (CandidateTreeEnsemble.h), written from an
# xgboost model, and the reference evaluation of that file.
#
# File format (text, one record per line, '#' lines are comments):
#   format candidateTreeEnsemble 1
#   learner <free text: package, version, objective>
#   baseMargin <float64, %.17g>            # margin added to the sum of the leaves (logit of xgboost's base_score)
#   inputs <n>
#   input <i> <expression>                 # input i = expression over candidate-tree variables (candidate_expression.py)
#   trees <nTrees>
#   tree <t> <nNodes>
#   node <id> <input> <threshold> <yes> <no> <missing> <leaf>
#     split node: input >= 0; go to <yes> if value < threshold (float32 comparison, as xgboost), else <no>;
#                 a NaN value goes to <missing>; <leaf> is 0
#     leaf node:  input = -1, threshold 0, yes = no = missing = -1, <leaf> = leaf value (float32, %.9g)
#   end
# Score (the Maker contract): margin = baseMargin + sum of the leaves (float64), p = 1 / (1 + exp(-margin)), score = 2 p - 1.
# Inputs are float32 (the value of the expression, rounded to float32, as the Maker's mInputs).
import json
import numpy as np

FORMAT_LINE = 'format candidateTreeEnsemble 1'


def from_xgboost(booster, expressions, learner_text):
    """Build the ensemble dict from an xgboost Booster (binary:logistic, gbtree) and the input expressions."""
    model = json.loads(booster.save_raw(raw_format='json'))
    learner = model['learner']
    objective = learner['objective']['name']
    if objective != 'binary:logistic':
        raise ValueError(f'objective must be binary:logistic, not {objective}')
    if learner['gradient_booster']['name'] != 'gbtree':
        raise ValueError('only the gbtree booster is supported')
    base_score = float(learner['learner_model_param']['base_score'])
    base_margin = float(np.log(base_score / (1. - base_score)))
    trees = []
    for t in learner['gradient_booster']['model']['trees']:
        left, right = t['left_children'], t['right_children']
        idx, cond, dleft = t['split_indices'], t['split_conditions'], t['default_left']
        nodes = []
        for i in range(len(left)):
            if left[i] == -1:
                nodes.append((i, -1, 0., -1, -1, -1, float(np.float32(cond[i]))))
            else:
                miss = left[i] if dleft[i] else right[i]
                nodes.append((i, int(idx[i]), float(np.float32(cond[i])), int(left[i]), int(right[i]), int(miss), 0.))
        trees.append(nodes)
    return dict(learner=learner_text, base_margin=base_margin, inputs=list(expressions), trees=trees)


def write(ens, path):
    with open(path, 'w') as o:
        o.write('# tree ensemble for StCandidateFemtoMaker (CandidateTreeEnsemble.h), written by common/macro/datadriven_bdt\n')
        o.write(FORMAT_LINE + '\n')
        o.write(f"learner {ens['learner']}\n")
        o.write(f"baseMargin {ens['base_margin']:.17g}\n")
        o.write(f"inputs {len(ens['inputs'])}\n")
        for i, e in enumerate(ens['inputs']):
            o.write(f'input {i} {e}\n')
        o.write(f"trees {len(ens['trees'])}\n")
        for t, nodes in enumerate(ens['trees']):
            o.write(f'tree {t} {len(nodes)}\n')
            for (i, f, thr, yes, no, miss, leaf) in nodes:
                o.write(f'node {i} {f} {np.float32(thr):.9g} {yes} {no} {miss} {np.float32(leaf):.9g}\n')
        o.write('end\n')


def read(path):
    """Parse a tree-ensemble file (strict: every record in order, as CandidateTreeEnsemble.h)."""
    lines = [l.rstrip('\n') for l in open(path) if l.strip() and not l.startswith('#')]
    it = iter(lines)
    def expect(prefix):
        line = next(it)
        if not line.startswith(prefix):
            raise ValueError(f'{path}: expected "{prefix}", got "{line}"')
        return line[len(prefix):].strip()
    if next(it) != FORMAT_LINE:
        raise ValueError(f'{path}: not a candidateTreeEnsemble 1 file')
    learner = expect('learner ')
    base_margin = float(expect('baseMargin '))
    n_in = int(expect('inputs '))
    inputs = []
    for i in range(n_in):
        rest = expect('input ')
        k, e = rest.split(' ', 1)
        if int(k) != i:
            raise ValueError(f'{path}: input {k} out of order')
        inputs.append(e.strip())
    n_trees = int(expect('trees '))
    trees = []
    for t in range(n_trees):
        k, n_nodes = expect('tree ').split()
        if int(k) != t:
            raise ValueError(f'{path}: tree {k} out of order')
        nodes = []
        for j in range(int(n_nodes)):
            p = expect('node ').split()
            nodes.append((int(p[0]), int(p[1]), float(np.float32(p[2])), int(p[3]), int(p[4]), int(p[5]), float(np.float32(p[6]))))
            if nodes[-1][0] != j:
                raise ValueError(f'{path}: node {p[0]} out of order in tree {t}')
        trees.append(nodes)
    if next(it) != 'end':
        raise ValueError(f'{path}: missing "end"')
    return dict(learner=learner, base_margin=base_margin, inputs=inputs, trees=trees)


def margin(ens, x):
    """x: float32 array (n, n_inputs). Returns float64 margins (leaf sum in float64 + baseMargin)."""
    x = np.asarray(x, dtype=np.float32)
    n = x.shape[0]
    total = np.zeros(n, dtype=np.float64)
    rows = np.arange(n)
    for nodes in ens['trees']:
        a = np.array(nodes, dtype=np.float64)
        feat = a[:, 1].astype(np.int64); thr = a[:, 2].astype(np.float32)
        yes = a[:, 3].astype(np.int64); no = a[:, 4].astype(np.int64); miss = a[:, 5].astype(np.int64); leaf = a[:, 6]
        node = np.zeros(n, dtype=np.int64)
        while True:
            f = feat[node]
            active = f >= 0
            if not active.any():
                break
            v = x[rows[active], f[active]]
            nxt = np.where(np.isnan(v), miss[node[active]], np.where(v < thr[node[active]], yes[node[active]], no[node[active]]))
            node[active] = nxt
        total += leaf[node]
    return total + ens['base_margin']


def score(ens, x):
    """Maker score in [-1, 1]: 2 / (1 + exp(-margin)) - 1."""
    m = margin(ens, x)
    return 2. / (1. + np.exp(-m)) - 1.
