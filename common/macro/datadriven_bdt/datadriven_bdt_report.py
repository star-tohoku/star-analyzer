# datadriven_bdt_report.py -- one PDF (with bookmarks) with every evaluation item of a datadriven_bdt recipe.
# Called by datadriven_bdt.py (action report); also provides Spectra for the training-size study.
#
# Inputs: the recipe's apply output (<output>/<name>_apply.root), its training records (<name>_f<fold>.training.json,
# <name>_f<fold>.shap.npz) and, if evaluation.checks has trainingSize, the study output (<name>_training_size.json).
# Further apply outputs can be compared: `report label=file.root ...` (any file with a mass_vs_score TH2, e.g. the
# apply output of another BDT); ratios are taken to this recipe's BDT.
#
# S and B at each score threshold come from a fit of the mass spectrum: a Student-t peak + pol2
# (B = pol2 integral in the window). Sculpting checks: background efficiency eps(m) in the side bands and the gaps next
# to the window, histogram ratios at equal purity, background-model variation (pol2 / pol4 / fit-free side-band count).
import json, os, subprocess, tempfile
import numpy as np
import uproot
from scipy.optimize import curve_fit
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402


def peak(x, a, nu, mu, s):
    # Student-t peak: the K0S and Xi- peaks have non-Gaussian tails; two Gaussians counted them as
    # background (K0S fit deviance/ndf 4.1 with two Gaussians, 1.2 with Student-t). mu stays parameter 2.
    return a * (1. + ((x - mu) / s) ** 2 / nu) ** (-(nu + 1.) / 2.)


def readable(expression):
    """Short display name of an input expression: log10(max(x,1e-6)) -> log10(x)."""
    import re
    return re.sub(r'log10\(max\(([A-Za-z0-9_]+),1e-6\)\)', r'log10(\1)', expression)


class Spectra:
    """Mass spectra versus score threshold (a mass x score 2D histogram)."""

    def __init__(self, c, xe, ye, r, cand=None):
        self.c, self.xe, self.ye = c, xe, ye
        self.x = 0.5 * (xe[1:] + xe[:-1])
        w = [float(v) for v in r['labels']['window']]
        self.wl, self.wh = w
        self.win = (self.x > w[0]) & (self.x < w[1])
        self.sb = np.zeros(len(self.x), bool)
        for lo, hi in r['labels']['sidebands']:
            self.sb |= (self.x > float(lo)) & (self.x < float(hi))
        self.mid = 0.5 * (w[0] + w[1])
        self.cand = cand
        self.curve = None

    @classmethod
    def from_root(cls, fn, r):
        f = uproot.open(fn)
        c, xe, ye = f['mass_vs_score'].to_numpy()
        cand = None
        if 'candidates' in f:
            cand = f['candidates'].arrays(library='np')
        return cls(c, xe, ye, r, cand)

    @classmethod
    def from_arrays(cls, mass, score, r):
        lo, hi = [float(v) for v in r['evaluation']['massRange']]
        nb = int(round((hi - lo) / float(r['evaluation']['binWidth'])))
        c, xe, ye = np.histogram2d(mass, score, bins=[nb, 4000], range=[[lo, hi], [-1., 1.]])
        return cls(c, xe, ye, r)

    def y(self, thr):
        return self.c[:, np.searchsorted(self.ye, thr - 1e-12):].sum(axis=1)

    def fit(self, y, npol=2):
        w = np.sqrt(np.maximum(y, 1))
        top = max(y.max() - np.median(y), 1)
        f = lambda xx, *p: peak(xx, *p[:4]) + np.polyval(p[4:][::-1], xx - self.mid)  # noqa: E731
        width = self.wh - self.wl
        p0 = [top * .9, 5., self.mid, width / 9, max(np.median(y[self.sb]) if self.sb.any() else 1, .1)] + [0] * npol
        lo = [0., .5, self.wl, width / 100] + [-np.inf] * (npol + 1)
        hi = [np.inf, 200., self.wh, width] + [np.inf] * (npol + 1)
        try:
            p, _ = curve_fit(f, self.x, y, p0=p0, sigma=w, bounds=(lo, hi), maxfev=40000)
        except Exception:
            return None
        bkg = np.polyval(p[4:][::-1], self.x - self.mid)
        chi2 = (((y - f(self.x, *p)) / w) ** 2).sum() / max(len(self.x) - len(p), 1)
        return dict(p=p, bkg=bkg, chi2=chi2)

    def count(self, y, fitres, half=None):
        win = self.win if half is None else (np.abs(self.x - fitres['p'][2]) < half)
        return (y - fitres['bkg'])[win].sum(), fitres['bkg'][win].sum()

    def scan(self, n=40):
        """Thresholds keeping 100% ... 1% of the window candidates; columns thr, S, B, purity, S/sqrt(S+B), chi2."""
        wy = self.c[self.win].sum(axis=0)
        cum = np.cumsum(wy[::-1])[::-1]
        tot = cum[0]
        rows = []
        for q in np.exp(np.linspace(0, np.log(0.01), n)):
            k = min(max(np.searchsorted(-cum, -q * tot), 0), len(self.ye) - 2)
            thr = self.ye[k]
            yy = self.y(thr)
            if yy[self.win].sum() < 50:
                continue
            res = self.fit(yy)
            if res is None:
                continue
            S, B = self.count(yy, res)
            if S <= 0:
                continue
            rows.append((thr, S, B, S / (S + B), S / np.sqrt(S + B), res['chi2']))
        a = np.array(rows) if rows else np.zeros((0, 6))
        a = a[np.argsort(a[:, 3])]
        _, iu = np.unique(a[:, 3], return_index=True)
        self.curve = a[iu]
        return self.curve

    def at(self, pur, col):
        a = self.curve
        if a is None or len(a) < 2 or not (a[:, 3].min() <= pur <= a[:, 3].max()):
            return float('nan')
        return float(np.interp(pur, a[:, 3], a[:, col]))


class Book:
    """PdfPages with one bookmark per page; bookmarks are added with ghostscript at the end."""

    def __init__(self, out):
        self.out = out
        self.tmp = out + '.tmp.pdf'
        self.pdf = PdfPages(self.tmp)
        self.titles = []

    def page(self, fig, title):
        self.pdf.savefig(fig)
        plt.close(fig)
        self.titles.append(title)

    def text(self, title, lines, size=8):
        fig = plt.figure(figsize=(11.7, 8.3))
        fig.text(0.03, 0.96, title, fontsize=13, weight='bold', va='top')
        fig.text(0.03, 0.91, '\n'.join(lines), family='monospace', fontsize=size, va='top')
        self.page(fig, title)

    def close(self):
        self.pdf.close()
        with tempfile.NamedTemporaryFile('w', suffix='.ps', delete=False) as m:
            for i, t in enumerate(self.titles):
                t = t.replace('(', '[').replace(')', ']')
                m.write(f'[/Title ({t}) /Page {i + 1} /OUT pdfmark\n')
        rc = subprocess.run(['gs', '-q', '-o', self.out, '-sDEVICE=pdfwrite', self.tmp, m.name], capture_output=True, text=True)
        os.unlink(m.name)
        if rc.returncode != 0:
            print(f'DATADRIVEN_BDT note: bookmarks not added (gs failed: {rc.stderr.strip()[:200]}); PDF without bookmarks')
            os.replace(self.tmp, self.out)
        else:
            os.unlink(self.tmp)


def make_report(r, extra):
    import datadriven_bdt as tool
    P = tool.paths(r)
    ev = r['evaluation']
    purs = [float(p) for p in ev['purities']]
    if not os.path.exists(P['apply']):
        tool.fatal(f"run apply first: {P['apply']}")
    if 'trainingSize' in ev['checks'] and not os.path.exists(P['study']):     # checked before the PDF is opened
        tool.fatal(f"evaluation.checks has trainingSize: run the study first ({P['study']})")
    labs = [r['name']] + [a.split('=', 1)[0] for a in extra]
    files = [P['apply']] + [a.split('=', 1)[1] for a in extra if '=' in a]
    if len(files) != len(labs):
        tool.fatal('report arguments: label=apply.root')
    sp = {l: Spectra.from_root(f, r) for l, f in zip(labs, files)}
    for l in labs:
        sp[l].scan()
    ref = labs[0]
    infos = [json.load(open(P['training'](f))) for f in (0, 1)]
    book = Book(P['report'])
    sp_ref = sp[ref]

    # 1. summary
    L = [f"recipe: {r['_path_in_repo']}", f"particle: {r['particle']}", f"trees: {r['data']['trees']} ({r['data']['tree']})",
         f"eventSelection: {r['data']['eventSelection']}   preselection: {r['data']['preselection']}",
         f"window {r['labels']['window']}   side bands {r['labels']['sidebands']}   fit: Student-t + pol2 in {ev['massRange']}",
         f"learner: {infos[0]['learner']} {json.dumps(r['learner']['parameters'], sort_keys=True)} seed {r['learner']['seed']}",
         f"inputs ({len(infos[0]['inputs'])}): " + ', '.join(infos[0]['inputs'][:12]) + (' ...' if len(infos[0]['inputs']) > 12 else ''), '']
    for i in infos:
        L.append(f"fold {i['fold']}: window {i['candidates']['window']}, side bands {i['candidates']['sidebands']}; "
                 f"AUC train {i['auc']['train']:.4f}, test {i['auc']['test']:.4f}; exported trees vs xgboost {i['exportedVsXgboostMaxDifference']:.1e}")
    L += ['', f"{'BDT':16s} {'max S/sqrt(S+B)':>16s} {'@purity':>8s} | S at purity " + ' '.join(f'{p:>8.2f}' for p in purs)]
    for l in labs:
        a = sp[l].curve
        if len(a) == 0:
            L.append(f'{l:16s} no valid fit'); continue
        i = np.argmax(a[:, 4])
        L.append(f"{l:16s} {a[i, 4]:16.1f} {a[i, 3]:8.3f} | {'':12s}" + ' '.join(f'{sp[l].at(p, 1):8.0f}' for p in purs))
    L += ['', f'S/sqrt(S+B) at equal purity, ratio to {ref}:', f"{'BDT':16s} " + ' '.join(f'{p:>8.2f}' for p in purs)]
    for l in labs:
        L.append(f'{l:16s} ' + ' '.join(f'{sp[l].at(p, 4) / sp[ref].at(p, 4):8.3f}' for p in purs))
    L += ['', 'Score thresholds at equal purity:', f"{'BDT':16s} " + ' '.join(f'{p:>8.2f}' for p in purs)]
    for l in labs:
        L.append(f'{l:16s} ' + ' '.join(f'{sp[l].at(p, 0):+8.3f}' for p in purs))
    book.text('Summary', L)

    # 2. curves
    fig, ax = plt.subplots(1, 3, figsize=(16, 5))
    pg = np.linspace(min(purs), max(purs), 80)
    for l in labs:
        a = sp[l].curve
        ax[0].plot(a[:, 3], a[:, 4], '-o', ms=2, label=l); ax[1].plot(a[:, 3], a[:, 1], '-o', ms=2, label=l)
        ax[2].plot(pg, [sp[l].at(p, 4) / sp[ref].at(p, 4) for p in pg], label=f'{l} / {ref}')
    ax[0].set_ylabel('S/sqrt(S+B)'); ax[1].set_ylabel('S'); ax[1].set_yscale('log'); ax[2].set_ylabel('ratio of S/sqrt(S+B)')
    ax[2].axhline(1, color='k', lw=.5)
    for a_ in ax:
        a_.set_xlabel('purity in the window'); a_.grid(alpha=.3); a_.legend(fontsize=8)
    plt.tight_layout(); book.page(fig, 'Significance and signal versus purity')

    # 3. spectra
    fig, ax = plt.subplots(2, len(labs), figsize=(5 * len(labs), 8), squeeze=False)
    for j, l in enumerate(labs):
        s = sp[l]
        for k, pur in enumerate(purs):
            t = s.at(pur, 0)
            if np.isnan(t):
                continue
            y = s.y(t)
            for i in (0, 1):
                ax[i, j].step(s.xe[:-1], y, where='post', color=plt.cm.viridis(k / max(len(purs) - 1, 1)),
                              label=f'purity {pur}: score>{t:+.3f}, N(window) {y[s.win].sum():.0f}')
        for i in (0, 1):
            for e in (s.wl, s.wh):
                ax[i, j].axvline(e, color='m', lw=.6)
            ax[i, j].set_xlabel('mass')
        ax[1, j].set_yscale('log'); ax[0, j].set_title(l); ax[0, j].legend(fontsize=6)
    plt.tight_layout(); book.page(fig, 'Mass spectra at fixed purity')

    pref = float(ev['sculptingPurity'])
    if 'sculpting' in ev['checks']:
        checks = [f'Background efficiency eps(m) = N(score > thr) / N(all) per 2 bins at purity {pref}; line fitted in the side bands;',
                  'pulls of the bins in the gaps between side bands and window (signal tails can raise the innermost ones).', '']
        fig, ax = plt.subplots(1, len(labs), figsize=(5 * len(labs), 4), squeeze=False)
        for j, l in enumerate(labs):
            s = sp[l]; t = s.at(pref, 0)
            if np.isnan(t):
                continue
            k = s.c[:, np.searchsorted(s.ye, t - 1e-12):].sum(axis=1); n = s.c.sum(axis=1)
            m = len(k) // 2 * 2
            k2 = k[:m].reshape(-1, 2).sum(1); n2 = n[:m].reshape(-1, 2).sum(1); x2 = s.x[:m].reshape(-1, 2).mean(1)
            sb2 = s.sb[:m].reshape(-1, 2).all(1); win2 = s.win[:m].reshape(-1, 2).any(1)
            e = k2 / np.maximum(n2, 1); er = np.sqrt(np.maximum(e * (1 - e), 1e-12) / np.maximum(n2, 1))
            p = np.polyfit(x2[sb2], e[sb2], 1, w=1 / er[sb2]); line = np.polyval(p, x2)
            chi = (((e - line) / er)[sb2] ** 2).sum() / max(sb2.sum() - 2, 1)
            gap = ~sb2 & ~win2 & (x2 > x2[sb2].min()) & (x2 < x2[sb2].max())
            checks.append(f'{l:16s} thr {t:+.3f}: side-band line chi2/ndf {chi:.2f}; gap-bin pulls ' +
                          ' '.join(f'{v:+.1f}' for v in ((e - line) / er)[gap]))
            ax[0, j].errorbar(x2, e, er, fmt='o', ms=3); ax[0, j].plot(x2, line, 'r-'); ax[0, j].set_title(f'{l}: eps(m), purity {pref}')
            for b in (s.wl, s.wh):
                ax[0, j].axvline(b, color='m', lw=.6)
        plt.tight_layout(); book.page(fig, 'Sculpting: background efficiency versus mass')
        fig, ax = plt.subplots(2, 1, figsize=(9, 7), sharex=True, gridspec_kw={'height_ratios': [2, 1]})
        b0 = None
        drawn = []
        for l in labs:
            s = sp[l]; t = s.at(pref, 0)
            if np.isnan(t):
                continue
            y = s.y(t); ax[0].step(s.xe[:-1], y, where='post', label=f'{l} (score > {t:+.3f})')
            drawn.append((s, y))
            if b0 is None:
                b0 = y
            elif len(y) == len(b0):
                q = y / np.maximum(b0, 1)
                ax[1].errorbar(s.x, q, q * np.sqrt(1 / np.maximum(y, 1) + 1 / np.maximum(b0, 1)), fmt='o', ms=2, label=f'{l} / {ref}')
        if len(drawn) == 1:
            # One model: no second BDT to divide by; show selected / before the score cut, scaled to 1 in the
            # side bands, so that the lower panel is not empty (a window-shaped step would again indicate sculpting).
            s, y = drawn[0]
            yall = s.y(-np.inf); q = y / np.maximum(yall, 1)
            qe = np.sqrt(q * np.maximum(1 - q, 0) / np.maximum(yall, 1))  # binomial
            norm = q[s.sb].mean() if s.sb.any() and q[s.sb].mean() > 0 else 1.
            ax[0].step(s.xe[:-1], yall, where='post', color='0.5', label='before the score cut')
            ax[1].errorbar(s.x, q / norm, qe / norm, fmt='o', ms=2, label='selected / before the cut (side bands = 1)')
            ax[1].axhline(1, color='0.6', lw=.6)
        ax[1].legend(fontsize=8)
        for a_ in ax:
            for b in (sp_ref.wl, sp_ref.wh):
                a_.axvline(b, color='m', lw=.6)
        ax[0].set_yscale('log'); ax[0].legend(fontsize=8)
        ax[0].set_title(f'equal purity {pref}: a window-shaped step in the ratio would indicate sculpting')
        ax[1].set_xlabel('mass'); ax[1].set_ylabel('ratio'); plt.tight_layout(); book.page(fig, 'Sculpting: ratios at equal purity')
    else:
        checks = []
    if 'backgroundModels' in ev['checks']:
        checks += ['', 'Background model: S and purity with Student-t + pol2 / + pol4, and fit-free (window - side bands x width ratio):']
        for l in labs:
            s = sp[l]
            for pur in purs:
                t = s.at(pur, 0)
                if np.isnan(t):
                    continue
                y = s.y(t); r2 = s.fit(y, 2); r4 = s.fit(y, 4)
                nw = y[s.win].sum(); Ssb = nw - y[s.sb].sum() * s.win.sum() / max(s.sb.sum(), 1)
                if r2 is None or r4 is None:
                    continue
                S2, B2 = s.count(y, r2); S4, B4 = s.count(y, r4)
                checks.append(f"{l:16s} purity {pur:.2f}: pol2 S {S2:9.0f} p {S2/(S2+B2):.3f} chi2 {r2['chi2']:.2f} | "
                              f"pol4 S {S4:9.0f} p {S4/(S4+B4):.3f} chi2 {r4['chi2']:.2f} | side bands S {Ssb:9.0f} p {Ssb/nw:.3f}")
        checks += ['', 'Window half-width scan (around the fitted mean, MeV if the mass is in GeV): S at purity ' +
                   ' / '.join(f'{p}' for p in purs)]
        for l in labs:
            s = sp[l]; full = (s.wh - s.wl) / 2
            for frac in ev['windowScan']:
                half = full * float(frac)
                rows = []
                for thr in s.curve[:, 0]:
                    y = s.y(thr); res = s.fit(y)
                    if res is None:
                        continue
                    S, B = s.count(y, res, half); rows.append((S / (S + B), S))
                a = np.array(sorted(rows))
                v = [np.interp(p, a[:, 0], a[:, 1]) if len(a) and a[:, 0].min() <= p <= a[:, 0].max() else np.nan for p in purs]
                checks.append(f'{l:16s} +-{half * 1e3:6.2f}: ' + ' '.join(f'{q:9.0f}' for q in v))
    if checks:
        book.text('Sculpting and background-model checks, window width', checks, 7)

    # purity in bins (this recipe's candidates only)
    if ev['purityBins'] and sp_ref.cand is not None:
        d = sp_ref.cand
        pt = np.hypot(d['px'], d['py']); pz = d['pz'].astype(np.float64); m = d['mass'].astype(np.float64)
        E = np.sqrt(pt ** 2 + pz ** 2 + m ** 2)
        var = dict(pt=pt, y=0.5 * np.log((E + pz) / (E - pz)), cent9=d['cent9'])
        wr = sp_ref.win.sum() / max(sp_ref.sb.sum(), 1)
        wm = (m > sp_ref.wl) & (m < sp_ref.wh)
        sbm = np.zeros(len(m), bool)
        for lo, hi in r['labels']['sidebands']:
            sbm |= (m > float(lo)) & (m < float(hi))
        lines = ['Purity from side bands (1 - N(side bands) x width ratio / N(window)) in bins (y: rapidity in the laboratory frame),',
                 'at the thresholds of overall purity ' + ' / '.join(str(p) for p in purs), '']
        for pur in purs:
            t = sp_ref.at(pur, 0)
            if np.isnan(t):
                continue
            sel = d['score'] > t
            for k, edges in ev['purityBins'].items():
                nw, _ = np.histogram(var[k][sel & wm], edges); ns, _ = np.histogram(var[k][sel & sbm], edges)
                lines.append(f'p{pur:.2f} {k:5s}: ' + '  '.join(f'[{edges[i]:g},{edges[i + 1]:g}) {nw[i]} {1 - ns[i] * wr / max(nw[i], 1):.2f}'
                                                              for i in range(len(nw)) if nw[i] > 0))
        book.text('Purity in kinematic and centrality bins', lines, 6.5)

    if 'overtraining' in ev['checks']:
        fig, ax = plt.subplots(1, 2, figsize=(14, 5))
        e = np.linspace(-1, 1, 51); xc = 0.5 * (e[1:] + e[:-1])
        for f, i in enumerate(infos):
            h = {k: np.array(v, float) for k, v in i['scoreHistograms'].items()}
            for k, st, col in (('trainWindow', '-', 'b'), ('testWindow', 'o', 'b'), ('trainSidebands', '-', 'r'), ('testSidebands', 'o', 'r')):
                y = h[k] / max(h[k].sum(), 1e-300)
                if st == '-':
                    ax[f].step(xc, y, where='mid', color=col, label=k)
                else:
                    ax[f].plot(xc, y, st, color=col, ms=3, label=k)
            ax[f].set_yscale('log'); ax[f].set_xlabel('score'); ax[f].legend(fontsize=8)
            ax[f].set_title(f"fold {f}: AUC train {i['auc']['train']:.4f}, test {i['auc']['test']:.4f} (weighted, window vs side bands)")
        plt.tight_layout(); book.page(fig, 'Overtraining: training versus test scores')

    if 'importance' in ev['checks']:
        import shap
        names = infos[0]['inputs']
        pretty = [readable(n) for n in names]
        # mean |SHAP| per input (both folds, sorted by their average): the main importance page
        ms = np.array([[i.get('meanAbsShap', {}).get(n, 0.) for n in names] for i in infos])
        order = np.argsort(-ms.mean(axis=0))
        fig, ax = plt.subplots(figsize=(11, max(5, 0.32 * len(names) + 1.5)))
        for f in (0, 1):
            ax.barh(np.arange(len(names)) + 0.2 * (f - 0.5), ms[f, order], height=0.4, label=f'fold {f}')
        ax.set_yticks(np.arange(len(names))); ax.set_yticklabels([pretty[k] for k in order], fontsize=8); ax.invert_yaxis()
        ax.set_xlabel('mean |SHAP value| (log-odds of window vs side bands)'); ax.legend(fontsize=8); ax.grid(axis='x', alpha=.3)
        ax.set_title(f"{r['name']}: input importance (SHAP, {int(ev['shapCandidates'])} test candidates per fold)")
        plt.tight_layout(); book.page(fig, 'SHAP: mean |SHAP| per input')
        for f in (0, 1):
            if not os.path.exists(P['shap'](f)):
                continue
            z = np.load(P['shap'](f))
            plt.figure(figsize=(11, max(5, 0.32 * len(names) + 1.5)))
            shap.summary_plot(z['shap'], z['x'], feature_names=[readable(n) for n in z['inputs']], show=False, max_display=len(names))
            fig = plt.gcf()
            fig.suptitle(f"{r['name']}, fold {f}: SHAP value of each input for {len(z['x'])} test candidates (window and side bands);\n"
                         'colour = input value; positive SHAP pushes towards the window (signal-like)', fontsize=10)
            plt.tight_layout(); book.page(fig, f'SHAP: summary plot, fold {f}')
        gain = np.array([[i.get('gain', {}).get(n, 0.) for n in names] for i in infos])
        gain = gain / np.maximum(gain.sum(axis=1, keepdims=True), 1e-300)
        L = ['Inputs ranked by mean |SHAP| (average of the two folds). SHAP values are in log-odds of "window" versus "side bands";',
             'the gain is the share of the total split gain. Correlation = Pearson correlation with the mass in the side bands.', '',
             f"{'rank':>4s}  {'input':44s} {'mean|SHAP| f0':>13s} {'f1':>8s} {'gain f0':>8s} {'f1':>8s} {'corr(mass)':>10s}"]
        corr = dict((e, v) for e, v in infos[0]['massCorrelation'])
        for k, j in enumerate(order):
            L.append(f'{k + 1:4d}  {pretty[j][:44]:44s} {ms[0, j]:13.4f} {ms[1, j]:8.4f} {gain[0, j]:8.3f} {gain[1, j]:8.3f} {corr.get(names[j], 0.):+10.3f}')
        book.text('SHAP: table of the inputs', L, 7)
        L = ['Pearson correlation of each input with the mass in the side bands (training sample, before training);',
             f"inputs above {r['inputs']['massCorrelationWarning']} can sculpt the background (see the sculpting pages):", '']
        for e, v in sorted(infos[0]['massCorrelation'], key=lambda t: -abs(t[1])):
            L.append(f'{v:+.3f}  {readable(e)}' + ('   <-- above the warning limit' if abs(v) > float(r['inputs']['massCorrelationWarning']) else ''))
        book.text('Correlation of the inputs with the mass', L, 7)

    if 'trainingSize' in ev['checks']:
        if not os.path.exists(P['study']):
            tool.fatal(f"evaluation.checks has trainingSize: run the study first ({P['study']})")
        st = json.load(open(P['study']))
        rows = st['rows']
        fr = np.array([q['fraction'] for q in rows])
        fig, ax = plt.subplots(1, 2, figsize=(14, 5))
        ax[0].plot(fr, [q['signalAtPurity'] for q in rows], 'o-')
        ax[0].set_ylabel(f"S at purity {st['purity']} (both folds)")
        for f in (0, 1):
            ax[1].plot(fr, [q['aucOtherFold'][f] for q in rows], 'o-', label=f'trained on fold {f}, AUC on the other fold')
        ax[1].set_ylabel('AUC (window vs side bands, weighted)'); ax[1].legend(fontsize=8)
        for a_ in ax:
            a_.set_xscale('log'); a_.set_xlabel('fraction of the training candidates used'); a_.grid(alpha=.3)
        plt.tight_layout(); book.page(fig, 'Training size: learning curve')
        L = [f"{'fraction':>9s} {'training candidates (fold 0, 1)':>32s} {'S at purity':>12s} {'S/sqrt(S+B)':>12s} {'AUC fold 0':>10s} {'AUC fold 1':>10s}"]
        for q in rows:
            L.append(f"{q['fraction']:9.3f} {str(q['trainingCandidates']):>32s} {q['signalAtPurity']:12.0f} {q['significanceAtPurity']:12.1f} "
                     f"{q['aucOtherFold'][0]:10.4f} {q['aucOtherFold'][1]:10.4f}")
        book.text('Training size: table', L)
    book.close()
    print(f"DATADRIVEN_BDT REPORT -> {P['report']}")
