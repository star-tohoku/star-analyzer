// Step 10 T5c Phase A -- mass-dependent background transfer, as a standalone diagnostic.
//
// The nominal subtracts a background whose normalisation is one number per k* bin:
//     S_X(m,k*) = F_X(m,k*) - alpha_X(k*) B_X(m,k*),   alpha from a single control window.
// That assumes F/B is flat in mass. It is not: T2 found alpha drifting with the window it is
// measured in, by 1.7-2.7% between 1.07-1.09 and 1.12-1.16, in the same direction for ROT and for
// MIX. This tool tests whether
//     log T_X(m,k*) = a_X(k*) + b_X(k*) (m - m0)   [+ c (m-m0)^2]
// fitted in a signal-free control region and extrapolated into the signal window closes better
// than the constant, and whether that extrapolation is stable.
//
// Phase A only. Nothing here changes the analysis: checkHistAnaFemtoPhi.C, the maker, the tree and
// the downstream output are untouched, and the nominal is not promoted. The exit of this tool is a
// PASS / FAIL / INCONCLUSIVE against the gates registered in the request, not a new nominal.
//
// Both event classes are corrected. Correcting the SE numerator and leaving the ME denominator on
// the old constant would manufacture a k*-dependent bias out of nothing.
//
// The prefix kmdt is deliberate: the existing "sbleak transfer" graphs in checkHistAnaFemtoPhi.C
// are t = B_window / B_sideband, a different quantity. Nothing here reuses or overwrites them.
//
// Arguments. Every analysis choice is mandatory and is given as "key=value;key=value". A missing
// key is an error that names it and stops -- there is no default window, binning or tolerance
// anywhere in this file. The spec-string form is used rather than 29 positional arguments because
// CINT argument lists of that length are where silent mistakes live; the requirement it implements
// is that nothing analysis-deciding is defaulted, and that is enforced in requireKey().
//
// Usage:
//   root -b -q 'tools/kmfMassDependentTransferClosure.C+("<spec>")'
// Keys: input, sidecar, outRoot, outPointCsv, outSummaryCsv, outPdf, centMin, centMax, dk,
//   lowMerge, kMin, kMax, massRebin, sigLo, sigHi, alphaLo, alphaHi, ctrlLo, ctrlHi, ctrlALo,
//   ctrlAHi, ctrlBLo, ctrlBHi, rsbLo, rsbHi, lsbLo, lsbHi, normLo, normHi, nFolds, nToys, seed,
//   toyTol, toyBiasTol, provenance, syntheticOnly

#include <TAxis.h>
#include <TCanvas.h>
#include <TDecompChol.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TLatex.h>
#include <TMath.h>
#include <TMatrixDSym.h>
#include <TNamed.h>
#include <TParameter.h>
#include <TRandom3.h>
#include <TString.h>
#include <TVectorD.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace kmdt {

// ---------------------------------------------------------------------------------------------
// argument handling: every key is required, and an unknown key is an error too
// ---------------------------------------------------------------------------------------------

typedef std::map<std::string, std::string> Spec;

bool gFatal = false;
bool gQuiet = false;   // set while the failure-mode tests deliberately trip the guards

void fail(const std::string& what) {
  if (gQuiet) return;  // the caller is testing that the guard fires, not reporting a defect
  std::cerr << "[kmdt] ERROR: " << what << std::endl;
  gFatal = true;
}

Spec parseSpec(const std::string& text) {
  Spec out;
  std::stringstream ss(text);
  std::string item;
  while (std::getline(ss, item, ';')) {
    if (item.empty()) continue;
    const size_t eq = item.find('=');
    if (eq == std::string::npos || eq == 0) {
      fail("bad spec item '" + item + "' (expected key=value)");
      continue;
    }
    const std::string k = item.substr(0, eq);
    const std::string v = item.substr(eq + 1);
    if (out.count(k)) fail("duplicate key '" + k + "'");
    out[k] = v;
  }
  return out;
}

std::string requireStr(const Spec& s, const char* key) {
  Spec::const_iterator it = s.find(key);
  if (it == s.end()) { fail(std::string("missing required key '") + key + "'"); return ""; }
  if (it->second.empty()) { fail(std::string("empty value for key '") + key + "'"); return ""; }
  return it->second;
}

double requireNum(const Spec& s, const char* key) {
  const std::string v = requireStr(s, key);
  if (v.empty()) return 0.0;
  return TString(v.c_str()).Atof();
}

int requireInt(const Spec& s, const char* key) {
  const std::string v = requireStr(s, key);
  if (v.empty()) return 0;
  return TString(v.c_str()).Atoi();
}

// ---------------------------------------------------------------------------------------------
// histogram helpers (same projection convention as tools/kmfCrossSourceClosure.C)
// ---------------------------------------------------------------------------------------------

int gUnique = 0;
std::string uniq(const char* stem) { return Form("kmdt_%s_%d", stem, ++gUnique); }

TH2* projectMkkKstar(TH3* h3, int centMin, int centMax, double dk, const char* stem) {
  if (!h3) return 0;
  const int z0 = h3->GetZaxis()->FindBin((double)centMin - 0.01);
  const int z1 = h3->GetZaxis()->FindBin((double)centMax + 0.01);
  TH2D* h2 = new TH2D(uniq(stem).c_str(), h3->GetTitle(),
                      h3->GetNbinsX(), h3->GetXaxis()->GetXmin(), h3->GetXaxis()->GetXmax(),
                      h3->GetNbinsY(), h3->GetYaxis()->GetXmin(), h3->GetYaxis()->GetXmax());
  h2->SetDirectory(0);
  h2->Sumw2();
  for (int ix = 1; ix <= h3->GetNbinsX(); ++ix) {
    for (int iy = 1; iy <= h3->GetNbinsY(); ++iy) {
      double v = 0.0, e2 = 0.0;
      for (int iz = z0; iz <= z1; ++iz) {
        v += h3->GetBinContent(ix, iy, iz);
        const double e = h3->GetBinError(ix, iy, iz);
        e2 += e * e;
      }
      h2->SetBinContent(ix, iy, v);
      h2->SetBinError(ix, iy, TMath::Sqrt(e2));
    }
  }
  const double native = h2->GetYaxis()->GetBinWidth(1);
  const int factor = TMath::Nint(dk / native);
  if (factor < 1 || TMath::Abs(factor * native - dk) > 1e-8) {
    fail(Form("dk=%g is not an integer multiple of the native k* bin %g", dk, native));
    delete h2;
    return 0;
  }
  if (factor > 1) h2->RebinY(factor);
  return h2;
}

TH1* projectMass(TH2* h2, int iy0, int iy1, int massRebin, const char* stem) {
  if (!h2) return 0;
  TH1* h = h2->ProjectionX(uniq(stem).c_str(), iy0, iy1);
  if (!h) return 0;
  h->SetDirectory(0);
  if (massRebin > 1) h->Rebin(massRebin);
  return h;
}

// A window must land on bin edges: FindBin would otherwise widen it silently.
bool windowBins(TAxis* ax, double lo, double hi, int& b0, int& b1, const char* tag) {
  if (!ax) return false;
  b0 = ax->FindBin(lo + 1e-9);
  b1 = ax->FindBin(hi - 1e-9);
  if (b1 < b0) { fail(Form("window %s [%g,%g] holds no bin", tag, lo, hi)); return false; }
  const double eLo = ax->GetBinLowEdge(b0);
  const double eHi = ax->GetBinUpEdge(b1);
  if (TMath::Abs(eLo - lo) > 1e-6 || TMath::Abs(eHi - hi) > 1e-6) {
    fail(Form("window %s [%g,%g] does not land on bin edges (got [%g,%g]); "
              "change the window or the rebin rather than letting it widen",
              tag, lo, hi, eLo, eHi));
    return false;
  }
  return true;
}

// ---------------------------------------------------------------------------------------------
// the transfer fit
// ---------------------------------------------------------------------------------------------

struct Fit {
  bool ok;
  int npar;          // 1 = M0, 2 = M1, 3 = M2
  double p[3];
  double cov[3][3];
  double chi2;
  int ndf;
  double prob;
  double aicc;
  int nUsed;
  int nNonPositive;
  std::string status;
  Fit() : ok(false), npar(0), chi2(0.0), ndf(0), prob(-1.0), aicc(0.0), nUsed(0),
          nNonPositive(0), status("uninitialised") {
    for (int i = 0; i < 3; ++i) { p[i] = 0.0; for (int j = 0; j < 3; ++j) cov[i][j] = 0.0; }
  }
};

double transferAt(const Fit& f, double m, double m0) {
  const double d = m - m0;
  double e = f.p[0];
  if (f.npar > 1) e += f.p[1] * d;
  if (f.npar > 2) e += f.p[2] * d * d;
  return TMath::Exp(e);
}

// Weighted least squares on r_i = log(F_i/B_i) with Var(r_i) = (sF/F)^2 + (sB/B)^2.
// skipLo..skipHi (inclusive, in bin numbers) is held out; pass 0,-1 to use every bin.
bool fitTransfer(TH1* F, TH1* B, double lo, double hi, int order, double m0,
                 int skipLo, int skipHi, Fit& out) {
  out = Fit();
  out.npar = order + 1;
  if (!F || !B) { out.status = "null histogram"; return false; }
  if (F->GetNbinsX() != B->GetNbinsX()) { out.status = "binning mismatch"; return false; }
  int b0 = 0, b1 = 0;
  if (!windowBins(F->GetXaxis(), lo, hi, b0, b1, "control")) { out.status = "bad window"; return false; }

  std::vector<double> x, y, w;
  for (int i = b0; i <= b1; ++i) {
    if (i >= skipLo && i <= skipHi) continue;
    const double fv = F->GetBinContent(i);
    const double bv = B->GetBinContent(i);
    if (!(fv > 0.0) || !(bv > 0.0)) { ++out.nNonPositive; continue; }
    const double ef = F->GetBinError(i);
    const double eb = B->GetBinError(i);
    const double var = (ef / fv) * (ef / fv) + (eb / bv) * (eb / bv);
    if (!(var > 0.0) || !TMath::Finite(var)) { ++out.nNonPositive; continue; }
    x.push_back(F->GetXaxis()->GetBinCenter(i) - m0);
    y.push_back(TMath::Log(fv / bv));
    w.push_back(1.0 / var);
  }
  const int npar = out.npar;
  out.nUsed = (int)x.size();
  if (out.nUsed < npar + 1) { out.status = "too few usable bins"; return false; }

  TMatrixDSym A(npar);
  TVectorD rhs(npar);
  A.Zero();
  rhs.Zero();
  for (size_t k = 0; k < x.size(); ++k) {
    double basis[3] = {1.0, x[k], x[k] * x[k]};
    for (int i = 0; i < npar; ++i) {
      rhs(i) += w[k] * basis[i] * y[k];
      for (int j = 0; j < npar; ++j) A(i, j) += w[k] * basis[i] * basis[j];
    }
  }
  TDecompChol chol(A);
  if (!chol.Decompose()) { out.status = "singular normal matrix"; return false; }
  TMatrixDSym cov(A);
  cov.Invert();
  TVectorD par = rhs;
  { TMatrixD c(cov); par = c * rhs; }
  for (int i = 0; i < npar; ++i) {
    out.p[i] = par(i);
    for (int j = 0; j < npar; ++j) out.cov[i][j] = cov(i, j);
  }
  out.chi2 = 0.0;
  for (size_t k = 0; k < x.size(); ++k) {
    double model = out.p[0];
    if (npar > 1) model += out.p[1] * x[k];
    if (npar > 2) model += out.p[2] * x[k] * x[k];
    const double d = y[k] - model;
    out.chi2 += w[k] * d * d;
  }
  out.ndf = out.nUsed - npar;
  out.prob = (out.ndf > 0) ? TMath::Prob(out.chi2, out.ndf) : -1.0;
  const double n = (double)out.nUsed;
  const double kpar = (double)npar;
  out.aicc = out.chi2 + 2.0 * kpar + ((n - kpar - 1.0) > 0.0 ? 2.0 * kpar * (kpar + 1.0) / (n - kpar - 1.0) : 0.0);
  for (int i = 0; i < npar; ++i)
    if (!TMath::Finite(out.p[i]) || !(out.cov[i][i] >= 0.0)) { out.status = "non-finite parameters"; return false; }
  out.ok = true;
  out.status = "ok";
  return true;
}

// Y = sum_w [F - T B] over bins b0..b1, with the transfer-parameter covariance carried.
bool yieldWithTransfer(TH1* F, TH1* B, int b0, int b1, const Fit& f, double m0,
                       double& y, double& eStat, double& eTransfer) {
  y = 0.0; eStat = 0.0; eTransfer = 0.0;
  if (!F || !B || !f.ok) return false;
  double varStat = 0.0;
  double g[3] = {0.0, 0.0, 0.0};
  for (int i = b0; i <= b1; ++i) {
    const double m = F->GetXaxis()->GetBinCenter(i);
    const double fv = F->GetBinContent(i);
    const double bv = B->GetBinContent(i);
    const double ef = F->GetBinError(i);
    const double eb = B->GetBinError(i);
    const double T = transferAt(f, m, m0);
    y += fv - T * bv;                       // negative bins are kept, not clipped
    varStat += ef * ef + T * T * eb * eb;
    const double d = m - m0;
    g[0] += -T * bv;
    if (f.npar > 1) g[1] += -d * T * bv;
    if (f.npar > 2) g[2] += -d * d * T * bv;
  }
  double varPar = 0.0;
  for (int i = 0; i < f.npar; ++i)
    for (int j = 0; j < f.npar; ++j) varPar += g[i] * f.cov[i][j] * g[j];
  if (!(varStat >= 0.0) || !(varPar >= 0.0)) return false;
  eStat = TMath::Sqrt(varStat);
  eTransfer = TMath::Sqrt(varPar);
  return TMath::Finite(y) && TMath::Finite(eStat) && TMath::Finite(eTransfer);
}

// N0: the production nominal. One constant per k* bin from the alpha window, counted yield.
bool constantAlphaYield(TH1* F, TH1* B, int aLo, int aHi, int sLo, int sHi,
                        double& alpha, double& alphaErr, double& y, double& eStat, double& eAlpha) {
  alpha = 0.0; alphaErr = 0.0; y = 0.0; eStat = 0.0; eAlpha = 0.0;
  if (!F || !B) return false;
  double eFa = 0.0, eBa = 0.0, eFs = 0.0, eBs = 0.0;
  const double fa = F->IntegralAndError(aLo, aHi, eFa);
  const double ba = B->IntegralAndError(aLo, aHi, eBa);
  if (!(ba > 0.0) || !(fa > 0.0)) return false;
  alpha = fa / ba;
  alphaErr = alpha * TMath::Sqrt((eFa / fa) * (eFa / fa) + (eBa / ba) * (eBa / ba));
  const double fs = F->IntegralAndError(sLo, sHi, eFs);
  const double bs = B->IntegralAndError(sLo, sHi, eBs);
  y = fs - alpha * bs;
  eStat = TMath::Sqrt(eFs * eFs + alpha * alpha * eBs * eBs);
  eAlpha = alphaErr * bs;   // coherent: alpha moves the whole background level together
  return TMath::Finite(y);
}

// local slope of log[F/(T B)] on native mass bins, for the external validation windows
bool correctedSlope(TH1* F, TH1* B, double lo, double hi, const Fit& f, double m0,
                    double& slope, double& err, int& nUsed) {
  slope = 0.0; err = 0.0; nUsed = 0;
  if (!F || !B || !f.ok) return false;
  int b0 = 0, b1 = 0;
  if (!windowBins(F->GetXaxis(), lo, hi, b0, b1, "external")) return false;
  double sw = 0, swx = 0, swy = 0, swxx = 0, swxy = 0;
  for (int i = b0; i <= b1; ++i) {
    const double m = F->GetXaxis()->GetBinCenter(i);
    const double fv = F->GetBinContent(i);
    const double bv = B->GetBinContent(i);
    if (!(fv > 0.0) || !(bv > 0.0)) continue;
    const double ef = F->GetBinError(i);
    const double eb = B->GetBinError(i);
    const double var = (ef / fv) * (ef / fv) + (eb / bv) * (eb / bv);
    if (!(var > 0.0)) continue;
    const double T = transferAt(f, m, m0);
    const double yv = TMath::Log(fv / (T * bv));
    const double w = 1.0 / var;
    const double xx = m - m0;
    sw += w; swx += w * xx; swy += w * yv; swxx += w * xx * xx; swxy += w * xx * yv;
    ++nUsed;
  }
  if (nUsed < 3) return false;
  const double det = sw * swxx - swx * swx;
  if (!(TMath::Abs(det) > 0.0)) return false;
  slope = (sw * swxy - swx * swy) / det;
  err = TMath::Sqrt(sw / det);
  return TMath::Finite(slope) && TMath::Finite(err) && err > 0.0;
}

}  // namespace kmdt

namespace kmdt {

// ---------------------------------------------------------------------------------------------
// synthetic tests: these run before any real histogram is opened
// ---------------------------------------------------------------------------------------------

struct Toy {
  TH1D* F;
  TH1D* B;
  Toy() : F(0), B(0) {}
};

// Background falling smoothly across the mass range, transfer exp(a+b dm), optional peak.
Toy makeToy(TRandom3& rng, double a, double b, double m0, double sigAmp, double sigMean,
            double sigSigma, int nbins, double mLo, double mHi, double scale) {
  Toy t;
  t.F = new TH1D(uniq("toyF").c_str(), "toy F", nbins, mLo, mHi);
  t.B = new TH1D(uniq("toyB").c_str(), "toy B", nbins, mLo, mHi);
  t.F->SetDirectory(0);
  t.B->SetDirectory(0);
  t.F->Sumw2();
  t.B->Sumw2();
  for (int i = 1; i <= nbins; ++i) {
    const double m = t.B->GetBinCenter(i);
    const double shape = scale * TMath::Exp(-3.0 * (m - mLo));
    const double bTrue = shape;
    const double T = TMath::Exp(a + b * (m - m0));
    double fTrue = T * bTrue;
    if (sigAmp > 0.0) fTrue += sigAmp * TMath::Gaus(m, sigMean, sigSigma, kTRUE) * t.B->GetBinWidth(i);
    const double bObs = rng.Poisson(bTrue);
    const double fObs = rng.Poisson(fTrue);
    t.B->SetBinContent(i, bObs);
    t.B->SetBinError(i, TMath::Sqrt(TMath::Max(bObs, 1.0)));
    t.F->SetBinContent(i, fObs);
    t.F->SetBinError(i, TMath::Sqrt(TMath::Max(fObs, 1.0)));
  }
  return t;
}

struct SyntheticResult {
  bool pass;
  std::string detail;
  std::map<std::string, double> values;
  SyntheticResult() : pass(false) {}
};

SyntheticResult runSynthetic(int seed, int nToys, double toyTol, double toyBiasTol,
                             double ctrlLo, double ctrlHi, double sigLo, double sigHi) {
  SyntheticResult r;
  TRandom3 rng(seed);
  const double m0 = 0.5 * (sigLo + sigHi);
  const int nb = 100;
  const double mLo = 0.98, mHi = 1.18;
  std::string why;
  bool ok = true;

  // 1. constant truth: b consistent with zero
  {
    Toy t = makeToy(rng, TMath::Log(0.1), 0.0, m0, 0.0, 0.0, 0.0, nb, mLo, mHi, 2.0e5);
    Fit f0, f1;
    fitTransfer(t.F, t.B, ctrlLo, ctrlHi, 0, m0, 0, -1, f0);
    fitTransfer(t.F, t.B, ctrlLo, ctrlHi, 1, m0, 0, -1, f1);
    const double pull = (f1.cov[1][1] > 0.0) ? f1.p[1] / TMath::Sqrt(f1.cov[1][1]) : 1e9;
    r.values["const_b_pull"] = pull;
    r.values["const_a_recovered"] = TMath::Exp(f0.p[0]);
    if (!f0.ok || !f1.ok || TMath::Abs(pull) > 3.0) { ok = false; why += "constant-truth b not zero; "; }
    delete t.F; delete t.B;
  }

  // 2. log-linear truth: recover a and b, and the signal-window background
  {
    const double aT = TMath::Log(0.1), bT = 1.5;
    Toy t = makeToy(rng, aT, bT, m0, 0.0, 0.0, 0.0, nb, mLo, mHi, 2.0e5);
    Fit f1;
    fitTransfer(t.F, t.B, ctrlLo, ctrlHi, 1, m0, 0, -1, f1);
    const double pa = (f1.cov[0][0] > 0.0) ? (f1.p[0] - aT) / TMath::Sqrt(f1.cov[0][0]) : 1e9;
    const double pb = (f1.cov[1][1] > 0.0) ? (f1.p[1] - bT) / TMath::Sqrt(f1.cov[1][1]) : 1e9;
    r.values["loglin_a_pull"] = pa;
    r.values["loglin_b_pull"] = pb;
    int s0 = 0, s1 = 0;
    windowBins(t.F->GetXaxis(), sigLo, sigHi, s0, s1, "toy signal");
    double y = 0, es = 0, et = 0;
    yieldWithTransfer(t.F, t.B, s0, s1, f1, m0, y, es, et);
    const double pullY = (es > 0 || et > 0) ? y / TMath::Sqrt(es * es + et * et) : 1e9;
    r.values["loglin_zero_yield_pull"] = pullY;
    if (!f1.ok || TMath::Abs(pa) > 3.0 || TMath::Abs(pb) > 3.0 || TMath::Abs(pullY) > 3.0) {
      ok = false; why += "log-linear truth not recovered; ";
    }
    delete t.F; delete t.B;
  }

  // 3. signal injection: control fit unmoved, injected yield returned.
  //    Averaged over independent toys: one toy only says the answer is within its own error, and
  //    a few-percent bias would hide there.
  {
    const double aT = TMath::Log(0.1), bT = 1.0;
    const double inject = 5.0e4;
    const double injSigma = 0.004;
    const double injInWindow =
        inject * 0.5 * (TMath::Erf((sigHi - m0) / (TMath::Sqrt2() * injSigma)) -
                        TMath::Erf((sigLo - m0) / (TMath::Sqrt2() * injSigma)));
    const int nRep = 20;
    double sumPull = 0.0, sumRel = 0.0, sumBpull = 0.0;
    int nOk = 0;
    double lastY = 0.0;
    for (int rep = 0; rep < nRep; ++rep) {
      Toy t = makeToy(rng, aT, bT, m0, inject, m0, injSigma, nb, mLo, mHi, 2.0e5);
      Fit f1;
      if (fitTransfer(t.F, t.B, ctrlLo, ctrlHi, 1, m0, 0, -1, f1)) {
        int s0 = 0, s1 = 0;
        windowBins(t.F->GetXaxis(), sigLo, sigHi, s0, s1, "toy signal");
        double yv = 0, es = 0, et = 0;
        if (yieldWithTransfer(t.F, t.B, s0, s1, f1, m0, yv, es, et)) {
          const double sig = TMath::Sqrt(es * es + et * et);
          if (sig > 0.0) {
            sumPull += (yv - injInWindow) / sig;
            sumRel += (yv - injInWindow) / injInWindow;
            if (f1.cov[1][1] > 0.0) sumBpull += (f1.p[1] - bT) / TMath::Sqrt(f1.cov[1][1]);
            lastY = yv;
            ++nOk;
          }
        }
      }
      delete t.F; delete t.B;
    }
    const double meanPull = (nOk > 0) ? sumPull / nOk : 1e9;
    const double meanRel = (nOk > 0) ? sumRel / nOk : 1e9;
    const double meanBpull = (nOk > 0) ? sumBpull / nOk : 1e9;
    r.values["inject_expected_in_window"] = injInWindow;
    r.values["inject_yield_last"] = lastY;
    r.values["inject_mean_pull"] = meanPull;
    r.values["inject_mean_rel_bias"] = meanRel;
    r.values["inject_mean_b_pull"] = meanBpull;
    r.values["inject_toys_ok"] = (double)nOk;
    if (nOk < nRep || TMath::Abs(meanPull) > 0.5 || TMath::Abs(meanRel) > 0.02 ||
        TMath::Abs(meanBpull) > 0.5) {
      ok = false; why += "signal injection biased over repeated toys; ";
    }
  }

  // 4. uncertainty toy: analytic error against the spread of refluctuated yields
  {
    const double aT = TMath::Log(0.1), bT = 1.0;
    Toy t = makeToy(rng, aT, bT, m0, 0.0, 0.0, 0.0, nb, mLo, mHi, 2.0e5);
    Fit f1;
    fitTransfer(t.F, t.B, ctrlLo, ctrlHi, 1, m0, 0, -1, f1);
    int s0 = 0, s1 = 0;
    windowBins(t.F->GetXaxis(), sigLo, sigHi, s0, s1, "toy signal");
    double y0 = 0, es = 0, et = 0;
    yieldWithTransfer(t.F, t.B, s0, s1, f1, m0, y0, es, et);
    const double analytic = TMath::Sqrt(es * es + et * et);
    // Cholesky of the 2x2 parameter covariance
    const double l11 = TMath::Sqrt(f1.cov[0][0]);
    const double l21 = (l11 > 0.0) ? f1.cov[1][0] / l11 : 0.0;
    const double l22sq = f1.cov[1][1] - l21 * l21;
    const double l22 = (l22sq > 0.0) ? TMath::Sqrt(l22sq) : 0.0;
    double sum = 0.0, sum2 = 0.0;
    int n = 0;
    for (int it = 0; it < nToys; ++it) {
      const double z1 = rng.Gaus(), z2 = rng.Gaus();
      Fit fd = f1;
      fd.p[0] = f1.p[0] + l11 * z1;
      fd.p[1] = f1.p[1] + l21 * z1 + l22 * z2;
      double yy = 0.0;
      for (int i = s0; i <= s1; ++i) {
        const double m = t.F->GetXaxis()->GetBinCenter(i);
        const double fv = rng.Gaus(t.F->GetBinContent(i), t.F->GetBinError(i));
        const double bv = rng.Gaus(t.B->GetBinContent(i), t.B->GetBinError(i));
        yy += fv - transferAt(fd, m, m0) * bv;
      }
      sum += yy; sum2 += yy * yy; ++n;
    }
    const double mean = sum / n;
    const double rms = TMath::Sqrt(TMath::Max(sum2 / n - mean * mean, 0.0));
    const double relDiff = (analytic > 0.0) ? TMath::Abs(rms - analytic) / analytic : 1e9;
    const double bias = (analytic > 0.0) ? TMath::Abs(mean - y0) / analytic : 1e9;
    r.values["toy_analytic"] = analytic;
    r.values["toy_rms"] = rms;
    r.values["toy_rel_diff"] = relDiff;
    r.values["toy_mean_bias_sigma"] = bias;
    if (relDiff > toyTol || bias > toyBiasTol) { ok = false; why += "toy vs analytic error mismatch; "; }
    delete t.F; delete t.B;
  }

  // 5. failure modes stop rather than return a number
  {
    gQuiet = true;
    Toy t = makeToy(rng, TMath::Log(0.1), 0.0, m0, 0.0, 0.0, 0.0, 8, mLo, mHi, 1.0e3);
    Fit f;
    const bool tooFew = fitTransfer(t.F, t.B, ctrlLo, ctrlHi, 2, m0, 0, -1, f);
    r.values["failure_toofew_rejected"] = tooFew ? 0.0 : 1.0;
    if (tooFew) { ok = false; why += "a fit with too few bins was accepted; "; }
    delete t.F; delete t.B;
    gQuiet = false;
  }

  r.pass = ok;
  r.detail = ok ? "all synthetic tests pass" : why;
  return r;
}

}  // namespace kmdt

namespace kmdt {

struct Point {
  std::string channel, tmpl, evt;
  double kLo, kHi, kC;
  // N0
  double alpha, alphaErr, yN0, eN0stat, eN0alpha;
  // per model
  Fit fit[3];                       // index 0 = M0, 1 = M1, 2 = M2
  double y[3], eStat[3], eTrans[3];
  // cross-checks
  double yRangeA[3], yRangeB[3];    // fit-range cross-checks, same models
  double slopeRsb[2], slopeRsbErr[2], slopeLsb[2], slopeLsbErr[2];  // M0, M1 corrected
  // CV
  double cvChi2[3];
  int cvNdf[3];
  double cvMaxPull[3];
  int nNegSignalBins;
  bool valid;
  Point() : kLo(0), kHi(0), kC(0), alpha(0), alphaErr(0), yN0(0), eN0stat(0), eN0alpha(0),
            nNegSignalBins(0), valid(false) {
    for (int i = 0; i < 3; ++i) {
      y[i] = eStat[i] = eTrans[i] = 0.0;
      yRangeA[i] = yRangeB[i] = 0.0;
      cvChi2[i] = 0.0; cvNdf[i] = 0; cvMaxPull[i] = 0.0;
    }
    for (int i = 0; i < 2; ++i) { slopeRsb[i] = slopeRsbErr[i] = slopeLsb[i] = slopeLsbErr[i] = 0.0; }
  }
};

// blocked cross-validation over contiguous mass blocks of the control region
void blockedCV(TH1* F, TH1* B, double lo, double hi, int order, double m0, int nFolds,
               double& chi2, int& ndf, double& maxPull) {
  chi2 = 0.0; ndf = 0; maxPull = 0.0;
  int b0 = 0, b1 = 0;
  if (!windowBins(F->GetXaxis(), lo, hi, b0, b1, "control")) return;
  const int n = b1 - b0 + 1;
  if (nFolds < 2 || n < nFolds * 2) return;
  const int per = n / nFolds;
  for (int fold = 0; fold < nFolds; ++fold) {
    const int sLo = b0 + fold * per;
    const int sHi = (fold == nFolds - 1) ? b1 : (sLo + per - 1);
    Fit f;
    if (!fitTransfer(F, B, lo, hi, order, m0, sLo, sHi, f)) continue;
    double d = 0.0, es = 0.0, et = 0.0;
    if (!yieldWithTransfer(F, B, sLo, sHi, f, m0, d, es, et)) continue;
    const double sig = TMath::Sqrt(es * es + et * et);
    if (!(sig > 0.0)) continue;
    const double pull = d / sig;
    chi2 += pull * pull;
    ++ndf;
    if (TMath::Abs(pull) > maxPull) maxPull = TMath::Abs(pull);
  }
}

}  // namespace kmdt

void kmfMassDependentTransferClosure(const char* specText) {
  using namespace kmdt;
  gFatal = false;
  const Spec S = parseSpec(specText ? specText : "");

  const std::string input   = requireStr(S, "input");
  const std::string sidecar = requireStr(S, "sidecar");
  const std::string outRoot = requireStr(S, "outRoot");
  const std::string outPt   = requireStr(S, "outPointCsv");
  const std::string outSum  = requireStr(S, "outSummaryCsv");
  const std::string outPdf  = requireStr(S, "outPdf");
  const std::string prov    = requireStr(S, "provenance");
  const int centMin   = requireInt(S, "centMin");
  const int centMax   = requireInt(S, "centMax");
  const double dk     = requireNum(S, "dk");
  const int lowMerge  = requireInt(S, "lowMerge");
  const double kMin   = requireNum(S, "kMin");
  const double kMax   = requireNum(S, "kMax");
  const int massRebin = requireInt(S, "massRebin");
  const double sigLo  = requireNum(S, "sigLo");
  const double sigHi  = requireNum(S, "sigHi");
  const double aLo    = requireNum(S, "alphaLo");
  const double aHi    = requireNum(S, "alphaHi");
  const double cLo    = requireNum(S, "ctrlLo");
  const double cHi    = requireNum(S, "ctrlHi");
  const double cALo   = requireNum(S, "ctrlALo");
  const double cAHi   = requireNum(S, "ctrlAHi");
  const double cBLo   = requireNum(S, "ctrlBLo");
  const double cBHi   = requireNum(S, "ctrlBHi");
  const double rLo    = requireNum(S, "rsbLo");
  const double rHi    = requireNum(S, "rsbHi");
  const double lLo    = requireNum(S, "lsbLo");
  const double lHi    = requireNum(S, "lsbHi");
  const double nLo    = requireNum(S, "normLo");
  const double nHi    = requireNum(S, "normHi");
  const int nFolds    = requireInt(S, "nFolds");
  const int nToys     = requireInt(S, "nToys");
  const int seed      = requireInt(S, "seed");
  const double toyTol = requireNum(S, "toyTol");
  const double toyBias= requireNum(S, "toyBiasTol");
  const int synthOnly = requireInt(S, "syntheticOnly");
  if (gFatal) { std::cerr << "[kmdt] stopping: the specification is incomplete" << std::endl; return; }

  // window sanity: the control regions must not touch the signal window
  if (!(sigHi <= cLo) || !(sigHi <= cALo) || !(sigHi <= cBLo)) {
    fail("a control window overlaps or precedes the signal window");
  }
  if (!(sigLo < sigHi) || !(cLo < cHi) || !(aLo < aHi) || !(nLo < nHi)) fail("a window has hi <= lo");
  if (gFatal) return;

  const double m0 = 0.5 * (sigLo + sigHi);
  std::cout << "[kmdt] provenance: " << prov << std::endl;
  std::cout << "[kmdt] m0 = " << m0 << " (from the signal window, not hardcoded)" << std::endl;

  // ---- synthetic gate first, before any real histogram is opened ----
  SyntheticResult syn = runSynthetic(seed, nToys, toyTol, toyBias, cLo, cHi, sigLo, sigHi);
  std::cout << "[kmdt] synthetic: " << (syn.pass ? "PASS" : "FAIL") << " -- " << syn.detail << std::endl;
  for (std::map<std::string, double>::const_iterator it = syn.values.begin(); it != syn.values.end(); ++it)
    std::cout << "[kmdt]   " << it->first << " = " << it->second << std::endl;
  if (synthOnly != 0) {
    TFile* fo = TFile::Open(outRoot.c_str(), "RECREATE");
    if (fo && !fo->IsZombie()) {
      TNamed("meta_provenance", prov.c_str()).Write();
      TNamed("meta_mode", "syntheticOnly").Write();
      TNamed("synthetic_detail", syn.detail.c_str()).Write();
      TParameter<double>("synthetic_pass", syn.pass ? 1.0 : 0.0).Write();
      for (std::map<std::string, double>::const_iterator it = syn.values.begin(); it != syn.values.end(); ++it)
        TParameter<double>(("synthetic_" + it->first).c_str(), it->second).Write();
      fo->Close();
    }
    std::cout << "[kmdt] synthetic-only run finished" << std::endl;
    return;
  }
  if (!syn.pass) {
    std::cerr << "[kmdt] synthetic tests failed; refusing to report real-data results" << std::endl;
    return;
  }

  TFile* fin = TFile::Open(input.c_str());
  if (!fin || fin->IsZombie()) { fail("cannot open " + input); return; }

  const char* channels[2] = {"phi_proton", "phi_deuteron"};
  const char* bach[2]     = {"proton", "deuteron"};
  const char* tmpls[2]    = {"rot", "mix"};
  const char* evts[2]     = {"SE", "ME"};

  std::vector<Point> points;
  for (int ic = 0; ic < 2; ++ic) {
    for (int it = 0; it < 2; ++it) {
      for (int ie = 0; ie < 2; ++ie) {
        const std::string fname = Form("hPhiMKK_vs_Kstar%s_%s_wide", evts[ie], channels[ic]);
        const std::string bname = Form("hPhiMKK_vs_Kstar%s_phi_%s_%s_wide", evts[ie], tmpls[it], bach[ic]);
        TH3* h3F = (TH3*)fin->Get(fname.c_str());
        TH3* h3B = (TH3*)fin->Get(bname.c_str());
        if (!h3F || !h3B) { fail("missing histogram " + (h3F ? bname : fname)); continue; }
        TH2* h2F = projectMkkKstar(h3F, centMin, centMax, dk, "F");
        TH2* h2B = projectMkkKstar(h3B, centMin, centMax, dk, "B");
        if (!h2F || !h2B) { fail("projection failed for " + fname); continue; }

        for (int iy = 1; iy <= h2F->GetNbinsY(); ++iy) {
          const int iyFirst = iy;
          const int iyLast = (iyFirst == 1) ? TMath::Min(h2F->GetNbinsY(), lowMerge) : iyFirst;
          iy = iyLast;
          const double kLo = h2F->GetYaxis()->GetBinLowEdge(iyFirst);
          const double kHi = h2F->GetYaxis()->GetBinUpEdge(iyLast);
          const double kC = 0.5 * (kLo + kHi);
          if (kC < kMin || kC > kMax) continue;

          Point P;
          P.channel = channels[ic]; P.tmpl = tmpls[it]; P.evt = evts[ie];
          P.kLo = kLo; P.kHi = kHi; P.kC = kC;

          TH1* F = projectMass(h2F, iyFirst, iyLast, massRebin, "mF");
          TH1* B = projectMass(h2B, iyFirst, iyLast, massRebin, "mB");
          TH1* Fn = projectMass(h2F, iyFirst, iyLast, 1, "nF");   // native bins for external slopes
          TH1* Bn = projectMass(h2B, iyFirst, iyLast, 1, "nB");
          if (!F || !B || !Fn || !Bn) { fail("mass projection failed"); continue; }

          int sb0 = 0, sb1 = 0, ab0 = 0, ab1 = 0;
          if (windowBins(F->GetXaxis(), sigLo, sigHi, sb0, sb1, "signal") &&
              windowBins(F->GetXaxis(), aLo, aHi, ab0, ab1, "alpha")) {
            constantAlphaYield(F, B, ab0, ab1, sb0, sb1, P.alpha, P.alphaErr, P.yN0, P.eN0stat, P.eN0alpha);
            for (int i = sb0; i <= sb1; ++i) if (F->GetBinContent(i) - P.alpha * B->GetBinContent(i) < 0.0) ++P.nNegSignalBins;
            for (int order = 0; order <= 2; ++order) {
              Fit f;
              if (fitTransfer(F, B, cLo, cHi, order, m0, 0, -1, f)) {
                P.fit[order] = f;
                yieldWithTransfer(F, B, sb0, sb1, f, m0, P.y[order], P.eStat[order], P.eTrans[order]);
              }
              Fit fa, fb;
              double dum1 = 0, dum2 = 0;
              if (fitTransfer(F, B, cALo, cAHi, order, m0, 0, -1, fa))
                yieldWithTransfer(F, B, sb0, sb1, fa, m0, P.yRangeA[order], dum1, dum2);
              if (fitTransfer(F, B, cBLo, cBHi, order, m0, 0, -1, fb))
                yieldWithTransfer(F, B, sb0, sb1, fb, m0, P.yRangeB[order], dum1, dum2);
              blockedCV(F, B, cLo, cHi, order, m0, nFolds, P.cvChi2[order], P.cvNdf[order], P.cvMaxPull[order]);
            }
            for (int m = 0; m < 2; ++m) {
              int nu = 0;
              correctedSlope(Fn, Bn, rLo, rHi, P.fit[m], m0, P.slopeRsb[m], P.slopeRsbErr[m], nu);
              correctedSlope(Fn, Bn, lLo, lHi, P.fit[m], m0, P.slopeLsb[m], P.slopeLsbErr[m], nu);
            }
            P.valid = P.fit[0].ok && P.fit[1].ok;
          }
          points.push_back(P);
          delete F; delete B; delete Fn; delete Bn;
        }
        delete h2F; delete h2B;
      }
    }
  }
  std::cout << "[kmdt] points: " << points.size() << std::endl;
  if (gFatal) {
    std::cerr << "[kmdt] a histogram or window failed while reading the real data; "
              << "refusing to report gates on an incomplete set" << std::endl;
    fin->Close();
    return;
  }

  // ---- correlation functions per model, normalised over the same k* range as the analysis ----
  // key: channel|template|model -> vectors over k*
  std::map<std::string, std::vector<double> > cfK, cfRaw, cfRawErr, cfNorm, cfNormErr;
  for (int ic = 0; ic < 2; ++ic) {
    for (int it = 0; it < 2; ++it) {
      for (int model = -1; model <= 2; ++model) {   // -1 = N0
        const std::string key = std::string(channels[ic]) + "|" + tmpls[it] + "|" +
                                (model < 0 ? "N0" : Form("M%d", model));
        std::vector<double> kk, rr, ee;
        for (size_t i = 0; i < points.size(); ++i) {
          const Point& p = points[i];
          if (p.channel != channels[ic] || p.tmpl != tmpls[it] || p.evt != "SE") continue;
          const Point* q = 0;
          for (size_t j = 0; j < points.size(); ++j)
            if (points[j].channel == p.channel && points[j].tmpl == p.tmpl && points[j].evt == "ME" &&
                TMath::Abs(points[j].kC - p.kC) < 1e-9) { q = &points[j]; break; }
          if (!q) continue;
          double ySE, eSE, yME, eME;
          if (model < 0) {
            ySE = p.yN0; eSE = TMath::Sqrt(p.eN0stat * p.eN0stat + p.eN0alpha * p.eN0alpha);
            yME = q->yN0; eME = TMath::Sqrt(q->eN0stat * q->eN0stat + q->eN0alpha * q->eN0alpha);
          } else {
            if (!p.fit[model].ok || !q->fit[model].ok) continue;
            ySE = p.y[model]; eSE = TMath::Sqrt(p.eStat[model] * p.eStat[model] + p.eTrans[model] * p.eTrans[model]);
            yME = q->y[model]; eME = TMath::Sqrt(q->eStat[model] * q->eStat[model] + q->eTrans[model] * q->eTrans[model]);
          }
          if (!(yME > 0.0) || !TMath::Finite(ySE) || !TMath::Finite(yME)) continue;
          const double c = ySE / yME;
          const double ec = TMath::Abs(c) * TMath::Sqrt((eSE / (ySE + 1e-12)) * (eSE / (ySE + 1e-12)) +
                                                        (eME / (yME + 1e-12)) * (eME / (yME + 1e-12)));
          kk.push_back(p.kC); rr.push_back(c); ee.push_back(ec);
        }
        double sum = 0.0; int nUsed = 0;
        for (size_t i = 0; i < kk.size(); ++i)
          if (kk[i] >= nLo && kk[i] <= nHi) { sum += rr[i]; ++nUsed; }
        const double scale = (nUsed > 0 && sum > 0.0) ? ((double)nUsed / sum) : 0.0;
        std::vector<double> nn, ne;
        for (size_t i = 0; i < kk.size(); ++i) { nn.push_back(rr[i] * scale); ne.push_back(ee[i] * scale); }
        cfK[key] = kk; cfRaw[key] = rr; cfRawErr[key] = ee; cfNorm[key] = nn; cfNormErr[key] = ne;
        if (nUsed == 0) std::cerr << "[kmdt] WARNING: no normalisation points for " << key << std::endl;
      }
    }
  }

  // ---- gates ----
  // G-A technical integrity
  int nExpected = 8 * 12;
  int nFinite = 0, nInvalid = 0;
  for (size_t i = 0; i < points.size(); ++i) {
    const Point& p = points[i];
    if (!p.valid) { ++nInvalid; continue; }
    bool fin = TMath::Finite(p.y[1]) && TMath::Finite(p.eStat[1]) && TMath::Finite(p.eTrans[1]);
    if (fin) ++nFinite;
  }
  // N0 against the production sidecar
  double maxN0Diff = -1.0;
  TFile* fsc = TFile::Open(sidecar.c_str());
  if (!fsc || fsc->IsZombie()) {
    fail("cannot open the sidecar " + sidecar + " (needed for the N0 regression in G-A)");
  } else {
    for (int ic = 0; ic < 2; ++ic) {
      for (int ie = 0; ie < 2; ++ie) {
        const std::string gn = Form("slice_pct_0_60_kmf_Y_%s_rot_%s_dk%.0f", evts[ie], channels[ic], dk * 1000.0);
        TGraphErrors* g = (TGraphErrors*)fsc->Get(gn.c_str());
        if (!g) { std::cerr << "[kmdt] WARNING: sidecar graph absent: " << gn << std::endl; continue; }
        for (size_t i = 0; i < points.size(); ++i) {
          const Point& p = points[i];
          if (p.channel != channels[ic] || p.tmpl != "rot" || p.evt != evts[ie]) continue;
          for (int q = 0; q < g->GetN(); ++q) {
            double x, y; g->GetPoint(q, x, y);
            if (TMath::Abs(x - p.kC) > 1e-6) continue;
            const double rel = (TMath::Abs(y) > 0.0) ? TMath::Abs(p.yN0 - y) / TMath::Abs(y) : TMath::Abs(p.yN0 - y);
            if (rel > maxN0Diff) maxN0Diff = rel;
          }
        }
      }
    }
  }
  const bool gA = (nInvalid == 0) && (nFinite == nExpected) && (maxN0Diff >= 0.0 && maxN0Diff < 1e-12) && syn.pass;

  // G-B closure improvement, ROT only
  double cvM0 = 0.0, cvM1 = 0.0; int ndfM0 = 0, ndfM1 = 0;
  double worstPull = 0.0, worstProb = 1.0, worstRsbSigma = 0.0;
  bool anyRotWorse = false;
  for (size_t i = 0; i < points.size(); ++i) {
    const Point& p = points[i];
    if (p.tmpl != "rot" || !p.valid) continue;
    if (p.cvNdf[0] != p.cvNdf[1]) continue;           // compare on the same block set only
    cvM0 += p.cvChi2[0]; ndfM0 += p.cvNdf[0];
    cvM1 += p.cvChi2[1]; ndfM1 += p.cvNdf[1];
    if (p.cvMaxPull[1] > worstPull) worstPull = p.cvMaxPull[1];
    if (p.cvChi2[1] > p.cvChi2[0]) anyRotWorse = true;
    if (p.slopeRsbErr[1] > 0.0) {
      const double s = TMath::Abs(p.slopeRsb[1] / p.slopeRsbErr[1]);
      if (s > worstRsbSigma) worstRsbSigma = s;
    }
  }
  const double probM1 = (ndfM1 > 0) ? TMath::Prob(cvM1, ndfM1) : -1.0;
  const bool gB = (ndfM1 > 0) && (cvM1 < cvM0) && (probM1 > 0.01) && (worstPull < 3.0) &&
                  (worstRsbSigma < 2.0) && !anyRotWorse;

  // G-C extrapolation stability below 0.4 GeV/c, ROT
  double medDiff = 0.0, maxDiff = 0.0, maxRangeDiff = 0.0;
  {
    std::vector<double> diffs;
    for (int ic = 0; ic < 2; ++ic) {
      const std::string k1 = std::string(channels[ic]) + "|rot|M1";
      const std::string k2 = std::string(channels[ic]) + "|rot|M2";
      if (!cfNorm.count(k1) || !cfNorm.count(k2)) continue;
      for (size_t i = 0; i < cfNorm[k1].size() && i < cfNorm[k2].size(); ++i) {
        if (cfK[k1][i] >= 0.4) continue;
        const double sig = cfNormErr[k1][i];
        if (!(sig > 0.0)) continue;
        const double d = TMath::Abs(cfNorm[k1][i] - cfNorm[k2][i]) / sig;
        diffs.push_back(d);
        if (d > maxDiff) maxDiff = d;
      }
    }
    if (!diffs.empty()) {
      std::sort(diffs.begin(), diffs.end());
      medDiff = diffs[diffs.size() / 2];
    }
    for (size_t i = 0; i < points.size(); ++i) {
      const Point& p = points[i];
      if (p.tmpl != "rot" || !p.valid || p.kC >= 0.4) continue;
      const double s = TMath::Sqrt(p.eStat[1] * p.eStat[1] + p.eTrans[1] * p.eTrans[1]);
      if (!(s > 0.0)) continue;
      const double dA = TMath::Abs(p.y[1] - p.yRangeA[1]) / s;
      const double dB = TMath::Abs(p.y[1] - p.yRangeB[1]) / s;
      if (dA > maxRangeDiff) maxRangeDiff = dA;
      if (dB > maxRangeDiff) maxRangeDiff = dB;
    }
  }
  const bool gC = (medDiff <= 0.25) && (maxDiff <= 1.0) && (maxRangeDiff <= 1.0);

  // G-D template robustness
  double rmsN0 = 0.0, rmsM1 = 0.0; int nRms = 0;
  for (int ic = 0; ic < 2; ++ic) {
    const std::string r0 = std::string(channels[ic]) + "|rot|N0";
    const std::string m0k = std::string(channels[ic]) + "|mix|N0";
    const std::string r1 = std::string(channels[ic]) + "|rot|M1";
    const std::string m1k = std::string(channels[ic]) + "|mix|M1";
    if (!cfNorm.count(r0) || !cfNorm.count(m0k) || !cfNorm.count(r1) || !cfNorm.count(m1k)) continue;
    for (size_t i = 0; i < cfNorm[r0].size(); ++i) {
      if (i >= cfNorm[m0k].size() || i >= cfNorm[r1].size() || i >= cfNorm[m1k].size()) break;
      if (cfK[r0][i] >= 0.4) continue;
      const double d0 = cfNorm[r0][i] - cfNorm[m0k][i];
      const double d1 = cfNorm[r1][i] - cfNorm[m1k][i];
      rmsN0 += d0 * d0; rmsM1 += d1 * d1; ++nRms;
    }
  }
  if (nRms > 0) { rmsN0 = TMath::Sqrt(rmsN0 / nRms); rmsM1 = TMath::Sqrt(rmsM1 / nRms); }
  const bool gD = (nRms > 0) && (rmsM1 <= rmsN0);

  std::string verdict = "INCONCLUSIVE";
  if (gA && gB && gC && gD) verdict = "PASS";
  else if (gA && (!gB || !gC || !gD)) verdict = "FAIL";

  std::cout << "[kmdt] G-A " << (gA ? "PASS" : "FAIL")
            << "  (points " << nFinite << "/" << nExpected << ", invalid " << nInvalid
            << ", maxN0RelDiff " << maxN0Diff << ")" << std::endl;
  std::cout << "[kmdt] G-B " << (gB ? "PASS" : "FAIL")
            << "  (CV chi2 M0 " << cvM0 << "/" << ndfM0 << " vs M1 " << cvM1 << "/" << ndfM1
            << ", p " << probM1 << ", maxPull " << worstPull << ", worstRSBsigma " << worstRsbSigma
            << ", anyRotWorse " << (anyRotWorse ? "yes" : "no") << ")" << std::endl;
  std::cout << "[kmdt] G-C " << (gC ? "PASS" : "FAIL")
            << "  (median |M1-M2| " << medDiff << " sigma, max " << maxDiff
            << ", max fit-range " << maxRangeDiff << ")" << std::endl;
  std::cout << "[kmdt] G-D " << (gD ? "PASS" : "FAIL")
            << "  (ROT-MIX RMS N0 " << rmsN0 << " -> M1 " << rmsM1 << ")" << std::endl;
  std::cout << "[kmdt] VERDICT: " << verdict << std::endl;

  // ---- outputs ----
  {
    std::ofstream csv(outPt.c_str());
    csv << "channel,template,evtclass,k_lo,k_hi,k_center,model,alpha,alpha_err,"
        << "a,b,c,cov_aa,cov_ab,cov_ac,cov_bb,cov_bc,cov_cc,chi2,ndf,prob,aicc,fit_status,"
        << "n_used_bins,n_nonpositive_bins,cv_chi2,cv_ndf,cv_maxpull,"
        << "yield,yield_stat_err,yield_transfer_err,yield_total_err,"
        << "yield_fitrangeA,yield_fitrangeB,slope_rsb,slope_rsb_err,slope_lsb,slope_lsb_err,"
        << "n_negative_signal_bins,valid\n";
    for (size_t i = 0; i < points.size(); ++i) {
      const Point& p = points[i];
      for (int model = -1; model <= 2; ++model) {
        const Fit& f = (model >= 0) ? p.fit[model] : p.fit[0];
        const double yy = (model < 0) ? p.yN0 : p.y[model];
        const double es = (model < 0) ? p.eN0stat : p.eStat[model];
        const double et = (model < 0) ? p.eN0alpha : p.eTrans[model];
        csv << p.channel << "," << p.tmpl << "," << p.evt << ","
            << p.kLo << "," << p.kHi << "," << p.kC << ","
            << (model < 0 ? "N0" : Form("M%d", model)) << ","
            << p.alpha << "," << p.alphaErr << ",";
        if (model < 0) csv << ",,,,,,,,,,,,,,,,";
        else {
          csv << f.p[0] << "," << (f.npar > 1 ? f.p[1] : 0.0) << "," << (f.npar > 2 ? f.p[2] : 0.0) << ","
              << f.cov[0][0] << "," << f.cov[0][1] << "," << f.cov[0][2] << ","
              << f.cov[1][1] << "," << f.cov[1][2] << "," << f.cov[2][2] << ","
              << f.chi2 << "," << f.ndf << "," << f.prob << "," << f.aicc << "," << f.status << ","
              << f.nUsed << "," << f.nNonPositive << ","
              << p.cvChi2[model] << "," << p.cvNdf[model] << "," << p.cvMaxPull[model] << ",";
        }
        csv << yy << "," << es << "," << et << "," << TMath::Sqrt(es * es + et * et) << ","
            << (model < 0 ? 0.0 : p.yRangeA[model]) << "," << (model < 0 ? 0.0 : p.yRangeB[model]) << ","
            << (model >= 0 && model < 2 ? p.slopeRsb[model] : 0.0) << ","
            << (model >= 0 && model < 2 ? p.slopeRsbErr[model] : 0.0) << ","
            << (model >= 0 && model < 2 ? p.slopeLsb[model] : 0.0) << ","
            << (model >= 0 && model < 2 ? p.slopeLsbErr[model] : 0.0) << ","
            << p.nNegSignalBins << "," << (p.valid ? 1 : 0) << "\n";
      }
    }
    csv.close();
    std::cout << "[kmdt] wrote " << outPt << std::endl;
  }
  {
    std::ofstream csv(outSum.c_str());
    csv << "quantity,value\n";
    csv << "points_expected," << nExpected << "\n";
    csv << "points_finite," << nFinite << "\n";
    csv << "points_invalid," << nInvalid << "\n";
    csv << "max_N0_rel_diff_vs_sidecar," << maxN0Diff << "\n";
    csv << "cv_chi2_M0," << cvM0 << "\ncv_ndf_M0," << ndfM0 << "\n";
    csv << "cv_chi2_M1," << cvM1 << "\ncv_ndf_M1," << ndfM1 << "\n";
    csv << "cv_prob_M1," << probM1 << "\n";
    csv << "cv_max_pull_M1," << worstPull << "\n";
    csv << "worst_rsb_slope_sigma_M1," << worstRsbSigma << "\n";
    csv << "any_rot_series_worse," << (anyRotWorse ? 1 : 0) << "\n";
    csv << "median_M1_M2_sigma," << medDiff << "\nmax_M1_M2_sigma," << maxDiff << "\n";
    csv << "max_fitrange_sigma," << maxRangeDiff << "\n";
    csv << "rot_mix_rms_N0," << rmsN0 << "\nrot_mix_rms_M1," << rmsM1 << "\n";
    csv << "gate_A," << (gA ? "PASS" : "FAIL") << "\n";
    csv << "gate_B," << (gB ? "PASS" : "FAIL") << "\n";
    csv << "gate_C," << (gC ? "PASS" : "FAIL") << "\n";
    csv << "gate_D," << (gD ? "PASS" : "FAIL") << "\n";
    csv << "verdict," << verdict << "\n";
    csv.close();
    std::cout << "[kmdt] wrote " << outSum << std::endl;
  }

  TFile* fo = TFile::Open(outRoot.c_str(), "RECREATE");
  if (fo && !fo->IsZombie()) {
    TNamed("meta_provenance", prov.c_str()).Write();
    TNamed("meta_input", input.c_str()).Write();
    TNamed("meta_sidecar", sidecar.c_str()).Write();
    TNamed("meta_spec", specText).Write();
    TNamed("meta_error_model",
           "diagonal bin errors; F and the ROT/MIX templates share candidates and events, so the "
           "histogram-to-histogram covariance is not stored and the p-values here are diagnostics "
           "on a diagonal-error model").Write();
    TNamed("meta_scope", "Phase A diagnostic; the nominal is unchanged and these errors are not a "
           "final statistical covariance matrix (Step 10 T8)").Write();
    TNamed("verdict", verdict.c_str()).Write();
    TNamed("synthetic_detail", syn.detail.c_str()).Write();
    TParameter<double>("gate_A", gA ? 1.0 : 0.0).Write();
    TParameter<double>("gate_B", gB ? 1.0 : 0.0).Write();
    TParameter<double>("gate_C", gC ? 1.0 : 0.0).Write();
    TParameter<double>("gate_D", gD ? 1.0 : 0.0).Write();
    TParameter<double>("cv_chi2_M0", cvM0).Write();
    TParameter<double>("cv_chi2_M1", cvM1).Write();
    TParameter<double>("cv_ndf_M1", ndfM1).Write();
    TParameter<double>("cv_prob_M1", probM1).Write();
    TParameter<double>("median_M1_M2_sigma", medDiff).Write();
    TParameter<double>("max_M1_M2_sigma", maxDiff).Write();
    TParameter<double>("max_fitrange_sigma", maxRangeDiff).Write();
    TParameter<double>("rot_mix_rms_N0", rmsN0).Write();
    TParameter<double>("rot_mix_rms_M1", rmsM1).Write();
    TParameter<double>("max_N0_rel_diff_vs_sidecar", maxN0Diff).Write();
    for (std::map<std::string, double>::const_iterator it = syn.values.begin(); it != syn.values.end(); ++it)
      TParameter<double>(("synthetic_" + it->first).c_str(), it->second).Write();
    // every required argument, so the run can be reproduced from the file alone
    for (Spec::const_iterator it = S.begin(); it != S.end(); ++it)
      TNamed(("arg_" + it->first).c_str(), it->second.c_str()).Write();
    for (std::map<std::string, std::vector<double> >::const_iterator it = cfNorm.begin();
         it != cfNorm.end(); ++it) {
      const std::string& key = it->first;
      std::string safe = key;
      for (size_t i = 0; i < safe.size(); ++i) if (safe[i] == '|') safe[i] = '_';
      if (cfK[key].empty()) continue;
      TGraphErrors g((int)cfK[key].size(), &cfK[key][0], &cfNorm[key][0], 0, &cfNormErr[key][0]);
      g.SetName(("CF_norm_" + safe).c_str());
      g.SetTitle(("normalised CF " + key).c_str());
      g.Write();
      TGraphErrors gr((int)cfK[key].size(), &cfK[key][0], &cfRaw[key][0], 0, &cfRawErr[key][0]);
      gr.SetName(("CF_raw_" + safe).c_str());
      gr.Write();
    }
    fo->Close();
    std::cout << "[kmdt] wrote " << outRoot << std::endl;
  }

  // PDF: guide, fit parameters, corrected slopes, model comparison, gate summary
  {
    TCanvas c("kmdtCanvas", "kmdt", 1400, 900);
    c.Print((outPdf + "[").c_str());
    TLatex tx;
    tx.SetNDC(kTRUE);
    tx.SetTextSize(0.024);
    c.Clear();
    double y = 0.94;
    tx.DrawLatex(0.05, y, "Step 10 T5c Phase A -- mass-dependent background transfer"); y -= 0.05;
    tx.DrawLatex(0.05, y, Form("provenance: %s", prov.c_str())); y -= 0.035;
    tx.DrawLatex(0.05, y, Form("input: %s", input.c_str())); y -= 0.035;
    tx.DrawLatex(0.05, y, "models: N0 single-window alpha, M0 log T=a, M1 log T=a+b(m-m0), M2 adds c(m-m0)^2"); y -= 0.035;
    tx.DrawLatex(0.05, y, Form("control %.3f-%.3f, signal %.3f-%.3f, alpha %.3f-%.3f", cLo, cHi, sigLo, sigHi, aLo, aHi)); y -= 0.035;
    tx.DrawLatex(0.05, y, "both SE and ME are corrected; negative bins are not clipped"); y -= 0.035;
    tx.DrawLatex(0.05, y, "errors are diagonal: F and the templates share candidates, and that covariance is not stored"); y -= 0.05;
    tx.DrawLatex(0.05, y, Form("G-A %s   G-B %s   G-C %s   G-D %s", gA ? "PASS" : "FAIL", gB ? "PASS" : "FAIL",
                               gC ? "PASS" : "FAIL", gD ? "PASS" : "FAIL")); y -= 0.04;
    tx.DrawLatex(0.05, y, Form("VERDICT: %s", verdict.c_str())); y -= 0.05;
    tx.DrawLatex(0.05, y, Form("CV chi2 M0 %.1f/%d, M1 %.1f/%d, p(M1) %.3f, max |pull| %.2f",
                               cvM0, ndfM0, cvM1, ndfM1, probM1, worstPull)); y -= 0.035;
    tx.DrawLatex(0.05, y, Form("M1 vs M2 below 0.4 GeV/c: median %.2f sigma, max %.2f sigma", medDiff, maxDiff)); y -= 0.035;
    tx.DrawLatex(0.05, y, Form("ROT-MIX RMS: N0 %.4f -> M1 %.4f", rmsN0, rmsM1)); y -= 0.035;
    tx.DrawLatex(0.05, y, Form("synthetic: %s", syn.detail.c_str()));
    c.Print(outPdf.c_str());

    for (int ic = 0; ic < 2; ++ic) {
      std::vector<double> kx, ba, bb;
      for (size_t i = 0; i < points.size(); ++i) {
        const Point& p = points[i];
        if (p.channel != channels[ic] || p.tmpl != "rot" || p.evt != "SE" || !p.fit[1].ok) continue;
        kx.push_back(p.kC); ba.push_back(p.fit[1].p[1]);
        bb.push_back(p.fit[1].cov[1][1] > 0 ? TMath::Sqrt(p.fit[1].cov[1][1]) : 0.0);
      }
      if (kx.empty()) continue;
      c.Clear();
      TGraphErrors g((int)kx.size(), &kx[0], &ba[0], 0, &bb[0]);
      g.SetTitle(Form("M1 slope b vs k*  %s ROT SE;k* (GeV/c);b (GeV/c^{2})^{-1}", channels[ic]));
      g.SetMarkerStyle(20);
      g.Draw("AP");
      c.Print(outPdf.c_str());
    }
    for (int ic = 0; ic < 2; ++ic) {
      const std::string kn0 = std::string(channels[ic]) + "|rot|N0";
      const std::string km1 = std::string(channels[ic]) + "|rot|M1";
      if (cfK[kn0].empty() || cfK[km1].empty()) continue;
      c.Clear();
      TGraphErrors g0((int)cfK[kn0].size(), &cfK[kn0][0], &cfNorm[kn0][0], 0, &cfNormErr[kn0][0]);
      TGraphErrors g1((int)cfK[km1].size(), &cfK[km1][0], &cfNorm[km1][0], 0, &cfNormErr[km1][0]);
      g0.SetTitle(Form("normalised CF  %s ROT;k* (GeV/c);C(k*)", channels[ic]));
      g0.SetMarkerStyle(20); g1.SetMarkerStyle(24); g1.SetMarkerColor(kRed + 1); g1.SetLineColor(kRed + 1);
      g0.Draw("AP"); g1.Draw("P SAME");
      c.Print(outPdf.c_str());
    }
    c.Print((outPdf + "]").c_str());
    std::cout << "[kmdt] wrote " << outPdf << std::endl;
  }
  fin->Close();
}
