// purityCorrectedCorrelation.C
//
// Correlation function of one StCandidateFemtoMaker channel (species A - species B) with the
// k*-differential purity correction, no sidebands:
//   CF(k*) = N_SE(k*) P_A^SE(k*) P_B^SE(k*) / [N_ME(k*) P_A^ME(k*) P_B^ME(k*)], normalised to 1
// in [cfNormMin, cfNormMax] (candidateFemto YAML). N = pairs with BOTH masses in the signal windows.
// P_A(k*) = purity of the A signal window from a fit of the A mass variable for pairs whose B is in
// the B signal window (same k* bin); P_B(k*) likewise with the roles swapped. The two purities
// are treated as independent (naive product).
// Per k* bin the peak shape (mean, width and, for Student-t / two Gaussians, the extra shape parameters) is fixed to
// the k*-integrated fit of the same sample (SE or ME); only the peak amplitude and the background float (binned likelihood). A species with
// purityModel tofMassSquaredKaonPion (TOF m^2; its peak width depends on p and so on k*) floats
// mean/sigma in every k*-resolved fit (per bin, sliding window, merged bins), started from the
// k*-integrated fit of the same sample; such a fit counts as valid only if it converged and its width
// is not at a limit. A species with purityByDefinition true gets purity 1 (its fits are QA only).
// Analysis settings come from the candidateFemto YAML: purityFitRange, purityFitSeed (per species),
// cfNormMin/Max, purityFitMinEntries, slidingWindowHalfWidth, purityCorrectionMaxMergedBins.
// Also: constant (k*-integrated) purity, the sliding-window purity correction
// and the purity correction with merged k* bins where a fit fails (the final result).
// Input: the THnSparse of StCandidateFemtoMaker, projected with CandidateFemtoProjection.h; configuration:
// CandidateFemtoConfig.h.
//
// Usage (via script/singularity_purityCorrectedCorrelation.sh):
//   purityCorrectedCorrelation("file.root", "mainconf.yaml", "k0s_xi_minus", "out/prefix", kstarRebin,
//                              "4,5,6,7,8", "")
//   cent9List "" = all centralities. BDT thresholds: species_<s>_scoreMin of the YAML (optional).
//   options (';' separated, for scans; a chosen value belongs in the YAML): purityA=config|fit|unit;
//   windowA=lo,hi (overrides the A signal window); scoreMinA=x, scoreMinB=x (override the BDT thresholds);
//   slideHalfWidth=N (overrides slidingWindowHalfWidth of the YAML).
#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH3.h>
#include <TF1.h>
#include <TVirtualFitter.h>
#include <TLine.h>
#include <TLegend.h>
#include <TMath.h>
#include <TString.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TPaveText.h>
#include <TDatabasePDG.h>
#include <TParticlePDG.h>
#include <iostream>
#include <map>
#include <vector>
#include "CandidateFemtoConfig.h"
#include "CandidateFemtoProjection.h"

namespace {

struct FitOut {
  Double_t purity, ePurity, nWin, sig;
  Bool_t ok;
  Int_t status;  // TFitResult status (0 = converged); -1 = not fitted (too few entries)
};

// Mass projection (axis 'y' = A or 'z' = B) for k* bins [kLo, kHi] with the OTHER particle's
// mass inside its signal window [oLo, oHi].
TH1D* massProj(TH3* h, Char_t axis, Int_t kLo, Int_t kHi, Double_t oLo, Double_t oHi, const char* name) {
  TH1D* p = 0;
  if (axis == 'y') {
    const Int_t z1 = h->GetZaxis()->FindBin(oLo + 1e-7), z2 = h->GetZaxis()->FindBin(oHi - 1e-7);
    p = h->ProjectionY(name, kLo, kHi, z1, z2);
  } else {
    const Int_t y1 = h->GetYaxis()->FindBin(oLo + 1e-7), y2 = h->GetYaxis()->FindBin(oHi - 1e-7);
    p = h->ProjectionZ(name, kLo, kHi, y1, y2);
  }
  p->SetDirectory(0);
  return p;
}

// gaus+pol2 fit; if fixMean/fixSigma > 0 they are fixed. Purity of [wLo, wHi] = S/N.
// model 0: gaus + pol2 (invariant mass); model 1: gaus (K) + gaus (pi, mean fixed to m_pi^2) + pol1
// (TOF m^2). The signal Gaussian is always parameters 0-2.
// Peak model of the invariant-mass fits (model 0 only), per species from purityModel of the candidateFemto YAML
//: invariantMassGausPol2 -> 0 (Gaussian), invariantMassStudentTPol2 -> 1 (Student-t: the K0S and Xi-
// peaks have non-Gaussian tails, which a Gaussian counts as background), invariantMassDoubleGausPol2 -> 2 (two
// Gaussians, common mean). The extra shape parameters (nu; or sigma2 and amplitude ratio) are taken from the
// k*-integrated fit and fixed in the per-k* fits, like mean and sigma. gPeakModel is set by fitPurityA/B.
Int_t gPeakModel = 0, gPeakModelA = 0, gPeakModelB = 0;
Int_t PeakModelIndex(const std::string& purityModel) {
  if (purityModel == "invariantMassStudentTPol2") return 1;
  if (purityModel == "invariantMassDoubleGausPol2") return 2;
  return 0;  // invariantMassGausPol2 (and the TOF model, which does not use it)
}
// purityFitMinEntries of the candidateFemto YAML: projections with fewer entries are not fitted.
Int_t gMinEntries = -1;

FitOut fitPurityShape(TH1D* h, Double_t seedMean, Double_t seedSigma, Double_t fixMean, Double_t fixSigma,
                      Double_t wLo, Double_t wHi, const char* fname, Double_t* outMean, Double_t* outSigma,
                      const Double_t* fixExtra, Double_t* outExtra) {
  FitOut r; r.purity = 0; r.ePurity = 0; r.nWin = 0; r.sig = 0; r.ok = kFALSE; r.status = -1;
  const Int_t b1 = h->FindBin(wLo + 1e-7), b2 = h->FindBin(wHi - 1e-7);
  r.nWin = h->Integral(b1, b2);
  const Double_t eLo = h->GetXaxis()->GetBinLowEdge(b1), eHi = h->GetXaxis()->GetBinUpEdge(b2);
  if (h->GetEntries() < gMinEntries || r.nWin <= 0) return r;
  const Double_t xMin = h->GetXaxis()->GetXmin(), xMax = h->GetXaxis()->GetXmax(), bw = h->GetBinWidth(1);
  const Double_t x0 = seedMean;  // background polynomial centred at the peak (well conditioned)
  TString sigExpr = (gPeakModel == 1)
      ? "[0]*TMath::Power(1.+((x-[1])/[2])*((x-[1])/[2])/[3],-([3]+1.)/2.)"
      : "[0]*(exp(-0.5*((x-[1])/[2])*((x-[1])/[2]))+[3]*exp(-0.5*((x-[1])/[4])*((x-[1])/[4])))";
  const Int_t nSig = (gPeakModel == 1) ? 4 : 5;
  TF1* f = new TF1(fname, sigExpr + TString::Format("+[%d]+[%d]*(x-%.6f)+[%d]*(x-%.6f)*(x-%.6f)", nSig, nSig + 1, x0, nSig + 2, x0, x0),
                   xMin, xMax);
  f->SetParameter(0, TMath::Max(1., h->GetMaximum() * 0.5));
  f->SetParLimits(0, 0., 10. * h->GetMaximum() + 10.);
  if (fixMean > 0) f->FixParameter(1, fixMean); else { f->SetParameter(1, seedMean); f->SetParLimits(1, seedMean - 3 * seedSigma, seedMean + 3 * seedSigma); }
  if (fixSigma > 0) f->FixParameter(2, fixSigma); else { f->SetParameter(2, 0.8 * seedSigma); f->SetParLimits(2, 0.2 * seedSigma, 3. * seedSigma); }
  if (gPeakModel == 1) {
    if (fixExtra) f->FixParameter(3, fixExtra[0]); else { f->SetParameter(3, 5.); f->SetParLimits(3, 0.5, 200.); }
  } else {
    // fixExtra = the k*-integrated fit's values: fixed even when the amplitude ratio fitted to 0
    if (fixExtra) f->FixParameter(3, fixExtra[0]); else { f->SetParameter(3, 0.2); f->SetParLimits(3, 0., 2.); }
    if (fixExtra) f->FixParameter(4, fixExtra[1]); else { f->SetParameter(4, 2.5 * seedSigma); f->SetParLimits(4, 1.2 * seedSigma, 10. * seedSigma); }
  }
  f->SetParameter(nSig, TMath::Max(0., h->GetMinimum()));
  f->SetParameter(nSig + 1, 0.); f->SetParameter(nSig + 2, 0.);
  Int_t st = h->Fit(f, "QLR0");
  if (st != 0) st = h->Fit(f, "QLR0");
  TF1 sig(TString(fname) + "_sig", sigExpr, xMin, xMax);
  for (Int_t i = 0; i < nSig; ++i) sig.SetParameter(i, f->GetParameter(i));
  const Double_t a = f->GetParameter(0), ea = f->GetParError(0);
  r.sig = sig.Integral(eLo, eHi) / bw;
  r.purity = TMath::Max(0., TMath::Min(1., r.sig / r.nWin));
  r.ePurity = (a > 0) ? r.purity * ea / a : 0.;
  r.ok = (st == 0); r.status = st;
  if (outMean) *outMean = f->GetParameter(1);
  if (outSigma) *outSigma = TMath::Abs(f->GetParameter(2));
  if (outExtra) { outExtra[0] = f->GetParameter(3); outExtra[1] = (gPeakModel == 2) ? f->GetParameter(4) : 0.; }
  h->GetListOfFunctions()->Add(f);
  return r;
}

FitOut fitPurity(TH1D* h, Int_t model, Double_t seedMean, Double_t seedSigma, Double_t fixMean, Double_t fixSigma,
                 Double_t wLo, Double_t wHi, const char* fname, Double_t* outMean = 0, Double_t* outSigma = 0,
                 const Double_t* fixExtra = 0, Double_t* outExtra = 0) {
  if (model == 0 && gPeakModel != 0)
    return fitPurityShape(h, seedMean, seedSigma, fixMean, fixSigma, wLo, wHi, fname, outMean, outSigma, fixExtra, outExtra);
  FitOut r; r.purity = 0; r.ePurity = 0; r.nWin = 0; r.sig = 0; r.ok = kFALSE; r.status = -1;
  // The window is applied as whole histogram bins (every bin overlapping [wLo, wHi]); the signal
  // integral below uses the SAME bin-aligned edges, so S and N refer to the same mass range (the
  // pair counts N in the CF use the identical bin range). Integrating the Gaussian over the
  // nominal [wLo, wHi] instead biased P low by ~0.015 in the toy closure test.
  const Int_t b1 = h->FindBin(wLo + 1e-7), b2 = h->FindBin(wHi - 1e-7);
  r.nWin = h->Integral(b1, b2);
  const Double_t eLo = h->GetXaxis()->GetBinLowEdge(b1), eHi = h->GetXaxis()->GetBinUpEdge(b2);
  if (h->GetEntries() < gMinEntries || r.nWin <= 0) return r;
  const Double_t xMin = h->GetXaxis()->GetXmin(), xMax = h->GetXaxis()->GetXmax();
  TF1* f = (model == 0) ? new TF1(fname, "gaus(0)+pol2(3)", xMin, xMax)
                        : new TF1(fname, "gaus(0)+gaus(3)+pol1(6)", xMin, xMax);
  const Double_t bw = h->GetBinWidth(1);
  if (model == 0) {
    f->SetParameters(TMath::Max(1., h->GetMaximum() * 0.5), seedMean, seedSigma, h->GetMinimum() + 0.5, 0., 0.);
  } else {
    // pion peak at the PDG pion mass squared; its width starts at 0.03 GeV^2/c^4 and stays in 0.005-0.15
    // (fit start values and limits, not analysis choices)
    const Double_t mPi = TDatabasePDG::Instance()->GetParticle(211)->Mass();
    const Double_t mPi2 = mPi * mPi;
    f->SetParameters(TMath::Max(1., h->GetMaximum() * 0.5), seedMean, seedSigma, h->GetBinContent(1) + 0.5, mPi2, 0.03, 0.5, 0.);
    f->FixParameter(4, mPi2);
    f->SetParLimits(3, 0., 10. * h->GetMaximum() + 10.);
    f->SetParLimits(5, 0.005, 0.15);
  }
  f->SetParLimits(0, 0., 10. * h->GetMaximum() + 10.);
  if (fixMean > 0) f->FixParameter(1, fixMean); else f->SetParLimits(1, seedMean - 3 * seedSigma, seedMean + 3 * seedSigma);
  if (fixSigma > 0) f->FixParameter(2, fixSigma); else f->SetParLimits(2, 0.3 * seedSigma, 3. * seedSigma);
  Int_t st = h->Fit(f, "QLR0");
  // The pol2 in absolute mass (~1.3 GeV) is ill-conditioned: MIGRAD often ends with status 4 at
  // high statistics although the minimum is reached. One refit from the result converges with an
  // unchanged amplitude (e.g. 7162.28 -> 7162.50); only then the fit counts as valid.
  if (st != 0) st = h->Fit(f, "QLR0");
  const Double_t a = f->GetParameter(0), ea = f->GetParError(0);
  const Double_t mu = f->GetParameter(1), sg = TMath::Abs(f->GetParameter(2));
  // Gaussian counts inside the window (integral / bin width).
  const Double_t gInt = a * sg * TMath::Sqrt(TMath::PiOver2()) *
                        (TMath::Erf((eHi - mu) / (TMath::Sqrt2() * sg)) - TMath::Erf((eLo - mu) / (TMath::Sqrt2() * sg))) / bw;
  r.sig = gInt;
  r.purity = TMath::Max(0., TMath::Min(1., gInt / r.nWin));
  r.ePurity = (a > 0) ? r.purity * ea / a : 0.;
  r.ok = (st == 0);
  r.status = st;
  // TOF model: a free width that ends at one of its limits (0.3 or 3 times the start value) is not a measurement
  if (r.ok && model == 1 && fixSigma <= 0 && (sg < 0.301 * seedSigma || sg > 2.99 * seedSigma)) { r.ok = kFALSE; r.status = -2; }
  if (outMean) *outMean = mu;
  if (outSigma) *outSigma = sg;
  h->GetListOfFunctions()->Add(f);
  return r;
}

// fitPurity for species A / B: selects that species' peak model.
FitOut fitPurityA(TH1D* h, Int_t model, Double_t seedMean, Double_t seedSigma, Double_t fixMean, Double_t fixSigma,
                  Double_t wLo, Double_t wHi, const char* fname, Double_t* outMean = 0, Double_t* outSigma = 0,
                  const Double_t* fixExtra = 0, Double_t* outExtra = 0) {
  gPeakModel = gPeakModelA;
  return fitPurity(h, model, seedMean, seedSigma, fixMean, fixSigma, wLo, wHi, fname, outMean, outSigma, fixExtra, outExtra);
}
FitOut fitPurityB(TH1D* h, Int_t model, Double_t seedMean, Double_t seedSigma, Double_t fixMean, Double_t fixSigma,
                  Double_t wLo, Double_t wHi, const char* fname, Double_t* outMean = 0, Double_t* outSigma = 0,
                  const Double_t* fixExtra = 0, Double_t* outExtra = 0) {
  gPeakModel = gPeakModelB;
  return fitPurity(h, model, seedMean, seedSigma, fixMean, fixSigma, wLo, wHi, fname, outMean, outSigma, fixExtra, outExtra);
}

// Draw the fit attached by fitPurity(), if any (low-statistics projections are not fitted).
void drawFit(TH1* h) {
  TObject* o = h->GetListOfFunctions()->Last();
  if (o && o->InheritsFrom(TF1::Class())) ((TF1*)o)->Draw("same");
}

void drawWin(Double_t lo, Double_t hi, Double_t ymax) {
  TLine* l1 = new TLine(lo, 0, lo, ymax); l1->SetLineColor(kBlue); l1->SetLineStyle(2); l1->Draw();
  TLine* l2 = new TLine(hi, 0, hi, ymax); l2->SetLineColor(kBlue); l2->SetLineStyle(2); l2->Draw();
}

// The k*-differential purity correction over the merged k* bins [k1, k2] (where a per-bin purity fit
// fails, merge it with the neighbouring bin and make only that CF point coarser).
// Pair counts are summed over [k1, k2] and the four purities (A and B, SE and ME) are fitted on
// the mass projections summed over the same range, with the same fixed shapes as the per-bin fits.
// Returns kFALSE if any required fit fails; then val/err are not set.
struct GroupFitCfg {
  Bool_t aTof, unitA; Int_t aModel; Double_t aLo, aHi, xiLo, xiHi, k0Mean0, k0Sig0, xiMean0, xiSig0;
  Double_t k0Mu[2], k0Sg[2], xiMu[2], xiSg[2], k0Ex[2][2], xiEx[2][2];
  Double_t aSeedMu[2], aSeedSg[2];  // start values of the free A fits (TOF model): the k*-integrated fit
};
Bool_t purityCorrectedGroup(TH3* hh[2], Int_t k1, Int_t k2, const GroupFitCfg& g, Double_t& val, Double_t& err,
               Double_t& nSE, TH1D** projOut) {
  Double_t n[2], en[2], pa[2], epa[2], px[2], epx[2];
  for (Int_t s = 0; s < 2; ++s) {
    const Int_t y1 = hh[s]->GetYaxis()->FindBin(g.aLo + 1e-7), y2 = hh[s]->GetYaxis()->FindBin(g.aHi - 1e-7);
    const Int_t z1 = hh[s]->GetZaxis()->FindBin(g.xiLo + 1e-7), z2 = hh[s]->GetZaxis()->FindBin(g.xiHi - 1e-7);
    n[s] = hh[s]->IntegralAndError(k1, k2, y1, y2, z1, z2, en[s]);
    if (n[s] <= 0) return kFALSE;
    TH1D* pA = massProj(hh[s], 'y', k1, k2, g.xiLo, g.xiHi, TString::Format("mergedMassA_%d_%d_%d", s, k1, k2));
    TH1D* pX = massProj(hh[s], 'z', k1, k2, g.aLo, g.aHi, TString::Format("mergedMassB_%d_%d_%d", s, k1, k2));
    FitOut ra = fitPurityA(pA, g.aModel, g.aTof ? g.aSeedMu[s] : g.k0Mean0, g.aTof ? g.aSeedSg[s] : g.k0Sig0,
                          g.aTof ? -1. : g.k0Mu[s], g.aTof ? -1. : g.k0Sg[s],
                          g.aLo, g.aHi, TString::Format("fitMergedA_%d_%d_%d", s, k1, k2), 0, 0, g.k0Ex[s]);
    FitOut rx = fitPurityB(pX, 0, g.xiMean0, g.xiSig0, g.xiMu[s], g.xiSg[s], g.xiLo, g.xiHi,
                          TString::Format("fitMergedB_%d_%d_%d", s, k1, k2), 0, 0, g.xiEx[s]);
    if (projOut) { projOut[2 * s] = pA; projOut[2 * s + 1] = pX; }
    if (!(g.unitA || ra.ok) || !rx.ok || (!g.unitA && ra.purity <= 0) || rx.purity <= 0) return kFALSE;
    pa[s] = g.unitA ? 1. : ra.purity; epa[s] = g.unitA ? 0. : ra.ePurity;
    px[s] = rx.purity; epx[s] = rx.ePurity;
  }
  val = n[0] * pa[0] * px[0] / (n[1] * pa[1] * px[1]);
  Double_t rel2 = TMath::Power(en[0] / n[0], 2) + TMath::Power(en[1] / n[1], 2);
  for (Int_t s = 0; s < 2; ++s) {
    if (pa[s] > 0) rel2 += TMath::Power(epa[s] / pa[s], 2);
    rel2 += TMath::Power(epx[s] / px[s], 2);
  }
  err = val * TMath::Sqrt(rel2);
  nSE = n[0];
  return kTRUE;
}


// LaTeX label of a species key for titles (falls back to the key).
const char* speciesLabel(const std::string& key) {
  if (key == "k0s") return "K^{0}_{S}";
  if (key == "xi_minus") return "#Xi^{-}";
  if (key == "xi_plus") return "#bar{#Xi}^{+}";
  if (key == "kaon_plus") return "K^{+}";
  if (key == "kaon_minus") return "K^{-}";
  return key.c_str();
}

// "a=b;c=d" -> map
std::map<std::string, std::string> parseOptions(const char* options) {
  std::map<std::string, std::string> out;
  TObjArray* tok = TString(options ? options : "").Tokenize(";");
  for (Int_t i = 0; i < tok->GetEntries(); ++i) {
    const TString item = ((TObjString*)tok->At(i))->GetString();
    const Int_t eq = item.Index("=");
    if (eq < 0) { out[item.Strip(TString::kBoth).Data()] = ""; continue; }
    out[TString(item(0, eq)).Strip(TString::kBoth).Data()] = TString(item(eq + 1, item.Length())).Strip(TString::kBoth).Data();
  }
  delete tok;
  return out;
}

}  // namespace

Int_t purityCorrectedCorrelation(const char* rootFile, const char* mainconf, const char* channel, const char* outPrefix,
                                Int_t kstarRebin, const char* cent9List, const char* options) {
  gStyle->SetOptStat(0);
  gStyle->SetPaperSize(29.7, 21.0);  // landscape A4 for the PDF pages
  // ROOT5 TMinuit default (5000 calls) ends many high-statistics likelihood fits with status 4
  // although the result is stable; allow more calls so that convergence is real.
  TVirtualFitter::SetMaxIterations(100000);
  std::map<std::string, std::string> opt = parseOptions(options);
  const char* known[] = {"purityA", "windowA", "slideHalfWidth", "scoreMinA", "scoreMinB"};
  for (std::map<std::string, std::string>::const_iterator it = opt.begin(); it != opt.end(); ++it) {
    Bool_t ok = kFALSE;
    for (Int_t i = 0; i < 5; ++i) ok = ok || it->first == known[i];
    if (!ok) { std::cerr << "unknown option '" << it->first << "'" << std::endl; return 1; }
  }
  if (!mainconf || !*mainconf) { std::cerr << "the mainconf is required" << std::endl; return 1; }
  if (kstarRebin < 1 || !cent9List || !options) { std::cerr << "kstarRebin >= 1, cent9List and options are required" << std::endl; return 1; }
  std::cout << "mainconf: " << mainconf << std::endl;
  CandidateFemtoConfig cfg;
  if (!cfg.LoadFromMainConfig(mainconf)) { std::cerr << "cannot load the candidateFemto YAML of " << mainconf << std::endl; return 1; }
  gMinEntries = cfg.purityFitMinEntries;
  const Int_t slideHalfWidth = opt.count("slideHalfWidth") ? TString(opt["slideHalfWidth"].c_str()).Atoi() : cfg.slidingWindowHalfWidth;
  std::cout << "k* rebin " << kstarRebin << ", sliding-window half width " << slideHalfWidth << ", fit minimum entries "
            << gMinEntries << ", at most " << cfg.purityCorrectionMaxMergedBins << " merged k* bins" << std::endl;
  Int_t ich = -1;
  for (size_t i = 0; i < cfg.channels.size(); ++i) if (cfg.channels[i].name == channel) ich = i;
  if (ich < 0) { std::cerr << "channel " << channel << " is not in the candidateFemto YAML" << std::endl; return 1; }
  const CandidateFemtoSpecies& spA = cfg.species[cfg.channels[ich].indexA];
  const CandidateFemtoSpecies& spB = cfg.species[cfg.channels[ich].indexB];
  if (spB.purityByDefinition) { std::cerr << "purityByDefinition is supported for species A only" << std::endl; return 1; }
  if (spB.purityModel == "tofMassSquaredKaonPion") { std::cerr << "species B must use an invariant-mass purity model" << std::endl; return 1; }
  gPeakModelA = PeakModelIndex(spA.purityModel); gPeakModelB = PeakModelIndex(spB.purityModel);
  std::cout << "purity peak models: A " << spA.purityModel << ", B " << spB.purityModel << std::endl;
  const char* aName = speciesLabel(spA.key);
  const char* bName = speciesLabel(spB.key);
  // A uses the TOF m^2 model: per-k* mean/sigma float, the integrated ME shape is not imposed.
  const Bool_t aTof = (spA.purityModel == "tofMassSquaredKaonPion");
  const Int_t aModel = aTof ? 1 : 0;
  Double_t aLo = spA.signal.min, aHi = spA.signal.max;
  const Double_t bLo = spB.signal.min, bHi = spB.signal.max;
  if (opt.count("windowA")) {  // cut-scan override of the A signal window
    const TString w = opt["windowA"].c_str();
    const Int_t c = w.Index(",");
    if (c < 0) { std::cerr << "windowA must be lo,hi" << std::endl; return 1; }
    aLo = TString(w(0, c)).Atof(); aHi = TString(w(c + 1, w.Length())).Atof();
    if (!(aHi > aLo)) { std::cerr << "windowA needs hi > lo" << std::endl; return 1; }
    std::cout << "A signal window overridden: [" << aLo << ", " << aHi << "]" << std::endl;
  }
  // Only the purity fit range is projected: a window reaching outside it would be cut silently.
  if (aLo < spA.purityFitRange.min || aHi > spA.purityFitRange.max) {
    std::cerr << "ERROR: A signal window [" << aLo << ", " << aHi << "] is not inside purityFitRange of " << spA.key << std::endl;
    return 1;
  }
  const TString purityA = opt.count("purityA") ? opt["purityA"].c_str() : "config";
  if (purityA != "config" && purityA != "fit" && purityA != "unit") { std::cerr << "purityA must be config, fit or unit" << std::endl; return 1; }
  const Bool_t unitA = (purityA == "unit") || (purityA == "config" && spA.purityByDefinition);
  std::cout << "purity of " << spA.key << ": " << (unitA ? "1 (by definition)" : "fitted") << std::endl;

  TFile* f = TFile::Open(rootFile);
  if (!f || f->IsZombie()) return 1;
  THnSparse* sparse[2] = {(THnSparse*)f->Get(TString::Format("hPairSE_%s", channel)),
                          (THnSparse*)f->Get(TString::Format("hPairME_%s", channel))};
  if (!sparse[0] || !sparse[1]) { std::cerr << "hPairSE/ME_" << channel << " missing in " << rootFile << std::endl; return 1; }
  CandidatePairSelection sel;
  sel.SetCentList(cent9List);
  sel.scoreMinA = spA.scoreMin;
  sel.scoreMinB = spB.scoreMin;
  // Threshold-scan override (like windowA): scoreMinA / scoreMinB replace species_<s>_scoreMin of the YAML; the value must be
  // a score-axis edge (checked below). For scans only: a chosen threshold belongs in the selection YAML.
  if (opt.count("scoreMinA")) { sel.scoreMinA = TString(opt["scoreMinA"].c_str()).Atof(); std::cout << "score A threshold overridden: " << sel.scoreMinA << std::endl; }
  if (opt.count("scoreMinB")) { sel.scoreMinB = TString(opt["scoreMinB"].c_str()).Atof(); std::cout << "score B threshold overridden: " << sel.scoreMinB << std::endl; }
  for (Int_t s = 0; s < 2; ++s)
    if (!CandidateScoreThresholdIsEdge(sparse[s], 6, sel.scoreMinA) || !CandidateScoreThresholdIsEdge(sparse[s], 7, sel.scoreMinB)) return 1;
  std::cout << "selection: " << sel.Describe() << std::endl;
  // (k*, M_A, M_B) over the purity-fit ranges, the input of the k*-differential purity correction.
  TH3* hS = CandidateProjectKstarMassMass(sparse[0], sel, spA.purityFitRange.min, spA.purityFitRange.max,
                                          spB.purityFitRange.min, spB.purityFitRange.max, "pairSE");
  TH3* hM = CandidateProjectKstarMassMass(sparse[1], sel, spA.purityFitRange.min, spA.purityFitRange.max,
                                          spB.purityFitRange.min, spB.purityFitRange.max, "pairME");
  // Fine-k* copies for the sliding-window correction, before the coarse k* rebin.
  TH3* hSf = (TH3*)hS->Clone("fineSE"); hSf->SetDirectory(0);
  TH3* hMf = (TH3*)hM->Clone("fineME"); hMf->SetDirectory(0);
  if (hS->GetNbinsX() % kstarRebin != 0) {
    std::cerr << "kstarRebin must divide the number of fine k* bins (" << hS->GetNbinsX() << ")" << std::endl;
    return 1;
  }
  hS->RebinX(kstarRebin); hM->RebinX(kstarRebin);
  const Int_t nK = hS->GetNbinsX();

  const Double_t k0Mean0 = spA.purityFitSeed.min, k0Sig0 = spA.purityFitSeed.max;
  const Double_t xiMean0 = spB.purityFitSeed.min, xiSig0 = spB.purityFitSeed.max;
  const TString pdf = TString(outPrefix) + ".pdf";
  {
    const TString dir = gSystem->DirName(outPrefix);
    if (dir != "" && dir != "." && gSystem->AccessPathName(dir)) gSystem->mkdir(dir, kTRUE);
    if (gSystem->AccessPathName(dir, kWritePermission)) { std::cerr << "output directory not writable: " << dir << std::endl; return 1; }
  }
  TCanvas* c = new TCanvas("c", "c", 1400, 900);
  c->Print(pdf + "[");
  // Everything drawn in the PDF is also written to the output ROOT file, so
  // figures can be restyled later without rerunning: mass projections with their fit functions
  // attached (fits/), sliding-window projections (fits_sliding/), merged-bin projections (fits_merged/).
  std::vector<TH1*> keepFit, keepSliding;

  // k*-integrated fits (all k* in the histogram) fix mean/sigma per sample.
  TH3* hh[2] = {hS, hM};
  const char* tag[2] = {"SE", "ME"};
  Double_t k0Mu[2], k0Sg[2], xiMu[2], xiSg[2];
  Double_t k0Ex[2][2] = {{0., 0.}, {0., 0.}}, xiEx[2][2] = {{0., 0.}, {0., 0.}};
  FitOut k0Int[2], xiInt[2];
  c->Clear(); c->Divide(2, 2);
  // ME first: with a non-Gaussian peak model the peak shape (mean, sigma, tail parameters) is
  // taken from the high-statistics ME sample and fixed in the SE fits -- SE and ME contain the
  // same candidates, so the shape is the same; free tails at SE statistics absorb background
  // (e.g. Student-t nu=0.8, a 16 MeV wide second Gaussian). With the Gaussian model SE and ME are fitted
  // independently.
  const Int_t order[2] = {1, 0};
  for (Int_t io = 0; io < 2; ++io) {
    const Int_t s = order[io];
    TH1D* pk = massProj(hh[s], 'y', 1, nK, bLo, bHi, TString::Format("integratedA_%s", tag[s]));
    TH1D* px = massProj(hh[s], 'z', 1, nK, aLo, aHi, TString::Format("integratedB_%s", tag[s]));
    keepFit.push_back(pk); keepFit.push_back(px);
    const Bool_t useMeShapeA = (gPeakModelA != 0 && s == 0), useMeShapeB = (gPeakModelB != 0 && s == 0);
    k0Int[s] = useMeShapeA && !aTof
        ? fitPurityA(pk, aModel, k0Mean0, k0Sig0, k0Mu[1], k0Sg[1], aLo, aHi, TString::Format("fitIntegratedA_%s", tag[s]), &k0Mu[s], &k0Sg[s], k0Ex[1], k0Ex[s])
        : fitPurityA(pk, aModel, k0Mean0, k0Sig0, -1, -1, aLo, aHi, TString::Format("fitIntegratedA_%s", tag[s]), &k0Mu[s], &k0Sg[s], 0, k0Ex[s]);
    xiInt[s] = useMeShapeB
        ? fitPurityB(px, 0, xiMean0, xiSig0, xiMu[1], xiSg[1], bLo, bHi, TString::Format("fitIntegratedB_%s", tag[s]), &xiMu[s], &xiSg[s], xiEx[1], xiEx[s])
        : fitPurityB(px, 0, xiMean0, xiSig0, -1, -1, bLo, bHi, TString::Format("fitIntegratedB_%s", tag[s]), &xiMu[s], &xiSg[s], 0, xiEx[s]);
    if (gPeakModelA != 0 || gPeakModelB != 0)
      std::cout << "integrated " << tag[s] << ": K0 shape extras " << k0Ex[s][0] << " " << k0Ex[s][1] << " status " << k0Int[s].status
                << " | Xi shape extras " << xiEx[s][0] << " " << xiEx[s][1] << " status " << xiInt[s].status << std::endl;
    c->cd(1 + s);
    pk->SetTitle(TString::Format("%s: %s, %s in window, all k*  P=%.3f#pm%.3f (#mu=%.4f #sigma=%.4f);M [GeV/c^{2}];pairs",
                                 tag[s], aName, bName, k0Int[s].purity, k0Int[s].ePurity, k0Mu[s], k0Sg[s]));
    pk->Draw("E"); drawFit(pk); drawWin(aLo, aHi, pk->GetMaximum());
    c->cd(3 + s);
    px->SetTitle(TString::Format("%s: M(%s), %s in window, all k*  P=%.3f#pm%.3f (#mu=%.4f #sigma=%.4f);M [GeV/c^{2}];pairs",
                                 tag[s], bName, aName, xiInt[s].purity, xiInt[s].ePurity, xiMu[s], xiSg[s]));
    px->Draw("E"); drawFit(px); drawWin(bLo, bHi, px->GetMaximum());
  }
  c->Print(pdf);
  // Start values of the free (TOF model) A fits: the k*-integrated fit of the same sample, if it converged.
  Double_t aSeedMu[2], aSeedSg[2];
  for (Int_t s = 0; s < 2; ++s) {
    const Bool_t useInt = aTof && k0Int[s].ok;
    aSeedMu[s] = useInt ? k0Mu[s] : k0Mean0;
    aSeedSg[s] = useInt ? k0Sg[s] : k0Sig0;
    if (aTof) std::cout << "A (TOF m^2) fits of " << tag[s] << " start at mean " << aSeedMu[s] << ", sigma " << aSeedSg[s]
                        << (useInt ? " (k*-integrated fit)" : " (purityFitSeed: integrated fit failed)") << std::endl;
  }

  // Per-k* fits. All k* histograms below are booked on [0, kMax]: the k* axis must start at 0.
  const Double_t kMax = hS->GetXaxis()->GetXmax();
  if (hS->GetXaxis()->GetXmin() != 0.) {
    std::cerr << "ERROR: kstarAxis must start at 0 (it starts at " << hS->GetXaxis()->GetXmin() << ")" << std::endl;
    return 1;
  }
  TH1D* hP[2][2]; TH1D* hN[2];
  Int_t nInvalidBins = 0, nSlidingFallbacks = 0, nSlidingInvalid = 0;
  for (Int_t s = 0; s < 2; ++s) {
    hP[s][0] = new TH1D(TString::Format("purityA_%s", tag[s]), ";k* [GeV/c];purity", nK, 0., kMax);
    hP[s][1] = new TH1D(TString::Format("purityB_%s", tag[s]), ";k* [GeV/c];purity", nK, 0., kMax);
    hN[s] = new TH1D(TString::Format("pairsInWindows_%s", tag[s]), ";k* [GeV/c];pairs (both in window)", nK, 0., kMax);
  }
  for (Int_t s = 0; s < 2; ++s) {
    for (Int_t sp = 0; sp < 2; ++sp) {
      const Int_t nCol = 5, nRow = (nK + nCol - 1) / nCol;
      c->Clear(); c->Divide(nCol, nRow);
      for (Int_t k = 1; k <= nK; ++k) {
        TH1D* p = (sp == 0)
            ? massProj(hh[s], 'y', k, k, bLo, bHi, TString::Format("massA_%s_%d", tag[s], k))
            : massProj(hh[s], 'z', k, k, aLo, aHi, TString::Format("massB_%s_%d", tag[s], k));
        FitOut r = (sp == 0)
            // K+ m^2 width depends on p, which varies with k*: let mean/sigma float per k* bin
            // (K+ statistics are large); K0S keeps them fixed to the k*-integrated fit.
            ? fitPurityA(p, aModel, aTof ? aSeedMu[s] : k0Mean0, aTof ? aSeedSg[s] : k0Sig0, aTof ? -1. : k0Mu[s], aTof ? -1. : k0Sg[s],
                        aLo, aHi, TString::Format("fitA_%s_%d", tag[s], k), 0, 0, k0Ex[s])
            : fitPurityB(p, 0, xiMean0, xiSig0, xiMu[s], xiSg[s], bLo, bHi, TString::Format("fitB_%s_%d", tag[s], k), 0, 0, xiEx[s]);
        keepFit.push_back(p);
        if (r.nWin > 0) {
          const Bool_t unitPurity = unitA && sp == 0;  // K+: purity 1 by definition (unless options purityA=fit)
          if (unitPurity || r.ok) {
            hP[s][sp]->SetBinContent(k, unitPurity ? 1. : r.purity);
            hP[s][sp]->SetBinError(k, unitPurity ? 0. : r.ePurity);
          } else {
            // Failed / too-few-entries fit: the purity is left at 0, which marks the
            // bin invalid -- the purity-corrected CF skips it instead of using an unconverged value.
            ++nInvalidBins;
            std::cout << "PURITY FIT INVALID: " << tag[s] << " " << (sp == 0 ? aName : bName) << " k* bin " << k
                      << " entries " << p->GetEntries() << " status " << r.status << " purity " << r.purity << std::endl;
          }
        }
        c->cd(k);
        p->SetTitle(TString::Format("%s %s k*=[%.2f,%.2f] P=%.2f", tag[s], sp == 0 ? spA.key.c_str() : spB.key.c_str(),
                                    hS->GetXaxis()->GetBinLowEdge(k), hS->GetXaxis()->GetBinUpEdge(k), r.purity));
        p->Draw("E");
        drawFit(p);
      }
      c->Print(pdf);
    }
    // Pairs with both masses in the signal windows, per k*.
    for (Int_t k = 1; k <= nK; ++k) {
      const Int_t y1 = hh[s]->GetYaxis()->FindBin(aLo + 1e-7), y2 = hh[s]->GetYaxis()->FindBin(aHi - 1e-7);
      const Int_t z1 = hh[s]->GetZaxis()->FindBin(bLo + 1e-7), z2 = hh[s]->GetZaxis()->FindBin(bHi - 1e-7);
      Double_t e = 0; const Double_t n = hh[s]->IntegralAndError(k, k, y1, y2, z1, z2, e);
      hN[s]->SetBinContent(k, n); hN[s]->SetBinError(k, e);
    }
  }

  // Sliding-window purity correction: purity of each FINE k* bin k is
  // taken from a fit of the fine bins [k - w, k + w] (clipped at the edges), so neighbouring bins
  // share statistics and the purity curve is locally smoothed. The CF is then built in the coarse
  // bins as sum_fine(N P_A P_B)_SE / sum_fine(N P_A P_B)_ME. Errors treat the fine-bin purities as
  // independent, which overstates the precision of the smoothing somewhat (adjacent windows overlap).
  TH3* hf[2] = {hSf, hMf};
  const Int_t nF = hSf->GetNbinsX();
  TH1D* hPuritySliding[2][2]; TH1D* hSumSliding[2];
  for (Int_t s = 0; s < 2; ++s) {
    hPuritySliding[s][0] = new TH1D(TString::Format("puritySlidingA_%s", tag[s]), ";k* [GeV/c];purity", nF, 0., kMax);
    hPuritySliding[s][1] = new TH1D(TString::Format("puritySlidingB_%s", tag[s]), ";k* [GeV/c];purity", nF, 0., kMax);
    hSumSliding[s] = new TH1D(TString::Format("sumSliding_%s", tag[s]), ";k* [GeV/c];N P_{A} P_{B}", nK, 0., kMax);
    const Int_t y1 = hf[s]->GetYaxis()->FindBin(aLo + 1e-7), y2 = hf[s]->GetYaxis()->FindBin(aHi - 1e-7);
    const Int_t z1 = hf[s]->GetZaxis()->FindBin(bLo + 1e-7), z2 = hf[s]->GetZaxis()->FindBin(bHi - 1e-7);
    std::vector<Double_t> sum(nK + 1, 0.), var(nK + 1, 0.);
    std::vector<Bool_t> badSliding(nK + 1, kFALSE);
    for (Int_t k = 1; k <= nF; ++k) {
      const Int_t lo = TMath::Max(1, k - slideHalfWidth), hi = TMath::Min(nF, k + slideHalfWidth);
      Double_t pr[2] = {0., 0.}, epr[2] = {0., 0.};
      for (Int_t sp = 0; sp < 2; ++sp) {
        if (unitA && sp == 0) { pr[0] = 1.; continue; }  // K+: purity 1 by definition
        TH1D* p = (sp == 0) ? massProj(hf[s], 'y', lo, hi, bLo, bHi, TString::Format("slidingA_%s_%d", tag[s], k))
                            : massProj(hf[s], 'z', lo, hi, aLo, aHi, TString::Format("slidingB_%s_%d", tag[s], k));
        keepSliding.push_back(p);
        // TOF model: free mean/sigma as in the per-bin fits (the width depends on k*)
        FitOut r = (sp == 0) ? fitPurityA(p, aModel, aTof ? aSeedMu[s] : k0Mean0, aTof ? aSeedSg[s] : k0Sig0, aTof ? -1. : k0Mu[s],
                                         aTof ? -1. : k0Sg[s], aLo, aHi, TString::Format("fSlidingA_%s_%d", tag[s], k), 0, 0, k0Ex[s])
                             : fitPurityB(p, 0, xiMean0, xiSig0, xiMu[s], xiSg[s], bLo, bHi, TString::Format("fSlidingB_%s_%d", tag[s], k), 0, 0, xiEx[s]);
        if (r.ok) { pr[sp] = r.purity; epr[sp] = r.ePurity; }
        else {
          // Invalid sliding-window fit: use the k*-differential purity of the enclosing
          // coarse k* bin instead of 0 (0 would silently drop these SE/ME pairs from the sums).
          const Int_t kcB = (k - 1) / kstarRebin + 1;
          pr[sp] = hP[s][sp]->GetBinContent(kcB); epr[sp] = hP[s][sp]->GetBinError(kcB);
          ++nSlidingFallbacks;
        }
        // (kept for the output file: keepSliding)
      }
      hPuritySliding[s][0]->SetBinContent(k, pr[0]); hPuritySliding[s][0]->SetBinError(k, epr[0]);
      hPuritySliding[s][1]->SetBinContent(k, pr[1]); hPuritySliding[s][1]->SetBinError(k, epr[1]);
      Double_t en = 0; const Double_t n = hf[s]->IntegralAndError(k, k, y1, y2, z1, z2, en);
      const Double_t w = pr[0] * pr[1];
      const Int_t kc = (k - 1) / kstarRebin + 1;
      if (n > 0 && w <= 0) badSliding[kc] = kTRUE;  // no valid purity even after the fallback
      sum[kc] += n * w;
      Double_t relw2 = 0.;
      if (pr[0] > 0) relw2 += TMath::Power(epr[0] / pr[0], 2);
      if (pr[1] > 0) relw2 += TMath::Power(epr[1] / pr[1], 2);
      var[kc] += w * w * en * en + n * n * w * w * relw2;
    }
    for (Int_t kc = 1; kc <= nK; ++kc) {
      if (badSliding[kc]) {  // coarse bin left empty: the sliding-window CF is not defined there
        ++nSlidingInvalid;
        std::cout << "SLIDING-WINDOW PURITY INVALID coarse bin " << tag[s] << " k* bin " << kc << std::endl;
        continue;
      }
      hSumSliding[s]->SetBinContent(kc, sum[kc]); hSumSliding[s]->SetBinError(kc, TMath::Sqrt(var[kc]));
    }
  }
  TH1D* cfSliding = (TH1D*)hSumSliding[0]->Clone("cfSlidingWindowPurity");
  cfSliding->Divide(hSumSliding[1]);

  // NOTE: the mass histogram bin edges need not coincide with the window edges; FindBin on the
  // window edges keeps every bin that contains part of the window (documented approximation).

  // CF: raw (window counts only), k*-differential purity correction, constant purity (integrated).
  TH1D* cfRaw = new TH1D("cfRaw", ";k* [GeV/c];C(k*)", nK, 0., kMax);
  TH1D* cfPurity = new TH1D("cfPurityCorrected", ";k* [GeV/c];C(k*)", nK, 0., kMax);
  TH1D* cfConstant = new TH1D("cfConstantPurity", ";k* [GeV/c];C(k*)", nK, 0., kMax);
  for (Int_t k = 1; k <= nK; ++k) {
    const Double_t ns = hN[0]->GetBinContent(k), nm = hN[1]->GetBinContent(k);
    if (ns <= 0 || nm <= 0) continue;
    const Double_t r = ns / nm, er = r * TMath::Sqrt(1. / ns + 1. / nm);
    cfRaw->SetBinContent(k, r); cfRaw->SetBinError(k, er);
    const Double_t p1s = hP[0][0]->GetBinContent(k), p2s = hP[0][1]->GetBinContent(k);
    const Double_t p1m = hP[1][0]->GetBinContent(k), p2m = hP[1][1]->GetBinContent(k);
    if (p1s > 0 && p2s > 0 && p1m > 0 && p2m > 0) {
      const Double_t v = r * p1s * p2s / (p1m * p2m);
      Double_t rel2 = 1. / ns + 1. / nm;
      rel2 += TMath::Power(hP[0][0]->GetBinError(k) / p1s, 2) + TMath::Power(hP[0][1]->GetBinError(k) / p2s, 2);
      rel2 += TMath::Power(hP[1][0]->GetBinError(k) / p1m, 2) + TMath::Power(hP[1][1]->GetBinError(k) / p2m, 2);
      cfPurity->SetBinContent(k, v); cfPurity->SetBinError(k, v * TMath::Sqrt(rel2));
    }
    const Double_t pa0 = unitA ? 1. : k0Int[0].purity, pa1 = unitA ? 1. : k0Int[1].purity;
    const Double_t pc = (pa0 * xiInt[0].purity) / TMath::Max(1e-9, pa1 * xiInt[1].purity);
    cfConstant->SetBinContent(k, r * pc); cfConstant->SetBinError(k, er * pc);
  }
  TH1D* cfs[4] = {cfRaw, cfConstant, cfPurity, cfSliding};
  for (Int_t i = 0; i < 4; ++i) {
    const Int_t b1 = cfs[i]->FindBin(cfg.cfNormMin + 1e-7), b2 = cfs[i]->FindBin(cfg.cfNormMax - 1e-7);
    Double_t sw = 0, s = 0;
    for (Int_t b = b1; b <= b2; ++b) {
      const Double_t e = cfs[i]->GetBinError(b);
      if (cfs[i]->GetBinContent(b) <= 0 || e <= 0) continue;
      s += cfs[i]->GetBinContent(b) / (e * e); sw += 1. / (e * e);
    }
    if (sw > 0) cfs[i]->Scale(1. / (s / sw));
  }

  // The k*-differential purity correction with merged bins where a fit fails (cfMerged, variable k* bins). A bin is kept alone when
  // all of its per-bin fits converged (then the point equals cfPurity); otherwise it is merged with the
  // following bin(s), at most kMaxMerge (purityCorrectionMaxMergedBins) bins in total, until the fits on the
  // merged range converge. A failing last group is merged into the previous group. Bins that stay invalid are dropped.
  const Int_t kMaxMerge = cfg.purityCorrectionMaxMergedBins;
  GroupFitCfg gc;
  gc.aTof = aTof; gc.unitA = unitA; gc.aModel = aModel; gc.aLo = aLo; gc.aHi = aHi; gc.xiLo = bLo; gc.xiHi = bHi;
  gc.k0Mean0 = k0Mean0; gc.k0Sig0 = k0Sig0; gc.xiMean0 = xiMean0; gc.xiSig0 = xiSig0;
  for (Int_t s = 0; s < 2; ++s) {
    gc.k0Mu[s] = k0Mu[s]; gc.k0Sg[s] = k0Sg[s]; gc.xiMu[s] = xiMu[s]; gc.xiSg[s] = xiSg[s];
    gc.aSeedMu[s] = aSeedMu[s]; gc.aSeedSg[s] = aSeedSg[s];
    for (Int_t j = 0; j < 2; ++j) { gc.k0Ex[s][j] = k0Ex[s][j]; gc.xiEx[s][j] = xiEx[s][j]; }
  }
  std::vector<Int_t> gLo, gHi; std::vector<Double_t> gVal, gErr;
  std::vector<TH1D*> mergedProj;  // mass projections of merged groups, for the QA page
  std::vector<TString> mergedTitle;
  Int_t nDropped = 0;
  for (Int_t k = 1; k <= nK;) {
    Double_t v = 0, e = 0, ns = 0;
    Int_t k2 = k;
    Bool_t ok = kFALSE;
    TH1D* proj[4] = {0, 0, 0, 0};
    for (; k2 <= TMath::Min(nK, k + kMaxMerge - 1); ++k2) {
      ok = purityCorrectedGroup(hh, k, k2, gc, v, e, ns, proj);
      if (ok) break;
    }
    if (!ok && !gLo.empty()) {  // merge the failing tail into the previous group
      const Int_t pk1 = gLo.back();
      const Int_t pk2 = TMath::Min(nK, k + kMaxMerge - 1);
      if (purityCorrectedGroup(hh, pk1, pk2, gc, v, e, ns, proj)) {
        gLo.back() = pk1; gHi.back() = pk2; gVal.back() = v; gErr.back() = e;
        std::cout << "merged into previous group: k* bins " << pk1 << "-" << pk2 << std::endl;
        for (Int_t j = 0; j < 4; ++j) if (proj[j]) { mergedProj.push_back(proj[j]); mergedTitle.push_back(TString::Format("bins %d-%d", pk1, pk2)); }
        k = pk2 + 1;
        continue;
      }
    }
    if (!ok) {
      ++nDropped;
      std::cout << "PURITY FIT INVALID after merging: k* bins " << k << "-" << TMath::Min(nK, k + kMaxMerge - 1) << " dropped" << std::endl;
      k = TMath::Min(nK, k + kMaxMerge - 1) + 1;
      continue;
    }
    if (k2 > k) {
      std::cout << "merged k* bins " << k << "-" << k2 << " ([" << hS->GetXaxis()->GetBinLowEdge(k) << ", "
                << hS->GetXaxis()->GetBinUpEdge(k2) << "] GeV/c)" << std::endl;
      for (Int_t j = 0; j < 4; ++j) if (proj[j]) { mergedProj.push_back(proj[j]); mergedTitle.push_back(TString::Format("bins %d-%d", k, k2)); }
    }
    gLo.push_back(k); gHi.push_back(k2); gVal.push_back(v); gErr.push_back(e);
    k = k2 + 1;
  }
  std::vector<Double_t> edges;
  for (size_t i = 0; i < gLo.size(); ++i) {
    if (i == 0 || gLo[i] != gHi[i - 1] + 1) edges.push_back(hS->GetXaxis()->GetBinLowEdge(gLo[i]));
    edges.push_back(hS->GetXaxis()->GetBinUpEdge(gHi[i]));
  }
  // Gaps (dropped bins) become empty bins of the variable-width histogram.
  TH1D* cfMerged = new TH1D("cfPurityCorrectedMergedBins", ";k* [GeV/c];C(k*)", TMath::Max(1, (Int_t)edges.size() - 1), edges.empty() ? 0 : &edges[0]);
  if (edges.empty()) cfMerged->SetBins(1, 0., kMax);
  for (size_t i = 0; i < gLo.size(); ++i) {
    const Int_t b = cfMerged->FindBin(0.5 * (hS->GetXaxis()->GetBinLowEdge(gLo[i]) + hS->GetXaxis()->GetBinUpEdge(gHi[i])));
    cfMerged->SetBinContent(b, gVal[i]); cfMerged->SetBinError(b, gErr[i]);
  }
  {
    Double_t sw = 0, s = 0;
    for (Int_t b = 1; b <= cfMerged->GetNbinsX(); ++b) {
      const Double_t lo = cfMerged->GetXaxis()->GetBinLowEdge(b), hi = cfMerged->GetXaxis()->GetBinUpEdge(b);
      if (lo < cfg.cfNormMin - 1e-7 || hi > cfg.cfNormMax + 1e-7) continue;
      const Double_t e = cfMerged->GetBinError(b);
      if (cfMerged->GetBinContent(b) <= 0 || e <= 0) continue;
      s += cfMerged->GetBinContent(b) / (e * e); sw += 1. / (e * e);
    }
    if (sw > 0) cfMerged->Scale(1. / (s / sw));
  }
  std::cout << "merged-bin groups " << gLo.size() << " merged-projection pads " << mergedProj.size() << " dropped " << nDropped << std::endl;

  c->Clear(); c->Divide(2, 2);
  for (Int_t sp = 0; sp < 2; ++sp) {
    c->cd(1 + sp);
    hP[0][sp]->SetTitle(TString::Format("%s purity vs k*: SE (black), ME (red);k* [GeV/c];purity", sp == 0 ? aName : bName));
    hP[0][sp]->SetMinimum(0); hP[0][sp]->SetMaximum(1.05); hP[0][sp]->SetMarkerStyle(20); hP[0][sp]->Draw("E");
    hP[1][sp]->SetLineColor(kRed); hP[1][sp]->SetMarkerColor(kRed); hP[1][sp]->SetMarkerStyle(24); hP[1][sp]->Draw("E SAME");
    // sliding-window fine-bin purities as lines.
    hPuritySliding[0][sp]->SetLineColor(kGray + 2); hPuritySliding[0][sp]->Draw("HIST SAME");
    hPuritySliding[1][sp]->SetLineColor(kOrange + 1); hPuritySliding[1][sp]->Draw("HIST SAME");
  }
  c->cd(3);
  hN[0]->SetTitle("SE pairs, both in signal window;k* [GeV/c];pairs"); hN[0]->Draw("E");
  c->cd(4);
  cfRaw->SetTitle(TString::Format("CF (norm. %.1f-%.1f): raw (gray), constant purity (blue), k*-differential purity (black), sliding window (green, w=#pm%d);k* [GeV/c];C(k*)", cfg.cfNormMin, cfg.cfNormMax, slideHalfWidth));
  {
    // y range from the drawn CFs (content +- error of the points measured better than 100%)
    Double_t lo = 0., hi = 2.;
    TH1D* drawn[4] = {cfRaw, cfConstant, cfPurity, cfSliding};
    for (Int_t i = 0; i < 4; ++i)
      for (Int_t b = 1; b <= drawn[i]->GetNbinsX(); ++b) {
        if (drawn[i]->GetBinError(b) <= 0 || drawn[i]->GetBinError(b) >= TMath::Abs(drawn[i]->GetBinContent(b))) continue;
        lo = TMath::Min(lo, drawn[i]->GetBinContent(b) - drawn[i]->GetBinError(b));
        hi = TMath::Max(hi, drawn[i]->GetBinContent(b) + drawn[i]->GetBinError(b));
      }
    cfRaw->SetMinimum(lo - 0.1 * (hi - lo)); cfRaw->SetMaximum(hi + 0.1 * (hi - lo));
  }
  cfRaw->SetLineColor(kGray + 1); cfRaw->SetMarkerColor(kGray + 1); cfRaw->Draw("E");
  cfConstant->SetLineColor(kBlue); cfConstant->SetMarkerColor(kBlue); cfConstant->Draw("E SAME");
  cfPurity->SetMarkerStyle(20); cfPurity->Draw("E SAME");
  cfSliding->SetLineColor(kGreen + 2); cfSliding->SetMarkerColor(kGreen + 2); cfSliding->SetMarkerStyle(25); cfSliding->Draw("E SAME");
  TLine* one = new TLine(0, 1, kMax, 1); one->SetLineStyle(2); one->Draw();
  c->Print(pdf);

  // QA page: mass fits of the merged k* ranges (only if some bins were merged).
  if (!mergedProj.empty()) {
    c->Clear();
    const Int_t nPad = (Int_t)mergedProj.size(), nCol = 4, nRow = (nPad + nCol - 1) / nCol;
    c->Divide(nCol, nRow);
    const char* lab[4] = {"SE A", "SE B", "ME A", "ME B"};
    for (Int_t i = 0; i < nPad; ++i) {
      c->cd(i + 1);
      mergedProj[i]->SetTitle(TString::Format("merged %s: %s", mergedTitle[i].Data(), lab[i % 4]));
      mergedProj[i]->Draw("E"); drawFit(mergedProj[i]);
    }
    c->Print(pdf);
  }

  // Final page: the k*-differential purity-corrected result alone (no overlays). Bins whose purity
  // fit failed are merged with their neighbour (cfMerged), so the result has variable bins.
  {
    c->Clear(); c->cd();
    gPad->SetLeftMargin(0.10); gPad->SetRightMargin(0.04);
    TH1D* fin = (TH1D*)cfMerged->Clone("cfFinal");
    fin->SetTitle(TString::Format("%s - %s correlation function (k*-differential purity correction);k* [GeV/c];C(k*)", aName, bName));
    fin->SetLineColor(kBlack); fin->SetMarkerColor(kBlack); fin->SetMarkerStyle(20); fin->SetMarkerSize(1.1);
    Double_t ymax = 1.5;
    for (Int_t b = 1; b <= fin->GetNbinsX(); ++b)  // points measured better than 100% set the range
      if (fin->GetBinError(b) > 0 && fin->GetBinError(b) < TMath::Abs(fin->GetBinContent(b)))
        ymax = TMath::Max(ymax, fin->GetBinContent(b) + fin->GetBinError(b));
    fin->SetMinimum(0.); fin->SetMaximum(1.6 * ymax);  // headroom for the text box (layout)
    fin->Draw("E1");
    TLine* one2 = new TLine(0, 1, kMax, 1); one2->SetLineStyle(2); one2->Draw();
    TPaveText* t = new TPaveText(0.40, 0.70, 0.95, 0.89, "NDC");
    t->SetFillStyle(0); t->SetBorderSize(0); t->SetTextAlign(12); t->SetTextSize(0.028);
    if (unitA) t->AddText(TString::Format("C(k*) = N #frac{N_{SE}(k*) P_{B}^{SE}(k*)}{N_{ME}(k*) P_{B}^{ME}(k*)}   (P_{%s} = 1 by definition)", aName));
    else t->AddText("C(k*) = N #frac{N_{SE}(k*) P_{A}^{SE}(k*) P_{B}^{SE}(k*)}{N_{ME}(k*) P_{A}^{ME}(k*) P_{B}^{ME}(k*)}");
    t->AddText(TString::Format("N: C = 1 in k* %.1f-%.1f GeV/c;  SE pairs in window: %.0f", cfg.cfNormMin, cfg.cfNormMax, hN[0]->Integral()));
    t->AddText(TString::Format("integrated SE purity: %s%s %.3f, %s %.3f", aName, unitA ? " (not used)" : "", k0Int[0].purity, bName, xiInt[0].purity));
    if (fin->GetNbinsX() != nK || nDropped > 0)
      t->AddText(TString::Format("k* bins with failed purity fits: merged (wide points), %d dropped", nDropped));
    t->Draw();
    c->Print(pdf);
  }
  c->Print(pdf + "]");

  TFile* fo = TFile::Open(TString(outPrefix) + ".root", "RECREATE");
  if (!fo || fo->IsZombie()) { std::cerr << "cannot create " << outPrefix << ".root" << std::endl; return 1; }
  for (Int_t s = 0; s < 2; ++s) { hP[s][0]->Write(); hP[s][1]->Write(); hN[s]->Write(); }
  cfRaw->Write(); cfConstant->Write(); cfPurity->Write(); cfSliding->Write(); cfMerged->Write();
  for (Int_t s = 0; s < 2; ++s) { hPuritySliding[s][0]->Write(); hPuritySliding[s][1]->Write(); hSumSliding[s]->Write(); }
  // inputs (rebinned k*) and all fitted projections, for redrawing without rerunning
  hS->Write("inputTH3_SE"); hM->Write("inputTH3_ME");
  fo->mkdir("fits")->cd();
  for (size_t i = 0; i < keepFit.size(); ++i) keepFit[i]->Write();
  fo->cd(); fo->mkdir("fits_sliding")->cd();
  for (size_t i = 0; i < keepSliding.size(); ++i) keepSliding[i]->Write();
  fo->cd(); fo->mkdir("fits_merged")->cd();
  for (size_t i = 0; i < mergedProj.size(); ++i) mergedProj[i]->Write();
  fo->cd();
  fo->Close();
  std::cout << "Wrote " << pdf << " and " << outPrefix << ".root" << std::endl;
  std::cout << "fit validity: invalid k* bins " << nInvalidBins << ", sliding-window fallbacks " << nSlidingFallbacks
            << ", sliding-window invalid coarse bins " << nSlidingInvalid << std::endl;
  std::cout << "integrated purities" << (unitA ? " (A not used: purity 1 by definition)" : "") << ": SE A " << k0Int[0].purity
            << " B " << xiInt[0].purity << " | ME A " << k0Int[1].purity << " B " << xiInt[1].purity << std::endl;
  return 0;
}
