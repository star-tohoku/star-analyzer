// checkHistCandidateFemto.C - QA / correlation-function PDF for StCandidateFemtoMaker output.
//
// Species and channels are read from the candidateFemto YAML of the mainconf
// (CandidateFemtoConfig.h), pair histograms are the THnSparseF hPairSE_<channel> / hPairME_<channel>
// (axes: 0 k*, 1 cent9, 2 mass class A, 3 mass class B, 4 mass A, 5 mass B, 6 score A, 7 score B;
// mass class 0 none, 1 signal, 2 left sideband, 3 right sideband), projected with
// CandidateFemtoProjection.h.
//
// Physics computations:
//   CF(bin) = scale * SE(bin) / ME(bin),  scale = Int(ME) / Int(SE) over k* in [cfNormMin, cfNormMax],
//   CF bin width qaCfKstarBinWidth of the YAML, Poisson errors;
//   individual purities from a fit of the mass spectrum of this species in the channel's SE pairs whose partner is
//   in its signal window (hPairSE_<channel>, both score thresholds; not hChannelMassVsCent9_*,
//   which hold every candidate of the events with both species without the partner's threshold); the species pages
//   show the per-species spectra
//   hCandidateMass_<s> / hCandidateMassVsCent9_<s> (every usable candidate) for information
//   (P = Gaussian / total integral over the signal window; local range and background order:
//   species_<s>_qaPurityFitSigmaRange / qaPurityBackgroundOrder of the YAML);
//   4-term purity-corrected CF
//     C_corr = [C_sig - (1-P1) P2 C_SB1 - P1 (1-P2) C_SB2 - (1-P1)(1-P2) C_SB12] / (P1 P2)
//   with C_SB1 = A sideband x B signal, C_SB2 = A signal x B sideband, C_SB12 = A sideband x B sideband
//   (left/right sideband CFs averaged bin by bin); SE/ME counts summed over the cent9 bins of a
//   centrality group before any ratio or fit.
// A species with purityByDefinition: true (e.g. K+) has P = 1; its sideband terms then have weight 0.
// BDT score thresholds: species_<s>_scoreMin of the candidateFemto YAML (no cut if absent). With a
// threshold, the channel purities (integrated and per group) use the pairs with score >= scoreMin of both species
// (thresholds must be score edges of hPairSE; otherwise, or without hPairSE, a "FALLBACK" to the per-species spectra).
//
// Usage (ROOT5, via analysis/run_checkHistCandidateFemto.C or script/singularity_checkHistCandidateFemto.sh):
//   checkHistCandidateFemto("out.root", "config/mainconf/main_....yaml", "qa.pdf", "8;6,7;4,5")  (all required)
// Writes qa.pdf (multi-page, landscape) and qa.root (all CFs, per channel directory).

#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TKey.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2.h>
#include <TH2D.h>
#include <TH3.h>
#include <TProfile.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TLine.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TLatex.h>
#include <TBox.h>
#include <TMath.h>
#include <TNamed.h>
#include <TString.h>
#include <TStyle.h>
#include <THnSparse.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <iostream>
#include <vector>
#include <string>
#include <cstring>

#include "PdfIOMan.h"
#include "CandidateFemtoConfig.h"
#include "CandidateFemtoProjection.h"

namespace {

// From the candidateFemto YAML (set at the start of checkHistCandidateFemto): CF normalisation range
// (cfNormMin / cfNormMax), CF k* bin width (qaCfKstarBinWidth), k* display range (kstarAxis max),
// close-pair |delta eta| (qaClosePairMaxDeltaEta).
Double_t gNormMin = -1.;
Double_t gNormMax = -1.;
Double_t gCfKstarBinWidth = -1.;
Double_t gDisplayMax = -1.;
Double_t gClosePairMaxDeltaEta = -1.;
Int_t gMissingHists = 0;

// y range of a CF plot from its points measured better than 100% (content +- error inside the k* display
// range), 10% margin; the lower edge stays at 0 unless negative values are allowed (purity-corrected CFs).
void setCfRange(TH1* h, Bool_t allowNegative) {
  Double_t lo = 1e30, hi = -1e30;
  for (Int_t b = 1; b <= h->GetNbinsX(); ++b) {
    if (h->GetBinError(b) <= 0. || h->GetBinError(b) >= TMath::Abs(h->GetBinContent(b)) ||
        h->GetXaxis()->GetBinLowEdge(b) >= gDisplayMax) continue;
    lo = TMath::Min(lo, h->GetBinContent(b) - h->GetBinError(b));
    hi = TMath::Max(hi, h->GetBinContent(b) + h->GetBinError(b));
  }
  if (lo > hi) { lo = 0.; hi = 2.; }
  hi = TMath::Max(hi, 1.);
  lo = allowNegative ? TMath::Min(lo, 0.) : 0.;
  const Double_t margin = 0.1 * (hi - lo);
  h->SetMinimum(allowNegative ? lo - margin : 0.);
  h->SetMaximum(hi + margin);
}

// cent9 -> centrality percentile (StRefMultCorr cent9 convention, larger = more central; a specification).
const char* kCent9Pct[9] = {"70-80%", "60-70%", "50-60%", "40-50%", "30-40%",
                            "20-30%", "10-20%", "5-10%", "0-5%"};
const Int_t kCent9PctLo[9] = {70, 60, 50, 40, 30, 20, 10, 5, 0};
const Int_t kCent9PctHi[9] = {80, 70, 60, 50, 40, 30, 20, 10, 5};

const char* kClassLabel[4] = {"none", "sig", "L-SB", "R-SB"};
const char* kClassTag[4] = {"none", "sig", "lsb", "rsb"};

// ---------------------------------------------------------------------------------------------
// Per-species information that is not in the candidateFemto YAML: histogram names of the
// reconstruction source (KfParticleCandidateSource QA, one charge per histogram), labels and daughter
// names (same table as script/generate_hist_anaCandidateFemto.py). Purity spectra always come from the
// Maker's per-species histograms (hCandidateMass*_<s>), never from the source QA histograms.
struct SpeciesInfo {
  const CandidateFemtoSpecies* cfg;
  TString key, label, massTitle;
  TString rawHist;                            // source finder-output mass histogram ("" = none)
  std::vector<TString> topology;              // source topology QA histograms
  std::vector<TString> daughters;
  Double_t localSigmaMultiplier;
  Int_t bgOrder;
  Bool_t hasScores;     // hCandidateScore_<s> has entries (BDT used)
  Bool_t scoreCut;      // scoreMin applied to the CFs
  // integrated purity
  TH1* hPurityMass;
  TF1* fit;
  Double_t purity;
  TString puritySource;
  // centrality groups
  std::vector<TH1D*> grpMass;
  std::vector<TF1*> grpFit;
  std::vector<Double_t> grpPurity;
};

void fillSpeciesInfo(const CandidateFemtoSpecies& s, SpeciesInfo& info) {
  info.cfg = &s;
  info.key = s.key.c_str();
  info.label = info.key;
  info.massTitle = "mass variable";
  info.localSigmaMultiplier = s.qaPurityFitSigmaRange;  // 0 / -1 if purityByDefinition (no fit)
  info.bgOrder = s.qaPurityBackgroundOrder;
  info.hasScores = kFALSE;
  info.scoreCut = kFALSE;
  info.hPurityMass = 0;
  info.fit = 0;
  info.purity = -1.;
  const std::string& pk = s.particleKey;
  if (pk == "k0s_kf") {
    info.label = "K^{0}_{S}";
    info.massTitle = "M(#pi^{+}#pi^{-}) [GeV/c^{2}]";
    info.rawHist = "hKfK0shortMassRaw";
    const char* t[4] = {"hK0short_DCA12", "hK0short_DCAV0", "hK0short_CosPointing", "hK0short_Pt"};
    for (Int_t i = 0; i < 4; ++i) info.topology.push_back(t[i]);
    info.daughters.push_back("pionPlus");
    info.daughters.push_back("pionMinus");
  } else if (pk == "xi_minus_kf") {
    info.label = "#Xi^{-}";
    info.massTitle = "M(#Lambda#pi^{-}) [GeV/c^{2}]";
    info.rawHist = "hKfXiMassRaw";
    const char* t[4] = {"hXiMinus_DCA12", "hXiMinus_DCAV0", "hXiMinus_CosPointing", "hXiMinus_Pt"};
    for (Int_t i = 0; i < 4; ++i) info.topology.push_back(t[i]);
    info.daughters.push_back("proton");
    info.daughters.push_back("pion");
    info.daughters.push_back("bachelor");
  } else if (pk == "xi_plus_kf") {
    info.label = "#bar{#Xi}^{+}";
    info.massTitle = "M(#bar{#Lambda}#pi^{+}) [GeV/c^{2}]";
    info.rawHist = "hKfAntiXiMassRaw";
    const char* t[4] = {"hXiPlus_DCA12", "hXiPlus_DCAV0", "hXiPlus_CosPointing", "hXiPlus_Pt"};
    for (Int_t i = 0; i < 4; ++i) info.topology.push_back(t[i]);
    info.daughters.push_back("proton");
    info.daughters.push_back("pion");
    info.daughters.push_back("bachelor");
  } else {
    if (pk == "kaon_plus_tof") info.label = "K^{+}";
    if (s.source == "track_pid") info.massTitle = "m^{2}_{TOF} [GeV^{2}/c^{4}]";
    info.daughters.push_back("track");
  }
}

// ---------------------------------------------------------------------------------------------
// Drawing helpers.

void drawWindowLine(Double_t x, Int_t color, Int_t style = 2) {
  if (!gPad || x < -1.0e8) return;
  gPad->Update();
  TLine* l = new TLine(x, gPad->GetUymin(), x, gPad->GetUymax());
  l->SetLineColor(color);
  l->SetLineStyle(style);
  l->SetLineWidth(1);
  l->Draw("same");
}

// Signal window = blue solid, sidebands = red dashed (sidebands first).
void drawSpeciesWindows(const CandidateFemtoSpecies& s) {
  if (s.leftSideband.defined) {
    drawWindowLine(s.leftSideband.min, kRed + 1);
    drawWindowLine(s.leftSideband.max, kRed + 1);
  }
  if (s.rightSideband.defined) {
    drawWindowLine(s.rightSideband.min, kRed + 1);
    drawWindowLine(s.rightSideband.max, kRed + 1);
  }
  drawWindowLine(s.signal.min, kBlue + 1, 1);
  drawWindowLine(s.signal.max, kBlue + 1, 1);
}

// Horizontal window lines (mass variable on the y axis), over [x0, x1].
void drawSpeciesWindowsY(const CandidateFemtoSpecies& s, Double_t x0, Double_t x1) {
  const Double_t y[2] = {s.signal.min, s.signal.max};
  for (Int_t i = 0; i < 2; ++i) {
    TLine* l = new TLine(x0, y[i], x1, y[i]);
    l->SetLineColor(kRed);
    l->Draw();
  }
}

void drawNormBand() {
  if (!gPad || gNormMax <= gNormMin) return;
  gPad->Update();
  TBox* b = new TBox(gNormMin, gPad->GetUymin(), gNormMax, gPad->GetUymax());
  b->SetFillColor(kGray + 1);
  b->SetFillStyle(3004);  // hatched: ROOT5 PDF has no alpha
  b->SetLineWidth(0);
  b->Draw("same");
  TLine* one = new TLine(0., 1., gDisplayMax, 1.);
  one->SetLineStyle(3);
  one->SetLineColor(kGray + 2);
  one->Draw("same");
  gPad->RedrawAxis();
}

void drawPadText(Double_t x, Double_t y, const char* text, Double_t size = 0.045) {
  TLatex* t = new TLatex();
  t->SetNDC();
  t->SetTextSize(size);
  t->DrawLatex(x, y, text);
}

// studentT (purityModel invariantMassStudentTPol2): pass 2 of FitMassSpectrum uses a Student-t peak
// [0]*(1+((x-[1])/[2])^2/[3])^(-([3]+1)/2) (parameter 3 named "nu") instead of the Gaussian; the K0S and Xi- peaks
// have non-Gaussian tails that a Gaussian counts as background. Parameters 0-2 keep their meaning.
Bool_t IsStudentT(const TF1* fit) { return fit && fit->GetNpar() > 3 && TString(fit->GetParName(3)) == "nu"; }
Int_t NPeakPar(const TF1* fit) { return IsStudentT(fit) ? 4 : 3; }
const char* kStudentTPeak = "[0]*TMath::Power(1.+((x-[1])/[2])*((x-[1])/[2])/[3],-([3]+1.)/2.)";

void drawMassLegend(TH1* h, TF1* fit) {
  TLegend* leg = new TLegend(0.60, 0.62, 0.89, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.035);
  if (h) leg->AddEntry(h, "data", "l");
  if (fit) leg->AddEntry(fit, IsStudentT(fit) ? "Student-t + pol fit" : "gaus + pol fit", "l");
  TLine* lSig = new TLine(); lSig->SetLineColor(kBlue + 1); lSig->SetLineStyle(1);
  TLine* lSb = new TLine();  lSb->SetLineColor(kRed + 1);   lSb->SetLineStyle(2);
  TLine* lFit = new TLine(); lFit->SetLineColor(kOrange + 7); lFit->SetLineStyle(1);
  leg->AddEntry(lSig, "signal window (config)", "l");
  leg->AddEntry(lSb, "sidebands (config)", "l");
  if (fit) leg->AddEntry(lFit, "fit mean #pm 3#sigma", "l");
  leg->Draw();
}

void drawFitSigmaLines(TF1* fit, Double_t nSigma) {
  if (!fit) return;
  const Double_t mean = fit->GetParameter(1);
  const Double_t sigma = TMath::Abs(fit->GetParameter(2));
  drawWindowLine(mean - nSigma * sigma, kOrange + 7, 1);
  drawWindowLine(mean + nSigma * sigma, kOrange + 7, 1);
}

TObject* getObjOrWarn(TFile* fin, const char* name) {
  TObject* o = fin->Get(name);
  if (!o) {
    ++gMissingHists;
    std::cerr << "[checkHistCandidateFemto] MISSING: " << name << std::endl;
  }
  return o;
}
TH1* getOrWarn(TFile* fin, const char* name) { return (TH1*)getObjOrWarn(fin, name); }

void zoomToFilled(TH1* h) {
  if (!h || h->GetEntries() <= 0) return;
  Int_t lo = -1, hi = -1;
  for (Int_t b = 1; b <= h->GetNbinsX(); ++b) {
    if (h->GetBinContent(b) <= 0) continue;
    if (lo < 0) lo = b;
    hi = b;
  }
  if (lo < 0) return;
  const Int_t m = TMath::Max(1, (hi - lo) / 20);
  h->GetXaxis()->SetRange(TMath::Max(1, lo - m), TMath::Min(h->GetNbinsX(), hi + m));
}

void acceptanceUnionRange(TFile* fin, const std::vector<TString>& names,
                          Double_t& xMin, Double_t& xMax, Double_t& yMin, Double_t& yMax) {
  xMin = yMin = 1.e30; xMax = yMax = -1.e30;
  for (size_t k = 0; k < names.size(); ++k) {
    TH2* h = (TH2*)fin->Get(names[k]);
    if (!h || h->GetEntries() <= 0) continue;
    for (Int_t ix = 1; ix <= h->GetNbinsX(); ++ix)
      for (Int_t iy = 1; iy <= h->GetNbinsY(); ++iy) {
        if (h->GetBinContent(ix, iy) <= 0) continue;
        const Int_t mx = 2;
        xMin = TMath::Min(xMin, h->GetXaxis()->GetBinLowEdge(TMath::Max(1, ix - mx)));
        xMax = TMath::Max(xMax, h->GetXaxis()->GetBinUpEdge(TMath::Min(h->GetNbinsX(), ix + mx)));
        yMin = TMath::Min(yMin, h->GetYaxis()->GetBinLowEdge(TMath::Max(1, iy - mx)));
        yMax = TMath::Max(yMax, h->GetYaxis()->GetBinUpEdge(TMath::Min(h->GetNbinsY(), iy + mx)));
      }
  }
}

// SE/ME ratio of a (dphi*, deta, k*) TH3 projected onto (dphi*, deta) for k* bins [zLo, zHi],
// each normalised to unit integral. Returns 0 if empty.
TH2D* closePairRatio2D(TH3* se, TH3* me, Int_t zLo, Int_t zHi, const char* name) {
  if (!se || !me) return 0;
  se->GetZaxis()->SetRange(zLo, zHi); me->GetZaxis()->SetRange(zLo, zHi);
  TH2D* s = (TH2D*)se->Project3D("yx")->Clone(TString::Format("%s_se", name));
  TH2D* m = (TH2D*)me->Project3D("yx")->Clone(TString::Format("%s_me", name));
  se->GetZaxis()->SetRange(0, 0); me->GetZaxis()->SetRange(0, 0);
  if (s->Integral() <= 0 || m->Integral() <= 0) return 0;
  s->Rebin2D(4, 4); m->Rebin2D(4, 4);
  s->Scale(1. / s->Integral()); m->Scale(1. / m->Integral());
  TH2D* r = (TH2D*)s->Clone(name);
  r->Divide(m);
  r->SetDirectory(0);
  return r;
}

// SE/ME ratio vs dphi* for |deta| < detaMax, k* bins [zLo, zHi], normalised to unit integral.
TH1D* closePairRatioDPhi(TH3* se, TH3* me, Double_t detaMax, Int_t zLo, Int_t zHi, const char* name) {
  if (!se || !me) return 0;
  const Int_t yLo = se->GetYaxis()->FindBin(-detaMax + 1e-6), yHi = se->GetYaxis()->FindBin(detaMax - 1e-6);
  TH1D* s = se->ProjectionX(TString::Format("%s_se", name), yLo, yHi, zLo, zHi);
  TH1D* m = me->ProjectionX(TString::Format("%s_me", name), yLo, yHi, zLo, zHi);
  if (s->Integral() <= 0 || m->Integral() <= 0) return 0;
  s->Rebin(4); m->Rebin(4);
  s->Scale(1. / s->Integral()); m->Scale(1. / m->Integral());
  TH1D* r = (TH1D*)s->Clone(name);
  r->Divide(m);
  r->SetDirectory(0);
  return r;
}

void drawAcceptanceMap(TFile* fin, const char* name, Double_t xMin, Double_t xMax, Double_t yMin, Double_t yMax) {
  TH2* h = (TH2*)getOrWarn(fin, name);
  if (!h) return;
  if (xMin < xMax && yMin < yMax) {
    h->GetXaxis()->SetRangeUser(xMin, xMax);
    h->GetYaxis()->SetRangeUser(yMin, yMax);
  }
  h->SetStats(kFALSE);
  h->SetTitle(TString::Format("%s (entries = %.0f)", h->GetTitle(), h->GetEntries()));
  gPad->SetRightMargin(0.15);
  h->Draw("COLZ");
}

void enlargeTitle(Double_t size = 0.07) {
  if (!gPad) return;
  gPad->Update();
  TPaveText* t = (TPaveText*)gPad->GetPrimitive("title");
  if (t) t->SetTextSize(size);
}

// ---------------------------------------------------------------------------------------------
// CF machinery.

Bool_t cfBinValid(const TH1* h, Int_t b) { return h && h->GetBinError(b) > 0.; }

// scale*SE/ME with simple Poisson errors; empty bins (SE or ME = 0) are marked invalid
// (content 0, error 0). SE/ME are rebinned to the CF bin width (qaCfKstarBinWidth) first (clones).
TH1D* buildCF(TH1* hSEIn, TH1* hMEIn, const char* name) {
  if (!hSEIn || !hMEIn) return 0;
  TH1D* hSE = (TH1D*)hSEIn->Clone(TString(name) + "_se_rebin");
  TH1D* hME = (TH1D*)hMEIn->Clone(TString(name) + "_me_rebin");
  hSE->SetDirectory(0);
  hME->SetDirectory(0);
  const Int_t rebin = TMath::Max(1, TMath::Nint(gCfKstarBinWidth / hSE->GetXaxis()->GetBinWidth(1)));
  if (rebin > 1) {
    hSE->Rebin(rebin);
    hME->Rebin(rebin);
  }
  const Int_t nb = hSE->GetNbinsX();
  TH1D* cf = new TH1D(name, name, nb, hSE->GetXaxis()->GetXmin(), hSE->GetXaxis()->GetXmax());
  cf->SetDirectory(0);
  const Int_t bLo = hSE->GetXaxis()->FindBin(gNormMin + 1.0e-6);
  const Int_t bHi = hSE->GetXaxis()->FindBin(gNormMax - 1.0e-6);
  const Double_t seNorm = hSE->Integral(bLo, bHi);
  const Double_t meNorm = hME->Integral(bLo, bHi);
  const Double_t scale = (seNorm > 0.) ? meNorm / seNorm : 0.;
  for (Int_t b = 1; b <= nb; ++b) {
    const Double_t se = hSE->GetBinContent(b);
    const Double_t me = hME->GetBinContent(b);
    if (se <= 0. || me <= 0.) { cf->SetBinContent(b, 0.); cf->SetBinError(b, 0.); continue; }
    const Double_t val = scale * se / me;
    cf->SetBinContent(b, val);
    cf->SetBinError(b, val * TMath::Sqrt(1. / se + 1. / me));
  }
  cf->GetXaxis()->SetRangeUser(0., gDisplayMax);
  setCfRange(cf, kFALSE);
  cf->SetTitle(TString(name) + ";k* [GeV/c];C(k*)");
  cf->SetStats(kFALSE);
  delete hSE;
  delete hME;
  return cf;
}

// Two-pass fit: pass 1 gaus+pol4 over the full
// histogram range, seeded at expectedMean/expectedSigma with bounded mean/sigma; pass 2
// gaus+pol[bgOrder] in [mean0 -+ localSigmaMultiplier*sigma0]. Returns the pass-2 TF1 (attached
// to the histogram), or 0 if the histogram is missing or empty.
TF1* FitMassSpectrum(TH1* hMass, const char* fitName, Double_t localSigmaMultiplier, Int_t bgOrder,
                     Double_t expectedMean, Double_t expectedSigma, Bool_t studentT) {
  if (!hMass) return 0;
  if (hMass->GetEntries() <= 0) {  // nothing to fit
    std::cout << "[FitMassSpectrum] " << fitName << ": empty histogram, no fit" << std::endl;
    return 0;
  }
  const Double_t fullMin = hMass->GetXaxis()->GetXmin();
  const Double_t fullMax = hMass->GetXaxis()->GetXmax();

  TF1 pre("pre_tmp_fit", "gaus(0)+pol4(3)", fullMin, fullMax);
  pre.SetParameter(0, hMass->GetMaximum());
  pre.SetParameter(1, expectedMean);
  pre.SetParameter(2, expectedSigma);
  pre.SetParLimits(1, expectedMean - 5. * expectedSigma, expectedMean + 5. * expectedSigma);
  pre.SetParLimits(2, 0.2 * expectedSigma, 3. * expectedSigma);
  hMass->Fit(&pre, "RQS0");
  const Double_t mean0 = pre.GetParameter(1);
  const Double_t sigma0 = TMath::Abs(pre.GetParameter(2));
  std::cout << "[FitMassSpectrum] " << fitName << " PASS1: mean0=" << mean0 << " sigma0=" << sigma0
            << " expectedMean=" << expectedMean << " expectedSigma=" << expectedSigma << std::endl;

  const Double_t lo = mean0 - localSigmaMultiplier * sigma0;
  const Double_t hi = mean0 + localSigmaMultiplier * sigma0;
  TString formula = studentT ? TString::Format("%s+pol%d(4)", kStudentTPeak, bgOrder)
                             : TString::Format("gaus(0)+pol%d(3)", bgOrder);
  TF1* fit = new TF1(fitName, formula.Data(), lo, hi);
  fit->SetParameter(0, pre.GetParameter(0));
  fit->SetParameter(1, mean0);
  fit->SetParameter(2, sigma0);
  if (studentT) {
    fit->SetParName(3, "nu");
    fit->SetParLimits(0, 0., 10. * hMass->GetMaximum() + 10.);
    fit->SetParLimits(1, mean0 - 2. * sigma0, mean0 + 2. * sigma0);
    fit->SetParameter(2, 0.8 * sigma0);  // Student-t scale < Gaussian sigma for the same width
    fit->SetParLimits(2, 0.2 * sigma0, 3. * sigma0);
    fit->SetParameter(3, 5.);
    fit->SetParLimits(3, 0.5, 200.);
    for (Int_t i = 4; i < fit->GetNpar(); ++i) fit->SetParameter(i, pre.GetParameter(i - 1));
  }
  hMass->Fit(fit, "RQS");
  std::cout << "[FitMassSpectrum] " << fitName << ": mean=" << fit->GetParameter(1)
            << " sigma=" << TMath::Abs(fit->GetParameter(2))
            << " chi2/ndf=" << (fit->GetNDF() > 0 ? fit->GetChisquare() / fit->GetNDF() : -1.)
            << " range=[" << lo << "," << hi << "]" << std::endl;
  return fit;
}

// P = Gaussian / full fit integral over [winMin, winMax]; -1 if unavailable or not in (0, 1].
Double_t PurityFromFit(TF1* fit, Double_t winMin, Double_t winMax) {
  if (!fit) return -1.;
  TF1 gausOnly("gausOnly_tmp_purity", IsStudentT(fit) ? kStudentTPeak : "gaus(0)", winMin, winMax);
  for (Int_t i = 0; i < NPeakPar(fit); ++i) gausOnly.SetParameter(i, fit->GetParameter(i));
  const Double_t sig = gausOnly.Integral(winMin, winMax);
  const Double_t total = fit->Integral(winMin, winMax);
  if (total <= 0.) return -1.;
  const Double_t purity = sig / total;
  if (purity <= 0. || purity > 1.) return -1.;
  return purity;
}

void drawFitComponents(TF1* fit) {
  if (!fit) return;
  const Double_t xMin = fit->GetXmin();
  const Double_t xMax = fit->GetXmax();
  const Int_t nPeak = NPeakPar(fit), nBgPar = fit->GetNpar() - nPeak;
  TF1* peakOnly = new TF1(TString(fit->GetName()) + "_peakOnly", IsStudentT(fit) ? kStudentTPeak : "gaus(0)", xMin, xMax);
  for (Int_t i = 0; i < nPeak; ++i) peakOnly->SetParameter(i, fit->GetParameter(i));
  peakOnly->SetLineColor(kBlue + 1);
  peakOnly->SetLineStyle(2);
  peakOnly->Draw("same");
  TF1* bgOnly = new TF1(TString(fit->GetName()) + "_bgOnly", TString::Format("pol%d(0)", nBgPar - 1).Data(), xMin, xMax);
  for (Int_t i = 0; i < nBgPar; ++i) bgOnly->SetParameter(i, fit->GetParameter(nPeak + i));
  bgOnly->SetLineColor(kGreen + 2);
  bgOnly->SetLineStyle(2);
  bgOnly->Draw("same");
}

// Purity at mean +/- {2, 2.5, 3} sigma (fit and sideband-integral method), diagnostic only.
void addSigmaWindowLines(TPaveText* note, TF1* fit, TH1* hMass, const CandidateFemtoSpecies& s) {
  if (!fit) return;
  const Double_t mean = fit->GetParameter(1);
  const Double_t sigma = TMath::Abs(fit->GetParameter(2));
  const Double_t nSigmas[3] = {2.0, 2.5, 3.0};
  const Bool_t haveSb = s.leftSideband.defined && s.rightSideband.defined;
  for (Int_t i = 0; i < 3; ++i) {
    const Double_t lo = mean - nSigmas[i] * sigma;
    const Double_t hi = mean + nSigmas[i] * sigma;
    const Double_t pFit = PurityFromFit(fit, lo, hi);
    TString line = TString::Format("%.1f#sigma [%.4f,%.4f]: fit P=%.3f", nSigmas[i], lo, hi, pFit);
    if (hMass && haveSb) {
      const Double_t sbLeftMin = s.leftSideband.min, sbLeftMax = s.leftSideband.max;
      const Double_t sbRightMin = s.rightSideband.min, sbRightMax = s.rightSideband.max;
      const Double_t sig = hMass->Integral(hMass->GetXaxis()->FindBin(lo + 1.0e-6), hMass->GetXaxis()->FindBin(hi - 1.0e-6));
      const Double_t bl = hMass->Integral(hMass->GetXaxis()->FindBin(sbLeftMin + 1.0e-6),
                                          hMass->GetXaxis()->FindBin(sbLeftMax - 1.0e-6));
      const Double_t br = hMass->Integral(hMass->GetXaxis()->FindBin(sbRightMin + 1.0e-6),
                                          hMass->GetXaxis()->FindBin(sbRightMax - 1.0e-6));
      const Double_t sigW = hi - lo;
      const Double_t sbW = (sbLeftMax - sbLeftMin) + (sbRightMax - sbRightMin);
      if (sig > 0. && sbW > 0.) {
        const Double_t bNorm = (bl + br) * (sigW / sbW);
        const Double_t pSb = (sig - bNorm) / sig;
        line += TString::Format("  sideband-method P=%.3f", pSb);
      }
    }
    note->AddText(line.Data());
  }
}

// CF subtraction with one purity (used only as the centrality-group fallback when
// the 4-term correction is unavailable): (CF_sig - (1-P) CF_SB) / P.
TH1D* subtractPurity(TH1D* cfSig, TH1D* cfSB, Double_t purity, const char* name) {
  if (!cfSig || !cfSB || purity <= 0.) return 0;
  TH1D* out = (TH1D*)cfSig->Clone(name);
  out->SetDirectory(0);
  const Int_t nb = out->GetNbinsX();
  for (Int_t b = 1; b <= nb; ++b) {
    const Double_t sig = cfSig->GetBinContent(b);
    const Double_t sb = cfSB->GetBinContent(b);
    if (!cfBinValid(cfSig, b) || !cfBinValid(cfSB, b)) { out->SetBinContent(b, 0.); out->SetBinError(b, 0.); continue; }
    const Double_t val = (sig - (1. - purity) * sb) / purity;
    const Double_t errSig = cfSig->GetBinError(b);
    const Double_t errSB = cfSB->GetBinError(b);
    const Double_t err = TMath::Sqrt(errSig * errSig + (1. - purity) * (1. - purity) * errSB * errSB) / purity;
    out->SetBinContent(b, val);
    out->SetBinError(b, err);
  }
  out->SetTitle(TString(name) + ";k* [GeV/c];C(k*) (purity-corrected)");
  return out;
}

// Full 4-term purity correction (side-band subtraction with two purities):
//   C_corr = [C_sig - (1-P1) P2 C_SB1 - P1 (1-P2) C_SB2 - (1-P1)(1-P2) C_SB12] / (P1 P2)
// cfSB1 = A SB x B sig, cfSB2 = A sig x B SB, cfSB12 = A SB x B SB. Generalisation: a background
// CF may be missing (no sideband class for that species) only if its weight is exactly 0 (P = 1);
// every CF that is given must be valid in a bin.
TH1D* subtractPurity4Term(TH1D* cfSig, TH1D* cfSB1, TH1D* cfSB2, TH1D* cfSB12, Double_t p1, Double_t p2,
                          const char* name) {
  if (!cfSig || p1 <= 0. || p2 <= 0.) return 0;
  const Double_t p1p2 = p1 * p2;
  const Double_t w1 = (1. - p1) * p2;
  const Double_t w2 = p1 * (1. - p2);
  const Double_t w12 = (1. - p1) * (1. - p2);
  if ((!cfSB1 && w1 != 0.) || (!cfSB2 && w2 != 0.) || (!cfSB12 && w12 != 0.)) return 0;
  TH1D* out = (TH1D*)cfSig->Clone(name);
  out->SetDirectory(0);
  const Int_t nb = out->GetNbinsX();
  for (Int_t b = 1; b <= nb; ++b) {
    if (!cfBinValid(cfSig, b) || (cfSB1 && !cfBinValid(cfSB1, b)) || (cfSB2 && !cfBinValid(cfSB2, b)) ||
        (cfSB12 && !cfBinValid(cfSB12, b))) {
      out->SetBinContent(b, 0.); out->SetBinError(b, 0.); continue;
    }
    const Double_t sig = cfSig->GetBinContent(b);
    const Double_t sb1 = cfSB1 ? cfSB1->GetBinContent(b) : 0.;
    const Double_t sb2 = cfSB2 ? cfSB2->GetBinContent(b) : 0.;
    const Double_t sb12 = cfSB12 ? cfSB12->GetBinContent(b) : 0.;
    const Double_t val = (sig - w1 * sb1 - w2 * sb2 - w12 * sb12) / p1p2;
    const Double_t eSig = cfSig->GetBinError(b);
    const Double_t eSB1 = cfSB1 ? cfSB1->GetBinError(b) : 0.;
    const Double_t eSB2 = cfSB2 ? cfSB2->GetBinError(b) : 0.;
    const Double_t eSB12 = cfSB12 ? cfSB12->GetBinError(b) : 0.;
    const Double_t err = TMath::Sqrt(eSig * eSig + w1 * w1 * eSB1 * eSB1 + w2 * w2 * eSB2 * eSB2 +
                                     w12 * w12 * eSB12 * eSB12) / p1p2;
    out->SetBinContent(b, val);
    out->SetBinError(b, err);
  }
  out->SetTitle(TString(name) + ";k* [GeV/c];C_{corr}(k*)");
  return out;
}

// Equal-weight average of up to 4 CFs, per bin over the valid inputs only; err = sqrt(sum e^2)/n.
TH1D* averageCFs(TH1D* a, TH1D* b, TH1D* c, TH1D* d, const char* name) {
  std::vector<TH1D*> in;
  if (a) in.push_back(a);
  if (b) in.push_back(b);
  if (c) in.push_back(c);
  if (d) in.push_back(d);
  if (in.empty()) return 0;
  TH1D* out = (TH1D*)in[0]->Clone(name);
  out->SetDirectory(0);
  out->SetTitle(TString(name) + ";k* [GeV/c];C(k*)");
  for (Int_t bin = 1; bin <= out->GetNbinsX(); ++bin) {
    Double_t sum = 0., err2 = 0.;
    Int_t n = 0;
    for (size_t i = 0; i < in.size(); ++i) {
      if (!cfBinValid(in[i], bin)) continue;
      sum += in[i]->GetBinContent(bin);
      err2 += in[i]->GetBinError(bin) * in[i]->GetBinError(bin);
      ++n;
    }
    out->SetBinContent(bin, n > 0 ? sum / n : 0.);
    out->SetBinError(bin, n > 0 ? TMath::Sqrt(err2) / n : 0.);
  }
  return out;
}

TH1D* cloneOrNull(TH1D* h, const char* name) {
  if (!h) return 0;
  TH1D* c = (TH1D*)h->Clone(name);
  c->SetDirectory(0);
  return c;
}

void drawConsolidatedCf(TH1D* cfSig, TH1D* cfSB1, TH1D* cfSB2, TH1D* cfSB12, const char* label,
                        const char* labA, const char* labB) {
  TH1D* frame = 0;
  if (cfSig)  { cfSig->SetLineColor(kBlack);       cfSig->SetMarkerColor(kBlack);       cfSig->SetMarkerStyle(20); if (!frame) frame = cfSig; }
  if (cfSB1)  { cfSB1->SetLineColor(kGreen + 2);    cfSB1->SetMarkerColor(kGreen + 2);   cfSB1->SetMarkerStyle(21); if (!frame) frame = cfSB1; }
  if (cfSB2)  { cfSB2->SetLineColor(kRed + 1);      cfSB2->SetMarkerColor(kRed + 1);     cfSB2->SetMarkerStyle(22); if (!frame) frame = cfSB2; }
  if (cfSB12) { cfSB12->SetLineColor(kMagenta + 1); cfSB12->SetMarkerColor(kMagenta + 1); cfSB12->SetMarkerStyle(23); if (!frame) frame = cfSB12; }
  if (!frame) {
    drawPadText(0.15, 0.5, TString::Format("%s: no raw CF", label).Data());
    return;
  }
  frame->SetTitle(TString::Format("raw CFs, %s;k* [GeV/c];C(k*)", label));
  {
    // common y range of all drawn CFs
    Double_t lo = 1e30, hi = -1e30;
    TH1D* all[4] = {cfSig, cfSB1, cfSB2, cfSB12};
    for (Int_t i = 0; i < 4; ++i) {
      if (!all[i]) continue;
      setCfRange(all[i], kFALSE);
      lo = TMath::Min(lo, all[i]->GetMinimum());
      hi = TMath::Max(hi, all[i]->GetMaximum());
    }
    frame->SetMinimum(lo);
    frame->SetMaximum(hi);
  }
  frame->Draw("e");
  enlargeTitle(0.06);
  drawNormBand();
  frame->Draw("e same");
  if (cfSig && cfSig != frame) cfSig->Draw("e same");
  if (cfSB1 && cfSB1 != frame) cfSB1->Draw("e same");
  if (cfSB2 && cfSB2 != frame) cfSB2->Draw("e same");
  if (cfSB12 && cfSB12 != frame) cfSB12->Draw("e same");
  TLegend* leg = new TLegend(0.50, 0.70, 0.89, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.035);
  if (cfSig) leg->AddEntry(cfSig, TString::Format("CF_{sig} (%s sig #times %s sig)", labA, labB), "lep");
  if (cfSB1) leg->AddEntry(cfSB1, TString::Format("CF_{SB1} (%s SB #times %s sig)", labA, labB), "lep");
  if (cfSB2) leg->AddEntry(cfSB2, TString::Format("CF_{SB2} (%s sig #times %s SB)", labA, labB), "lep");
  if (cfSB12) leg->AddEntry(cfSB12, "CF_{SB12} (both SB)", "lep");
  leg->Draw();
}

void drawFinalCf(TH1D* cfFinal, const char* label, Double_t p1, Double_t p2, const char* labA, const char* labB) {
  if (!cfFinal) {
    drawPadText(0.15, 0.5, TString::Format("%s: purity-corrected CF unavailable", label).Data());
    return;
  }
  cfFinal->SetLineColor(kBlue + 1);
  cfFinal->SetMarkerColor(kBlue + 1);
  cfFinal->SetMarkerStyle(34);
  cfFinal->SetTitle(TString::Format("purity-corrected CF, %s;k* [GeV/c];C_{corr}(k*)", label));
  cfFinal->GetXaxis()->SetRangeUser(0., gDisplayMax);
  setCfRange(cfFinal, kTRUE);
  cfFinal->Draw("e");
  enlargeTitle(0.06);
  drawNormBand();
  cfFinal->Draw("e same");
  drawPadText(0.40, 0.82, TString::Format("P_{1}(%s) = %.3f,  P_{2}(%s) = %.3f", labA, p1, labB, p2).Data(), 0.04);
}

// Sum of ProjectionX (k* vs cent9, y = cent9) over a group of cent9 bins before any ratio.
TH1D* projectCentGroupSumX(TH2* h2, const std::vector<Int_t>& bins, const char* name) {
  if (!h2) return 0;
  TH1D* sum = 0;
  for (size_t i = 0; i < bins.size(); ++i) {
    const Int_t yBin = h2->GetYaxis()->FindBin((Double_t)bins[i]);
    TH1D* p = h2->ProjectionX(TString::Format("%s_tmp%d", name, (Int_t)i), yBin, yBin);
    if (!sum) { sum = (TH1D*)p->Clone(name); sum->SetDirectory(0); delete p; }
    else { sum->Add(p); delete p; }
  }
  return sum;
}

// Same for mass vs cent9 histograms (x = cent9, y = mass): sum of ProjectionY.
// ---------------------------------------------------------------------------------------------
// Pair projections of one channel for one selection: k* and k* vs cent9 per mass-class pair.
struct ClassProjection {
  ClassProjection() { for (Int_t i = 0; i < 4; ++i) for (Int_t j = 0; j < 4; ++j) { se[i][j] = me[i][j] = 0; seCent[i][j] = meCent[i][j] = 0; } }
  TH1D* se[4][4];
  TH1D* me[4][4];
  TH2D* seCent[4][4];
  TH2D* meCent[4][4];
};

Bool_t classExists(const CandidateFemtoSpecies& s, Int_t cls) {
  if (cls == 1) return kTRUE;
  if (cls == 2) return s.leftSideband.defined;
  if (cls == 3) return s.rightSideband.defined;
  return kFALSE;
}

void projectClasses(THnSparse* hSE, THnSparse* hME, const CandidatePairSelection& sel, const CandidateFemtoSpecies& sa,
                    const CandidateFemtoSpecies& sb, const char* prefix, ClassProjection& out) {
  for (Int_t ca = 1; ca <= 3; ++ca)
    for (Int_t cb = 1; cb <= 3; ++cb) {
      if (!classExists(sa, ca) || !classExists(sb, cb)) continue;
      const TString n = TString::Format("%s_%s_%s", prefix, kClassTag[ca], kClassTag[cb]);
      out.se[ca][cb] = CandidateProjectKstar(hSE, sel, ca, cb, (n + "_se").Data(), &out.seCent[ca][cb]);
      out.me[ca][cb] = CandidateProjectKstar(hME, sel, ca, cb, (n + "_me").Data(), &out.meCent[ca][cb]);
    }
}

// The CF set of one channel for one cent9 selection (integrated or a group).
struct ChannelCFs {
  ChannelCFs() : bSbAvg(0), aSbAvg(0), sb12(0), final(0), bOnlyCorrected(0) {
    for (Int_t i = 0; i < 4; ++i) for (Int_t j = 0; j < 4; ++j) cf[i][j] = 0;
  }
  TH1D* cf[4][4];
  TH1D* bSbAvg;  // A sig x B SB (L+R avg)  = CF_SB2
  TH1D* aSbAvg;  // A SB (L+R avg) x B sig  = CF_SB1
  TH1D* sb12;    // A SB x B SB (corner avg)
  TH1D* final;
  TH1D* bOnlyCorrected;  // fallback: B-background-only correction
};

void buildChannelCFs(TH1D* se[4][4], TH1D* me[4][4], const char* prefix, Double_t pA, Double_t pB, ChannelCFs& out) {
  for (Int_t ca = 1; ca <= 3; ++ca)
    for (Int_t cb = 1; cb <= 3; ++cb)
      if (se[ca][cb] && me[ca][cb])
        out.cf[ca][cb] = buildCF(se[ca][cb], me[ca][cb],
                                 TString::Format("%s_cf_%s_%s", prefix, kClassTag[ca], kClassTag[cb]).Data());
  out.bSbAvg = averageCFs(out.cf[1][2], out.cf[1][3], 0, 0, TString::Format("%s_cf_sig_sbavg", prefix).Data());
  out.aSbAvg = averageCFs(out.cf[2][1], out.cf[3][1], 0, 0, TString::Format("%s_cf_sbavg_sig", prefix).Data());
  out.sb12 = averageCFs(out.cf[2][2], out.cf[2][3], out.cf[3][2], out.cf[3][3],
                        TString::Format("%s_cf_sb12avg", prefix).Data());
  out.final = subtractPurity4Term(out.cf[1][1], out.aSbAvg, out.bSbAvg, out.sb12, pA, pB,
                                  TString::Format("%s_cf_final", prefix).Data());
  if (out.cf[1][1] && out.bSbAvg && pB > 0.)
    out.bOnlyCorrected = subtractPurity(out.cf[1][1], out.bSbAvg, pB, TString::Format("%s_cf_bonlycorr", prefix).Data());
}

// Mass spectrum (x) of hCandidateMassVsCent9VsScore_<s> (mass, cent9, score) for score >= scoreMin,
// summed over the given cent9 bins (empty = all, including under/overflow). 0 if scoreMin is not a
// score bin edge.
TH1D* projectMassCentScore(TH3* h3, const std::vector<Int_t>& cent9, Double_t scoreMin, const char* name) {
  if (!h3) return 0;
  const Int_t zLo = h3->GetZaxis()->FindBin(scoreMin + 1e-9);
  if (TMath::Abs(h3->GetZaxis()->GetBinLowEdge(zLo) - scoreMin) > 1e-9) {
    std::cerr << "[checkHistCandidateFemto] ERROR: scoreMin " << scoreMin << " is not a score bin edge of "
              << h3->GetName() << std::endl;
    return 0;
  }
  const Int_t zHi = h3->GetNbinsZ() + 1;
  TH1D* sum = 0;
  if (cent9.empty()) {
    sum = h3->ProjectionX(name, 0, h3->GetNbinsY() + 1, zLo, zHi);
    sum->SetDirectory(0);
    return sum;
  }
  for (size_t i = 0; i < cent9.size(); ++i) {
    const Int_t yBin = h3->GetYaxis()->FindBin((Double_t)cent9[i]);
    TH1D* p = h3->ProjectionX(TString::Format("%s_tmp%d", name, (Int_t)i), yBin, yBin, zLo, zHi);
    if (!sum) { sum = (TH1D*)p->Clone(name); sum->SetDirectory(0); delete p; }
    else { sum->Add(p); delete p; }
  }
  return sum;
}

// Mass spectrum used for a species purity: hCandidateMass_<s> (the Maker's candidates of this species and
// charge, every usable candidate: the per-species histograms are filled before the per-channel shared-track rule); with a score cut, hCandidateMassVsCent9VsScore_<s> for
// score >= scoreMin (fallback if missing: hCandidateScoreVsMass_<s>, integrated only).
TH1* puritySpectrum(TFile* fin, SpeciesInfo& sp) {
  const CandidateFemtoSpecies& s = *sp.cfg;
  if (sp.scoreCut) {
    TH3* h3 = (TH3*)fin->Get(TString::Format("hCandidateMassVsCent9VsScore_%s", sp.key.Data()));
    if (h3) {
      TH1D* p = projectMassCentScore(h3, std::vector<Int_t>(), s.scoreMin, TString::Format("hPurityMass_%s", sp.key.Data()).Data());
      sp.puritySource = TString::Format("hCandidateMassVsCent9VsScore_%s, score >= %g", sp.key.Data(), s.scoreMin);
      return p;
    }
    std::cerr << "[checkHistCandidateFemto] WARNING: hCandidateMassVsCent9VsScore_" << sp.key
              << " missing: fallback to hCandidateScoreVsMass_" << sp.key << " (no centrality-group purities)" << std::endl;
    TH2* h2 = (TH2*)getOrWarn(fin, TString::Format("hCandidateScoreVsMass_%s", sp.key.Data()));
    if (!h2) return 0;
    const Int_t yLo = h2->GetYaxis()->FindBin(s.scoreMin + 1e-9);
    if (TMath::Abs(h2->GetYaxis()->GetBinLowEdge(yLo) - s.scoreMin) > 1e-9) {
      std::cerr << "[checkHistCandidateFemto] ERROR: scoreMin " << s.scoreMin << " of " << sp.key
                << " is not a bin edge of hCandidateScoreVsMass_" << sp.key << ": purity unavailable" << std::endl;
      return 0;
    }
    TH1D* p = h2->ProjectionX(TString::Format("hPurityMass_%s", sp.key.Data()), yLo, h2->GetNbinsY() + 1);
    p->SetDirectory(0);
    sp.puritySource = TString::Format("FALLBACK hCandidateScoreVsMass_%s, score >= %g (hCandidateMassVsCent9VsScore_%s missing)",
                                      sp.key.Data(), s.scoreMin, sp.key.Data());
    return p;
  }
  sp.puritySource = TString::Format("hCandidateMass_%s", sp.key.Data());
  return getOrWarn(fin, sp.puritySource);
}

Bool_t parseCentGroups(const char* spec, std::vector<std::vector<Int_t> >& groups, std::vector<TString>& labels,
                       std::vector<TString>& tags) {
  groups.clear(); labels.clear(); tags.clear();
  TObjArray* g = TString(spec ? spec : "").Tokenize(";");
  for (Int_t i = 0; i < g->GetEntries(); ++i) {
    TString item = ((TObjString*)g->At(i))->GetString();
    item = item.Strip(TString::kBoth);
    if (!item.Length()) continue;
    std::vector<Int_t> bins;
    TObjArray* b = item.Tokenize(",");
    for (Int_t j = 0; j < b->GetEntries(); ++j) {
      TString t = ((TObjString*)b->At(j))->GetString().Strip(TString::kBoth);
      if (!t.IsDigit() || t.Atoi() < 0 || t.Atoi() > 8) {
        std::cerr << "[checkHistCandidateFemto] ERROR: bad cent9 value '" << t << "' in centGroups '" << spec << "'" << std::endl;
        delete b; delete g;
        return kFALSE;
      }
      bins.push_back(t.Atoi());
    }
    delete b;
    if (bins.empty()) continue;
    Int_t lo = 100, hi = -1;
    TString list;
    for (size_t j = 0; j < bins.size(); ++j) {
      lo = TMath::Min(lo, kCent9PctLo[bins[j]]);
      hi = TMath::Max(hi, kCent9PctHi[bins[j]]);
      list += TString::Format("%s%d", j ? "," : "", bins[j]);
    }
    groups.push_back(bins);
    labels.push_back(TString::Format("cent9=%s (%d-%d%%)", list.Data(), lo, hi));
    tags.push_back(TString::Format("g_%d_%d", lo, hi));
  }
  delete g;
  return !groups.empty();
}

void printCF(const TH1D* h, const char* tag) {
  if (!h) { std::cout << "[CFDUMP] " << tag << " unavailable" << std::endl; return; }
  std::cout << "[CFDUMP] " << tag << ":";
  Int_t n = 0;
  for (Int_t b = 1; b <= h->GetNbinsX(); ++b) {
    if (h->GetBinError(b) <= 0.) continue;
    std::cout << TString::Format(" [%.2f]%.6g+-%.4g", h->GetXaxis()->GetBinCenter(b), h->GetBinContent(b), h->GetBinError(b)).Data();
    ++n;
  }
  if (!n) std::cout << " (no valid bin)";
  std::cout << std::endl;
}

void newPage(TCanvas* c, Int_t nx, Int_t ny) {
  c->Clear();
  if (nx * ny > 1) c->Divide(nx, ny);
  c->cd(1);
}

void textPage(TCanvas* c, const char* pdf, const char* text, Double_t size = 0.04) {
  c->Clear();
  TPaveText* t = new TPaveText(0.1, 0.35, 0.9, 0.65, "NDC");
  t->SetFillColor(0);
  t->SetBorderSize(1);
  t->SetTextSize(size);
  t->AddText(text);
  t->Draw();
  c->Print(pdf);
}

TPaveText* newNote(Double_t x1, Double_t y1, Double_t x2, Double_t y2, Double_t size) {
  TPaveText* t = new TPaveText(x1, y1, x2, y2, "NDC");
  t->SetFillColor(0);
  t->SetBorderSize(0);
  t->SetTextAlign(12);
  t->SetTextSize(size);
  return t;
}

}  // namespace

// Mass spectrum of one side of the SE pairs (mass axis 4 = A, 5 = B) whose partner is in its signal window
// (class axis 3 = B, 2 = A, class 1), with both score thresholds of sel and the given cent9 list (empty = all).
// Each candidate counts once per partner, as in the signal x signal pairs that the purity corrects.
TH1D* projectPairMass(THnSparse* hSE, Int_t side, const CandidatePairSelection& sel, const char* name) {
  const Int_t massAxis = side == 0 ? 4 : 5;
  const Int_t partnerClassAxis = side == 0 ? 3 : 2;
  const TAxis* m = hSE->GetAxis(massAxis);
  TH1D* out = m->GetXbins()->GetSize() > 0 ? new TH1D(name, "", m->GetNbins(), m->GetXbins()->GetArray())
                                           : new TH1D(name, "", m->GetNbins(), m->GetXmin(), m->GetXmax());
  out->SetDirectory(0);
  out->Sumw2();
  std::vector<Int_t> idx(hSE->GetNdimensions());
  const TAxis* centAxis = hSE->GetAxis(1);
  for (Long64_t i = 0; i < hSE->GetNbins(); ++i) {
    const Double_t w = hSE->GetBinContent(i, &idx[0]);
    if (w == 0) continue;
    if (TMath::Nint(hSE->GetAxis(partnerClassAxis)->GetBinCenter(idx[partnerClassAxis])) != 1) continue;
    if (!sel.AcceptCent(TMath::Nint(centAxis->GetBinCenter(idx[1])))) continue;
    if (hSE->GetAxis(6)->GetBinLowEdge(idx[6]) < sel.scoreMinA - 1e-9) continue;
    if (hSE->GetAxis(7)->GetBinLowEdge(idx[7]) < sel.scoreMinB - 1e-9) continue;
    const Int_t b = idx[massAxis];
    out->SetBinContent(b, out->GetBinContent(b) + w);
    out->GetSumw2()->fArray[b] += hSE->GetBinError2(i);
  }
  return out;
}

// =============================================================================================
// Purity of one species in one channel (integrated and per centrality group) from the SE pairs of the channel
// (hPairSE_<channel>): mass of this species in the pairs whose partner is in its signal window, with the score
// thresholds of both species (the channel candidate histograms hChannel* hold every candidate of
// the events with both species, without the partner's threshold, so their purity is not that of the pairs).
// S starts as a copy of the species information (whose purities remain the fallback when hPairSE is missing).
void channelPurity(THnSparse* hSE, SpeciesInfo& S, const CandidateFemtoSpecies& cs, const CandidateFemtoSpecies& sa,
                   const CandidateFemtoSpecies& sb, Int_t side, const char* channel,
                   const std::vector<std::vector<Int_t> >& groups, const std::vector<TString>& groupTag) {
  if (cs.purityByDefinition || cs.purityModel == "tofMassSquaredKaonPion" || cs.purityModel == "invariantMassDoubleGausPol2") return;  // as for the species
  const Int_t nG = groups.size();
  const TString tag = TString::Format("%s_%s", channel, S.key.Data());
  if (!hSE) {
    std::cerr << "[checkHistCandidateFemto] WARNING: hPairSE_" << channel << " missing: FALLBACK to the per-species purity of "
              << S.key << std::endl;
    S.puritySource = "FALLBACK per-species: " + S.puritySource;
    return;
  }
  CandidatePairSelection sel;
  if (sa.hasScoreMin) sel.scoreMinA = sa.scoreMin;
  if (sb.hasScoreMin) sel.scoreMinB = sb.scoreMin;
  if ((sa.hasScoreMin && !CandidateScoreThresholdIsEdge(hSE, 6, sel.scoreMinA)) ||
      (sb.hasScoreMin && !CandidateScoreThresholdIsEdge(hSE, 7, sel.scoreMinB))) {
    S.puritySource = "FALLBACK per-species (threshold not a pair score edge): " + S.puritySource;
    return;
  }
  S.hPurityMass = projectPairMass(hSE, side, sel, ("hPurityMass_" + tag).Data());
  S.puritySource = TString::Format("hPairSE_%s, partner in its signal window, %s", channel, sel.Describe().Data());
  S.fit = FitMassSpectrum(S.hPurityMass, (tag + "_mass_fit").Data(), S.localSigmaMultiplier, S.bgOrder,
                          cs.purityFitSeed.min, cs.purityFitSeed.max, cs.purityModel == "invariantMassStudentTPol2");
  S.purity = PurityFromFit(S.fit, cs.signal.min, cs.signal.max);
  S.grpMass.assign(nG, (TH1D*)0);
  S.grpFit.assign(nG, (TF1*)0);
  S.grpPurity.assign(nG, -1.);
  for (Int_t g = 0; g < nG; ++g) {
    const TString name = TString::Format("%s_mass_%s", tag.Data(), groupTag[g].Data());
    CandidatePairSelection gsel = sel;
    gsel.cent9 = groups[g];
    S.grpMass[g] = projectPairMass(hSE, side, gsel, name.Data());
    S.grpFit[g] = FitMassSpectrum(S.grpMass[g], TString::Format("%s_mass_fit_%s", tag.Data(), groupTag[g].Data()).Data(),
                                  S.localSigmaMultiplier, S.bgOrder, cs.purityFitSeed.min, cs.purityFitSeed.max, cs.purityModel == "invariantMassStudentTPol2");
    S.grpPurity[g] = PurityFromFit(S.grpFit[g], cs.signal.min, cs.signal.max);
  }
}

Int_t checkHistCandidateFemto(const char* rootFile, const char* mainconf, const char* outPdf,
                             const char* centGroups) {
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(10);
  gStyle->SetStatW(0.18);
  gStyle->SetStatH(0.08);
  gMissingHists = 0;
  TH1::AddDirectory(kTRUE);

  if (!rootFile || !*rootFile || !mainconf || !*mainconf || !outPdf || !*outPdf || !centGroups || !*centGroups) {
    std::cerr << "[checkHistCandidateFemto] ERROR: usage checkHistCandidateFemto(rootFile, mainconf, outPdf, centGroups)"
              << std::endl;
    return 1;
  }
  TString mainconfPath(mainconf);
  if (!mainconfPath.BeginsWith("/")) mainconfPath = TString(gSystem->WorkingDirectory()) + "/" + mainconfPath;
  std::cout << "[checkHistCandidateFemto] mainconf: " << mainconfPath.Data() << std::endl;
  CandidateFemtoConfig cfg;
  if (!cfg.LoadFromMainConfig(mainconfPath.Data())) {
    std::cerr << "[checkHistCandidateFemto] ERROR: candidateFemto config not loaded; no PDF written" << std::endl;
    return 1;
  }
  std::cout << "[checkHistCandidateFemto] candidateFemto config: " << cfg.path << std::endl;
  cfg.Print(std::cout);
  gNormMin = cfg.cfNormMin;
  gNormMax = cfg.cfNormMax;
  gCfKstarBinWidth = cfg.qaCfKstarBinWidth;
  gDisplayMax = cfg.kstarAxis.max;
  gClosePairMaxDeltaEta = cfg.qaClosePairMaxDeltaEta;

  std::vector<std::vector<Int_t> > groups;
  std::vector<TString> groupLabel, groupTag;
  if (!parseCentGroups(centGroups, groups, groupLabel, groupTag)) {
    std::cerr << "[checkHistCandidateFemto] ERROR: invalid centGroups '" << (centGroups ? centGroups : "") << "'" << std::endl;
    return 1;
  }
  const Int_t nG = groups.size();
  for (Int_t g = 0; g < nG; ++g) std::cout << "[checkHistCandidateFemto] centrality group " << g << ": " << groupLabel[g] << std::endl;

  TFile* fin = TFile::Open(rootFile);
  if (!fin || fin->IsZombie()) {
    std::cerr << "[checkHistCandidateFemto] ERROR: cannot open " << rootFile << std::endl;
    return 1;
  }
  gROOT->cd();

  // ---- species
  const Int_t nS = cfg.species.size();
  std::vector<SpeciesInfo> sp(nS);
  for (Int_t s = 0; s < nS; ++s) {
    fillSpeciesInfo(cfg.species[s], sp[s]);
    TH1* hs = (TH1*)fin->Get(TString::Format("hCandidateScore_%s", sp[s].key.Data()));
    sp[s].hasScores = hs && hs->GetEntries() > 0;
    sp[s].scoreCut = cfg.species[s].hasScoreMin;
  }

  // ---- channels: pair histograms and score thresholds
  const Int_t nC = cfg.channels.size();
  std::vector<THnSparse*> hSE(nC, (THnSparse*)0), hME(nC, (THnSparse*)0);
  std::vector<CandidatePairSelection> sel(nC);
  for (Int_t c = 0; c < nC; ++c) {
    const CandidateFemtoChannel& ch = cfg.channels[c];
    hSE[c] = (THnSparse*)getObjOrWarn(fin, TString::Format("hPairSE_%s", ch.name.c_str()));
    hME[c] = (THnSparse*)getObjOrWarn(fin, TString::Format("hPairME_%s", ch.name.c_str()));
    if (!hSE[c] || !hME[c]) {
      std::cerr << "[checkHistCandidateFemto] ERROR: pair histograms of channel " << ch.name << " missing; stop" << std::endl;
      fin->Close();
      return 1;
    }
    const CandidateFemtoSpecies& a = cfg.species[ch.indexA];
    const CandidateFemtoSpecies& b = cfg.species[ch.indexB];
    sel[c].scoreMinA = a.hasScoreMin ? a.scoreMin : -1e9;
    sel[c].scoreMinB = b.hasScoreMin ? b.scoreMin : -1e9;
    if (!CandidateScoreThresholdIsEdge(hSE[c], 6, sel[c].scoreMinA) || !CandidateScoreThresholdIsEdge(hSE[c], 7, sel[c].scoreMinB)) {
      std::cerr << "[checkHistCandidateFemto] ERROR: scoreMin of channel " << ch.name
                << " is not a score-axis bin edge (or the file has no BDT scores); stop" << std::endl;
      fin->Close();
      return 1;
    }
    std::cout << "[checkHistCandidateFemto] channel " << ch.name << " selection: " << sel[c].Describe() << std::endl;
  }

  // ---- output
  TString pdfName(outPdf);
  TString rootOutName(outPdf);
  if (rootOutName.EndsWith(".pdf")) rootOutName.Remove(rootOutName.Length() - 4);
  rootOutName += ".root";
  {
    TString dir = gSystem->DirName(pdfName);
    if (gSystem->AccessPathName(dir)) gSystem->mkdir(dir, kTRUE);
  }
  std::vector<TObject*> outObjects;     // written into the output ROOT file
  std::vector<TString> outObjectDir;    // channel directory of each object ("" = top)

  // Statistics for the header page.
  TH1* hNAcc = (TH1*)fin->Get("hN");
  TH1* hEvSel = (TH1*)fin->Get("hKfEventSelection");
  const Double_t nEvAcc = hNAcc ? hNAcc->GetEntries() : -1.;
  const Double_t nEvRead = hEvSel ? hEvSel->GetBinContent(1) : -1.;

  gStyle->SetPaperSize(29.7, 21.0);
  PdfHeader::OpenPdf(pdfName);
  {
    std::vector<std::string> inputs;
    inputs.push_back(rootFile);
    inputs.push_back(mainconfPath.Data());
    inputs.push_back(cfg.path);
    TString note = "StCandidateFemtoMaker: CFs of the channels of the candidateFemto config.\n";
    note += TString::Format("Events read = %.0f, accepted (hN) = %.0f\n", nEvRead, nEvAcc);
    for (Int_t c = 0; c < nC; ++c) {
      TH1* hs = (TH1*)fin->Get(TString::Format("hKstarSE_%s", cfg.channels[c].name.c_str()));
      TH1* hm = (TH1*)fin->Get(TString::Format("hKstarME_%s", cfg.channels[c].name.c_str()));
      note += TString::Format("%s (sig x sig): SE = %.0f, ME = %.0f; %s\n", cfg.channels[c].name.c_str(),
                              hs ? hs->Integral(0, hs->GetNbinsX() + 1) : -1., hm ? hm->Integral(0, hm->GetNbinsX() + 1) : -1.,
                              sel[c].Describe().Data());
    }
    note += TString::Format("CF = scale*SE/ME, scale = Int(ME)/Int(SE) over k* in [%.2f, %.2f] GeV/c, bin %.2f GeV/c\n",
                            gNormMin, gNormMax, gCfKstarBinWidth);
    note += TString::Format("centrality groups: %s\n", centGroups);
    PdfHeader::MakePdfHeaderPage(pdfName, "checkHistCandidateFemto.C", inputs, note.Data(), true, 0);
  }

  TCanvas* c1 = new TCanvas("c1", "canvas", 1200, 800);
  TH1* h1 = 0;

  // ---- species purities (integrated and per centrality group), computed once
  for (Int_t s = 0; s < nS; ++s) {
    SpeciesInfo& S = sp[s];
    const CandidateFemtoSpecies& cs = cfg.species[s];
    S.grpMass.assign(nG, (TH1D*)0);
    S.grpFit.assign(nG, (TF1*)0);
    S.grpPurity.assign(nG, -1.);
    if (cs.purityByDefinition) {
      S.purity = 1.;
      S.puritySource = "purity 1 by definition (config)";
      for (Int_t g = 0; g < nG; ++g) S.grpPurity[g] = 1.;
      continue;
    }
    if (cs.purityModel != "invariantMassGausPol2" && cs.purityModel != "invariantMassStudentTPol2") {
      std::cerr << "[checkHistCandidateFemto] WARNING: purity model " << cs.purityModel << " of " << S.key
                << " is not implemented: purity unavailable" << std::endl;
      S.puritySource = "purity model " + TString(cs.purityModel.c_str()) + " not implemented";
      continue;
    }
    std::cout << "[checkHistCandidateFemto] " << S.key << ": purity fit mean +- " << S.localSigmaMultiplier
              << " sigma, pol" << S.bgOrder << " background" << std::endl;
    S.hPurityMass = puritySpectrum(fin, S);
    S.fit = FitMassSpectrum(S.hPurityMass, TString::Format("%s_mass_fit", S.key.Data()).Data(), S.localSigmaMultiplier,
                            S.bgOrder, cs.purityFitSeed.min, cs.purityFitSeed.max, cs.purityModel == "invariantMassStudentTPol2");
    S.purity = PurityFromFit(S.fit, cs.signal.min, cs.signal.max);
    // Groups: without a score cut hCandidateMassVsCent9_<s> (mass x cent9); with a score cut
    // hCandidateMassVsCent9VsScore_<s> (unavailable if missing).
    TH2* hVsCent = (TH2*)fin->Get(TString::Format("hCandidateMassVsCent9_%s", S.key.Data()));
    if (S.scoreCut) {
      hVsCent = 0;
      TH3* h3 = (TH3*)fin->Get(TString::Format("hCandidateMassVsCent9VsScore_%s", S.key.Data()));
      if (!h3)
        std::cerr << "[checkHistCandidateFemto] WARNING: " << S.key << ": hCandidateMassVsCent9VsScore_" << S.key
                  << " missing; centrality-group purities unavailable" << std::endl;
      for (Int_t g = 0; g < nG && h3; ++g) {
        S.grpMass[g] = projectMassCentScore(h3, groups[g], cs.scoreMin,
                                            TString::Format("%s_mass_%s", S.key.Data(), groupTag[g].Data()).Data());
        S.grpFit[g] = FitMassSpectrum(S.grpMass[g], TString::Format("%s_mass_fit_%s", S.key.Data(), groupTag[g].Data()).Data(),
                                      S.localSigmaMultiplier, S.bgOrder, cs.purityFitSeed.min, cs.purityFitSeed.max, cs.purityModel == "invariantMassStudentTPol2");
        S.grpPurity[g] = PurityFromFit(S.grpFit[g], cs.signal.min, cs.signal.max);
      }
    } else if (!hVsCent) {
      std::cerr << "[checkHistCandidateFemto] WARNING: " << S.key << ": no hCandidateMassVsCent9_" << S.key
                << "; centrality-group purities unavailable" << std::endl;
    }
    for (Int_t g = 0; g < nG && hVsCent; ++g) {
      S.grpMass[g] = projectCentGroupSumX(hVsCent, groups[g], TString::Format("%s_mass_%s", S.key.Data(), groupTag[g].Data()).Data());
      S.grpFit[g] = FitMassSpectrum(S.grpMass[g], TString::Format("%s_mass_fit_%s", S.key.Data(), groupTag[g].Data()).Data(),
                                    S.localSigmaMultiplier, S.bgOrder, cs.purityFitSeed.min, cs.purityFitSeed.max, cs.purityModel == "invariantMassStudentTPol2");
      S.grpPurity[g] = PurityFromFit(S.grpFit[g], cs.signal.min, cs.signal.max);
    }
  }
  for (Int_t s = 0; s < nS; ++s) {
    std::cout << "[PURITY] species (information) " << sp[s].key << " integrated P = " << sp[s].purity << " (" << sp[s].puritySource << ")";
    for (Int_t g = 0; g < nG; ++g) std::cout << " | " << groupTag[g] << " P = " << sp[s].grpPurity[g];
    std::cout << std::endl;
  }

  // ================= Species pages =================
  // Page A: mass spectra (purity spectrum with fit, usable candidates, KF finder output or m^2 vs p).
  newPage(c1, 3, nS);
  for (Int_t s = 0; s < nS; ++s) {
    SpeciesInfo& S = sp[s];
    const CandidateFemtoSpecies& cs = cfg.species[s];
    c1->cd(3 * s + 1);
    if (S.hPurityMass && !cs.purityByDefinition) {
      TH1* h = S.hPurityMass;
      h->SetMinimum(0.);
      h->SetTitle(TString::Format("%s, purity spectrum (%s);%s;Counts", S.label.Data(), S.puritySource.Data(), S.massTitle.Data()));
      h->Draw();
      drawFitComponents(S.fit);
      drawFitSigmaLines(S.fit, 3.0);
      drawSpeciesWindows(cs);
      drawMassLegend(h, S.fit);
      drawPadText(0.15, 0.82, TString::Format("P = %.3f", S.purity).Data(), 0.05);
    } else {
      drawPadText(0.1, 0.6, TString::Format("%s: %s", S.label.Data(), S.puritySource.Data()).Data(), 0.05);
    }
    c1->cd(3 * s + 2);
    if ((h1 = getOrWarn(fin, TString::Format("hCandidateMass_%s", S.key.Data())))) {
      h1->SetMinimum(0.);
      h1->SetTitle(TString::Format("%s, usable candidates (mass axis);%s;Counts", S.label.Data(), S.massTitle.Data()));
      h1->Draw();
      drawSpeciesWindows(cs);
    }
    c1->cd(3 * s + 3);
    if (S.rawHist.Length()) {
      if ((h1 = getOrWarn(fin, S.rawHist))) {
        h1->SetMinimum(0.);
        h1->SetTitle(TString::Format("%s, KF finder output before candidate cuts;%s;Counts", S.label.Data(), S.massTitle.Data()));
        h1->Draw();
        drawSpeciesWindows(cs);
      }
    } else {
      TH2* h2 = (TH2*)getOrWarn(fin, TString::Format("hCandidateMassVsP_%s", S.key.Data()));
      if (h2) {
        h2->SetStats(kFALSE);
        gPad->SetRightMargin(0.15);
        gPad->SetLogz();
        h2->SetTitle(TString::Format("%s candidates;|p| [GeV/c];%s", S.label.Data(), S.massTitle.Data()));
        h2->Draw("COLZ");
        drawSpeciesWindowsY(cs, h2->GetXaxis()->GetXmin(), h2->GetXaxis()->GetXmax());
      }
    }
  }
  c1->Print(pdfName);

  // Page B: topology QA of the KF species (selected candidates).
  {
    std::vector<Int_t> kfSpecies;
    for (Int_t s = 0; s < nS; ++s) if (!sp[s].topology.empty()) kfSpecies.push_back(s);
    if (!kfSpecies.empty()) {
      newPage(c1, 4, kfSpecies.size());
      for (size_t k = 0; k < kfSpecies.size(); ++k) {
        const SpeciesInfo& S = sp[kfSpecies[k]];
        for (Int_t i = 0; i < 4 && i < (Int_t)S.topology.size(); ++i) {
          c1->cd(4 * k + i + 1);
          const Bool_t logy = (i < 3);
          gPad->SetLeftMargin(0.15);
          gPad->SetRightMargin(0.05);
          gPad->SetLogy(logy);
          if ((h1 = getOrWarn(fin, S.topology[i]))) {
            h1->SetTitle(TString::Format("%s: %s", S.label.Data(), h1->GetTitle()));
            h1->GetXaxis()->SetNdivisions(505);
            if (logy) h1->SetMinimum(0.5);
            zoomToFilled(h1);
            h1->Draw();
          }
        }
      }
      c1->Print(pdfName);
    }
  }

  // Page 1: event QA and statistics.
  newPage(c1, 2, 1);
  c1->cd(1);
  if ((h1 = getOrWarn(fin, "hVz"))) {
    h1->SetTitle("V_{z}, all read events;V_{z} [cm];Events");
    h1->Draw();
  }
  c1->cd(2);
  {
    TPaveText* st = newNote(0.05, 0.05, 0.95, 0.95, 0.03);
    st->AddText("Statistics used");
    st->AddText(TString::Format("events read: %.0f", nEvRead));
    st->AddText(TString::Format("events accepted (hN): %.0f", nEvAcc));
    TH1* hOff = (TH1*)fin->Get("hKfOfficialEventSelection");
    const Bool_t official = hOff && hOff->GetBinContent(1) > 0;
    if (official) {
      st->SetTextSize(0.024);
      st->AddText("official event selection (rejected / kept):");
      for (Int_t b = 2; b <= hOff->GetNbinsX(); ++b)
        st->AddText(TString::Format("  %s: %.0f", hOff->GetXaxis()->GetBinLabel(b), hOff->GetBinContent(b)));
      st->AddText("event-cut detail (hKfEventSelection):");
    }
    if (hEvSel) {
      for (Int_t b = 2; b <= hEvSel->GetNbinsX(); ++b) {
        const char* lab = hEvSel->GetXaxis()->GetBinLabel(b);
        if (!lab || !lab[0]) continue;
        if (official && hEvSel->GetBinContent(b) == 0) continue;
        st->AddText(TString::Format("  %s: %.0f", lab, hEvSel->GetBinContent(b)));
      }
    }
    for (Int_t c = 0; c < nC; ++c) {
      const char* n = cfg.channels[c].name.c_str();
      TH1* hs = (TH1*)fin->Get(TString::Format("hKstarSE_%s", n));
      TH1* hm = (TH1*)fin->Get(TString::Format("hKstarME_%s", n));
      TH1* hr = (TH1*)fin->Get(TString::Format("hShareRejectSE_%s", n));
      st->AddText(TString::Format("%s sig #times sig: SE %.0f, ME %.0f; SE rejected (shared track): %.0f", n,
                                  hs ? hs->Integral(0, hs->GetNbinsX() + 1) : -1., hm ? hm->Integral(0, hm->GetNbinsX() + 1) : -1.,
                                  hr ? hr->GetEntries() : -1.));
    }
    st->AddText("(ME pairs: not checked, different events)");
    st->Draw();
  }
  c1->Print(pdfName);

  // ================= Channel pages =================
  TString channelPurityLines;
  for (Int_t c = 0; c < nC; ++c) {
    const CandidateFemtoChannel& ch = cfg.channels[c];
    const char* cn = ch.name.c_str();
    // Copies: the purities of this channel come from its own candidates (channelPurity).
    SpeciesInfo A = sp[ch.indexA];
    SpeciesInfo B = sp[ch.indexB];
    const CandidateFemtoSpecies& sa = cfg.species[ch.indexA];
    const CandidateFemtoSpecies& sb = cfg.species[ch.indexB];
    channelPurity(hSE[c], A, sa, sa, sb, 0, cn, groups, groupTag);
    channelPurity(hSE[c], B, sb, sa, sb, 1, cn, groups, groupTag);
    for (Int_t k = 0; k < 2; ++k) {
      const SpeciesInfo& S = k == 0 ? A : B;
      TString line = TString::Format("%s %s integrated %.6f (%s)", cn, S.key.Data(), S.purity, S.puritySource.Data());
      for (Int_t g = 0; g < nG; ++g) line += TString::Format("; %s %.6f", groupLabel[g].Data(), S.grpPurity[g]);
      std::cout << "[PURITY] channel " << line << std::endl;
      channelPurityLines += line + "\n";
    }
    const char* lA = A.label.Data();
    const char* lB = B.label.Data();
    const Double_t pA = A.purity, pB = B.purity;

    textPage(c1, pdfName, TString::Format("Channel %s: %s #times %s", cn, lA, lB).Data());

    // Event QA maps of this channel.
    {
      const TString names[8] = {"hEvtQA_MultVsNTof_all", TString::Format("hEvtQA_MultVsNTof_%s_withB", cn),
                                TString::Format("hEvtQA_NAVsNB_%s", cn), TString::Format("hEvtQA_MassAVsNB_%s", cn),
                                TString::Format("hEvtQA_MassBVsNA_%s", cn), TString::Format("hEvtQA_MassBVsMult_%s", cn),
                                TString::Format("hEvtQA_MassAVsMult_%s_withB", cn), TString::Format("hEvtQA_MassAVsVr_%s_withB", cn)};
      newPage(c1, 4, 2);
      for (Int_t i = 0; i < 8; ++i) {
        c1->cd(i + 1);
        TH2* h = (TH2*)getOrWarn(fin, names[i]);
        if (!h) continue;
        h->SetStats(kFALSE);
        gPad->SetRightMargin(0.15);
        gPad->SetLogz();
        h->SetTitle(TString::Format("%s (%.0f)", h->GetTitle(), h->GetEntries()));
        h->Draw("COLZ");
      }
      c1->Print(pdfName);
    }

    // Projections (integrated selection of the config).
    ClassProjection proj;
    projectClasses(hSE[c], hME[c], sel[c], sa, sb, TString::Format("%s_int", cn).Data(), proj);
    ChannelCFs cfs;
    buildChannelCFs(proj.se, proj.me, TString::Format("%s", cn).Data(), pA, pB, cfs);
    std::cout << "[CHANNEL] " << cn << " SE/ME pairs per class (A,B), k* < " << hSE[c]->GetAxis(0)->GetXmax() << ":";
    for (Int_t ca = 1; ca <= 3; ++ca)
      for (Int_t cb = 1; cb <= 3; ++cb)
        if (proj.se[ca][cb])
          std::cout << TString::Format(" [%s,%s] %.0f/%.0f", kClassTag[ca], kClassTag[cb],
                                       proj.se[ca][cb]->Integral(1, proj.se[ca][cb]->GetNbinsX()),
                                       proj.me[ca][cb]->Integral(1, proj.me[ca][cb]->GetNbinsX())).Data();
    std::cout << std::endl;
    for (Int_t ca = 1; ca <= 3; ++ca)
      for (Int_t cb = 1; cb <= 3; ++cb)
        if (cfs.cf[ca][cb]) printCF(cfs.cf[ca][cb], TString::Format("%s cf %s x %s", cn, kClassTag[ca], kClassTag[cb]).Data());
    printCF(cfs.final, TString::Format("%s cf_final", cn).Data());

    // Page 2: raw SE/ME k*, A signal x B {sig, L-SB, R-SB}.
    {
      newPage(c1, 3, 2);
      for (Int_t i = 0; i < 6; ++i) {
        const Int_t cb = i % 3 + 1;
        c1->cd(i + 1);
        TH1D* h = (i < 3) ? proj.se[1][cb] : proj.me[1][cb];
        if (!h) { drawPadText(0.15, 0.5, TString::Format("no %s class for %s", kClassLabel[cb], lB).Data()); continue; }
        h->SetTitle(TString::Format("%s: %s sig #times %s %s;k* [GeV/c];Pairs", i < 3 ? "SE" : "ME", lA, lB, kClassLabel[cb]));
        h->GetXaxis()->SetRangeUser(0., gDisplayMax);
        h->Draw();
      }
      c1->Print(pdfName);
    }

    // Page 3: CF, A sig x B {sig, L-SB, R-SB}, and signal vs B sideband average.
    {
      newPage(c1, 2, 2);
      for (Int_t cb = 1; cb <= 3; ++cb) {
        c1->cd(cb);
        TH1D* h = cloneOrNull(cfs.cf[1][cb], TString::Format("%s_p3_%d", cn, cb).Data());
        if (!h) { drawPadText(0.15, 0.5, TString::Format("no %s class for %s", kClassLabel[cb], lB).Data()); continue; }
        h->SetTitle(TString::Format("CF: %s sig #times %s %s;k* [GeV/c];C(k*)", lA, lB, kClassLabel[cb]));
        h->Draw("e");
        drawNormBand();
        h->Draw("e same");
      }
      c1->cd(4);
      if (cfs.cf[1][1]) {
        TH1D* s = cloneOrNull(cfs.cf[1][1], TString::Format("%s_sig_overlay", cn).Data());
        s->SetTitle(TString::Format("CF: %s signal vs %s sideband average;k* [GeV/c];C(k*)", lB, lB));
        s->SetLineColor(kBlue + 1);
        s->SetMarkerColor(kBlue + 1);
        s->Draw("e");
        drawNormBand();
        s->Draw("e same");
        TLegend* leg = new TLegend(0.45, 0.72, 0.89, 0.88);
        leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.04);
        leg->AddEntry(s, TString::Format("%s sig #times %s sig", lA, lB), "lep");
        if (cfs.bSbAvg) {
          TH1D* b = cloneOrNull(cfs.bSbAvg, TString::Format("%s_sbavg_overlay", cn).Data());
          b->SetLineColor(kRed + 1);
          b->SetMarkerColor(kRed + 1);
          b->Draw("e same");
          leg->AddEntry(b, TString::Format("%s sig #times %s SB (L+R avg)", lA, lB), "lep");
        }
        leg->Draw();
      }
      c1->Print(pdfName);
    }

    // Page 3b: individual purities.
    {
      c1->Clear();
      TPaveText* note = newNote(0.04, 0.04, 0.96, 0.96, 0.028);
      note->AddText("Purity = Gaussian / (Gaussian + background) integral of the local gaus+pol fit over the config signal window.");
      note->AddText("Fit mean/sigma and the n#sigma rows are diagnostic only: SE/ME pairing uses the config windows.");
      SpeciesInfo* list[2] = {&B, &A};
      const CandidateFemtoSpecies* lcs[2] = {&sb, &sa};
      const char* pn[2] = {"P_{2}", "P_{1}"};
      for (Int_t k = 0; k < 2; ++k) {
        note->AddText(" ");
        SpeciesInfo& S = *list[k];
        if (lcs[k]->purityByDefinition) {
          note->AddText(TString::Format("%s: %s = 1 by definition (species_%s_purityByDefinition)", S.label.Data(), pn[k], S.key.Data()));
        } else if (S.purity > 0.) {
          note->AddText(TString::Format("%s: %s = %.3f, window [%.4f, %.4f] (%s)", S.label.Data(), pn[k], S.purity,
                                        lcs[k]->signal.min, lcs[k]->signal.max, S.puritySource.Data()));
          if (S.fit) {
            note->AddText(TString::Format("    fit mean = %.5f, #sigma = %.5f", S.fit->GetParameter(1), TMath::Abs(S.fit->GetParameter(2))));
            addSigmaWindowLines(note, S.fit, S.hPurityMass, *lcs[k]);
          }
        } else {
          note->AddText(TString::Format("%s purity could not be computed (%s: missing, empty or degenerate fit)", S.label.Data(),
                                        S.puritySource.Data()));
        }
      }
      note->Draw();
      c1->Print(pdfName);
    }

    // Page 3c: 4-term purity-corrected CF.
    {
      newPage(c1, 2, 2);
      c1->cd(1);
      TPaveText* note = newNote(0.05, 0.1, 0.95, 0.9, 0.045);
      note->SetTextAlign(22);
      note->AddText(TString::Format("P_{1}(%s) = %.3f   P_{2}(%s) = %.3f   P_{1}P_{2} = %.3f", lA, pA, lB, pB, pA * pB));
      note->AddText("C_{corr} = [C_{sig} - (1-P_{1})P_{2}C_{SB1} - P_{1}(1-P_{2})C_{SB2}");
      note->AddText("                - (1-P_{1})(1-P_{2})C_{SB12}] / (P_{1}P_{2})");
      note->AddText(TString::Format("SB1 = %s SB #times %s sig,  SB2 = %s sig #times %s SB,  SB12 = both SB", lA, lB, lA, lB));
      note->AddText("Purity uncertainty is not propagated into C_{corr}.");
      note->Draw();
      c1->cd(2);
      if (cfs.aSbAvg) {
        TH1D* h = cloneOrNull(cfs.aSbAvg, TString::Format("%s_p3c_sb1", cn).Data());
        h->SetTitle(TString::Format("C_{SB1}: %s SB (L+R avg) #times %s sig;k* [GeV/c];C(k*)", lA, lB));
        h->Draw("e"); drawNormBand(); h->Draw("e same");
      } else {
        drawPadText(0.1, 0.5, TString::Format("C_{SB1}: no %s sideband (weight (1-P_{1})P_{2} = %.3f)", lA, (1. - pA) * pB).Data(), 0.05);
      }
      c1->cd(3);
      if (cfs.sb12) {
        TH1D* h = cloneOrNull(cfs.sb12, TString::Format("%s_p3c_sb12", cn).Data());
        h->SetTitle(TString::Format("C_{SB12}: %s SB #times %s SB (corner avg);k* [GeV/c];C(k*)", lA, lB));
        h->Draw("e"); drawNormBand(); h->Draw("e same");
      } else {
        drawPadText(0.1, 0.5, TString::Format("C_{SB12}: no %s/%s sideband pair", lA, lB).Data(), 0.05);
      }
      c1->cd(4);
      if (cfs.final) drawFinalCf(cloneOrNull(cfs.final, TString::Format("%s_p3c_final", cn).Data()), "cent-integrated", pA, pB, lA, lB);
      else drawPadText(0.1, 0.5, "purity-corrected CF unavailable (purity or background CF missing)", 0.045);
      c1->Print(pdfName);
    }

    // Page 3d: full mass-class matrix (rows = A class, columns = B class), existing classes only.
    {
      std::vector<Int_t> rows, cols;
      for (Int_t k = 1; k <= 3; ++k) { if (classExists(sa, k)) rows.push_back(k); if (classExists(sb, k)) cols.push_back(k); }
      newPage(c1, cols.size(), rows.size());
      for (size_t r = 0; r < rows.size(); ++r)
        for (size_t q = 0; q < cols.size(); ++q) {
          c1->cd(r * cols.size() + q + 1);
          gPad->SetTopMargin(0.12);
          gPad->SetBottomMargin(0.14);
          TH1D* h = cloneOrNull(cfs.cf[rows[r]][cols[q]], TString::Format("%s_grid_%d_%d", cn, rows[r], cols[q]).Data());
          if (!h) continue;
          h->SetTitle(TString::Format("%s %s #times %s %s;k* [GeV/c];C(k*)", lA, kClassLabel[rows[r]], lB, kClassLabel[cols[q]]));
          h->GetXaxis()->SetLabelSize(0.06);
          h->GetYaxis()->SetLabelSize(0.06);
          h->GetXaxis()->SetTitleSize(0.06);
          h->GetYaxis()->SetTitleSize(0.06);
          h->Draw("e");
          enlargeTitle(0.08);
          drawNormBand();
          h->Draw("e same");
        }
      c1->Print(pdfName);
    }

    // Page 3e: all CFs of the correction and the corrected CF, cent-integrated.
    {
      newPage(c1, 2, 1);
      c1->cd(1);
      drawConsolidatedCf(cloneOrNull(cfs.cf[1][1], TString::Format("%s_p3e_sig", cn).Data()),
                         cloneOrNull(cfs.aSbAvg, TString::Format("%s_p3e_sb1", cn).Data()),
                         cloneOrNull(cfs.bSbAvg, TString::Format("%s_p3e_sb2", cn).Data()),
                         cloneOrNull(cfs.sb12, TString::Format("%s_p3e_sb12", cn).Data()), "cent-integrated", lA, lB);
      c1->cd(2);
      drawFinalCf(cloneOrNull(cfs.final, TString::Format("%s_p3e_final", cn).Data()), "cent-integrated", pA, pB, lA, lB);
      c1->Print(pdfName);
    }

    // Keep the integrated CFs for the output file.
    {
      TH1D* keep[13] = {cfs.cf[1][1], cfs.cf[1][2], cfs.cf[1][3], cfs.cf[2][1], cfs.cf[2][2], cfs.cf[2][3],
                        cfs.cf[3][1], cfs.cf[3][2], cfs.cf[3][3], cfs.bSbAvg, cfs.aSbAvg, cfs.sb12, cfs.final};
      for (Int_t i = 0; i < 13; ++i) if (keep[i]) { outObjects.push_back(keep[i]); outObjectDir.push_back(cn); }
      for (Int_t ca = 1; ca <= 3; ++ca)
        for (Int_t cb = 1; cb <= 3; ++cb)
          if (proj.se[ca][cb]) {
            outObjects.push_back(proj.se[ca][cb]); outObjectDir.push_back(cn);
            outObjects.push_back(proj.me[ca][cb]); outObjectDir.push_back(cn);
          }
    }

    // Pages E/F: centrality groups (SE/ME and mass counts summed over the group's cent9 bins first).
    {
      std::vector<ChannelCFs> grp(nG);
      for (Int_t g = 0; g < nG; ++g) {
        TH1D* gse[4][4];
        TH1D* gme[4][4];
        for (Int_t ca = 0; ca < 4; ++ca)
          for (Int_t cb = 0; cb < 4; ++cb) {
            gse[ca][cb] = gme[ca][cb] = 0;
            if (!proj.seCent[ca][cb]) continue;
            gse[ca][cb] = projectCentGroupSumX(proj.seCent[ca][cb], groups[g],
                                               TString::Format("%s_%s_se_%d_%d", cn, groupTag[g].Data(), ca, cb).Data());
            gme[ca][cb] = projectCentGroupSumX(proj.meCent[ca][cb], groups[g],
                                               TString::Format("%s_%s_me_%d_%d", cn, groupTag[g].Data(), ca, cb).Data());
          }
        buildChannelCFs(gse, gme, TString::Format("%s_%s", cn, groupTag[g].Data()).Data(), A.grpPurity[g], B.grpPurity[g], grp[g]);
        std::cout << "[PURITY] " << cn << " " << groupTag[g] << ": P1(" << A.key << ") = " << A.grpPurity[g] << ", P2(" << B.key
                  << ") = " << B.grpPurity[g] << std::endl;
        printCF(grp[g].cf[1][1], TString::Format("%s %s cf sig x sig", cn, groupTag[g].Data()).Data());
        printCF(grp[g].final, TString::Format("%s %s cf_final", cn, groupTag[g].Data()).Data());
        TH1D* keep[6] = {grp[g].cf[1][1], grp[g].aSbAvg, grp[g].bSbAvg, grp[g].sb12, grp[g].final, grp[g].bOnlyCorrected};
        for (Int_t i = 0; i < 6; ++i) if (keep[i]) { outObjects.push_back(keep[i]); outObjectDir.push_back(cn); }
        for (Int_t ca = 0; ca < 4; ++ca)
          for (Int_t cb = 0; cb < 4; ++cb) { delete gse[ca][cb]; delete gme[ca][cb]; }
      }
      newPage(c1, nG, 2);
      for (Int_t g = 0; g < nG; ++g) {
        c1->cd(g + 1);
        drawConsolidatedCf(cloneOrNull(grp[g].cf[1][1], TString::Format("%s_E_sig_%d", cn, g).Data()),
                           cloneOrNull(grp[g].aSbAvg, TString::Format("%s_E_sb1_%d", cn, g).Data()),
                           cloneOrNull(grp[g].bSbAvg, TString::Format("%s_E_sb2_%d", cn, g).Data()),
                           cloneOrNull(grp[g].sb12, TString::Format("%s_E_sb12_%d", cn, g).Data()), groupLabel[g].Data(), lA, lB);
      }
      for (Int_t g = 0; g < nG; ++g) {
        c1->cd(nG + g + 1);
        if (grp[g].final) {
          drawFinalCf(cloneOrNull(grp[g].final, TString::Format("%s_E_final_%d", cn, g).Data()), groupLabel[g].Data(),
                      A.grpPurity[g], B.grpPurity[g], lA, lB);
        } else if (grp[g].bOnlyCorrected) {
          // A purity unavailable for this group: fall back to the B-background-only correction.
          drawFinalCf(cloneOrNull(grp[g].bOnlyCorrected, TString::Format("%s_E_bonly_%d", cn, g).Data()),
                      TString::Format("%s (%s-only correction)", groupLabel[g].Data(), lB).Data(), A.grpPurity[g], B.grpPurity[g], lA, lB);
        } else {
          drawPadText(0.1, 0.5, TString::Format("%s: purity-corrected CF unavailable", groupLabel[g].Data()).Data(), 0.05);
        }
      }
      c1->Print(pdfName);

      // Page F: mass spectra of the groups (top row = B, bottom row = A).
      newPage(c1, nG, 2);
      SpeciesInfo* rowS[2] = {&B, &A};
      const CandidateFemtoSpecies* rowC[2] = {&sb, &sa};
      const char* pn[2] = {"P_{2}", "P_{1}"};
      for (Int_t r = 0; r < 2; ++r)
        for (Int_t g = 0; g < nG; ++g) {
          c1->cd(r * nG + g + 1);
          SpeciesInfo& S = *rowS[r];
          if (rowC[r]->purityByDefinition) {
            drawPadText(0.1, 0.5, TString::Format("%s, %s: %s = 1 by definition", S.label.Data(), groupLabel[g].Data(), pn[r]).Data(), 0.05);
            continue;
          }
          TH1D* m = S.grpMass[g];
          if (!m) {
            drawPadText(0.1, 0.5, TString::Format("%s, %s: no per-cent9 mass spectrum", S.label.Data(), groupLabel[g].Data()).Data(), 0.05);
            continue;
          }
          m->SetMinimum(0.);
          m->SetTitle(TString::Format("%s mass, %s;%s;Counts", S.label.Data(), groupLabel[g].Data(), S.massTitle.Data()));
          m->Draw();
          drawFitComponents(S.grpFit[g]);
          drawFitSigmaLines(S.grpFit[g], 3.0);
          drawSpeciesWindows(*rowC[r]);
          drawPadText(0.15, 0.82, TString::Format("%s = %.3f", pn[r], S.grpPurity[g]).Data(), 0.05);
          if (g == 0 && r == 0) drawMassLegend(m, S.grpFit[g]);
        }
      c1->Print(pdfName);
    }

    // BDT score page (only if a BDT was used for A or B): score distributions, score vs mass, and the
    // raw signal CF for every score-axis bin edge of that species as threshold.
    if (A.hasScores || B.hasScores) {
      newPage(c1, 3, 2);
      SpeciesInfo* rowS[2] = {&A, &B};
      const CandidateFemtoSpecies* rowC[2] = {&sa, &sb};
      for (Int_t r = 0; r < 2; ++r) {
        SpeciesInfo& S = *rowS[r];
        c1->cd(3 * r + 1);
        gPad->SetLogy();
        TH1* hs = getOrWarn(fin, TString::Format("hCandidateScore_%s", S.key.Data()));
        if (hs) {
          hs->SetTitle(TString::Format("%s BDT score (%s);score;Candidates", S.label.Data(), S.hasScores ? "all usable" : "no BDT"));
          hs->SetMinimum(0.5);
          hs->Draw();
          if (rowC[r]->hasScoreMin) drawWindowLine(rowC[r]->scoreMin, kRed, 1);
        }
        c1->cd(3 * r + 2);
        TH2* hsm = (TH2*)getOrWarn(fin, TString::Format("hCandidateScoreVsMass_%s", S.key.Data()));
        if (hsm) {
          hsm->SetStats(kFALSE);
          gPad->SetRightMargin(0.15);
          hsm->SetTitle(TString::Format("%s;%s;BDT score", S.label.Data(), S.massTitle.Data()));
          hsm->Draw("COLZ");
          if (rowC[r]->hasScoreMin) {
            TLine* l = new TLine(hsm->GetXaxis()->GetXmin(), rowC[r]->scoreMin, hsm->GetXaxis()->GetXmax(), rowC[r]->scoreMin);
            l->SetLineColor(kRed);
            l->Draw();
          }
          drawSpeciesWindows(*rowC[r]);
        }
        c1->cd(3 * r + 3);
        const TAxis* ax = hSE[c]->GetAxis(6 + r);
        if (ax->GetNbins() <= 1) {
          drawPadText(0.1, 0.5, TString::Format("%s: one score bin, no threshold scan", S.label.Data()).Data(), 0.05);
          continue;
        }
        TLegend* leg = new TLegend(0.40, 0.60, 0.89, 0.88);
        leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.035);
        const Int_t colors[8] = {kBlack, kBlue + 1, kRed + 1, kGreen + 2, kMagenta + 1, kOrange + 7, kCyan + 2, kGray + 2};
        Bool_t first = kTRUE;
        for (Int_t b = 1; b <= ax->GetNbins(); ++b) {
          CandidatePairSelection s2 = sel[c];
          const Double_t thr = ax->GetBinLowEdge(b);
          if (r == 0) s2.scoreMinA = thr; else s2.scoreMinB = thr;
          const TString n = TString::Format("%s_score%d_%d", cn, r, b);
          TH1D* se = CandidateProjectKstar(hSE[c], s2, 1, 1, (n + "_se").Data());
          TH1D* me = CandidateProjectKstar(hME[c], s2, 1, 1, (n + "_me").Data());
          TH1D* cf = buildCF(se, me, (n + "_cf").Data());
          printCF(cf, TString::Format("%s sig x sig, score %s >= %g", cn, S.key.Data(), thr).Data());
          outObjects.push_back(cf); outObjectDir.push_back(cn);
          cf->SetLineColor(colors[(b - 1) % 8]);
          cf->SetMarkerColor(colors[(b - 1) % 8]);
          cf->SetMarkerStyle(20);
          cf->SetMarkerSize(0.6);
          if (first) {
            cf->SetTitle(TString::Format("CF %s sig #times %s sig vs %s score threshold;k* [GeV/c];C(k*)", lA, lB, S.label.Data()));
            cf->Draw("e");
            drawNormBand();
            cf->Draw("e same");
            first = kFALSE;
          } else {
            cf->Draw("e same");
          }
          leg->AddEntry(cf, TString::Format("score #geq %g (SE %.0f)", thr, se->Integral()), "lep");
          delete se;
          delete me;
        }
        leg->Draw();
      }
      c1->Print(pdfName);
    } else {
      std::cout << "[checkHistCandidateFemto] " << cn << ": no BDT scores in the input, score page skipped" << std::endl;
    }

    // Backup: CF (signal x signal) per cent9.
    textPage(c1, pdfName, TString::Format("Backup (%s): full 9-bin (cent9=0..8) breakdown", cn).Data());
    if (proj.seCent[1][1] && proj.meCent[1][1]) {
      newPage(c1, 3, 3);
      for (Int_t cent9 = 0; cent9 < 9; ++cent9) {
        c1->cd(cent9 + 1);
        const Int_t yBin = proj.seCent[1][1]->GetYaxis()->FindBin((Double_t)cent9);
        TH1D* se = proj.seCent[1][1]->ProjectionX(TString::Format("%s_c%d_se", cn, cent9), yBin, yBin);
        TH1D* me = proj.meCent[1][1]->ProjectionX(TString::Format("%s_c%d_me", cn, cent9), yBin, yBin);
        TH1D* cf = buildCF(se, me, TString::Format("%s_cf_sig_sig_cent9_%d", cn, cent9).Data());
        if (cf) {
          cf->SetTitle(TString::Format("cent9=%d (%s), SE=%.0f;k* [GeV/c];C(k*)", cent9, kCent9Pct[cent9], se->Integral()));
          cf->Draw("e");
          enlargeTitle(0.07);
          drawNormBand();
          cf->Draw("e same");
          outObjects.push_back(cf); outObjectDir.push_back(cn);
        }
        delete se;
        delete me;
      }
      c1->Print(pdfName);
    }

    // Lab acceptance: candidates and daughters (signal window), common axis range.
    {
      std::vector<TString> names;
      names.push_back(TString::Format("hCandidatePtVsYLab_%s_signal", A.key.Data()));
      for (size_t d = 0; d < A.daughters.size(); ++d)
        names.push_back(TString::Format("hDaughterPtVsYLab_%s_%s_signal", A.key.Data(), A.daughters[d].Data()));
      const Int_t nARow = names.size();
      names.push_back(TString::Format("hCandidatePtVsYLab_%s_signal", B.key.Data()));
      for (size_t d = 0; d < B.daughters.size(); ++d)
        names.push_back(TString::Format("hDaughterPtVsYLab_%s_%s_signal", B.key.Data(), B.daughters[d].Data()));
      const Int_t nBRow = names.size() - nARow;
      const Int_t nx = TMath::Max(nARow, nBRow) + 1;
      Double_t axMin, axMax, ayMin, ayMax;
      acceptanceUnionRange(fin, names, axMin, axMax, ayMin, ayMax);
      newPage(c1, nx, 2);
      for (Int_t i = 0; i < (Int_t)names.size(); ++i) {
        c1->cd(i < nARow ? i + 1 : nx + (i - nARow) + 1);
        drawAcceptanceMap(fin, names[i], axMin, axMax, ayMin, ayMax);
      }
      c1->cd(nx);
      drawPadText(0.05, 0.70, "p_{T} vs y_{lab}, signal-window", 0.07);
      drawPadText(0.05, 0.60, "candidates and their daughters", 0.07);
      drawPadText(0.05, 0.50, "(occupancy, not efficiency;", 0.07);
      drawPadText(0.05, 0.40, " common axis range)", 0.07);
      c1->Print(pdfName);
    }

    // Shared tracks: SE pairs that still share a track (both outside the signal windows) per daughter combination;
    // candidates removed by the shared-track rule vs kept, per species.
    {
      newPage(c1, 3, 1);
      c1->cd(1);
      TH1* hCombo = getOrWarn(fin, TString::Format("hShareRejectSE_%s", cn));
      if (hCombo) {
        const Int_t nA = A.daughters.size(), nB = B.daughters.size();
        if (hCombo->GetNbinsX() == nA * nB)
          for (Int_t i = 0; i < nA; ++i)
            for (Int_t j = 0; j < nB; ++j)
              hCombo->GetXaxis()->SetBinLabel(i * nB + j + 1, TString::Format("%s,%s", A.daughters[i].Data(), B.daughters[j].Data()));
        hCombo->SetTitle(TString::Format("SE pairs sharing a track, both outside the signal windows (%s dau, %s dau);;Pairs", lA, lB));
        hCombo->SetStats(kFALSE); hCombo->SetMinimum(0); hCombo->Draw("HIST TEXT0");
      }
      // Candidates removed by this channel's shared-track rule (StCandidateFemtoMaker): mass spectra of every usable
      // candidate of the species and of the candidates this channel removed.
      SpeciesInfo* both[2] = {&A, &B};
      for (Int_t i = 0; i < 2; ++i) {
        c1->cd(2 + i);
        TH1* kept = getOrWarn(fin, TString::Format("hCandidateMass_%s", both[i]->key.Data()));
        TH1* rem = getOrWarn(fin, TString::Format("hChannelSharedTrackRemovedMass_%s_%s", cn, both[i]->key.Data()));
        if (!kept || !rem) continue;
        TH1* k2 = (TH1*)kept->Clone(TString::Format("%s_qaKept%d", cn, i)); k2->SetDirectory(0);
        TH1* r2 = (TH1*)rem->Clone(TString::Format("%s_qaRemoved%d", cn, i)); r2->SetDirectory(0);
        const Int_t s1 = r2->FindBin(both[i]->cfg->signal.min + 1e-9), s2 = r2->FindBin(both[i]->cfg->signal.max - 1e-9);
        k2->SetTitle(TString::Format("%s: all usable (black), removed in this channel by the shared-track rule (red);%s;Counts",
                                     both[i]->label.Data(), both[i]->massTitle.Data()));
        k2->SetStats(kFALSE); k2->SetMinimum(0.5);
        gPad->SetLogy();
        k2->Draw("HIST"); r2->SetLineColor(kRed); r2->Draw("HIST SAME");
        drawSpeciesWindows(*both[i]->cfg);
        TLegend* lg = new TLegend(0.15, 0.75, 0.88, 0.88);
        lg->SetBorderSize(0); lg->SetFillStyle(0);
        lg->AddEntry(k2, TString::Format("usable %.0f (signal window %.0f)", k2->Integral(), k2->Integral(s1, s2)), "l");
        lg->AddEntry(r2, TString::Format("removed %.0f (signal window %.0f)", r2->Integral(), r2->Integral(s1, s2)), "l");
        lg->Draw();
      }
      c1->Print(pdfName);
    }

    // Close pairs: SE/ME in (dphi*, deta) per daughter combination, k* < 0.5 and k* < 1; and vs dphi*.
    {
      std::vector<TString> combo, comboLabel;
      for (size_t i = 0; i < A.daughters.size(); ++i)
        for (size_t j = 0; j < B.daughters.size(); ++j) {
          combo.push_back(TString::Format("%s_%s_%s", cn, A.daughters[i].Data(), B.daughters[j].Data()));
          comboLabel.push_back(TString::Format("%s %s - %s %s", lA, A.daughters[i].Data(), lB, B.daughters[j].Data()));
        }
      const Int_t nCombo = combo.size();
      const Int_t nx = TMath::Min(nCombo, 3), ny = (nCombo + nx - 1) / nx;
      for (Int_t pass = 0; pass < 2; ++pass) {
        const Int_t zHi = (pass == 0) ? 2 : 4;
        newPage(c1, nx, ny);
        for (Int_t i = 0; i < nCombo; ++i) {
          c1->cd(i + 1);
          TH3* se = (TH3*)getOrWarn(fin, TString::Format("hDEtaDPhiStarSE_%s", combo[i].Data()));
          TH3* me = (TH3*)getOrWarn(fin, TString::Format("hDEtaDPhiStarME_%s", combo[i].Data()));
          TH2D* r = closePairRatio2D(se, me, 1, zHi, TString::Format("%s_close2D_%d_%d", cn, pass, i));
          if (!r) { drawPadText(0.1, 0.5, TString::Format("%s: empty SE or ME", comboLabel[i].Data()).Data(), 0.05); continue; }
          r->SetTitle(TString::Format("SE/ME %s, k* < %.2f (SE %.0f);#Delta#phi*_{min};#Delta#eta", comboLabel[i].Data(),
                                      se->GetZaxis()->GetBinUpEdge(zHi), se->Integral(0, -1, 0, -1, 1, zHi)));
          r->SetStats(kFALSE); r->SetMinimum(0.); r->SetMaximum(2.);
          gPad->SetRightMargin(0.15); r->Draw("COLZ");
        }
        c1->Print(pdfName);
      }
      newPage(c1, nx, ny);
      for (Int_t i = 0; i < nCombo; ++i) {
        c1->cd(i + 1);
        TH3* se = (TH3*)getOrWarn(fin, TString::Format("hDEtaDPhiStarSE_%s", combo[i].Data()));
        TH3* me = (TH3*)getOrWarn(fin, TString::Format("hDEtaDPhiStarME_%s", combo[i].Data()));
        TH1D* lo = closePairRatioDPhi(se, me, gClosePairMaxDeltaEta, 1, 2, TString::Format("%s_closeLo_%d", cn, i));
        TH1D* hi = closePairRatioDPhi(se, me, gClosePairMaxDeltaEta, 3, 4, TString::Format("%s_closeHi_%d", cn, i));
        if (!lo && !hi) { drawPadText(0.1, 0.5, TString::Format("%s: empty SE or ME", comboLabel[i].Data()).Data(), 0.05); continue; }
        TH1D* f = lo ? lo : hi;
        f->SetTitle(TString::Format("SE/ME %s, |#Delta#eta|<%g;#Delta#phi*_{min};SE/ME (norm.)", comboLabel[i].Data(),
                                    gClosePairMaxDeltaEta));
        f->SetStats(kFALSE); f->SetMinimum(0.); f->SetMaximum(2.5);
        f->Draw("E");
        if (lo && hi) { hi->SetLineColor(kRed); hi->SetMarkerColor(kRed); hi->Draw("E SAME"); }
        TLegend* lg = new TLegend(0.55, 0.75, 0.88, 0.88);
        const TAxis* kz = se->GetZaxis();
        if (lo) lg->AddEntry(lo, TString::Format("k* < %.2f", kz->GetBinUpEdge(2)), "l");
        if (hi) lg->AddEntry(hi, TString::Format("%.2f < k* < %.2f", kz->GetBinLowEdge(3), kz->GetBinUpEdge(4)), "l");
        lg->Draw();
        TLine* one = new TLine(f->GetXaxis()->GetXmin(), 1., f->GetXaxis()->GetXmax(), 1.);
        one->SetLineStyle(2); one->Draw();
      }
      c1->Print(pdfName);
    }

    // Pair-level acceptance (sig x sig): pair y_cm and m_T vs k*, SE vs ME.
    {
      const char* nm[2] = {"YCm", "Mt"};
      const char* yl[2] = {"y_{cm}^{pair}", "m_{T} [GeV/c^{2}]"};
      newPage(c1, 3, 2);
      for (Int_t v = 0; v < 2; ++v) {
        TH2* se = (TH2*)getOrWarn(fin, TString::Format("hPair%sVsKstarSE_%s", nm[v], cn));
        TH2* me = (TH2*)getOrWarn(fin, TString::Format("hPair%sVsKstarME_%s", nm[v], cn));
        TH2* hh[2] = {se, me};
        for (Int_t k = 0; k < 2; ++k) {
          c1->cd(3 * v + k + 1);
          if (!hh[k]) continue;
          hh[k]->SetStats(kFALSE); hh[k]->GetXaxis()->SetRangeUser(0., gDisplayMax);
          hh[k]->SetTitle(TString::Format("%s (%.0f);k* [GeV/c];%s", k == 0 ? "SE" : "ME", hh[k]->GetEntries(), yl[v]));
          gPad->SetRightMargin(0.15); hh[k]->Draw("COLZ");
        }
        c1->cd(3 * v + 3);
        if (!se || !me) continue;
        TH2* se2 = (TH2*)se->Clone(TString::Format("%s_pairSe%d", cn, v)); se2->RebinX(10);
        TH2* me2 = (TH2*)me->Clone(TString::Format("%s_pairMe%d", cn, v)); me2->RebinX(10);
        TProfile* ps = se2->ProfileX(TString::Format("%s_pairSeProf%d", cn, v));
        TProfile* pm = me2->ProfileX(TString::Format("%s_pairMeProf%d", cn, v));
        ps->SetTitle(TString::Format("mean %s vs k*: SE (black) / ME (red);k* [GeV/c];<%s>", yl[v], yl[v]));
        ps->SetStats(kFALSE); ps->GetXaxis()->SetRangeUser(0., gDisplayMax);
        ps->Draw("E"); pm->SetLineColor(kRed); pm->SetMarkerColor(kRed); pm->Draw("E SAME");
      }
      c1->Print(pdfName);
    }

    // CM-rapidity acceptance of A and B and daughter eta_lab.
    {
      newPage(c1, 3, 1);
      SpeciesInfo* both[2] = {&A, &B};
      for (Int_t i = 0; i < 2; ++i) {
        c1->cd(i + 1);
        std::vector<TString> one(1, TString::Format("hCandidatePtVsYCm_%s_signal", both[i]->key.Data()));
        Double_t xMin, xMax, yMin, yMax;
        acceptanceUnionRange(fin, one, xMin, xMax, yMin, yMax);
        drawAcceptanceMap(fin, one[0], xMin, xMax, yMin, yMax);
        gPad->Update();
        TLine* l0 = new TLine(0., gPad->GetUymin(), 0., gPad->GetUymax()); l0->SetLineStyle(2); l0->Draw();
      }
      c1->cd(3);
      const Int_t col[8] = {kBlack, kGray + 1, kRed, kBlue, kGreen + 2, kMagenta + 1, kOrange + 7, kCyan + 2};
      TLegend* lg = new TLegend(0.15, 0.6, 0.45, 0.88);
      Bool_t first = kTRUE;
      Int_t k = 0;
      Double_t ymax = 0.;
      std::vector<TH1*> etas;
      for (Int_t i = 0; i < 2; ++i)
        for (size_t d = 0; d < both[i]->daughters.size(); ++d, ++k) {
          TH1* h = getOrWarn(fin, TString::Format("hDaughterEtaLab_%s_%s_signal", both[i]->key.Data(), both[i]->daughters[d].Data()));
          if (!h || h->Integral() <= 0) continue;
          TH1* e = (TH1*)h->Clone(TString::Format("%s_eta%d", cn, k));
          e->SetDirectory(0);
          e->Scale(1. / e->Integral()); e->SetLineColor(col[k % 8]); e->SetStats(kFALSE);
          e->SetTitle("daughter #eta_{lab} (global, unit area);#eta_{lab};a.u.");
          ymax = TMath::Max(ymax, e->GetMaximum());
          etas.push_back(e);
          lg->AddEntry(e, TString::Format("%s %s", both[i]->label.Data(), both[i]->daughters[d].Data()), "l");
        }
      for (size_t i = 0; i < etas.size(); ++i) {
        // y maximum: 0.05 or 1.1 x the largest bin, whichever is larger.
        etas[i]->SetMaximum(TMath::Max(0.05, 1.1 * ymax));
        etas[i]->Draw(first ? "HIST" : "HIST SAME");
        first = kFALSE;
      }
      if (!etas.empty()) lg->Draw();
      c1->Print(pdfName);
    }
  }  // channels

  // ================= Global backup / QA pages =================
  // Mass vs cent9 of every species (the Maker's hCandidateMassVsCent9_<s>), config windows overlaid.
  for (Int_t s = 0; s < nS; ++s) {
    TH2* hVsCent = (TH2*)getOrWarn(fin, TString::Format("hCandidateMassVsCent9_%s", sp[s].key.Data()));
    if (!hVsCent) continue;
    newPage(c1, 3, 3);
    for (Int_t cent9 = 0; cent9 < 9; ++cent9) {
      c1->cd(cent9 + 1);
      const Int_t yBin = hVsCent->GetYaxis()->FindBin((Double_t)cent9);
      TH1D* p = hVsCent->ProjectionX(TString::Format("%s_mass_cent%d", sp[s].key.Data(), cent9), yBin, yBin);
      p->SetTitle(TString::Format("%s mass, cent9=%d (%s);%s;Counts", sp[s].label.Data(), cent9, kCent9Pct[cent9], sp[s].massTitle.Data()));
      p->Draw();
      enlargeTitle(0.07);
      drawSpeciesWindows(cfg.species[s]);
    }
    c1->Print(pdfName);
  }

  // Track PID QA: all tracks (summed over cent9) and the track-source species m^2 vs p with the window.
  {
    newPage(c1, 3, 1);
    const char* qa[2] = {"hTrackQA_dEdxVsQP_cent", "hTrackQA_m2VsQP_cent"};
    for (Int_t k = 0; k < 2; ++k) {
      c1->cd(k + 1);
      TH2* sum = 0;
      for (Int_t cent9 = 0; cent9 < 9; ++cent9) {
        TH2* h = (TH2*)getOrWarn(fin, TString::Format("%s%d", qa[k], cent9));
        if (!h) continue;
        if (!sum) { sum = (TH2*)h->Clone(TString::Format("%ssum", qa[k])); sum->SetDirectory(0); }
        else sum->Add(h);
      }
      if (!sum) continue;
      sum->SetStats(kFALSE);
      sum->SetTitle(TString::Format("%s (all cent9, %.0f)", sum->GetTitle(), sum->GetEntries()));
      gPad->SetRightMargin(0.15);
      gPad->SetLogz();
      sum->Draw("COLZ");
    }
    c1->cd(3);
    Bool_t drawn = kFALSE;
    for (Int_t s = 0; s < nS && !drawn; ++s) {
      if (cfg.species[s].source != "track_pid") continue;
      TH2* h2 = (TH2*)getOrWarn(fin, TString::Format("hCandidateMassVsP_%s", sp[s].key.Data()));
      if (!h2) continue;
      h2->SetStats(kFALSE);
      gPad->SetRightMargin(0.15);
      gPad->SetLogz();
      h2->Draw("COLZ");
      drawSpeciesWindowsY(cfg.species[s], h2->GetXaxis()->GetXmin(), h2->GetXaxis()->GetXmax());
      drawn = kTRUE;
    }
    c1->Print(pdfName);
  }

  // ---- output ROOT file with the CFs (one directory per channel) and the purities.
  {
    TFile fout(rootOutName, "RECREATE");
    if (fout.IsZombie()) {
      std::cerr << "[checkHistCandidateFemto] ERROR: cannot write " << rootOutName << std::endl;
      return 1;
    } else {
      for (Int_t c = 0; c < nC; ++c) fout.mkdir(cfg.channels[c].name.c_str());
      for (size_t i = 0; i < outObjects.size(); ++i) {
        fout.cd(outObjectDir[i]);
        outObjects[i]->Write();
      }
      fout.cd();
      TString pur = "per channel (used for the purity-corrected CFs):\n" + channelPurityLines + "per species (information):\n";
      for (Int_t s = 0; s < nS; ++s) {
        pur += TString::Format("%s integrated %.6f (%s)", sp[s].key.Data(), sp[s].purity, sp[s].puritySource.Data());
        for (Int_t g = 0; g < nG; ++g) pur += TString::Format("; %s %.6f", groupLabel[g].Data(), sp[s].grpPurity[g]);
        pur += "\n";
      }
      TNamed("purities", pur.Data()).Write();
      TNamed("centGroups", centGroups).Write();
      TNamed("mainconf", mainconfPath.Data()).Write();
      fout.Close();
      std::cout << "Wrote " << rootOutName.Data() << std::endl;
    }
  }

  PdfHeader::ClosePdf(pdfName);
  fin->Close();
  std::cout << "Wrote " << pdfName.Data() << std::endl;
  if (gMissingHists > 0)
    std::cerr << "[checkHistCandidateFemto] WARNING: " << gMissingHists << " expected histogram(s) missing from ROOT" << std::endl;
  else
    std::cout << "[checkHistCandidateFemto] All drawn histograms present." << std::endl;
  return 0;
}
