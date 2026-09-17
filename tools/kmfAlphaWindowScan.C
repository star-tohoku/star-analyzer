// Step 10 T2: alpha as a function of where it is measured.
//
// S = F - alpha B takes its shape from the rotated template B and its normalisation from
// alpha = sum_win F / sum_win B, measured in a mass window that is supposed to hold no signal.
// Two things can bias it. A window too close to the phi keeps part of the peak, which inflates
// alpha and over-subtracts. A window far away is signal-free but only useful if the template has
// the right shape there -- and if it does, alpha does not depend on the window at all.
//
// So the measurement is alpha versus window position. A plateau is a template that works; a drift
// is a template-shape error, and no choice of window fixes it.
//
// The leak column is the phi yield predicted inside the window by a Voigt profile of the fitted
// mean and resolution with Gamma fixed to the PDG width, normalised so that its integral over the
// signal window equals the counted yield there. It says how much of any rise at small offsets is
// the peak's own tail.
//
// Usage: kmfAlphaWindowScan.C("<downstream merged>.root", "phi_deuteron", "rot", cent9Min, cent9Max)
#include "TFile.h"
#include "TH3.h"
#include "TH1.h"
#include "TMath.h"
#include <cstdio>
#include <string>

namespace {
Double_t integralIn(TH1* h, Double_t lo, Double_t hi, Double_t& err) {
  err = 0.0;
  if (!h) return 0.0;
  const Int_t b0 = h->GetXaxis()->FindBin(lo + 1e-9);
  const Int_t b1 = h->GetXaxis()->FindBin(hi - 1e-9);
  if (b1 < b0) return 0.0;
  return h->IntegralAndError(b0, b1, err);
}

Double_t voigtIntegral(Double_t lo, Double_t hi, Double_t mean, Double_t sigma, Double_t gamma) {
  const Int_t n = 2000;
  const Double_t dx = (hi - lo) / n;
  Double_t s = 0.0;
  for (Int_t i = 0; i < n; ++i) {
    const Double_t x = lo + (i + 0.5) * dx;
    s += TMath::Voigt(x - mean, sigma, gamma);
  }
  return s * dx;
}
}  // namespace

void kmfAlphaWindowScan(const char* mergedFile, const char* base = "phi_deuteron",
                        const char* tmpl = "rot", Int_t cent9Min = 2, Int_t cent9Max = 8,
                        Double_t kLo = 0.0, Double_t kHi = 0.2, Double_t mean = 1.0191,
                        Double_t sigma = 0.0034, Double_t gamma = 0.004249,
                        Double_t sigWinLo = 1.01, Double_t sigWinHi = 1.03) {
  TFile* f = TFile::Open(mergedFile);
  if (!f || f->IsZombie()) { printf("cannot open %s\n", mergedFile); return; }
  const std::string bach = (std::string(base).find("deuteron") != std::string::npos) ? "deuteron" : "proton";
  TH3* h3F = (TH3*)f->Get(Form("hPhiMKK_vs_KstarSE_%s_wide", base));
  TH3* h3B = (TH3*)f->Get(Form("hPhiMKK_vs_KstarSE_phi_%s_%s_wide", tmpl, bach.c_str()));
  if (!h3F || !h3B) {
    printf("missing TH3: F=%p B=%p (looked for hPhiMKK_vs_KstarSE_%s_wide and ..._phi_%s_%s_wide)\n",
           (void*)h3F, (void*)h3B, base, tmpl, bach.c_str());
    return;
  }
  const Int_t zLo = h3F->GetZaxis()->FindBin(cent9Min - 0.01);
  const Int_t zHi = h3F->GetZaxis()->FindBin(cent9Max + 0.01);
  const Int_t yLo = h3F->GetYaxis()->FindBin(kLo + 1e-9);
  const Int_t yHi = h3F->GetYaxis()->FindBin(kHi - 1e-9);
  h3F->GetZaxis()->SetRange(zLo, zHi);
  h3B->GetZaxis()->SetRange(zLo, zHi);
  h3F->GetYaxis()->SetRange(yLo, yHi);
  h3B->GetYaxis()->SetRange(yLo, yHi);
  TH1* hF = (TH1*)h3F->Project3D("x");
  TH1* hB = (TH1*)h3B->Project3D("x");
  if (!hF || !hB) { printf("projection failed\n"); return; }

  Double_t eSig = 0.0;
  const Double_t fSig = integralIn(hF, sigWinLo, sigWinHi, eSig);
  const Double_t vSig = voigtIntegral(sigWinLo, sigWinHi, mean, sigma, gamma);

  printf("\n%s  %s template  cent9 %d-%d  k* %.2f-%.2f\n", base, tmpl, cent9Min, cent9Max, kLo, kHi);
  printf("F in signal window %.1f\n", fSig);
  printf("%-16s %12s %12s %10s %10s %12s\n", "window", "sumF", "sumB", "alpha", "dalpha", "leak/F [%]");
  const Double_t wins[][2] = {{1.030, 1.050}, {1.035, 1.055}, {1.040, 1.060}, {1.045, 1.065},
                              {1.050, 1.070}, {1.060, 1.080}, {1.070, 1.090}, {1.080, 1.100},
                              {1.100, 1.120}, {1.120, 1.160}, {0.990, 1.005}, {0.988, 1.000}};
  const Int_t nW = sizeof(wins) / sizeof(wins[0]);
  // Scale factor that turns the Voigt into counts: counted phi yield in the signal window is
  // approximated by F_sig minus a flat-ish background, so use the fitted-yield-free estimate
  // S_sig = F_sig - alpha_far * B_sig with alpha from the farthest window.
  Double_t eFfar = 0, eBfar = 0;
  const Double_t fFar = integralIn(hF, 1.120, 1.160, eFfar);
  const Double_t bFar = integralIn(hB, 1.120, 1.160, eBfar);
  const Double_t aFar = (bFar > 0.0) ? fFar / bFar : 0.0;
  Double_t eBsig = 0.0;
  const Double_t bSig = integralIn(hB, sigWinLo, sigWinHi, eBsig);
  const Double_t sSig = fSig - aFar * bSig;
  const Double_t scale = (vSig > 0.0) ? sSig / vSig : 0.0;
  printf("(reference: alpha from 1.120-1.160 = %.5f, S in signal window = %.1f)\n", aFar, sSig);

  for (Int_t i = 0; i < nW; ++i) {
    Double_t eF = 0, eB = 0;
    const Double_t sF = integralIn(hF, wins[i][0], wins[i][1], eF);
    const Double_t sB = integralIn(hB, wins[i][0], wins[i][1], eB);
    if (sB <= 0.0) { printf("%-16s  (empty template)\n", Form("%.3f-%.3f", wins[i][0], wins[i][1])); continue; }
    const Double_t a = sF / sB;
    const Double_t da = a * TMath::Sqrt((eF / (sF + 1e-12)) * (eF / (sF + 1e-12)) +
                                        (eB / (sB + 1e-12)) * (eB / (sB + 1e-12)));
    const Double_t leak = scale * voigtIntegral(wins[i][0], wins[i][1], mean, sigma, gamma);
    printf("%-16s %12.1f %12.1f %10.5f %10.5f %12.3f\n", Form("%.3f-%.3f", wins[i][0], wins[i][1]),
           sF, sB, a, da, (sF > 0.0) ? 100.0 * leak / sF : 0.0);
  }
}
