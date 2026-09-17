// Proton purity from a fit to m2, rather than by counting inside a window.
//
// The window count in ProtonPurity.C answers "what fraction of TOF-confirmed tracks sits inside
// the proton window". It cannot see contamination that sits UNDER the proton peak: a kaon whose
// m2 fluctuates into the window, or the flat background left by mismatched TOF hits. Those are
// counted as protons. This decomposes the distribution instead, so the proton yield inside the
// window is separated from whatever else is there.
//
// Model: three Gaussians (proton, kaon, pion, at m2 = 0.880, 0.244, 0.019) plus a linear term for
// the mismatch background. The kaon and pion means are held at their physical values because the
// sample is already |nSigma_p| < 2 and their peaks are small -- letting them float lets them wander
// into the proton peak and absorb signal. Their widths float, since the m2 resolution grows with
// momentum.
//
// The quantity wanted is the purity of a TPC-ONLY selection, which applies no m2 window at all:
//
//   purity = (proton yield over the whole m2 range) / (all tracks)
//
// A first version of this file divided the proton yield inside the window by the model inside the
// window, which answers a different question and came out at 0.999 everywhere -- it was measuring
// how proton-dominated the window is, not how pure the selection is.
//
// Two estimates are reported because the fit alone cannot be trusted here. At a million entries
// per slice a single Gaussian does not describe a real m2 peak (chi2/ndf runs to 100), and the
// fitted normalisation misses the data by a couple of percent. So the fit is shown next to a
// sideband estimate that does not depend on the peak shape at all: everything below m2 = 0.45 is
// unambiguously kaon or pion, and the mismatch background is taken as flat at the level measured
// above the proton peak, where nothing physical sits (the deuteron is at m2 = 3.5).
//
// Usage: ProtonPurityFit.C("<proton_purity_m2vsp.root>", "<out.pdf>")
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TMath.h"
#include "TString.h"
#include <cstdio>

namespace {
const Double_t kM2Proton = 0.8803;   // (0.938272)^2
const Double_t kM2Kaon = 0.2437;     // (0.493677)^2
const Double_t kM2Pion = 0.0195;     // (0.139570)^2
const Double_t kWinLo = 0.6, kWinHi = 1.2;
}  // namespace

void ProtonPurityFit(const char* inFile, const char* outPdf = "") {
  TFile* f = TFile::Open(inFile);
  if (!f || f->IsZombie()) { printf("cannot open %s\n", inFile); return; }
  TH2D* h2 = (TH2D*)f->Get("hM2vsP");
  if (!h2) { printf("no hM2vsP in %s\n", inFile); return; }

  const Int_t kNS = 5;
  const Double_t edge[kNS + 1] = {0.4, 0.8, 1.2, 1.6, 2.0, 2.4};

  TCanvas* c = 0;
  if (TString(outPdf).Length()) {
    c = new TCanvas("cPur", "", 1100, 750);
    c->Divide(3, 2);
  }

  printf("\n  %-12s %10s %10s %10s %10s %10s %8s\n", "p (GeV/c)", "tracks", "K/pi<0.45",
         "flat bkg", "purity(sb)", "purity(fit)", "chi2/ndf");
  Double_t sumAll = 0, sumP = 0;
  for (Int_t s = 0; s < kNS; ++s) {
    const Int_t b1 = h2->GetXaxis()->FindBin(edge[s] + 1e-6);
    const Int_t b2 = h2->GetXaxis()->FindBin(edge[s + 1] - 1e-6);
    TH1D* h = h2->ProjectionY(Form("m2_%d", s), b1, b2);
    if (h->Integral() < 100) continue;
    const Double_t bw = h->GetBinWidth(1);

    TF1* fit = new TF1(Form("fit%d", s), "gaus(0)+gaus(3)+gaus(6)+pol1(9)", -0.2, 1.8);
    const Double_t peak = h->GetBinContent(h->GetMaximumBin());
    fit->SetParameters(peak, kM2Proton, 0.04, peak * 0.01, kM2Kaon, 0.03, peak * 0.01, kM2Pion,
                       0.02, peak * 1e-4, 0.0);
    fit->SetParLimits(0, 0, peak * 10);
    fit->SetParLimits(1, kM2Proton - 0.06, kM2Proton + 0.06);
    fit->SetParLimits(2, 0.008, 0.30);
    fit->FixParameter(4, kM2Kaon);
    fit->SetParLimits(3, 0, peak);
    fit->SetParLimits(5, 0.008, 0.25);
    fit->FixParameter(7, kM2Pion);
    fit->SetParLimits(6, 0, peak);
    fit->SetParLimits(8, 0.008, 0.25);
    h->Fit(fit, "QNR");

    // Fit estimate: the proton Gaussian over the whole range, against every track in the slice.
    TF1* pOnly = new TF1(Form("p%d", s), "gaus", -0.5, 2.5);
    pOnly->SetParameters(fit->GetParameter(0), fit->GetParameter(1), fit->GetParameter(2));
    const Double_t nTotal = h->Integral();
    const Double_t nProtonFit = pOnly->Integral(-0.5, 2.5) / bw;
    const Double_t purFit = nTotal > 0 ? nProtonFit / nTotal : 0.0;

    // Sideband estimate, independent of the peak shape.
    const Double_t nKpi = h->Integral(1, h->FindBin(0.45));
    const Int_t sb1 = h->FindBin(1.40), sb2 = h->FindBin(2.20);
    const Double_t perBin = (sb2 > sb1) ? h->Integral(sb1, sb2) / (Double_t)(sb2 - sb1 + 1) : 0.0;
    // The flat term under the proton region only; the part below 0.45 is already inside nKpi.
    const Double_t nFlat = perBin * (Double_t)(h->GetNbinsX() - h->FindBin(0.45));
    const Double_t nProtonSb = nTotal - nKpi - nFlat;
    const Double_t purSb = nTotal > 0 ? nProtonSb / nTotal : 0.0;

    const Double_t ndf = fit->GetNDF();
    printf("  %4.1f - %4.1f  %10.0f %10.0f %10.0f %10.4f %10.4f %8.1f\n", edge[s], edge[s + 1],
           nTotal, nKpi, nFlat, purSb, purFit, ndf > 0 ? fit->GetChisquare() / ndf : 0.0);
    sumAll += nTotal;
    sumP += nProtonSb;

    if (c) {
      c->cd(s + 1);
      gPad->SetLogy();
      h->SetTitle(Form("%.1f < p < %.1f GeV/c;m^{2} (GeV^{2});tracks", edge[s], edge[s + 1]));
      h->GetXaxis()->SetRangeUser(-0.2, 1.8);
      h->SetLineColor(kBlack);
      h->Draw("hist");
      fit->SetLineColor(kRed);
      fit->Draw("same");
      pOnly->SetLineColor(kBlue);
      pOnly->SetLineStyle(2);
      pOnly->Draw("same");
      TLatex* t = new TLatex();
      t->SetNDC();
      t->SetTextSize(0.055);
      t->DrawLatex(0.15, 0.85, Form("purity(sb) = %.4f", purSb));
      t->DrawLatex(0.15, 0.79, Form("#sigma_{p} = %.3f", fit->GetParameter(2)));
    }
  }
  printf("  %-12s %10.0f %31.4f\n", "all", sumAll, sumAll > 0 ? sumP / sumAll : 0.0);
  printf("\npurity(sb): every track, minus the kaons and pions below m2 = 0.45, minus a flat\n"
         "mismatch background at the level measured in m2 = 1.4-2.2. No peak shape assumed.\n"
         "purity(fit): the proton Gaussian over the whole range, against every track. Read it only\n"
         "as a cross-check -- chi2/ndf shows the single Gaussian does not describe a peak this\n"
         "well populated.\n");

  if (c) {
    c->SaveAs(outPdf);
    printf("figure: %s\n", outPdf);
  }
  f->Close();
}
