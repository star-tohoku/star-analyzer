// One figure per bachelor species: m2 against rigidity for the TPC-only selection, with the
// configured m2 window drawn on top, plus a summary page of purity against momentum.
//
// The m2 axis is m2/Z^2 for 3He and 4He, which is the convention the YAML windows use and what
// TrackRowV2::Mass2() returns; the momentum axis is rigidity p/Z for the same reason.
//
// Usage: PurityPlots.C("<purity dir>", "<out.pdf>")
#include "TFile.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TLatex.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TMultiGraph.h"
#include "TString.h"
#include "TStyle.h"
#include <cstdio>

void PurityPlots(const char* dir, const char* outPdf) {
  const Int_t kN = 5;
  const char* sp[kN] = {"proton", "deuteron", "triton", "he3", "he4"};
  const char* file[kN] = {"proton_purity_m2vsp.root", "purity_deuteron.root", "purity_triton.root",
                          "purity_he3.root", "purity_he4.root"};
  const char* hist[kN] = {"hM2vsP", "hM2vsP_deuteron", "hM2vsP_triton", "hM2vsP_he3",
                          "hM2vsP_he4"};
  const Double_t lo[kN] = {0.6, 2.5, 5.8, 0.8, 2.5};
  const Double_t hi[kN] = {1.2, 4.5, 9.8, 2.8, 4.5};
  // Overall purity from NuclearPurity.C / ProtonPurity.C, for the label.
  const Double_t pur[kN] = {0.9873, 0.9865, 0.9778, 0.9869, 0.9275};

  gStyle->SetOptStat(0);
  TCanvas* c = new TCanvas("cPurAll", "", 1200, 800);
  c->Print(TString(outPdf) + "[");

  for (Int_t i = 0; i < kN; ++i) {
    TFile* f = TFile::Open(TString(dir) + "/" + file[i]);
    if (!f || f->IsZombie()) { printf("cannot open %s\n", file[i]); continue; }
    TH2D* h = (TH2D*)f->Get(hist[i]);
    if (!h) { printf("no %s in %s\n", hist[i], file[i]); f->Close(); continue; }
    c->Clear();
    c->SetLogz();
    c->SetRightMargin(0.13);
    h->SetTitle(Form("%s: TPC-only selection (|n#sigma| < 2, no TOF rule);p/Z (GeV/c);m^{2}/Z^{2} (GeV^{2})",
                     sp[i]));
    h->Draw("colz");
    TLine* l1 = new TLine(h->GetXaxis()->GetXmin(), lo[i], h->GetXaxis()->GetXmax(), lo[i]);
    TLine* l2 = new TLine(h->GetXaxis()->GetXmin(), hi[i], h->GetXaxis()->GetXmax(), hi[i]);
    l1->SetLineColor(kRed); l1->SetLineWidth(2); l1->Draw();
    l2->SetLineColor(kRed); l2->SetLineWidth(2); l2->Draw();
    TLatex* t = new TLatex();
    t->SetNDC(); t->SetTextSize(0.030);
    t->DrawLatex(0.13, 0.94, Form("red: configured m^{2} window [%.1f, %.1f]", lo[i], hi[i]));
    t->DrawLatex(0.13, 0.90, Form("overall purity (sideband) = %.4f", pur[i]));
    c->Print(outPdf);
    f->Close();
  }

  c->Print(TString(outPdf) + "]");
  printf("wrote %s\n", outPdf);
}
