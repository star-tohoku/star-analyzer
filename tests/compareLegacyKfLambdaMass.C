// Read-only final-selection comparison; deliberately no purity fit/subtraction.
#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TPad.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TString.h"
#include "TNamed.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace LegacyKfMassOverlay {
void Require(bool ok, const char* text) { if (!ok) throw std::runtime_error(text); }
double Window(TH1& h, double low, double high) {
  const int first = h.GetXaxis()->FindFixBin(low + 1e-12);
  const int last = h.GetXaxis()->FindFixBin(high - 1e-12);
  Require(first >= 1 && last <= h.GetNbinsX() && first <= last, "window outside histogram axis");
  Require(std::fabs(h.GetXaxis()->GetBinLowEdge(first) - low) < 1e-8 &&
          std::fabs(h.GetXaxis()->GetBinUpEdge(last) - high) < 1e-8,
          "display window must align with histogram bins; no partial-bin assumption is made");
  return h.Integral(first, last);
}
void Format(TH1& h, int color, int style, double low, double high) {
  h.SetDirectory(0); h.SetStats(kFALSE); h.SetLineColor(color);
  h.SetMarkerColor(color); h.SetLineStyle(style); h.SetLineWidth(2);
  h.GetXaxis()->SetRangeUser(low, high - 1e-12);
  h.GetXaxis()->SetTitle("M(p#pi^{-}) (GeV/c^{2})");
  h.GetXaxis()->SetTitleSize(0.048); h.GetXaxis()->SetLabelSize(0.042);
  h.GetYaxis()->SetTitleSize(0.047); h.GetYaxis()->SetLabelSize(0.042);
  h.GetYaxis()->SetTitleOffset(1.25);
  h.SetMinimum(0); h.SetMaximum();
}
}

void compareLegacyKfLambdaMass(const char* legacyFile, const char* kfFile,
    const char* outputPdf, const char* energyLabel, double massLow = 1.10, double massHigh = 1.14) {
  try {
    using namespace LegacyKfMassOverlay;
    Require(legacyFile && *legacyFile && kfFile && *kfFile && outputPdf && *outputPdf &&
            energyLabel && *energyLabel && massLow < massHigh, "explicit files/label and valid range required");
    Require(gSystem->AccessPathName(outputPdf), "refusing existing comparison PDF");
    TFile legacy(legacyFile, "READ"), kf(kfFile, "READ");
    Require(!legacy.IsZombie() && !kf.IsZombie(), "cannot read comparison ROOT files");
    TNamed* status = dynamic_cast<TNamed*>(kf.Get("FemtoLambdaRunStatus"));
    Require(status && TString(status->GetTitle()) == "completed", "KF input run is incomplete");
    TH1* original = dynamic_cast<TH1*>(legacy.Get("hLambda_InvMass"));
    TH1* selected = dynamic_cast<TH1*>(kf.Get("hLambda_InvMass"));
    Require(original && selected, "missing inclusive hLambda_InvMass");
    Require(original->GetDimension() == 1 && selected->GetDimension() == 1 &&
            original->GetNbinsX() == selected->GetNbinsX(), "incompatible mass axes");
    for (int b = 1; b <= original->GetNbinsX()+1; ++b)
      Require(std::fabs(original->GetXaxis()->GetBinLowEdge(b) - selected->GetXaxis()->GetBinLowEdge(b)) < 1e-10,
              "different mass bin edges");
    const double nLegacy = Window(*original, massLow, massHigh);
    const double nKf = Window(*selected, massLow, massHigh);
    Require(nLegacy > 0 && nKf > 0, "cannot normalize an empty display window");
    const TString directory = gSystem->DirName(outputPdf);
    if (gSystem->mkdir(directory, kTRUE) != 0)
      Require(!gSystem->AccessPathName(directory), "cannot create comparison output directory");

    TH1* rawLegacy = static_cast<TH1*>(original->Clone("overlayRawLegacy"));
    TH1* rawKf = static_cast<TH1*>(selected->Clone("overlayRawKf"));
    Format(*rawLegacy, kBlue+1, 2, massLow, massHigh);
    Format(*rawKf, kRed+1, 1, massLow, massHigh);
    TH1* normLegacy = static_cast<TH1*>(rawLegacy->Clone("overlayNormLegacy"));
    TH1* normKf = static_cast<TH1*>(rawKf->Clone("overlayNormKf"));
    normLegacy->SetDirectory(0); normKf->SetDirectory(0);
    normLegacy->Scale(1.0/nLegacy); normKf->Scale(1.0/nKf);

    gStyle->SetOptStat(0);
    gStyle->SetPaperSize(26.0, 15.6);
    TCanvas canvas("legacyKfLambdaMass", "Original Helix vs KF high-purity Lambda", 1400, 840);
    canvas.SetCanvasSize(1400, 840);
    TPad left("massRaw", "Raw", 0.01, 0.23, 0.50, 0.98);
    TPad right("massNorm", "Normalized", 0.50, 0.23, 0.99, 0.98);
    left.Draw(); right.Draw();
    left.cd(); left.SetLeftMargin(0.14); left.SetRightMargin(0.03); left.SetBottomMargin(0.13);
    rawLegacy->SetTitle(TString::Format("%s;M(p#pi^{-}) (GeV/c^{2});Raw candidates / bin", energyLabel));
    rawLegacy->SetMaximum(1.36*std::max(rawLegacy->GetMaximum(),rawKf->GetMaximum()));
    rawLegacy->Draw("HIST"); rawKf->Draw("HIST SAME");
    TLegend legend(0.16,0.73,0.97,0.91); legend.SetFillStyle(0); legend.SetBorderSize(0); legend.SetTextSize(0.031);
    legend.AddEntry(rawLegacy, "anaLambdaNuclearId (Helix Imp5)", "l");
    legend.AddEntry(rawKf, "StFemtoMaker + KF high purity", "l"); legend.Draw();
    right.cd(); right.SetLeftMargin(0.14); right.SetRightMargin(0.03); right.SetBottomMargin(0.13);
    normLegacy->SetTitle("Equal area in display window;M(p#pi^{-}) (GeV/c^{2});Fraction / bin");
    normLegacy->SetMaximum(1.36*std::max(normLegacy->GetMaximum(),normKf->GetMaximum()));
    normLegacy->Draw("HIST"); normKf->Draw("HIST SAME");
    TLegend normLegend(0.16,0.73,0.97,0.91); normLegend.SetFillStyle(0); normLegend.SetBorderSize(0); normLegend.SetTextSize(0.031);
    normLegend.AddEntry(normLegacy, "Helix: scaled by 1 / N_{window}", "l");
    normLegend.AddEntry(normKf, "KF: scaled by 1 / N_{window}", "l"); normLegend.Draw();
    canvas.cd();
    TLatex note; note.SetNDC(); note.SetTextSize(0.022);
    note.DrawLatex(0.045,0.185,TString::Format("All-mass raw entries (GetEntries, includes flow): Helix %.0f; KF %.0f.",original->GetEntries(),selected->GetEntries()));
    note.DrawLatex(0.045,0.142,TString::Format("Display-window bin sum [%.2f, %.2f) GeV/c^{2}: Helix %.0f; KF %.0f. Right: each window integral = 1.",massLow,massHigh,nLegacy,nKf));
    note.DrawLatex(0.045,0.099,"Same first 10,000 input events; different reconstruction AND final cuts. Not an isolated KF-algorithm comparison.");
    note.DrawLatex(0.045,0.056,"Candidate counts include background; their ratio is not efficiency. No background subtraction or purity fit.");
    canvas.Modified(); canvas.Update(); canvas.Print(outputPdf);
    FileStat_t info;
    Require(gSystem->GetPathInfo(outputPdf,info)==0 && info.fSize>0, "comparison PDF was not written");
    std::cout << "[mass-overlay] energy=" << energyLabel
              << " legacy_all_entries=" << original->GetEntries() << " kf_all_entries=" << selected->GetEntries()
              << " window=[" << massLow << "," << massHigh << ") legacy_window=" << nLegacy
              << " kf_window=" << nKf << " normalized_window_integral=1,1 output=" << outputPdf << std::endl;
    delete rawLegacy; delete rawKf; delete normLegacy; delete normKf;
  } catch (const std::exception& error) {
    std::cerr << "ERROR: " << error.what() << std::endl;
    gSystem->Exit(1);
  }
}
