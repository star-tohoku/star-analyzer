// checkHistAnaPhi.C - Draw histograms from run_anaPhi.C output and write PDF.
// Usage: root4star -l -b -q 'common/macro/checkHistAnaPhi.C("rootfile/auau19_anaPhi/xxx.root","auau19_anaPhi")'
// anaName is typically obtained by: script/analysis_info_helper.py --ana-name --mainconf config/mainconf/main_auau19_anaPhi.yaml
// PDF: share/figure/anaName/anaName_checkHistAnaPhi.pdf

#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TString.h>
#include <TStyle.h>
#include <iostream>
#include <vector>

#include "../../include/PdfIOMan.h"

void checkHistAnaPhi(const Char_t* inputRootFile,
                     const Char_t* anaName = "auau19_anaPhi")
{
  gROOT->SetBatch(kTRUE);

  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
  gStyle->SetTitleOffset(1.2, "x");
  gStyle->SetTitleOffset(1.4, "y");
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.12);
  gStyle->SetPadBottomMargin(0.12);

  TFile* fin = TFile::Open(inputRootFile);
  if (!fin || fin->IsZombie()) {
    std::cerr << "Error: Cannot open file " << inputRootFile << std::endl;
    return;
  }

  TString outDir = TString("share/figure/") + anaName + "/";
  if (gSystem->AccessPathName(outDir)) {
    gSystem->mkdir(outDir, kTRUE);
  }

  TString pdfName = TString(outDir) + anaName + "_checkHistAnaPhi.pdf";
  std::cout << "Output PDF: " << pdfName.Data() << std::endl;

  PdfHeader::OpenPdf(pdfName);

  std::vector<std::string> inputs;
  inputs.push_back((const char*)inputRootFile);

  TString note = "Check histograms from run_anaPhi.C (StPhiMaker output).\n";
  note += "Phi KK reconstruction, opening angle / pair rapidity cuts.\n";

  PdfHeader::MakePdfHeaderPage(pdfName, "checkHistAnaPhi.C", inputs, note.Data(), true, anaName);

  TCanvas* c1 = new TCanvas("c1", "canvas", 1200, 800);
  TH1* h1 = 0;
  TH2* h2 = 0;

  // Page 1: Event Level
  c1->Clear();
  c1->Divide(3, 2);
  c1->cd(1); h1 = (TH1*)fin->Get("hVz"); if (h1) h1->Draw();
  c1->cd(2); h1 = (TH1*)fin->Get("hVzDiff"); if (h1) h1->Draw();
  c1->cd(3); h2 = (TH2*)fin->Get("hVxVy"); if (h2) h2->Draw("colz");
  c1->cd(4); h1 = (TH1*)fin->Get("hRefMult"); if (h1) h1->Draw();
  c1->cd(5); h2 = (TH2*)fin->Get("hVzVsRun"); if (h2) h2->Draw("colz");
  c1->cd(6); h2 = (TH2*)fin->Get("hRefMultVsVz"); if (h2) h2->Draw("colz");
  c1->Print(pdfName);

  // Page 2: Track Kinematics & Quality
  c1->Clear();
  c1->Divide(4, 2);
  c1->cd(1); gPad->SetLogy(); h1 = (TH1*)fin->Get("hPt"); if (h1) h1->Draw();
  c1->cd(2); gPad->SetLogy(0); h1 = (TH1*)fin->Get("hEta"); if (h1) h1->Draw();
  c1->cd(3); h1 = (TH1*)fin->Get("hPhi"); if (h1) h1->Draw();
  c1->cd(4); h1 = (TH1*)fin->Get("hCharge"); if (h1) h1->Draw("hist");
  c1->cd(5); h1 = (TH1*)fin->Get("hNHitsFit"); if (h1) h1->Draw();
  c1->cd(6); h1 = (TH1*)fin->Get("hNHitsRatio"); if (h1) h1->Draw();
  c1->cd(7); gPad->SetLogy(); h1 = (TH1*)fin->Get("hDCA"); if (h1) h1->Draw();
  c1->cd(8); gPad->SetLogy(); h1 = (TH1*)fin->Get("hChi2"); if (h1) h1->Draw();
  c1->Print(pdfName);

  // Page 3: PID (TPC & TOF)
  c1->Clear();
  c1->Divide(3, 2);
  double pmax = 5.0;
  double dedxmax = 1e-6;
  c1->cd(1); gPad->SetLogz(); h2 = (TH2*)fin->Get("hDedxVsP"); if (h2) { h2->GetXaxis()->SetRangeUser(0, pmax); h2->GetYaxis()->SetRangeUser(0, dedxmax); h2->Draw("colz"); }
  c1->cd(2); gPad->SetLogz(); h2 = (TH2*)fin->Get("hBetaVsP"); if (h2) { h2->GetXaxis()->SetRangeUser(0, pmax); h2->Draw("colz"); }
  c1->cd(3); gPad->SetLogz(); h2 = (TH2*)fin->Get("hMass2VsP"); if (h2) { h2->GetXaxis()->SetRangeUser(0, pmax); h2->Draw("colz"); }
  c1->cd(4); gPad->SetLogz(0); h2 = (TH2*)fin->Get("hNSigmaPionVsP"); if (h2) { h2->GetXaxis()->SetRangeUser(0, pmax); h2->Draw("colz"); }
  c1->cd(5); h2 = (TH2*)fin->Get("hNSigmaKaonVsP"); if (h2) { h2->GetXaxis()->SetRangeUser(0, pmax); h2->Draw("colz"); }
  c1->cd(6); h2 = (TH2*)fin->Get("hNSigmaProtonVsP"); if (h2) { h2->GetXaxis()->SetRangeUser(0, pmax); h2->Draw("colz"); }
  c1->Print(pdfName);

  // Page 4: Event Plane & Misc
  c1->Clear();
  c1->Divide(2, 2);
  c1->cd(1); gPad->SetLogy(); h1 = (TH1*)fin->Get("hTriggerIds"); if (h1) { h1->SetFillColor(17); h1->Draw(); }
  c1->cd(2); gPad->SetLogy(0); h1 = (TH1*)fin->Get("hTofMatchMult"); if (h1) h1->Draw();
  c1->cd(3); h2 = (TH2*)fin->Get("hQxQy"); if (h2) h2->Draw("colz");
  c1->cd(4); h1 = (TH1*)fin->Get("hPsi2"); if (h1) { h1->SetMinimum(0); h1->Draw(); }
  c1->Print(pdfName);

  // Page 5: KK invariant mass (main phi analysis)
  c1->Clear();
  c1->Divide(2, 2);
  c1->cd(1); h1 = (TH1*)fin->Get("hMKK_SameEvent"); if (h1) h1->Draw();
  c1->cd(2); h1 = (TH1*)fin->Get("hMKK_AllCombinations"); if (h1) h1->Draw();
  c1->cd(3); h1 = (TH1*)fin->Get("hMKK_OpeningAngleCut"); if (h1) h1->Draw();
  c1->cd(4); h1 = (TH1*)fin->Get("hMKK_RapidityCut"); if (h1) h1->Draw();
  c1->Print(pdfName);

  // Page 6: KK mass (both cuts, mixed, background sub)
  c1->Clear();
  c1->Divide(2, 2);
  c1->cd(1); h1 = (TH1*)fin->Get("hMKK_BothCuts"); if (h1) h1->Draw();
  c1->cd(2); h1 = (TH1*)fin->Get("hMKK_MixedEvent"); if (h1) h1->Draw();
  c1->cd(3); h1 = (TH1*)fin->Get("hMKK_BackgroundSubtracted"); if (h1) h1->Draw();
  c1->cd(4); h2 = (TH2*)fin->Get("hMKK_vs_Pt"); if (h2) h2->Draw("colz");
  c1->Print(pdfName);

  // Page 7: Opening angle & pair rapidity (raw and 2D)
  c1->Clear();
  c1->Divide(3, 2);
  c1->cd(1); h1 = (TH1*)fin->Get("hOpeningAngle_Raw"); if (h1) h1->Draw();
  c1->cd(2); h1 = (TH1*)fin->Get("hPairRapidity_Raw"); if (h1) h1->Draw();
  c1->cd(3); h2 = (TH2*)fin->Get("hOpeningAngle_vs_MKK"); if (h2) h2->Draw("colz");
  c1->cd(4); h2 = (TH2*)fin->Get("hPairRapidity_vs_MKK"); if (h2) h2->Draw("colz");
  c1->cd(5); h2 = (TH2*)fin->Get("hOpeningAngle_vs_Pt"); if (h2) h2->Draw("colz");
  c1->cd(6); h2 = (TH2*)fin->Get("hOpeningAngle_vs_Rapidity"); if (h2) h2->Draw("colz");
  c1->Print(pdfName);

  // Page 8: KK pair pT (raw), rapidity vs pT, and after-cuts 1D
  c1->Clear();
  c1->Divide(3, 2);
  c1->cd(1); h1 = (TH1*)fin->Get("hPairPt_Raw"); if (h1) h1->Draw();
  c1->cd(2); h2 = (TH2*)fin->Get("hPairRapidity_vs_Pt"); if (h2) h2->Draw("colz");
  c1->cd(3); h1 = (TH1*)fin->Get("hOpeningAngle_AfterCuts"); if (h1) h1->Draw();
  c1->cd(4); h1 = (TH1*)fin->Get("hPairRapidity_AfterCuts"); if (h1) h1->Draw();
  c1->cd(5); h1 = (TH1*)fin->Get("hPairPt_AfterCuts"); if (h1) h1->Draw();
  c1->cd(6); h1 = (TH1*)fin->Get("hMKK_BothCuts"); if (h1) h1->Draw();
  c1->Print(pdfName);

  // Page 9: Kaon QA (K from phi candidates)
  c1->Clear();
  c1->Divide(3, 1);
  c1->cd(1); h1 = (TH1*)fin->Get("hK_Pt"); if (h1) h1->Draw();
  c1->cd(2); h1 = (TH1*)fin->Get("hK_Eta"); if (h1) h1->Draw();
  c1->cd(3); h1 = (TH1*)fin->Get("hK_NSigma"); if (h1) h1->Draw();
  c1->Print(pdfName);

  // Page 10: Rapidity window variations (|y| < 0.4, 0.3, 0.2, 0.1)
  c1->Clear();
  c1->Divide(2, 2);
  c1->cd(1); h1 = (TH1*)fin->Get("hMKK_RapidityCut_0p4"); if (h1) h1->Draw();
  c1->cd(2); h1 = (TH1*)fin->Get("hMKK_RapidityCut_0p3"); if (h1) h1->Draw();
  c1->cd(3); h1 = (TH1*)fin->Get("hMKK_RapidityCut_0p2"); if (h1) h1->Draw();
  c1->cd(4); h1 = (TH1*)fin->Get("hMKK_RapidityCut_0p1"); if (h1) h1->Draw();
  c1->Print(pdfName);

  // Page 11: M_KK in K+ multiplicity bins (angle+rapidity cut)
  c1->Clear();
  c1->Divide(2, 2);
  c1->cd(1); h1 = (TH1*)fin->Get("hMKK_Mult0to5_KPlus"); if (h1) h1->Draw();
  c1->cd(2); h1 = (TH1*)fin->Get("hMKK_Mult5to10_KPlus"); if (h1) h1->Draw();
  c1->cd(3); h1 = (TH1*)fin->Get("hMKK_Mult10to20_KPlus"); if (h1) h1->Draw();
  c1->cd(4); h1 = (TH1*)fin->Get("hMKK_Mult20up_KPlus"); if (h1) h1->Draw();
  c1->Print(pdfName);

  // Page 12: M_KK in K- multiplicity bins (angle+rapidity cut)
  c1->Clear();
  c1->Divide(2, 2);
  c1->cd(1); h1 = (TH1*)fin->Get("hMKK_Mult0to5_KMinus"); if (h1) h1->Draw();
  c1->cd(2); h1 = (TH1*)fin->Get("hMKK_Mult5to10_KMinus"); if (h1) h1->Draw();
  c1->cd(3); h1 = (TH1*)fin->Get("hMKK_Mult10to20_KMinus"); if (h1) h1->Draw();
  c1->cd(4); h1 = (TH1*)fin->Get("hMKK_Mult20up_KMinus"); if (h1) h1->Draw();
  c1->Print(pdfName);

  // Page 13: K+ and K- multiplicity (per event) and correlation
  c1->Clear();
  c1->Divide(2, 2);
  c1->cd(1); h2 = (TH2*)fin->Get("hNKaonPlusVsNKaonMinus"); if (h2) h2->Draw("colz");
  c1->cd(2); h1 = (TH1*)fin->Get("hNKaonPlus"); if (h1) h1->Draw();
  c1->cd(3); h1 = (TH1*)fin->Get("hNKaonMinus"); if (h1) h1->Draw();
  c1->cd(4); /* spare */;
  c1->Print(pdfName);

  PdfHeader::ClosePdf(pdfName);

  delete c1;
  fin->Close();

  std::cout << "Done. PDF saved to: " << pdfName.Data() << std::endl;
}
