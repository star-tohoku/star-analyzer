// Read-only QA for merged/new Lambda-nucleus ROOT output. No Phi mass model is reused.
#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TNamed.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TString.h"
#include "TSystem.h"
#include "ConfigManager.h"
#include "cuts/FemtoConfig.h"
#include <iostream>
#include <string>

namespace FemtoLambdaQa {
TH2* SumMass(TFile& file, const TString& pattern, const char* name, Int_t first, Int_t last) {
  TH2* sum = 0;
  for (Int_t cent = first; cent <= last; ++cent) {
    TH2* histogram = dynamic_cast<TH2*>(file.Get(TString::Format(pattern, cent)));
    if (!histogram) { delete sum; return 0; }
    if (!sum) { sum = static_cast<TH2*>(histogram->Clone(name)); sum->SetDirectory(0); }
    else sum->Add(histogram);
  }
  return sum;
}
void Message(const char* text) {
  TLatex label; label.SetNDC(); label.SetTextSize(0.035); label.DrawLatex(0.12, 0.5, text);
}

Int_t CheckAcceptanceAndSplit(TFile& file, const TString& channel,
                            const TString& species, TH2* acceptance[4]) {
  // Old ROOT output has none of the additions. A partial contract must not look successful.
  const TString names[] = {
    "hLambda_PtVsYLab_signal", "hLambdaProton_PtVsYLab_signal",
    "hLambdaPion_PtVsYLab_signal", "hNucleus_PtVsYLab_" + species,
    "hKstarSE_" + channel + "_signalLow", "hKstarSEVsCent_" + channel + "_signalLow",
    "hKstarME_" + channel + "_signalLow", "hKstarMEVsCent_" + channel + "_signalLow",
    "hKstarSE_" + channel + "_signalHigh", "hKstarSEVsCent_" + channel + "_signalHigh",
    "hKstarME_" + channel + "_signalHigh", "hKstarMEVsCent_" + channel + "_signalHigh"
  };
  const Int_t dimensions[] = {2, 2, 2, 2, 1, 2, 1, 2, 1, 2, 1, 2};
  Int_t present = 0;
  for (Int_t i = 0; i < 12; ++i) {
    TObject* object = file.Get(names[i]);
    if (!object) continue;
    ++present;
    TH1* histogram = dynamic_cast<TH1*>(object);
    if (!histogram || histogram->GetDimension() != dimensions[i]) {
      std::cerr << "ERROR: wrong histogram type/dimension for " << names[i] << std::endl;
      return -1;
    }
    if (i < 4) acceptance[i] = dynamic_cast<TH2*>(histogram);
  }
  if (present == 0) {
    std::cerr << "WARNING: new acceptance/split-signal histograms absent; rerun analysis "
              << "to create them. Drawing the original four QA pages only." << std::endl;
    return 0;
  }
  if (present != 12) {
    std::cerr << "ERROR: incomplete acceptance/split-signal histogram contract ("
              << present << "/12); rerun analysis with matching Maker and YAML."
              << std::endl;
    return -1;
  }
  return 1;
}
void DrawAcceptance(TCanvas& canvas, TH2* acceptance[4], const char* outputPdf) {
  canvas.Clear(); canvas.Divide(2, 2);
  for (Int_t i = 0; i < 4; ++i) {
    canvas.cd(i + 1);
    acceptance[i]->SetStats(kFALSE);
    acceptance[i]->Draw("COLZ");
  }
  canvas.Print(outputPdf);
}
}

Int_t checkHistAnaFemtoLambda(const char* rootFile, const char* mainconf, const char* outputPdf) {
  if (!rootFile || !*rootFile || !mainconf || !*mainconf || !outputPdf || !*outputPdf) return 1;
  if (!gSystem->AccessPathName(outputPdf)) {
    std::cerr << "ERROR: refusing existing QA PDF " << outputPdf << std::endl; return 1;
  }
  TString config(mainconf);
  if (!config.BeginsWith("/")) config = TString(gSystem->WorkingDirectory()) + "/" + config;
  std::cout << "[mainconf] (argument) " << config << std::endl;
  if (!ConfigManager::GetInstance().LoadConfig(config)) return 1;
  const FemtoConfig& settings = ConfigManager::GetInstance().GetFemtoConfig();
  if (settings.channels.empty() || settings.channels[0].partA != "lambda") return 1;
  const FemtoConfig::ChannelDef& base = settings.channels[0];
  const std::string species = base.partB;
  TString shortName;
  if (species == "deuteron") shortName = "d";
  else if (species == "triton") shortName = "t";
  else if (species == "he3") shortName = "3He";
  else if (species == "he4") shortName = "4He";
  else return 1;
  TFile file(rootFile, "READ");
  if (file.IsZombie()) return 1;
  TNamed* status = dynamic_cast<TNamed*>(file.Get("FemtoLambdaRunStatus"));
  if (status && TString(status->GetTitle()) != "completed") {
    std::cerr << "ERROR: input run was not completed" << std::endl; return 1;
  }
  const Int_t first = settings.cfCent9Min, last = settings.cfCent9Max;
  if (first < 0 || last > 8 || first > last || settings.cfRebinFactor < 1) return 1;
  TH1* inclusive = dynamic_cast<TH1*>(file.Get("hLambda_InvMass"));
  TH2* massSE = FemtoLambdaQa::SumMass(file, "true/hKstarMass_" + shortName + "_CentBin%d", "qaMassSE", first, last);
  TH2* massME = FemtoLambdaQa::SumMass(file, "mix/hKstarMass_Mixed_" + shortName + "_CentBin%d", "qaMassME", first, last);
  if (!inclusive || !massSE || !massME) {
    std::cerr << "ERROR: missing Lambda mass or full-mass legacy histogram contract" << std::endl;
    delete massSE; delete massME; return 1;
  }
  TH2* acceptance[4] = {0, 0, 0, 0};
  const Int_t extensions = FemtoLambdaQa::CheckAcceptanceAndSplit(
      file, base.name.c_str(), species.c_str(), acceptance);
  if (extensions < 0) { delete massSE; delete massME; return 1; }
  const TString outputDirectory = gSystem->DirName(outputPdf);
  if (gSystem->mkdir(outputDirectory.Data(), kTRUE) != 0
      && gSystem->AccessPathName(outputDirectory.Data())) {
    std::cerr << "ERROR: cannot create QA output directory " << outputDirectory << std::endl;
    delete massSE; delete massME; return 1;
  }
  TCanvas canvas("lambdaFemtoQa", "Lambda-nucleus QA", 1200, 900);
  canvas.Print(TString(outputPdf) + "[");
  canvas.Divide(2, 2);
  canvas.cd(1); inclusive->SetTitle("Inclusive selected KF Lambda;M(p#pi) (GeV/c^{2});Candidates"); inclusive->Draw("E");
  canvas.cd(2); massSE->SetTitle("SE pair Lambda mass vs k*;k* (GeV/c);M(p#pi) (GeV/c^{2})"); massSE->Draw("COLZ");
  canvas.cd(3); massME->SetTitle("ME pair Lambda mass vs k*;k* (GeV/c);M(p#pi) (GeV/c^{2})"); massME->Draw("COLZ");
  canvas.cd(4);
  TH1* pairMass = massSE->ProjectionY("qaPairLambdaMass", 0, massSE->GetNbinsX()+1);
  pairMass->SetDirectory(0);
  pairMass->SetTitle("Pair-associated Lambda mass (SE, all k* including flow);M(p#pi) (GeV/c^{2});Pairs");
  pairMass->Draw("E");
  canvas.Print(outputPdf);
  delete pairMass;

  Int_t result = 0;
  const char* suffixes[] = {"signal", "leftSB", "rightSB", "signalLow", "signalHigh"};
  const Int_t regions = extensions ? 5 : 3;
  for (Int_t region = 0; region < regions; ++region) {
    // Preserve the original four pages first, then add acceptance and the two split CFs.
    if (region == 3) FemtoLambdaQa::DrawAcceptance(canvas, acceptance, outputPdf);
    canvas.Clear(); canvas.Divide(2, 1);
    const TString channel = TString(base.name.c_str()) + "_" + suffixes[region];
    TH2* seCent = dynamic_cast<TH2*>(file.Get("hKstarSEVsCent_" + channel));
    TH2* meCent = dynamic_cast<TH2*>(file.Get("hKstarMEVsCent_" + channel));
    if (!seCent || !meCent) {
      std::cerr << "ERROR: missing SE/ME centrality histogram for " << channel << std::endl;
      result = 1; break;
    }
    TH1* se = seCent->ProjectionX("qaSE", seCent->GetYaxis()->FindBin(first), seCent->GetYaxis()->FindBin(last));
    TH1* me = meCent->ProjectionX("qaME", meCent->GetYaxis()->FindBin(first), meCent->GetYaxis()->FindBin(last));
    se->SetDirectory(0); me->SetDirectory(0);
    se->Rebin(settings.cfRebinFactor); me->Rebin(settings.cfRebinFactor);
    const Int_t normFirst = se->GetXaxis()->FindFixBin(base.normQMin);
    const Int_t normLast = se->GetXaxis()->FindFixBin(base.normQMax - 1e-9);
    const Double_t seNorm = se->Integral(normFirst, normLast), meNorm = me->Integral(normFirst, normLast);
    canvas.cd(1);
    se->SetTitle(TString::Format("%s, cent9 %d-%d;k* (GeV/c);Pairs", channel.Data(), first, last));
    se->SetLineColor(kBlack); me->SetLineColor(kRed+1);
    se->SetStats(kFALSE); me->SetStats(kFALSE);
    if (seNorm > 0 && meNorm > 0) me->Scale(seNorm / meNorm);
    se->SetMaximum(TMath::Max(1.0, 1.15 * TMath::Max(se->GetMaximum(), me->GetMaximum())));
    se->Draw("E"); me->Draw("HIST SAME");
    TLegend legend(0.53, 0.73, 0.89, 0.89); legend.AddEntry(se, "SE", "l");
    legend.AddEntry(me, seNorm > 0 && meNorm > 0 ? "ME (norm region scaled)" : "ME (unscaled)", "l"); legend.Draw();
    canvas.cd(2);
    TH1* cf = 0;
    TGraphErrors validCf;
    validCf.SetMarkerStyle(20);
    if (seNorm > 0 && meNorm > 0) {
      cf = static_cast<TH1*>(se->Clone("qaRawCF")); cf->SetDirectory(0); cf->Divide(me);
      cf->SetTitle("Raw CF (zero-ME bins omitted); k* (GeV/c);C(k*)");
      // Clone inherits SE display bounds; reset before querying ratio data maximum.
      cf->SetMaximum(); cf->SetStats(kFALSE);
      cf->SetMinimum(0); cf->SetMaximum(TMath::Max(1.5, 1.15*cf->GetMaximum()));
      cf->Draw("AXIS");
      for (Int_t bin = 1; bin <= cf->GetNbinsX(); ++bin) {
        if (me->GetBinContent(bin) <= 0) continue;
        const Int_t point = validCf.GetN();
        validCf.SetPoint(point, cf->GetBinCenter(bin), cf->GetBinContent(bin));
        validCf.SetPointError(point, 0, cf->GetBinError(bin));
      }
      validCf.Draw("P SAME");
      std::cout << "[QA] " << channel << " normalization SE=" << seNorm << " ME=" << meNorm << std::endl;
    } else {
      FemtoLambdaQa::Message("Insufficient normalization counts; CF is not defined.");
      std::cout << "[QA] " << channel << " CF undefined: zero normalization counts" << std::endl;
    }
    canvas.Print(outputPdf);
    delete cf; delete se; delete me;
  }
  canvas.Print(TString(outputPdf) + "]");
  delete massSE; delete massME;
  FileStat_t pdfInfo;
  if (gSystem->GetPathInfo(outputPdf, pdfInfo) != 0 || pdfInfo.fSize <= 0) {
    std::cerr << "ERROR: QA PDF was not written" << std::endl; return 1;
  }
  return result;
}
