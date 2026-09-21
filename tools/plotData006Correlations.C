// Diagnostic figures for a finalized DATA-006 package (ROOT 5 / C++98 compatible).
// Display rebinning always rebuilds C(k*) from rebinned SE and ME; native objects are untouched.
#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TParameter.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include "Data006Config.h"

#include <algorithm>
#include <cmath>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

#include <vector>

namespace {

const Int_t kColors[7] = {kBlack, kRed + 1, kBlue + 1, kGreen + 2,
                          kMagenta + 1, kOrange + 7, kCyan + 2};

TString ChannelLabel(const std::string& name) {
  if (name == "proton_proton") return "p-p";
  if (name == "proton_deuteron") return "p-d";
  if (name == "deuteron_deuteron") return "d-d";
  return name.c_str();
}

TString CentLabel(const Data006CentralityConfig& c) {
  return TString::Format("%d-%d%%", c.lowPercent, c.highPercent);
}
Bool_t CentLess(const Data006CentralityConfig* a, const Data006CentralityConfig* b) {
  if (a->lowPercent != b->lowPercent) return a->lowPercent < b->lowPercent;
  return a->highPercent < b->highPercent;
}

TH1D* DisplayCF(TFile* f, const std::string& channel, const std::string& cent,
                Int_t rebin, Double_t xmax) {
  const TString base = TString::Format("%s/%s", channel.c_str(), cent.c_str());
  TH1D* srcSE = (TH1D*)f->Get(base + "/hSE");
  TH1D* srcME = (TH1D*)f->Get(base + "/hME");
  TParameter<Double_t>* pAlpha =
      (TParameter<Double_t>*)f->Get(base + "/normalizationFactor");
  TParameter<Double_t>* pAlphaErr =
      (TParameter<Double_t>*)f->Get(base + "/normalizationError");
  if (!srcSE || !srcME || !pAlpha || !pAlphaErr) return 0;
  if (rebin < 1) rebin = 1;
  TH1D* se = (TH1D*)srcSE->Clone(TString::Format("displaySE_%s_%s", channel.c_str(), cent.c_str()));
  TH1D* me = (TH1D*)srcME->Clone(TString::Format("displayME_%s_%s", channel.c_str(), cent.c_str()));
  se->SetDirectory(0);
  me->SetDirectory(0);
  if (rebin > 1) { se->Rebin(rebin); me->Rebin(rebin); }
  TH1D* cf = (TH1D*)se->Clone(TString::Format("displayCF_%s_%s", channel.c_str(), cent.c_str()));
  cf->Reset();
  const Double_t alpha = pAlpha->GetVal();
  const Double_t alphaErr = pAlphaErr->GetVal();
  for (Int_t ib = 1; ib <= cf->GetNbinsX(); ++ib) {
    const Double_t s = se->GetBinContent(ib), es = se->GetBinError(ib);
    const Double_t m = me->GetBinContent(ib), em = me->GetBinError(ib);
    if (!(alpha > 0.0) || !(m > 0.0)) continue;
    const Double_t value = s / (alpha * m);
    const Double_t dS = 1.0 / (alpha * m);
    const Double_t dM = -s / (alpha * m * m);
    const Double_t dA = -s / (alpha * alpha * m);
    const Double_t error = std::sqrt(dS*dS*es*es + dM*dM*em*em +
                                     dA*dA*alphaErr*alphaErr);
    cf->SetBinContent(ib, value);
    cf->SetBinError(ib, error);
  }
  cf->GetXaxis()->SetRangeUser(0.0, xmax);
  delete se;
  delete me;
  return cf;
}

TString gRunLabel = "DATA-006 TEST";
TString gFilePrefix = "data006_test";

void DrawTestLabel(Int_t rebin, Long64_t eventCount) {
  TLatex label;
  label.SetNDC();
  label.SetTextSize(0.035);
  label.SetTextColor(kGray + 2);
  label.DrawLatex(0.13, 0.92, TString::Format(
      "%s; %lld events; display rebin = %d #times 5 MeV/c",
      gRunLabel.Data(), eventCount, rebin));
}

void DrawUnity(Double_t xmax) {
  TLine* line = new TLine(0.0, 1.0, xmax, 1.0);
  line->SetLineStyle(2);
  line->SetLineColor(kGray + 2);
  line->Draw();
}

}  // namespace

void plotData006Correlations(const Char_t* packageDirectory, const Char_t* data006Config,
                             const Char_t* outputDirectory, Int_t displayRebin = 4,
                             Double_t cfYMin = 0.0, Double_t cfYMax = 2.5,
                             const Char_t* runLabel = "DATA-006 TEST",
                             const Char_t* filePrefix = "data006_test") {
  gRunLabel = runLabel;
  gFilePrefix = filePrefix;
  Data006Config cfg;
  if (!cfg.Load(data006Config)) {
    std::cerr << "ERROR: invalid DATA-006 config: " << data006Config << std::endl;
    gSystem->Exit(1);
  }
  if (displayRebin < 1 || cfg.kstarBins % displayRebin != 0 || !(cfYMax > cfYMin)) {
    std::cerr << "ERROR: invalid display rebin or CF y range" << std::endl;
    gSystem->Exit(1);
  }
  const TString rootPath = TString::Format("%s/data006_correlations.root", packageDirectory);
  TFile* f = TFile::Open(rootPath, "READ");
  if (!f || f->IsZombie()) {
    std::cerr << "ERROR: cannot open " << rootPath << std::endl;
    gSystem->Exit(1);
  }
  if (gSystem->mkdir(outputDirectory, kTRUE) != 0 && gSystem->AccessPathName(outputDirectory)) {
    std::cerr << "ERROR: cannot create " << outputDirectory << std::endl;
    gSystem->Exit(1);
  }
  gStyle->SetOptStat(0);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetPadGridX(kTRUE);
  gStyle->SetPadGridY(kTRUE);
  std::ofstream acceptanceCsv(
      TString::Format("%s/acceptance_shape_diagnostics.csv", outputDirectory).Data());
  if (!acceptanceCsv) {
    std::cerr << "ERROR: cannot create acceptance diagnostics CSV" << std::endl;
    gSystem->Exit(1);
  }
  acceptanceCsv << std::setprecision(17)
                << "channel,centrality,se_pairs,me_pairs,ks_max_distance,total_variation,status\n";

  std::vector<const Data006CentralityConfig*> comparison;
  std::vector<const Data006CentralityConfig*> native;
  for (size_t i = 0; i < cfg.centralities.size(); ++i) {
    const Data006CentralityConfig& c = cfg.centralities[i];
    if (c.nativeClass) native.push_back(&c);
    if (!c.nativeClass || (c.lowPercent == 10 && c.highPercent == 20)) comparison.push_back(&c);
  }
  std::sort(comparison.begin(), comparison.end(), CentLess);
  TParameter<Long64_t>* pEventCount =
      (TParameter<Long64_t>*)f->Get("acceptedEventCount");
  const Long64_t eventCount = pEventCount ? pEventCount->GetVal() : 0;
  const Double_t xmax = cfg.deliveryKstarMax;
  const TString pdfPath = TString::Format("%s/%s_correlation_qa.pdf", outputDirectory, gFilePrefix.Data());

  // Page 1: paper-comparison classes overlaid for each channel.
  TCanvas overview("data006Overview", "DATA-006 comparison centralities", 1800, 600);
  overview.Divide(3, 1);
  for (size_t ich = 0; ich < cfg.channels.size(); ++ich) {
    if (!cfg.channels[ich].enabled) continue;
    overview.cd((Int_t)ich + 1);
    TLegend* leg = new TLegend(0.57, 0.66, 0.88, 0.88);
    leg->SetBorderSize(0);
    Bool_t first = kTRUE;
    for (size_t iz = 0; iz < comparison.size(); ++iz) {
      TH1D* cf = DisplayCF(f, cfg.channels[ich].name, comparison[iz]->id,
                           displayRebin, xmax);
      if (!cf) continue;
      cf->SetLineColor(kColors[iz % 7]);
      cf->SetMarkerColor(kColors[iz % 7]);
      cf->SetMarkerStyle(20 + iz);
      cf->SetMarkerSize(0.65);
      cf->SetMinimum(cfYMin);
      cf->SetMaximum(cfYMax);
      cf->SetTitle(TString::Format("%s comparison classes;k* (GeV/c);C(k*)",
                                   ChannelLabel(cfg.channels[ich].name).Data()));
      cf->Draw(first ? "E1" : "E1 SAME");
      leg->AddEntry(cf, CentLabel(*comparison[iz]), "lep");
      first = kFALSE;
    }
    DrawUnity(xmax);
    leg->Draw();
    DrawTestLabel(displayRebin, eventCount);
  }
  overview.Print(pdfPath + "(");
  overview.SaveAs(TString::Format("%s/%s_cf_overview.png", outputDirectory, gFilePrefix.Data()));

  // Pages 2-4: every native centrality, one page per channel.
  for (size_t ich = 0; ich < cfg.channels.size(); ++ich) {
    if (!cfg.channels[ich].enabled) continue;
    TCanvas nativeCanvas(TString::Format("native_%s", cfg.channels[ich].name.c_str()),
                         cfg.channels[ich].name.c_str(), 1600, 900);
    nativeCanvas.Divide(4, 2);
    for (size_t iz = 0; iz < native.size() && iz < 8; ++iz) {
      nativeCanvas.cd((Int_t)iz + 1);
      TH1D* cf = DisplayCF(f, cfg.channels[ich].name, native[iz]->id,
                           displayRebin, xmax);
      if (!cf) continue;
      cf->SetMarkerStyle(20);
      cf->SetMarkerSize(0.55);
      cf->SetMinimum(cfYMin);
      cf->SetMaximum(cfYMax);
      cf->SetTitle(TString::Format("%s %s;k* (GeV/c);C(k*)",
                                   ChannelLabel(cfg.channels[ich].name).Data(),
                                   CentLabel(*native[iz]).Data()));
      cf->Draw("E1");
      DrawUnity(xmax);
      DrawTestLabel(displayRebin, eventCount);
    }
    nativeCanvas.Print(pdfPath);
    nativeCanvas.SaveAs(TString::Format("%s/%s_cf_native_%s.png", outputDirectory, gFilePrefix.Data(),
                                        cfg.channels[ich].name.c_str()));
  }

  // Page 5: pair-y acceptance in the normalization interval, SE vs shape-scaled ME.
  TCanvas acceptance("data006Acceptance", "normalization-region pair-y acceptance", 1800, 900);
  acceptance.Divide(4, 3);
  Int_t pad = 1;
  for (size_t ich = 0; ich < cfg.channels.size(); ++ich) {
    if (!cfg.channels[ich].enabled) continue;
    for (size_t iz = 0; iz < comparison.size() && pad <= 12; ++iz, ++pad) {
      acceptance.cd(pad);
      const TString base = TString::Format("%s/%s", cfg.channels[ich].name.c_str(),
                                           comparison[iz]->id.c_str());
      TH1D* se = (TH1D*)f->Get(base + "/hPairRapidityNormSE");
      TH1D* me = (TH1D*)f->Get(base + "/hPairRapidityNormME");
      if (!se || !me) continue;
      TH1D* seDraw = (TH1D*)se->Clone(TString::Format("accSE_%d", pad));
      TH1D* meDraw = (TH1D*)me->Clone(TString::Format("accME_%d", pad));
      seDraw->SetDirectory(0); meDraw->SetDirectory(0);
      seDraw->Rebin(4); meDraw->Rebin(4);
      const Double_t seIntegral = seDraw->Integral();
      const Double_t meIntegral = meDraw->Integral();
      Double_t totalVariation = 0.0;
      if (seIntegral > 0.0 && meIntegral > 0.0) {
        for (Int_t ib = 1; ib <= seDraw->GetNbinsX(); ++ib)
          totalVariation += TMath::Abs(seDraw->GetBinContent(ib) / seIntegral -
                                      meDraw->GetBinContent(ib) / meIntegral);
        totalVariation *= 0.5;
      }
      const Double_t ksDistance = (seIntegral > 0.0 && meIntegral > 0.0)
          ? seDraw->KolmogorovTest(meDraw, "M") : -1.0;
      if (meIntegral > 0.0) meDraw->Scale(seIntegral / meIntegral);
      seDraw->SetLineColor(kBlack); meDraw->SetLineColor(kRed + 1);
      acceptanceCsv << cfg.channels[ich].name << ',' << comparison[iz]->id << ','
                    << seIntegral << ',' << meIntegral << ',' << ksDistance << ','
                    << totalVariation << ','
                    << (seIntegral > 0.0 && meIntegral > 0.0 ? "OK" : "EMPTY") << '\n';
      seDraw->SetTitle(TString::Format("%s %s;pair y_{CM};pairs",
                                      ChannelLabel(cfg.channels[ich].name).Data(),
                                      CentLabel(*comparison[iz]).Data()));
      seDraw->Draw("HIST"); meDraw->Draw("HIST SAME");
      TLatex label; label.SetNDC(); label.SetTextSize(0.04);
      label.DrawLatex(0.15, 0.86, "SE black; shape-scaled ME red");
    }
  }
  acceptance.Print(pdfPath);
  acceptance.SaveAs(TString::Format("%s/%s_norm_acceptance.png", outputDirectory, gFilePrefix.Data()));

  // Final page: SE pair-mT/2 distributions for comparison centralities.
  TCanvas mtCanvas("data006PairMt", "pair mT", 1800, 600);
  mtCanvas.Divide(3, 1);
  for (size_t ich = 0; ich < cfg.channels.size(); ++ich) {
    if (!cfg.channels[ich].enabled) continue;
    mtCanvas.cd((Int_t)ich + 1);
    TLegend* leg = new TLegend(0.58, 0.66, 0.88, 0.88);
    leg->SetBorderSize(0);
    Bool_t first = kTRUE;
    for (size_t iz = 0; iz < comparison.size(); ++iz) {
      const TString path = TString::Format("%s/%s/hPairMtHalfSE",
                                           cfg.channels[ich].name.c_str(),
                                           comparison[iz]->id.c_str());
      TH1D* src = (TH1D*)f->Get(path);
      if (!src || src->Integral() <= 0.0) continue;
      TH1D* h = (TH1D*)src->Clone(TString::Format("mt_%d_%d", (Int_t)ich, (Int_t)iz));
      h->SetDirectory(0); h->Rebin(3); h->Scale(1.0 / h->Integral("width"));
      h->SetLineColor(kColors[iz % 7]); h->SetLineWidth(2);
      h->SetTitle(TString::Format("%s SE pair m_{T}/2;pair m_{T}/2 (GeV/c^{2});density",
                                  ChannelLabel(cfg.channels[ich].name).Data()));
      h->Draw(first ? "HIST" : "HIST SAME");
      leg->AddEntry(h, CentLabel(*comparison[iz]), "l");
      first = kFALSE;
    }
    leg->Draw();
    TLatex label; label.SetNDC(); label.SetTextSize(0.035); label.SetTextColor(kGray + 2);
    label.DrawLatex(0.13, 0.92, TString::Format("%s; k* < 500 MeV/c", gRunLabel.Data()));
  }
  mtCanvas.Print(pdfPath + ")");
  mtCanvas.SaveAs(TString::Format("%s/%s_pair_mt.png", outputDirectory, gFilePrefix.Data()));

  acceptanceCsv.close();
  std::cout << "[plotData006Correlations] PDF=" << pdfPath << std::endl;
  f->Close();
}
