// checkHistAnaFemtoPhiProton.C - Draw histograms from run_anaFemtoPhiProton.C output and write PDF.
// Invoke via: ./script/singularity_checkHistAnaFemtoPhiProton.sh <root_file> <mainconf_path>
// Or: root4star -b -q 'analysis/run_checkHistAnaFemtoPhiProton.C("rootfile/...","anaName","config/mainconf/...")'
// If input is anaName_jobid_merge.root, jobid (32 hex) is parsed and PDF becomes anaName_checkHistAnaFemtoPhiProton_jobid.pdf.
// With config loaded, cut regions are overlaid on pre-cut histograms.

#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TGraphErrors.h>
#include <TLine.h>
#include <TMath.h>
#include <TObject.h>
#include <TString.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TEllipse.h>
#include <iostream>
#include <vector>
#include <limits.h>
#include <stdlib.h>

#include "../../include/PdfIOMan.h"
#include "ConfigManager.h"
#include "cuts/EventCutConfig.h"
#include "cuts/TrackCutConfig.h"
#include "cuts/PhiCutConfig.h"
#include "cuts/CentralityCutConfig.h"
#include "cuts/PIDCutConfig.h"
#include "cuts/FemtoConfig.h"
#include <map>

static Bool_t gConfigLoaded = kFALSE;

static TString resolveFigureRoot(const char* pwd) {
  const char* envFigureRoot = gSystem->Getenv("STAR_QA_FIGURE_ROOT");
  if (envFigureRoot && envFigureRoot[0] != '\0') {
    return TString(envFigureRoot);
  }
  TString figureRoot = TString(pwd ? pwd : ".") + "/share/figure";
  char resolved[PATH_MAX];
  if (realpath(figureRoot.Data(), resolved)) {
    return TString(resolved);
  }
  return figureRoot;
}

static void drawCutLines1D(TH1* h, Double_t x1, Double_t x2, Int_t color = kRed, Int_t style = 2) {
  if (!h || !gPad) return;
  Double_t ylo = gPad->GetUymin();
  Double_t yhi = gPad->GetUymax();
  const Double_t kBig = 1e30;
  if (x1 > -kBig && x1 < kBig) {
    TLine* l1 = new TLine(x1, ylo, x1, yhi);
    l1->SetLineColor(color);
    l1->SetLineStyle(style);
    l1->Draw("same");
  }
  if (x2 > -kBig && x2 < kBig && TMath::Abs(x2 - x1) > 1e-9) {
    TLine* l2 = new TLine(x2, ylo, x2, yhi);
    l2->SetLineColor(color);
    l2->SetLineStyle(style);
    l2->Draw("same");
  }
}

static void drawCutLine1D(TH1* h, Double_t x, Int_t color = kRed, Int_t style = 2) {
  if (!h || !gPad) return;
  Double_t ylo = gPad->GetUymin();
  Double_t yhi = gPad->GetUymax();
  TLine* l = new TLine(x, ylo, x, yhi);
  l->SetLineColor(color);
  l->SetLineStyle(style);
  l->Draw("same");
}

static void drawCutLine2DH(TH2* h, Double_t yVal, Int_t color = kRed, Int_t style = 2) {
  if (!h || !gPad) return;
  Double_t xlo = gPad->GetUxmin();
  Double_t xhi = gPad->GetUxmax();
  TLine* l = new TLine(xlo, yVal, xhi, yVal);
  l->SetLineColor(color);
  l->SetLineStyle(style);
  l->Draw("same");
}

static void drawVtxCutCircle(Double_t centerX, Double_t centerY, Double_t radius,
                             Int_t color = kRed, Int_t style = 2) {
  if (!gPad || radius <= 0.0) return;
  TEllipse* circle = new TEllipse(centerX, centerY, radius, radius);
  circle->SetFillStyle(0);
  circle->SetLineColor(color);
  circle->SetLineStyle(style);
  circle->Draw("same");
}

static void drawCent9ConventionNote() {
  if (!gPad) return;
  TLatex* note = new TLatex();
  note->SetNDC(kTRUE);
  note->SetTextSize(0.028);
  note->SetTextColor(kBlue + 1);
  note->DrawLatex(0.12, 0.96, "cent9: StRefMultCorr (0=peripheral, 8=central)");
}

static Bool_t isHex32(const TString& s) {
  if (s.Length() != 32) return kFALSE;
  for (Int_t i = 0; i < 32; i++) {
    Char_t c = s[i];
    if (!((c >= '0' && c <= '9') || (c >= 'A' && c <= 'F') || (c >= 'a' && c <= 'f')))
      return kFALSE;
  }
  return kTRUE;
}

static Double_t getHistEntries(TFile* fin, const char* key) {
  if (!fin || !key || key[0] == '\0') return -1.0;
  TObject* obj = fin->Get(key);
  if (!obj || !obj->InheritsFrom("TH1")) return -1.0;
  return ((TH1*)obj)->GetEntries();
}

static std::string cfHistKey(const std::string& channel) {
  return std::string("hCF_") + channel;
}

static std::string kstarSeHistKey(const std::string& channel) {
  return std::string("hKstarSE_") + channel;
}

static std::string kstarMeHistKey(const std::string& channel) {
  return std::string("hKstarME_") + channel;
}

static std::string kstarSeVsCentHistKey(const std::string& channel) {
  return std::string("hKstarSEVsCent_") + channel;
}

static std::string kstarMeVsCentHistKey(const std::string& channel) {
  return std::string("hKstarMEVsCent_") + channel;
}

static std::string cfCentCacheKey(const std::string& channel) {
  return std::string("cent:") + channel;
}

static void getCfCent9Range(Int_t& cent9Min, Int_t& cent9Max) {
  cent9Min = 2;
  cent9Max = 8;
  if (gConfigLoaded) {
    const FemtoConfig& femtoCfg = ConfigManager::GetInstance().GetFemtoConfig();
    cent9Min = femtoCfg.cfCent9Min;
    cent9Max = femtoCfg.cfCent9Max;
  }
}

static TH1* projectKstarVsCent9(TH2* h2, Int_t cent9Min, Int_t cent9Max, const char* hname) {
  if (!h2) return 0;
  Int_t yLo = h2->GetYaxis()->FindBin((Double_t)cent9Min - 0.01);
  Int_t yHi = h2->GetYaxis()->FindBin((Double_t)cent9Max + 0.01);
  TH1* h1 = h2->ProjectionX(hname, yLo, yHi);
  if (h1) {
    h1->SetDirectory(0);
    h1->SetTitle(Form("%s (cent9 %d-%d);k^{*} [GeV/c];Counts", h2->GetTitle(), cent9Min, cent9Max));
  }
  return h1;
}

static Int_t getCfRebinFactor() {
  const Int_t kFallbackCfRebinFactor = 5;
  if (!gConfigLoaded) return kFallbackCfRebinFactor;
  const FemtoConfig& femtoCfg = ConfigManager::GetInstance().GetFemtoConfig();
  return (femtoCfg.cfRebinFactor >= 1) ? femtoCfg.cfRebinFactor : 1;
}

static TH1* rebinHistCopy(TH1* h, Int_t factor, const char* cloneSuffix) {
  if (!h || factor <= 1) return 0;
  const Int_t nBins = h->GetNbinsX();
  if (nBins % factor != 0) {
    std::cout << "[checkHistAnaFemtoPhiProton] WARNING: cannot rebin " << h->GetName() << " with factor "
              << factor << " (nbins=" << nBins << ")\n";
    return 0;
  }
  TH1* hRb = (TH1*)h->Clone(cloneSuffix ? cloneSuffix : "_rebinned");
  hRb->SetDirectory(0);
  hRb->Rebin(factor);
  return hRb;
}

static TGraphErrors* computeCfGraphFromSeMe(TH1* hSE, TH1* hME, Double_t normQMin, Double_t normQMax,
                                            const char* graphTitle) {
  if (!hSE || !hME) return 0;
  Int_t binLo = hSE->FindBin(normQMin + 1e-9);
  Int_t binHi = hSE->FindBin(normQMax - 1e-9);
  Double_t seNorm = hSE->Integral(binLo, binHi);
  Double_t meNorm = hME->Integral(binLo, binHi);
  if (seNorm <= 0 || meNorm <= 0) return 0;

  const Double_t scale = meNorm / seNorm;
  std::vector<Double_t> x;
  std::vector<Double_t> y;
  std::vector<Double_t> ey;
  x.reserve(hSE->GetNbinsX());
  y.reserve(hSE->GetNbinsX());
  ey.reserve(hSE->GetNbinsX());

  for (Int_t ib = 1; ib <= hSE->GetNbinsX(); ++ib) {
    Double_t se = hSE->GetBinContent(ib);
    Double_t me = hME->GetBinContent(ib);
    if (se <= 0 || me <= 0) continue;
    Double_t cf = scale * se / me;
    Double_t err = cf * TMath::Sqrt(1.0 / se + 1.0 / me);
    x.push_back(hSE->GetBinCenter(ib));
    y.push_back(cf);
    ey.push_back(err);
  }
  if (x.empty()) return 0;

  TGraphErrors* gCF = new TGraphErrors((Int_t)x.size(), &x[0], &y[0], 0, &ey[0]);
  gCF->SetTitle(graphTitle ? graphTitle : "C(k^{*})");
  gCF->SetMarkerStyle(20);
  gCF->SetMarkerSize(0.8);
  gCF->SetLineColor(kBlack);
  gCF->SetMarkerColor(kBlack);
  return gCF;
}

static TGraphErrors* computeCfAndCache(const std::string& channel, const std::string& cacheKey, TH1* hSE, TH1* hME,
                                       Double_t normQMin, Double_t normQMax, std::map<std::string, TGraphErrors*>& cfCache,
                                       const char* logTag) {
  std::map<std::string, TGraphErrors*>::const_iterator cached = cfCache.find(cacheKey);
  if (cached != cfCache.end()) return cached->second;

  if (!hSE || !hME) {
    std::cout << "[checkHistAnaFemtoPhiProton] CF " << channel.c_str();
    if (logTag) std::cout << " (" << logTag << ")";
    std::cout << ": missing SE/ME histograms\n";
    cfCache[cacheKey] = 0;
    return 0;
  }

  const Int_t cfRebinFactor = getCfRebinFactor();
  TH1* hSEForCf = hSE;
  TH1* hMEForCf = hME;
  TH1* hSERebinned = rebinHistCopy(hSE, cfRebinFactor, "_se_rebin");
  TH1* hMERebinned = rebinHistCopy(hME, cfRebinFactor, "_me_rebin");
  if (hSERebinned) hSEForCf = hSERebinned;
  if (hMERebinned) hMEForCf = hMERebinned;

  Int_t binLo = hSEForCf->FindBin(normQMin + 1e-9);
  Int_t binHi = hSEForCf->FindBin(normQMax - 1e-9);
  Double_t seNorm = hSEForCf->Integral(binLo, binHi);
  Double_t meNorm = hMEForCf->Integral(binLo, binHi);
  Double_t aNorm = (seNorm > 0) ? meNorm / seNorm : 0.0;
  std::cout << "[checkHistAnaFemtoPhiProton] CF " << channel.c_str();
  if (logTag) std::cout << " (" << logTag << ")";
  std::cout << ": normQ=[" << normQMin << ", " << normQMax << "] GeV/c";
  if (cfRebinFactor > 1 && hSERebinned && hMERebinned) {
    std::cout << ", cfRebinFactor=" << cfRebinFactor << " (" << hSE->GetNbinsX() << " -> "
              << hSEForCf->GetNbinsX() << " bins)";
  }
  std::cout << ", seNorm=" << seNorm << ", meNorm=" << meNorm << ", a=meNorm/seNorm=" << aNorm << std::endl;

  TString cfTitle = Form("CF %s (checkHist);k^{*} [GeV/c];C(k^{*})", channel.c_str());
  if (logTag) cfTitle = Form("CF %s %s;k^{*} [GeV/c];C(k^{*})", channel.c_str(), logTag);
  TGraphErrors* gCF = computeCfGraphFromSeMe(hSEForCf, hMEForCf, normQMin, normQMax, cfTitle.Data());
  delete hSERebinned;
  delete hMERebinned;
  if (!gCF) {
    std::cout << "[checkHistAnaFemtoPhiProton] CF " << channel.c_str() << ": failed (empty norm integrals)\n";
  } else {
    std::cout << "[checkHistAnaFemtoPhiProton] CF " << channel.c_str() << ": " << gCF->GetN()
              << " points with Poisson stat errors\n";
  }
  cfCache[cacheKey] = gCF;
  return gCF;
}

static TGraphErrors* getOrComputeCf(TFile* fin, const std::string& channel, Double_t normQMin, Double_t normQMax,
                                    std::map<std::string, TGraphErrors*>& cfCache) {
  TH1* hSE = (TH1*)fin->Get(kstarSeHistKey(channel).c_str());
  TH1* hME = (TH1*)fin->Get(kstarMeHistKey(channel).c_str());
  return computeCfAndCache(channel, channel, hSE, hME, normQMin, normQMax, cfCache, 0);
}

static TH1* getProjectedSeMeFromCent(TFile* fin, const std::string& channel, Bool_t isSE, Int_t cent9Min,
                                   Int_t cent9Max) {
  const std::string key = isSE ? kstarSeVsCentHistKey(channel) : kstarMeVsCentHistKey(channel);
  TH2* h2 = (TH2*)fin->Get(key.c_str());
  if (!h2) return 0;
  TString hname = Form("_proj_%s_%s_cent%d_%d", isSE ? "se" : "me", channel.c_str(), cent9Min, cent9Max);
  return projectKstarVsCent9(h2, cent9Min, cent9Max, hname.Data());
}

static TGraphErrors* getOrComputeCfCentSlice(TFile* fin, const std::string& channel, Double_t normQMin,
                                             Double_t normQMax, Int_t cent9Min, Int_t cent9Max,
                                             std::map<std::string, TGraphErrors*>& cfCache) {
  TH1* hSE = getProjectedSeMeFromCent(fin, channel, kTRUE, cent9Min, cent9Max);
  TH1* hME = getProjectedSeMeFromCent(fin, channel, kFALSE, cent9Min, cent9Max);
  TString logTag = Form("cent9 %d-%d from hKstar*VsCent", cent9Min, cent9Max);
  TGraphErrors* gCF =
      computeCfAndCache(channel, cfCentCacheKey(channel), hSE, hME, normQMin, normQMax, cfCache, logTag.Data());
  delete hSE;
  delete hME;
  return gCF;
}

static Double_t channelNormQMin(const std::string& channel) {
  if (gConfigLoaded) {
    const FemtoConfig& femtoCfg = ConfigManager::GetInstance().GetFemtoConfig();
    const FemtoConfig::ChannelDef* ch = femtoCfg.FindChannel(channel);
    if (ch) return ch->normQMin;
  }
  return 0.5;
}

static Double_t channelNormQMax(const std::string& channel) {
  if (gConfigLoaded) {
    const FemtoConfig& femtoCfg = ConfigManager::GetInstance().GetFemtoConfig();
    const FemtoConfig::ChannelDef* ch = femtoCfg.FindChannel(channel);
    if (ch) return ch->normQMax;
  }
  return 1.0;
}

static void populateCfCache(TFile* fin, std::map<std::string, TGraphErrors*>& cfCache) {
  const Double_t kFallbackNormQMin = 0.5;
  const Double_t kFallbackNormQMax = 1.0;
  if (gConfigLoaded) {
    const FemtoConfig& femtoCfg = ConfigManager::GetInstance().GetFemtoConfig();
    for (size_t ic = 0; ic < femtoCfg.channels.size(); ic++) {
      const FemtoConfig::ChannelDef& ch = femtoCfg.channels[ic];
      if (!ch.enabled) continue;
      getOrComputeCf(fin, ch.name, ch.normQMin, ch.normQMax, cfCache);
    }
    return;
  }
  const char* fallbackChannels[] = {"phi_proton",         "phi_proton_signal", "phi_proton_leftSB",
                                    "phi_proton_rightSB", "phi_rot_proton",    0};
  for (Int_t i = 0; fallbackChannels[i]; ++i) {
    getOrComputeCf(fin, fallbackChannels[i], kFallbackNormQMin, kFallbackNormQMax, cfCache);
  }
}

static void populateCfCentCache(TFile* fin, std::map<std::string, TGraphErrors*>& cfCache) {
  Int_t cent9Min = 0;
  Int_t cent9Max = 0;
  getCfCent9Range(cent9Min, cent9Max);
  const char* centChannels[] = {"phi_proton_signal", "phi_proton_leftSB", "phi_proton_rightSB", "phi_rot_proton", 0};
  for (Int_t i = 0; centChannels[i]; ++i) {
    const std::string channel(centChannels[i]);
    getOrComputeCfCentSlice(fin, channel, channelNormQMin(channel), channelNormQMax(channel), cent9Min, cent9Max,
                            cfCache);
  }
}

static void drawCentProjectedSeMe(TCanvas* canvas, Int_t pad, TFile* fin, const std::string& channel, Bool_t isSE,
                                  Int_t cent9Min, Int_t cent9Max) {
  if (!canvas) return;
  canvas->cd(pad);
  TH1* h = getProjectedSeMeFromCent(fin, channel, isSE, cent9Min, cent9Max);
  if (h) {
    h->Draw();
    delete h;
  }
}

static void drawCentSliceCf(TCanvas* canvas, Int_t pad, const std::string& channel,
                            std::map<std::string, TGraphErrors*>& cfCache) {
  if (!canvas) return;
  canvas->cd(pad);
  std::map<std::string, TGraphErrors*>::const_iterator it = cfCache.find(cfCentCacheKey(channel));
  if (it != cfCache.end() && it->second) it->second->Draw("AP");
}

static void drawComputedCf(TCanvas* canvas, Int_t pad, TFile* fin, const std::string& channel, Double_t normQMin,
                           Double_t normQMax, std::map<std::string, TGraphErrors*>& cfCache) {
  if (!canvas) return;
  canvas->cd(pad);
  TGraphErrors* gCF = getOrComputeCf(fin, channel, normQMin, normQMax, cfCache);
  if (gCF) gCF->Draw("AP");
}

static Int_t getCachedCfPointCount(const std::map<std::string, TGraphErrors*>& cfCache, const char* channel) {
  std::map<std::string, TGraphErrors*>::const_iterator it = cfCache.find(channel);
  if (it == cfCache.end() || !it->second) return -1;
  return it->second->GetN();
}

static void freeCfCache(std::map<std::string, TGraphErrors*>& cfCache) {
  for (std::map<std::string, TGraphErrors*>::iterator it = cfCache.begin(); it != cfCache.end(); ++it) {
    if (it->second) delete it->second;
  }
  cfCache.clear();
}

void checkHistAnaFemtoPhiProton(const Char_t* inputRootFile,
                                const Char_t* anaNameArg = "auau3p85fxt_anaFemtoPhiProton",
                                const Char_t* mainconfPath = 0)
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

  // Parse input basename for anaName_jobid_merge.root -> anaName, jobid (32 hex)
  TString anaName(anaNameArg);
  TString jobid;
  TString base = gSystem->BaseName(inputRootFile);
  base.ReplaceAll(".root", "");
  std::vector<TString> tokens;
  for (Int_t i = 0; i < base.Length(); ) {
    Int_t j = base.Index("_", i);
    if (j < 0) {
      tokens.push_back(TString(base(i, base.Length() - i)));
      break;
    }
    tokens.push_back(TString(base(i, j - i)));
    i = j + 1;
  }
  for (size_t k = 0; k < tokens.size(); k++) {
    if (isHex32(tokens[k])) {
      jobid = tokens[k];
      anaName = tokens[0];
      for (size_t m = 1; m < k; m++) anaName += "_" + tokens[m];
      break;
    }
  }

  gConfigLoaded = kFALSE;
  const char* pwd = gSystem->Getenv("PWD");
  if (!pwd) pwd = ".";
  if (gSystem->Load(TString(pwd) + "/lib/libStarAnaConfig.so") >= 0) {
    TString mainconf;
    if (mainconfPath && strlen(mainconfPath) > 0) {
      mainconf = mainconfPath;
      if (mainconf[0] != '/') mainconf = TString(pwd) + "/" + mainconf;
    } else {
      mainconf = TString(pwd) + "/config/mainconf/main_" + anaName + ".yaml";
    }
    if (ConfigManager::GetInstance().LoadConfig(mainconf.Data())) {
      gConfigLoaded = kTRUE;
      if (!ConfigManager::GetInstance().GetPhiCuts().FinalizeRapidityFrame(
              ConfigManager::GetInstance().GetCentralityCuts())) {
        std::cerr << "[checkHistAnaFemtoPhiProton] WARNING: FinalizeRapidityFrame failed; PDF note may be incomplete."
                  << std::endl;
      }
    } else if (mainconfPath && strlen(mainconfPath) > 0) {
      std::cerr << "[checkHistAnaFemtoPhiProton] WARNING: Failed to load config " << mainconf.Data() << "; cut lines skipped." << std::endl;
    }
  }

  // Resolve symlinked share/figure to a physical path so ROOT PDF output
  // works the same in host and singularity-based runs.
  TString outDir = resolveFigureRoot(pwd) + "/" + anaName + "/";
  if (gSystem->AccessPathName(outDir)) {
    gSystem->mkdir(outDir, kTRUE);
  }

  TString pdfName = TString(outDir) + anaName + "_checkHistAnaFemtoPhiProton";
  if (jobid.Length()) pdfName += "_" + jobid;
  pdfName += ".pdf";
  std::cout << "Output PDF: " << pdfName.Data() << std::endl;

  PdfHeader::OpenPdf(pdfName);

  std::vector<std::string> inputs;
  inputs.push_back((const char*)inputRootFile);

  TString note = "Check histograms from run_anaFemtoPhiProton.C (StFemtoMaker output).\n";
  note += "Phi candidates built via ResonanceBuilder (KK); proton candidates via TrackPidBuilder.\n";
  note += "Pair QA stage0 / strict TOF histograms mirror anaPhi hPhiPair_* naming.\n";
  const Double_t nEvtAll = getHistEntries(fin, "hVz");
  const Double_t nEvtAfter = getHistEntries(fin, "hVz_After");
  const Double_t nPhiCand = getHistEntries(fin, "hPhi_NCand");
  const Double_t nProtonCand = getHistEntries(fin, "hP_NCand");
  const Double_t nKstarSE = getHistEntries(fin, "hKstarSE_phi_proton");
  if (nEvtAll >= 0.0) {
    note += Form("Total statistics (input ROOT): events(all) = %.0f", nEvtAll);
    if (nEvtAfter >= 0.0) {
      note += Form(", events(after event cuts) = %.0f", nEvtAfter);
    }
    if (nPhiCand >= 0.0) {
      note += Form(", phi cand fills = %.0f", nPhiCand);
    }
    if (nProtonCand >= 0.0) {
      note += Form(", proton cand fills = %.0f", nProtonCand);
    }
    if (nKstarSE >= 0.0) {
      note += Form(", k* SE pairs = %.0f", nKstarSE);
    }
    note += "\n";
    TH1* hVzCheck = (TH1*)fin->Get("hVz");
    if (hVzCheck && nEvtAll > 0.0 && hVzCheck->GetMaximum() < 1e-6) {
      note += "WARNING: hVz has entries but visible bins are empty — likely Vz axis mismatch (FXT ~200 cm). Re-run with updated hist YAML.\n";
    }
  } else {
    note += "Total statistics (input ROOT): hVz not found.\n";
  }
  if (gConfigLoaded) {
    note += ConfigManager::GetInstance().GetPhiCuts().GetRapidityFrameSummary().c_str();
    note += "\n";
  }
  note += "QA: Event/centrality/track/PID pages follow anaPhi layout; femto pages at end (proton, phi cand, k*).\n";
  note += Form("CF computed in checkHist from merged SE/ME (TGraphErrors, Poisson stat errors); cfRebinFactor=%d; norm region from maker YAML.\n",
               getCfRebinFactor());
  Int_t cfCent9MinNote = 2;
  Int_t cfCent9MaxNote = 8;
  getCfCent9Range(cfCent9MinNote, cfCent9MaxNote);
  note += Form("CF cent slice: cent9 [%d, %d] projected from hKstarSEVsCent/hKstarMEVsCent (0-60%% for default 2-8).\n",
               cfCent9MinNote, cfCent9MaxNote);
  if (gConfigLoaded) {
    const CentralityCutConfig& centCfg = ConfigManager::GetInstance().GetCentralityCuts();
    if (centCfg.cent9MaxRefMultCorrBin >= 0 && centCfg.cent9MaxRefMultCorr > 0.0) {
      note += Form("Centrality cut (YAML): cent9=%d requires refMult_corr <= %.1f (Maker; re-run needed in ROOT).\n",
                   centCfg.cent9MaxRefMultCorrBin, centCfg.cent9MaxRefMultCorr);
    }
  }
  note += "phi_proton legacy channel has no hKstar*VsCent 2D; cent-slice CF pages skip it.\n";
  note += "Re-run analysis after hist/Maker changes so new keys exist in the ROOT file.\n";

  std::map<std::string, TGraphErrors*> cfCache;
  populateCfCache(fin, cfCache);
  populateCfCentCache(fin, cfCache);

  PdfHeader::MakePdfHeaderPage(pdfName, "checkHistAnaFemtoPhiProton.C", inputs, note.Data(), true, anaName);

  TCanvas* c1 = new TCanvas("c1", "canvas", 1200, 800);
  TH1* h1 = 0;
  TH2* h2 = 0;

  // Page 1: Event Level
  c1->Clear();
  c1->Divide(3, 3);
  c1->cd(1); h1 = (TH1*)fin->Get("hVz"); if (h1) { h1->Draw(); if (gConfigLoaded) { EventCutConfig& ev = ConfigManager::GetInstance().GetEventCuts(); drawCutLines1D(h1, ev.minVz, ev.maxVz); } }
  c1->cd(2); h1 = (TH1*)fin->Get("hVzDiff"); if (h1) { h1->Draw(); if (gConfigLoaded) { EventCutConfig& ev = ConfigManager::GetInstance().GetEventCuts(); drawCutLines1D(h1, -ev.maxVzDiff, ev.maxVzDiff); } }
  c1->cd(3); gPad->SetLogz(); h2 = (TH2*)fin->Get("hVxVy"); if (h2) {
    h2->Draw("colz");
    if (gConfigLoaded) {
      EventCutConfig& ev = ConfigManager::GetInstance().GetEventCuts();
      drawVtxCutCircle(ev.vtxCenterX, ev.vtxCenterY, ev.maxVr);
    }
  }
  c1->cd(4); h1 = (TH1*)fin->Get("hRefMult"); if (h1) { h1->Draw(); if (gConfigLoaded) { EventCutConfig& ev = ConfigManager::GetInstance().GetEventCuts(); drawCutLines1D(h1, ev.minRefMult, ev.maxRefMult); } }
  c1->cd(5); gPad->SetLogz(); h2 = (TH2*)fin->Get("hVzVsRun"); if (h2) h2->Draw("colz");
  c1->cd(6); gPad->SetLogz(); h2 = (TH2*)fin->Get("hRefMultVsVz"); if (h2) h2->Draw("colz");
  c1->cd(7); h1 = (TH1*)fin->Get("hNTracks"); if (h1) { h1->Draw(); if (gConfigLoaded) { PhiCutConfig& phi = ConfigManager::GetInstance().GetPhiCuts(); if (phi.maxNTr > 0) drawCutLine1D(h1, (Double_t)phi.maxNTr); } }
  c1->cd(8); h1 = (TH1*)fin->Get("hVr"); if (h1) { h1->Draw(); if (gConfigLoaded) { EventCutConfig& ev = ConfigManager::GetInstance().GetEventCuts(); drawCutLine1D(h1, ev.maxVr); } }
  c1->cd(9); /* spare */;
  c1->Print(pdfName);

  // Page 1b: Centrality QA (event-level)
  c1->Clear();
  c1->Divide(3, 3);
  c1->cd(1); h1 = (TH1*)fin->Get("hCentrality"); if (h1) h1->Draw();
  c1->cd(2); h1 = (TH1*)fin->Get("hCentralityRaw"); if (h1) h1->Draw();
  c1->cd(3); h1 = (TH1*)fin->Get("hCentrality16"); if (h1) h1->Draw();
  c1->cd(4); h1 = (TH1*)fin->Get("hRefMultCorr"); if (h1) h1->Draw();
  c1->cd(5); h1 = (TH1*)fin->Get("hRefMultWeight"); if (h1) h1->Draw();
  c1->cd(6); gPad->SetLogz(); h2 = (TH2*)fin->Get("hRefMultVsNTOFMatch"); if (h2) h2->Draw("colz");
  c1->cd(7); gPad->SetLogz(); h2 = (TH2*)fin->Get("hRefMultVsNTOFMatchAfter"); if (h2) h2->Draw("colz");
  c1->cd(8); gPad->SetLogz(); h2 = (TH2*)fin->Get("hCentralityVsVz"); if (h2) h2->Draw("colz");
  c1->cd(9); h1 = (TH1*)fin->Get("hRawMult"); if (h1) h1->Draw();
  c1->Print(pdfName);

  // Page 1c: Centrality correlations (event observables vs cent9)
  c1->Clear();
  c1->Divide(3, 3);
  c1->cd(1); gPad->SetLogz(); h2 = (TH2*)fin->Get("hRawMult_vs_Cent9"); if (h2) h2->Draw("colz");
  c1->cd(2); gPad->SetLogz(); h2 = (TH2*)fin->Get("hRefMultCorr_vs_Cent9"); if (h2) {
    h2->Draw("colz");
    if (gConfigLoaded) {
      const CentralityCutConfig& centCfg = ConfigManager::GetInstance().GetCentralityCuts();
      if (centCfg.cent9MaxRefMultCorrBin >= 0 && centCfg.cent9MaxRefMultCorr > 0.0) {
        drawCutLine2DH(h2, centCfg.cent9MaxRefMultCorr);
      }
    }
  }
  c1->cd(3); gPad->SetLogz(); h2 = (TH2*)fin->Get("hNTracks_vs_Cent9"); if (h2) h2->Draw("colz");
  c1->cd(4); gPad->SetLogz(); h2 = (TH2*)fin->Get("hTofMatchMult_vs_Cent9"); if (h2) h2->Draw("colz");
  c1->cd(5); gPad->SetLogz(); h2 = (TH2*)fin->Get("hNKaonPlus_vs_Cent9"); if (h2) h2->Draw("colz");
  c1->cd(6); gPad->SetLogz(); h2 = (TH2*)fin->Get("hNKaonMinus_vs_Cent9"); if (h2) h2->Draw("colz");
  c1->cd(7); gPad->SetLogz(); h2 = (TH2*)fin->Get("hNPhiPairs_vs_Cent9"); if (h2) h2->Draw("colz");
  c1->cd(8); gPad->SetLogz(); h2 = (TH2*)fin->Get("hNProton_vs_Cent9"); if (h2) h2->Draw("colz");
  c1->cd(9); gPad->SetLogz(); h2 = (TH2*)fin->Get("hRawMult_vs_RefMultCorr"); if (h2) h2->Draw("colz");
  drawCent9ConventionNote();
  c1->Print(pdfName);

  // Page 1e: Event Level (Post-Cut)
  c1->Clear();
  c1->Divide(3, 3);
  c1->cd(1); h1 = (TH1*)fin->Get("hVz_After"); if (h1) { h1->Draw(); if (gConfigLoaded) { EventCutConfig& ev = ConfigManager::GetInstance().GetEventCuts(); drawCutLines1D(h1, ev.minVz, ev.maxVz); } }
  c1->cd(2); h1 = (TH1*)fin->Get("hVzDiff_After"); if (h1) { h1->Draw(); if (gConfigLoaded) { EventCutConfig& ev = ConfigManager::GetInstance().GetEventCuts(); drawCutLines1D(h1, -ev.maxVzDiff, ev.maxVzDiff); } }
  c1->cd(3); gPad->SetLogz(); h2 = (TH2*)fin->Get("hVxVy_After"); if (h2) {
    h2->Draw("colz");
    if (gConfigLoaded) {
      EventCutConfig& ev = ConfigManager::GetInstance().GetEventCuts();
      drawVtxCutCircle(ev.vtxCenterX, ev.vtxCenterY, ev.maxVr);
    }
  }
  c1->cd(4); h1 = (TH1*)fin->Get("hRefMult_After"); if (h1) { h1->Draw(); if (gConfigLoaded) { EventCutConfig& ev = ConfigManager::GetInstance().GetEventCuts(); drawCutLines1D(h1, ev.minRefMult, ev.maxRefMult); } }
  c1->cd(5); /* spare */;
  c1->cd(6); /* spare */;
  c1->cd(7); /* spare */;
  c1->cd(8); h1 = (TH1*)fin->Get("hVr_After"); if (h1) { h1->Draw(); if (gConfigLoaded) { EventCutConfig& ev = ConfigManager::GetInstance().GetEventCuts(); drawCutLine1D(h1, ev.maxVr); } }
  c1->cd(9); /* spare */;
  c1->Print(pdfName);

  // Page 2: Track Kinematics & Quality (Post-Cut)
  c1->Clear();
  c1->Divide(3, 3);
  c1->cd(1); gPad->SetLogy(); h1 = (TH1*)fin->Get("hPt"); if (h1) h1->Draw();
  c1->cd(2); gPad->SetLogy(0); h1 = (TH1*)fin->Get("hEta"); if (h1) h1->Draw();
  c1->cd(3); h1 = (TH1*)fin->Get("hPhi"); if (h1) h1->Draw();
  c1->cd(4); h1 = (TH1*)fin->Get("hCharge"); if (h1) h1->Draw("hist");
  c1->cd(5); h1 = (TH1*)fin->Get("hNHitsFit"); if (h1) h1->Draw();
  c1->cd(6); h1 = (TH1*)fin->Get("hNHitsRatio"); if (h1) h1->Draw();
  c1->cd(7); h1 = (TH1*)fin->Get("hNHitsDedx"); if (h1) h1->Draw();
  c1->cd(8); gPad->SetLogy(); h1 = (TH1*)fin->Get("hDCA"); if (h1) h1->Draw();
  c1->cd(9); gPad->SetLogy(); h1 = (TH1*)fin->Get("hChi2"); if (h1) h1->Draw();
  c1->Print(pdfName);

  // Page 2b: Track Kinematics & Quality (Pre-Cut)
  c1->Clear();
  c1->Divide(3, 3);
  if (gConfigLoaded) {
    TrackCutConfig& tr = ConfigManager::GetInstance().GetTrackCuts();
    c1->cd(1); gPad->SetLogy(); h1 = (TH1*)fin->Get("hPt_Raw"); if (h1) { h1->Draw(); drawCutLines1D(h1, tr.minPt, tr.maxPt); }
    c1->cd(2); gPad->SetLogy(0); h1 = (TH1*)fin->Get("hEta_Raw"); if (h1) { h1->Draw(); drawCutLines1D(h1, tr.minEta, tr.maxEta); }
    c1->cd(3); h1 = (TH1*)fin->Get("hPhi_Raw"); if (h1) h1->Draw();
    c1->cd(4); h1 = (TH1*)fin->Get("hCharge"); if (h1) h1->Draw("hist");
    c1->cd(5); h1 = (TH1*)fin->Get("hNHitsFit_Raw"); if (h1) { h1->Draw(); drawCutLine1D(h1, (Double_t)tr.minNHitsFit); }
    c1->cd(6); h1 = (TH1*)fin->Get("hNHitsRatio_Raw"); if (h1) { h1->Draw(); drawCutLine1D(h1, tr.minNHitsRatio); }
    c1->cd(7); h1 = (TH1*)fin->Get("hNHitsDedx_Raw"); if (h1) { h1->Draw(); drawCutLine1D(h1, (Double_t)tr.minNHitsDedx); }
    c1->cd(8); gPad->SetLogy(); h1 = (TH1*)fin->Get("hDCA_Raw"); if (h1) { h1->Draw(); drawCutLine1D(h1, tr.maxDCA); }
    c1->cd(9); gPad->SetLogy(); h1 = (TH1*)fin->Get("hChi2_Raw"); if (h1) { h1->Draw(); drawCutLine1D(h1, tr.maxChi2); }
  } else {
    c1->cd(1); gPad->SetLogy(); h1 = (TH1*)fin->Get("hPt_Raw"); if (h1) h1->Draw();
    c1->cd(2); gPad->SetLogy(0); h1 = (TH1*)fin->Get("hEta_Raw"); if (h1) h1->Draw();
    c1->cd(3); h1 = (TH1*)fin->Get("hPhi_Raw"); if (h1) h1->Draw();
    c1->cd(4); h1 = (TH1*)fin->Get("hCharge"); if (h1) h1->Draw("hist");
    c1->cd(5); h1 = (TH1*)fin->Get("hNHitsFit_Raw"); if (h1) h1->Draw();
    c1->cd(6); h1 = (TH1*)fin->Get("hNHitsRatio_Raw"); if (h1) h1->Draw();
    c1->cd(7); h1 = (TH1*)fin->Get("hNHitsDedx_Raw"); if (h1) h1->Draw();
    c1->cd(8); gPad->SetLogy(); h1 = (TH1*)fin->Get("hDCA_Raw"); if (h1) h1->Draw();
    c1->cd(9); gPad->SetLogy(); h1 = (TH1*)fin->Get("hChi2_Raw"); if (h1) h1->Draw();
  }
  c1->Print(pdfName);

  // Page 3: PID (TPC & TOF)
  c1->Clear();
  c1->Divide(3, 2);
  double pmax = 5.0;
  double dedxmax = 1e-6;
  c1->cd(1); gPad->SetLogz(); h2 = (TH2*)fin->Get("hDedxVsP"); if (h2) { h2->GetXaxis()->SetRangeUser(0, pmax); h2->GetYaxis()->SetRangeUser(0, dedxmax); h2->Draw("colz"); }
  c1->cd(2); gPad->SetLogz(); h2 = (TH2*)fin->Get("hBetaVsP"); if (h2) { h2->GetXaxis()->SetRangeUser(0, pmax); h2->Draw("colz"); }
  c1->cd(3); gPad->SetLogz(); h2 = (TH2*)fin->Get("hMass2VsP"); if (h2) {
    h2->GetXaxis()->SetRangeUser(0, pmax);
    h2->Draw("colz");
    if (gConfigLoaded) {
      PIDCutConfig& pid = ConfigManager::GetInstance().GetPIDCuts();
      drawCutLine2DH(h2, pid.minMass2Kaon);
      drawCutLine2DH(h2, pid.maxMass2Kaon);
    }
  }
  c1->cd(4); gPad->SetLogz(); h2 = (TH2*)fin->Get("hNSigmaPionVsP"); if (h2) { h2->GetXaxis()->SetRangeUser(0, pmax); h2->Draw("colz"); }
  c1->cd(5); gPad->SetLogz(); h2 = (TH2*)fin->Get("hNSigmaKaonVsP"); if (h2) { h2->GetXaxis()->SetRangeUser(0, pmax); h2->Draw("colz"); if (gConfigLoaded) { PhiCutConfig& phi = ConfigManager::GetInstance().GetPhiCuts(); drawCutLine2DH(h2, phi.nSigmaKaon); drawCutLine2DH(h2, -phi.nSigmaKaon); } }
  c1->cd(6); gPad->SetLogz(); h2 = (TH2*)fin->Get("hNSigmaProtonVsP"); if (h2) { h2->GetXaxis()->SetRangeUser(0, pmax); h2->Draw("colz"); }
  c1->Print(pdfName);

  // Page 3d: PID vs pT (mirror of Page 3)
  c1->Clear();
  c1->Divide(3, 2);
  double ptmax = 5.0;
  c1->cd(1); gPad->SetLogz(); h2 = (TH2*)fin->Get("hDedxVsPt"); if (h2) { h2->GetXaxis()->SetRangeUser(0, ptmax); h2->GetYaxis()->SetRangeUser(0, dedxmax); h2->Draw("colz"); }
  c1->cd(2); gPad->SetLogz(); h2 = (TH2*)fin->Get("hBetaVsPt"); if (h2) { h2->GetXaxis()->SetRangeUser(0, ptmax); h2->Draw("colz"); }
  c1->cd(3); gPad->SetLogz(); h2 = (TH2*)fin->Get("hMass2VsPt"); if (h2) {
    h2->GetXaxis()->SetRangeUser(0, ptmax);
    h2->Draw("colz");
    if (gConfigLoaded) {
      PIDCutConfig& pid = ConfigManager::GetInstance().GetPIDCuts();
      drawCutLine2DH(h2, pid.minMass2Kaon);
      drawCutLine2DH(h2, pid.maxMass2Kaon);
    }
  }
  c1->cd(4); gPad->SetLogz(); h2 = (TH2*)fin->Get("hNSigmaPionVsPt"); if (h2) { h2->GetXaxis()->SetRangeUser(0, ptmax); h2->Draw("colz"); }
  c1->cd(5); gPad->SetLogz(); h2 = (TH2*)fin->Get("hNSigmaKaonVsPt"); if (h2) { h2->GetXaxis()->SetRangeUser(0, ptmax); h2->Draw("colz"); if (gConfigLoaded) { PhiCutConfig& phi = ConfigManager::GetInstance().GetPhiCuts(); drawCutLine2DH(h2, phi.nSigmaKaon); drawCutLine2DH(h2, -phi.nSigmaKaon); } }
  c1->cd(6); gPad->SetLogz(); h2 = (TH2*)fin->Get("hNSigmaProtonVsPt"); if (h2) { h2->GetXaxis()->SetRangeUser(0, ptmax); h2->Draw("colz"); }
  c1->Print(pdfName);

  // Page 3bp: TOF m2 vs pT (TPC K vs final K)
  c1->Clear();
  c1->Divide(2, 1);
  c1->cd(1); gPad->SetLogz(); h2 = (TH2*)fin->Get("hMass2VsPt_TpcKaon"); if (h2) {
    h2->GetXaxis()->SetRangeUser(0, ptmax);
    h2->Draw("colz");
    if (gConfigLoaded) {
      PIDCutConfig& pid = ConfigManager::GetInstance().GetPIDCuts();
      drawCutLine2DH(h2, pid.minMass2Kaon);
      drawCutLine2DH(h2, pid.maxMass2Kaon);
    }
  }
  c1->cd(2); gPad->SetLogz(); h2 = (TH2*)fin->Get("hMass2VsPt_IsKaon"); if (h2) {
    h2->GetXaxis()->SetRangeUser(0, ptmax);
    h2->Draw("colz");
    if (gConfigLoaded) {
      PIDCutConfig& pid = ConfigManager::GetInstance().GetPIDCuts();
      drawCutLine2DH(h2, pid.minMass2Kaon);
      drawCutLine2DH(h2, pid.maxMass2Kaon);
    }
  }
  c1->Print(pdfName);

  // Page 3b: TOF m2 (TPC K vs final K) + pair DCA_KK
  c1->Clear();
  c1->Divide(2, 2);
  c1->cd(1); gPad->SetLogz(); h2 = (TH2*)fin->Get("hMass2VsP_TpcKaon"); if (h2) {
    h2->GetXaxis()->SetRangeUser(0, pmax);
    h2->Draw("colz");
    if (gConfigLoaded) {
      PIDCutConfig& pid = ConfigManager::GetInstance().GetPIDCuts();
      drawCutLine2DH(h2, pid.minMass2Kaon);
      drawCutLine2DH(h2, pid.maxMass2Kaon);
    }
  }
  c1->cd(2); gPad->SetLogz(); h2 = (TH2*)fin->Get("hMass2VsP_IsKaon"); if (h2) {
    h2->GetXaxis()->SetRangeUser(0, pmax);
    h2->Draw("colz");
    if (gConfigLoaded) {
      PIDCutConfig& pid = ConfigManager::GetInstance().GetPIDCuts();
      drawCutLine2DH(h2, pid.minMass2Kaon);
      drawCutLine2DH(h2, pid.maxMass2Kaon);
    }
  }
  c1->cd(3); h1 = (TH1*)fin->Get("hDCAKK_All"); if (h1) {
    h1->Draw();
    if (gConfigLoaded) {
      PhiCutConfig& phi = ConfigManager::GetInstance().GetPhiCuts();
      drawCutLine1D(h1, phi.maxDCAKK);
    }
  }
  c1->cd(4); h1 = (TH1*)fin->Get("hDCAKK_Pass"); if (h1) {
    h1->Draw();
    if (gConfigLoaded) {
      PhiCutConfig& phi = ConfigManager::GetInstance().GetPhiCuts();
      drawCutLine1D(h1, phi.maxDCAKK);
    }
  }
  c1->Print(pdfName);

  // Page 3c: TOF diagnostics (m2/q2 vs p/q, delta(1/beta) vs p and 1D, beta vs p)
  c1->Clear();
  c1->Divide(2, 2);
  c1->cd(1); gPad->SetLogz(); h2 = (TH2*)fin->Get("hM2q2VsPq"); if (h2) h2->Draw("colz");
  c1->cd(2); gPad->SetLogz(); h2 = (TH2*)fin->Get("hDeltaOneOverBetaVsP"); if (h2) {
    h2->GetXaxis()->SetRangeUser(0, pmax);
    h2->Draw("colz");
    if (gConfigLoaded) {
      PIDCutConfig& pid = ConfigManager::GetInstance().GetPIDCuts();
      drawCutLine2DH(h2, pid.maxAbsDeltaOneOverBetaKaon);
      drawCutLine2DH(h2, -pid.maxAbsDeltaOneOverBetaKaon);
    }
  }
  c1->cd(3); h1 = (TH1*)fin->Get("hDeltaOneOverBetaKaon"); if (h1) {
    h1->Draw();
    if (gConfigLoaded) {
      PIDCutConfig& pid = ConfigManager::GetInstance().GetPIDCuts();
      drawCutLines1D(h1, -pid.maxAbsDeltaOneOverBetaKaon, pid.maxAbsDeltaOneOverBetaKaon);
    }
  }
  c1->cd(4); gPad->SetLogz(); h2 = (TH2*)fin->Get("hBetaVsP"); if (h2) {
    h2->GetXaxis()->SetRangeUser(0, pmax);
    h2->Draw("colz");
  }
  c1->Print(pdfName);

  // Page 4: Event Plane & Misc
  c1->Clear();
  c1->Divide(2, 2);
  c1->cd(1); gPad->SetLogy(); h1 = (TH1*)fin->Get("hTriggerIds"); if (h1) { h1->SetFillColor(17); h1->Draw(); }
  c1->cd(2); gPad->SetLogy(0); h1 = (TH1*)fin->Get("hTofMatchMult"); if (h1) h1->Draw();
  c1->cd(3); gPad->SetLogz(); h2 = (TH2*)fin->Get("hQxQy"); if (h2) h2->Draw("colz");
  c1->cd(4); h1 = (TH1*)fin->Get("hPsi2"); if (h1) { h1->SetMinimum(0); h1->Draw(); }
  c1->Print(pdfName);

  // Page 8b: staged pair QA mass and overlay
  c1->Clear();
  c1->Divide(2, 2);
  c1->cd(1); h1 = (TH1*)fin->Get("hPhiPair_Mass_stage0"); if (h1) h1->Draw();
  c1->cd(2); h1 = (TH1*)fin->Get("hPhiPair_Mass_tofStrict"); if (h1) h1->Draw();
  c1->cd(3); h1 = (TH1*)fin->Get("hPhiPair_NPairs_stage0"); if (h1) h1->Draw();
  c1->cd(4); h1 = (TH1*)fin->Get("hPhiPair_NPairs_tofStrict"); if (h1) h1->Draw();
  c1->Print(pdfName);

  c1->Clear();
  c1->Divide(1, 1);
  c1->cd(1);
  TH1* hStage0 = (TH1*)fin->Get("hPhiPair_Mass_stage0");
  TH1* hTofStrict = (TH1*)fin->Get("hPhiPair_Mass_tofStrict");
  if (hStage0) {
    hStage0->SetLineColor(kBlue + 1);
    hStage0->SetTitle("Pair QA stage0 / strict TOF overlay;M_{KK} [GeV/c^{2}];Counts");
    hStage0->Draw("hist");
    if (hTofStrict) {
      hTofStrict->SetLineColor(kRed + 1);
      hTofStrict->Draw("hist same");
    }
    TLatex* tag = new TLatex();
    tag->SetNDC(kTRUE);
    tag->SetTextSize(0.035);
    tag->DrawLatex(0.15, 0.86, "Blue: stage0");
    if (hTofStrict) tag->DrawLatex(0.15, 0.80, "Red: strict TOF pair");
  }
  c1->Print(pdfName);

  // Page 8c: strict TOF pair QA mass projections
  c1->Clear();
  c1->Divide(2, 2);
  c1->cd(1); gPad->SetLogz(); h2 = (TH2*)fin->Get("hPhiPair_MassVsPt_tofStrict"); if (h2) h2->Draw("colz");
  c1->cd(2); gPad->SetLogz(); h2 = (TH2*)fin->Get("hPhiPair_MassVsCent_tofStrict"); if (h2) { h2->Draw("colz"); drawCent9ConventionNote(); }
  c1->cd(3);
  h2 = (TH2*)fin->Get("hPhiPair_MassVsPt_tofStrict");
  if (h2) {
    Int_t b1 = h2->GetXaxis()->FindBin(0.2 + 1e-6);
    Int_t b2 = h2->GetXaxis()->FindBin(1.6 - 1e-6);
    TH1D* hProjPt = h2->ProjectionY("hPhiPair_Mass_tofStrict_pt0p2to1p6", b1, b2);
    if (hProjPt) {
      hProjPt->SetTitle("Strict TOF M_{KK} (0.2<p_{T}<1.6 GeV/c);M_{KK} [GeV/c^{2}];Counts");
      hProjPt->Draw();
    }
  }
  c1->cd(4);
  h2 = (TH2*)fin->Get("hPhiPair_MassVsCent_tofStrict");
  if (h2) {
    Int_t b1 = h2->GetXaxis()->FindBin(6.0 + 1e-6);
    Int_t b2 = h2->GetXaxis()->FindBin(8.0 - 1e-6);
    TH1D* hProjCent = h2->ProjectionY("hPhiPair_Mass_tofStrict_cent6to8", b1, b2);
    if (hProjCent) {
      hProjCent->SetTitle("Strict TOF M_{KK} (cent9: 6-8);M_{KK} [GeV/c^{2}];Counts");
      hProjCent->Draw();
    }
  }
  c1->Print(pdfName);

  // Page 8d: phi pair kinematics (stage0 vs strict TOF)
  c1->Clear();
  c1->Divide(3, 2);
  c1->cd(1); h1 = (TH1*)fin->Get("hPhiPair_Pt_stage0"); if (h1) h1->Draw();
  c1->cd(2); h1 = (TH1*)fin->Get("hPhiPair_Pt_tofStrict"); if (h1) h1->Draw();
  c1->cd(3); h1 = (TH1*)fin->Get("hPhiPair_Rapidity_stage0"); if (h1) {
    h1->Draw();
    if (gConfigLoaded) {
      PhiCutConfig& phi = ConfigManager::GetInstance().GetPhiCuts();
      drawCutLines1D(h1, phi.minPairRapidity, phi.maxPairRapidity);
    }
  }
  c1->cd(4); h1 = (TH1*)fin->Get("hPhiPair_Rapidity_tofStrict"); if (h1) {
    h1->Draw();
    if (gConfigLoaded) {
      PhiCutConfig& phi = ConfigManager::GetInstance().GetPhiCuts();
      drawCutLines1D(h1, phi.minPairRapidity, phi.maxPairRapidity);
    }
  }
  c1->cd(5); h1 = (TH1*)fin->Get("hPhiPair_OpeningAngle_stage0"); if (h1) {
    h1->Draw();
    if (gConfigLoaded) {
      PhiCutConfig& phi = ConfigManager::GetInstance().GetPhiCuts();
      drawCutLines1D(h1, phi.minOpeningAngle, phi.maxOpeningAngle);
    }
  }
  c1->cd(6); h1 = (TH1*)fin->Get("hPhiPair_OpeningAngle_tofStrict"); if (h1) {
    h1->Draw();
    if (gConfigLoaded) {
      PhiCutConfig& phi = ConfigManager::GetInstance().GetPhiCuts();
      drawCutLines1D(h1, phi.minOpeningAngle, phi.maxOpeningAngle);
    }
  }
  c1->Print(pdfName);

  // Page 9: Kaon QA (before vs after PID cut)
  c1->Clear();
  c1->Divide(3, 2);
  // Row 1: Before PID cut (all accepted tracks)
  c1->cd(1); h1 = (TH1*)fin->Get("hPt"); if (h1) {
    h1->Draw();
    if (gConfigLoaded) {
      TrackCutConfig& tr = ConfigManager::GetInstance().GetTrackCuts();
      drawCutLines1D(h1, tr.minPt, tr.maxPt);
    }
  }
  c1->cd(2); h1 = (TH1*)fin->Get("hEta"); if (h1) {
    h1->Draw();
    if (gConfigLoaded) {
      TrackCutConfig& tr = ConfigManager::GetInstance().GetTrackCuts();
      drawCutLines1D(h1, tr.minEta, tr.maxEta);
    }
  }
  c1->cd(3); h1 = (TH1*)fin->Get("hNSigmaKaon_Raw"); if (h1) {
    h1->Draw();
    if (gConfigLoaded) {
      PhiCutConfig& phi = ConfigManager::GetInstance().GetPhiCuts();
      drawCutLines1D(h1, -phi.nSigmaKaon, phi.nSigmaKaon);
    }
  }
  // Row 2: After PID cut (final Kaons)
  c1->cd(4); h1 = (TH1*)fin->Get("hK_Pt"); if (h1) {
    h1->Draw();
    if (gConfigLoaded) {
      TrackCutConfig& tr = ConfigManager::GetInstance().GetTrackCuts();
      drawCutLines1D(h1, tr.minPt, tr.maxPt);
    }
  }
  c1->cd(5); h1 = (TH1*)fin->Get("hK_Eta"); if (h1) {
    h1->Draw();
    if (gConfigLoaded) {
      TrackCutConfig& tr = ConfigManager::GetInstance().GetTrackCuts();
      drawCutLines1D(h1, tr.minEta, tr.maxEta);
    }
  }
  c1->cd(6); h1 = (TH1*)fin->Get("hK_NSigma"); if (h1) {
    h1->Draw();
    if (gConfigLoaded) {
      PhiCutConfig& phi = ConfigManager::GetInstance().GetPhiCuts();
      drawCutLines1D(h1, -phi.nSigmaKaon, phi.nSigmaKaon);
    }
  }
  c1->Print(pdfName);

  // Page 10: K+ and K- multiplicity (per event) and correlation
  c1->Clear();
  c1->Divide(2, 2);
  c1->cd(1); gPad->SetLogz(); h2 = (TH2*)fin->Get("hNKaonPlusVsNKaonMinus"); if (h2) h2->Draw("colz");
  c1->cd(2); h1 = (TH1*)fin->Get("hNKaonPlus"); if (h1) h1->Draw();
  c1->cd(3); h1 = (TH1*)fin->Get("hNKaonMinus"); if (h1) h1->Draw();
  c1->cd(4); gPad->SetLogz(); h2 = (TH2*)fin->Get("hNSigmaProtonVsPt"); if (h2) h2->Draw("colz");
  c1->Print(pdfName);

  // Page 11: Proton candidate QA
  c1->Clear();
  c1->Divide(3, 2);
  c1->cd(1); gPad->SetLogy(); h1 = (TH1*)fin->Get("hP_Pt"); if (h1) h1->Draw();
  c1->cd(2); h1 = (TH1*)fin->Get("hP_Eta"); if (h1) h1->Draw();
  c1->cd(3); h1 = (TH1*)fin->Get("hP_Phi"); if (h1) h1->Draw();
  c1->cd(4); h1 = (TH1*)fin->Get("hP_NSigmaProton"); if (h1) h1->Draw();
  c1->cd(5); gPad->SetLogy(); h1 = (TH1*)fin->Get("hP_Mass2"); if (h1) h1->Draw();
  c1->cd(6); gPad->SetLogy(); h1 = (TH1*)fin->Get("hP_DCA"); if (h1) h1->Draw();
  c1->Print(pdfName);

  // Page 12: Phi candidate QA (final femto resonance list)
  c1->Clear();
  c1->Divide(2, 2);
  c1->cd(1); h1 = (TH1*)fin->Get("hPhi_MKK"); if (h1) h1->Draw();
  c1->cd(2); h1 = (TH1*)fin->Get("hPhi_Pt"); if (h1) h1->Draw();
  c1->cd(3); h1 = (TH1*)fin->Get("hPhi_Rapidity"); if (h1) h1->Draw();
  c1->cd(4); h1 = (TH1*)fin->Get("hPhi_NCand"); if (h1) h1->Draw();
  c1->Print(pdfName);

  // Page 13: Femto k* (legacy phi-proton channel)
  c1->Clear();
  c1->Divide(2, 2);
  c1->cd(1); h1 = (TH1*)fin->Get("hKstarSE_phi_proton"); if (h1) h1->Draw();
  c1->cd(2); h1 = (TH1*)fin->Get("hKstarME_phi_proton"); if (h1) h1->Draw();
  c1->cd(3);
  drawComputedCf(c1, 3, fin, "phi_proton", channelNormQMin("phi_proton"), channelNormQMax("phi_proton"), cfCache);
  c1->cd(4); h1 = (TH1*)fin->Get("hP_NCand"); if (h1) h1->Draw();
  c1->Print(pdfName);

  // Page 14: Proton femto PID (Zhangwei-like bachelor cuts)
  c1->Clear();
  c1->Divide(3, 2);
  c1->cd(1); h1 = (TH1*)fin->Get("hP_Y_PreFemtoCut"); if (h1) h1->Draw();
  c1->cd(2); gPad->SetLogz(); h2 = (TH2*)fin->Get("hP_PtVsY_PreFemtoCut"); if (h2) h2->Draw("colz");
  c1->cd(3); h1 = (TH1*)fin->Get("hP_Y_FemtoCut"); if (h1) h1->Draw();
  c1->cd(4); gPad->SetLogz(); h2 = (TH2*)fin->Get("hP_PtVsY_FemtoCut"); if (h2) h2->Draw("colz");
  c1->cd(5); gPad->SetLogz(); h2 = (TH2*)fin->Get("hP_Mass2VsP"); if (h2) h2->Draw("colz");
  c1->cd(6); h1 = (TH1*)fin->Get("hP_NHitsRatio_FemtoCut"); if (h1) h1->Draw();
  c1->Print(pdfName);

  // Page 15: Phi mass windows / sidebands / rotation
  c1->Clear();
  c1->Divide(3, 2);
  c1->cd(1); h1 = (TH1*)fin->Get("hPhi_MKK_signal"); if (h1) h1->Draw();
  c1->cd(2); h1 = (TH1*)fin->Get("hPhi_MKK_leftSB"); if (h1) h1->Draw();
  c1->cd(3); h1 = (TH1*)fin->Get("hPhi_MKK_rightSB"); if (h1) h1->Draw();
  c1->cd(4); h1 = (TH1*)fin->Get("hPhi_MKK_rot"); if (h1) h1->Draw();
  c1->cd(5); h1 = (TH1*)fin->Get("hPhiRot_MKK"); if (h1) h1->Draw();
  c1->cd(6); h1 = (TH1*)fin->Get("hPhiRot_NCand"); if (h1) h1->Draw();
  c1->Print(pdfName);

  // Page 16: k* SE/ME/CF — signal channel
  c1->Clear();
  c1->Divide(2, 2);
  c1->cd(1); h1 = (TH1*)fin->Get("hKstarSE_phi_proton_signal"); if (h1) h1->Draw();
  c1->cd(2); h1 = (TH1*)fin->Get("hKstarME_phi_proton_signal"); if (h1) h1->Draw();
  c1->cd(3);
  drawComputedCf(c1, 3, fin, "phi_proton_signal", channelNormQMin("phi_proton_signal"),
                 channelNormQMax("phi_proton_signal"), cfCache);
  c1->cd(4); gPad->SetLogz(); h2 = (TH2*)fin->Get("hKstarSEVsCent_phi_proton_signal"); if (h2) h2->Draw("colz");
  c1->Print(pdfName);

  // Page 17: k* SE/ME/CF — sideband channels
  c1->Clear();
  c1->Divide(2, 3);
  c1->cd(1); h1 = (TH1*)fin->Get("hKstarSE_phi_proton_leftSB"); if (h1) h1->Draw();
  c1->cd(2); h1 = (TH1*)fin->Get("hKstarME_phi_proton_leftSB"); if (h1) h1->Draw();
  c1->cd(3);
  drawComputedCf(c1, 3, fin, "phi_proton_leftSB", channelNormQMin("phi_proton_leftSB"),
                 channelNormQMax("phi_proton_leftSB"), cfCache);
  c1->cd(4); h1 = (TH1*)fin->Get("hKstarSE_phi_proton_rightSB"); if (h1) h1->Draw();
  c1->cd(5); h1 = (TH1*)fin->Get("hKstarME_phi_proton_rightSB"); if (h1) h1->Draw();
  c1->cd(6);
  drawComputedCf(c1, 6, fin, "phi_proton_rightSB", channelNormQMin("phi_proton_rightSB"),
                 channelNormQMax("phi_proton_rightSB"), cfCache);
  c1->Print(pdfName);

  // Page 18: k* SE/ME/CF — rotation background channel
  c1->Clear();
  c1->Divide(2, 2);
  c1->cd(1); h1 = (TH1*)fin->Get("hKstarSE_phi_rot_proton"); if (h1) h1->Draw();
  c1->cd(2); h1 = (TH1*)fin->Get("hKstarME_phi_rot_proton"); if (h1) h1->Draw();
  c1->cd(3);
  drawComputedCf(c1, 3, fin, "phi_rot_proton", channelNormQMin("phi_rot_proton"),
                 channelNormQMax("phi_rot_proton"), cfCache);
  c1->cd(4); gPad->SetLogz(); h2 = (TH2*)fin->Get("hKstarSEVsCent_phi_rot_proton"); if (h2) h2->Draw("colz");
  c1->Print(pdfName);

  // Page 19-20: k* SE/ME/CF — cent9 slice (default 0-60% = cent9 2-8), projected from hKstar*VsCent
  Int_t cfCent9Min = 0;
  Int_t cfCent9Max = 0;
  getCfCent9Range(cfCent9Min, cfCent9Max);
  c1->Clear();
  c1->Divide(2, 3);
  drawCentProjectedSeMe(c1, 1, fin, "phi_proton_signal", kTRUE, cfCent9Min, cfCent9Max);
  drawCentProjectedSeMe(c1, 2, fin, "phi_proton_signal", kFALSE, cfCent9Min, cfCent9Max);
  drawCentSliceCf(c1, 3, "phi_proton_signal", cfCache);
  drawCentProjectedSeMe(c1, 4, fin, "phi_proton_leftSB", kTRUE, cfCent9Min, cfCent9Max);
  drawCentProjectedSeMe(c1, 5, fin, "phi_proton_leftSB", kFALSE, cfCent9Min, cfCent9Max);
  drawCentSliceCf(c1, 6, "phi_proton_leftSB", cfCache);
  c1->Print(pdfName);

  c1->Clear();
  c1->Divide(2, 3);
  drawCentProjectedSeMe(c1, 1, fin, "phi_proton_rightSB", kTRUE, cfCent9Min, cfCent9Max);
  drawCentProjectedSeMe(c1, 2, fin, "phi_proton_rightSB", kFALSE, cfCent9Min, cfCent9Max);
  drawCentSliceCf(c1, 3, "phi_proton_rightSB", cfCache);
  drawCentProjectedSeMe(c1, 4, fin, "phi_rot_proton", kTRUE, cfCent9Min, cfCent9Max);
  drawCentProjectedSeMe(c1, 5, fin, "phi_rot_proton", kFALSE, cfCent9Min, cfCent9Max);
  drawCentSliceCf(c1, 6, "phi_rot_proton", cfCache);
  c1->Print(pdfName);

  // Console: verify key histograms exist and have entries
  {
    std::cout << "\n=== checkHistAnaFemtoPhiProton: key histogram entries ===\n";
    const char* keys[] = {"hVz",
                           "hVz_After",
                           "hCentrality",
                           "hNSigmaKaon_Raw",
                           "hK_Pt",
                           "hPhiPair_Mass_stage0",
                           "hPhiPair_Mass_tofStrict",
                           "hPhi_MKK",
                           "hPhi_MKK_signal",
                           "hPhi_MKK_leftSB",
                           "hPhi_MKK_rightSB",
                           "hPhi_MKK_rot",
                           "hPhiRot_MKK",
                           "hPhi_NCand",
                           "hP_Pt",
                           "hP_Y_FemtoCut",
                           "hP_NCand",
                           "hKstarSE_phi_proton",
                           "hKstarME_phi_proton",
                           "hKstarSE_phi_proton_signal",
                           "hKstarME_phi_proton_signal",
                           "hKstarSE_phi_proton_leftSB",
                           "hKstarME_phi_proton_leftSB",
                           "hKstarSE_phi_proton_rightSB",
                           "hKstarME_phi_proton_rightSB",
                           "hKstarSE_phi_rot_proton",
                           "hKstarME_phi_rot_proton",
                           "hNProton_vs_Cent9",
                           0};
    for (Int_t i = 0; keys[i]; ++i) {
      TObject* o = fin->Get(keys[i]);
      if (!o) {
        std::cout << "  " << keys[i] << ": NOT IN FILE\n";
        continue;
      }
      if (o->InheritsFrom("TH1")) {
        TH1* hh = (TH1*)o;
        Double_t n = hh->GetEntries();
        std::cout << "  " << keys[i] << ": entries=" << n;
        if (n < 1.0)
          std::cout << "  [empty — use ROOT from a run after new hist fills]";
        std::cout << "\n";
      } else {
        std::cout << "  " << keys[i] << ": unexpected class " << o->ClassName() << "\n";
      }
    }
    const char* cfKeys[] = {"hCF_phi_proton",         "hCF_phi_proton_signal", "hCF_phi_proton_leftSB",
                            "hCF_phi_proton_rightSB", "hCF_phi_rot_proton",    0};
    const char* cfChannels[] = {"phi_proton", "phi_proton_signal", "phi_proton_leftSB", "phi_proton_rightSB",
                                "phi_rot_proton", 0};
    for (Int_t i = 0; cfKeys[i]; ++i) {
      Int_t nPts = getCachedCfPointCount(cfCache, cfChannels[i]);
      std::cout << "  " << cfKeys[i] << " (computed): nPoints=" << nPts;
      if (nPts < 1) std::cout << "  [empty — check SE/ME or norm region]";
      std::cout << "\n";
    }
    const char* cfCentChannels[] = {"phi_proton_signal", "phi_proton_leftSB", "phi_proton_rightSB", "phi_rot_proton",
                                     0};
    std::cout << Form("  CF cent slice (cent9 %d-%d, projected):\n", cfCent9Min, cfCent9Max);
    for (Int_t i = 0; cfCentChannels[i]; ++i) {
      Int_t nPts = getCachedCfPointCount(cfCache, cfCentCacheKey(cfCentChannels[i]).c_str());
      std::cout << "    " << cfCentChannels[i] << ": nPoints=" << nPts;
      if (nPts < 1) std::cout << "  [empty — check hKstar*VsCent in ROOT]";
      std::cout << "\n";
    }
    std::cout << "=============================================================\n\n";
  }

  PdfHeader::ClosePdf(pdfName);

  freeCfCache(cfCache);
  delete c1;
  fin->Close();

  std::cout << "Done. PDF saved to: " << pdfName.Data() << std::endl;
}
