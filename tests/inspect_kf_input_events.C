// Read-only ROOT5 diagnostic of stored PicoEvent leaves and current event cuts.
// Does not reconstruct particles, fit a vertex, change cuts, or write ROOT data.
// Counts event-cut failures before centrality; compare the final surviving count
// with hKfStages to distinguish event acceptance from later centrality rejection.
#include "TEnv.h"
#include "TFile.h"
#include "TLeaf.h"
#include "TMath.h"
#include "TROOT.h"
#include "TString.h"
#include "TSystem.h"
#include "TTree.h"
#include <cfloat>
#include <fstream>
#include <iostream>

TLeaf* KfInputEventLeaf(TTree* tree, const char* name)
{
  TLeaf* leaf = tree->GetLeaf(TString("Event.") + name);
  if (!leaf) leaf = tree->GetLeaf(name);
  return leaf;
}

void inspect_kf_input_events(
    const char* filename, Long64_t maximumEvents = 1000,
    const char* eventConfig = "config/cuts/event/event_auau13p5_anaLambda_KFParticle.yaml")
{
  // TEnv is NOT a YAML parser and silently ignores nesting. Refuse mode-aware
  // input before reading any events instead of reporting a misleading (0,0)
  // center. The production Maker's stored metadata/event QA are authoritative.
  if (!eventConfig || !*eventConfig) {
    std::cerr << "ERROR: an event YAML path is required" << std::endl;
    gSystem->Exit(1); return;
  }
  std::ifstream configStream(eventConfig);
  if (!configStream.good()) {
    std::cerr << "ERROR: cannot read event cuts " << eventConfig << std::endl;
    gSystem->Exit(1); return;
  }
  TString configLine;
  while (configLine.ReadLine(configStream, kFALSE)) {
    const Ssiz_t comment = configLine.First('#');
    if (comment != kNPOS) configLine.Remove(comment);
    if (configLine.Contains("vertexByMode")) {
      std::cerr << "ERROR: event YAML contains vertexByMode; this legacy flat-TEnv "
                   "diagnostic cannot resolve analysis.mode. Use the mode-aware Maker's "
                   "KFEventSelectionConfiguration and hKfEventSelection/vertex QA via "
                   "tests/check_kfparticle_output.C instead." << std::endl;
      configStream.close();
      gSystem->Exit(1); return;
    }
  }
  configStream.close();
  TEnv cuts;
  // This diagnostic accepts the current flat numeric key: value YAML only.
  // TEnv's colon-delimited reader handles these numeric entries and # comments.
  if (cuts.ReadFile(eventConfig, kEnvLocal) != 0) {
    std::cerr << "ERROR: cannot read event cuts " << eventConfig << std::endl;
    gSystem->Exit(1); return;
  }
  const char* cutNames[8] = {"minVz", "maxVz", "maxVr", "minRefMult",
      "maxRefMult", "maxVzDiff", "maxAbsVzVpd", "maxNTr"};
  Int_t j = 0;
  for (j = 0; j < 8; ++j) {
    if (!cuts.Defined(cutNames[j])) {
      std::cerr << "ERROR: missing numeric event-cut key " << cutNames[j] << std::endl;
      gSystem->Exit(1); return;
    }
  }
  const Double_t minVz = cuts.GetValue("minVz", 0.);
  const Double_t maxVz = cuts.GetValue("maxVz", 0.);
  const Double_t maxVr = cuts.GetValue("maxVr", 0.);
  // These are EventCutConfig's documented defaults when absent in the YAML.
  const Double_t centerX = cuts.GetValue("vtxCenterX", 0.);
  const Double_t centerY = cuts.GetValue("vtxCenterY", 0.);
  const Double_t minRef = cuts.GetValue("minRefMult", 0.);
  const Double_t maxRef = cuts.GetValue("maxRefMult", 0.);
  const Double_t maxVzDiff = cuts.GetValue("maxVzDiff", 0.);
  const Double_t maxAbsVpd = cuts.GetValue("maxAbsVzVpd", 0.);
  const Int_t maxNTr = cuts.GetValue("maxNTr", 0);
  if (maximumEvents <= 0 || minVz >= maxVz || maxVr <= 0. || maxRef <= minRef || maxNTr > 0) {
    std::cerr << "ERROR: invalid limits; this Event-only diagnostic requires maxNTr disabled "
                 "(stored numberOfGlobalTracks is not the Pico Track array count)" << std::endl;
    gSystem->Exit(1); return;
  }
  std::cout << "Input: " << filename << "\nEvent cuts: " << eventConfig
            << "\nVz=[" << minVz << "," << maxVz << "]; Vr<=" << maxVr
            << " about (" << centerX << "," << centerY << "); refMult=["
            << minRef << "," << maxRef << "]; valid |VzVPD|<" << maxAbsVpd
            << " then |Vz-VzVPD|<=" << maxVzDiff << "; maxNTr=" << maxNTr << std::endl;
  if (gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C") < 0) {
    gSystem->Exit(1); return;
  }
  Int_t loadError = 0;
  gROOT->ProcessLine("loadSharedLibraries();", &loadError);
  if (loadError || gSystem->Load("StPicoEvent") < 0) { gSystem->Exit(1); return; }
  TFile* file = TFile::Open(filename, "READ");
  if (!file || file->IsZombie()) { gSystem->Exit(1); return; }
  TTree* tree = dynamic_cast<TTree*>(file->Get("PicoDst"));
  if (!tree) { file->Close(); delete file; gSystem->Exit(1); return; }
  tree->SetBranchStatus("*", 0);
  tree->SetBranchStatus("Event*", 1);
  const char* names[10] = {"mPrimaryVertexX", "mPrimaryVertexY", "mPrimaryVertexZ",
      "mVzVpd", "mRefMultPos", "mRefMultNeg", "mNumberOfPrimaryTracks",
      "mNBTOFMatch", "mRunId", "mNumberOfGlobalTracks"};
  TLeaf* leaves[10];
  for (j = 0; j < 10; ++j) {
    leaves[j] = KfInputEventLeaf(tree, names[j]);
    if (!leaves[j]) {
      std::cerr << "ERROR: missing Event leaf " << names[j] << std::endl;
      file->Close(); delete file; gSystem->Exit(1); return;
    }
  }
  Long64_t total = tree->GetEntries();
  Long64_t requested = total < maximumEvents ? total : maximumEvents;
  Long64_t failures[6] = {0,0,0,0,0,0};
  Long64_t sequential[6] = {0,0,0,0,0,0};
  Long64_t zeroVpd = 0, validVpd = 0, invalidVpd = 0, vpdRejectBeforeOther = 0;
  Double_t sums[4] = {0.,0.,0.,0.};
  Double_t minimum[4] = {DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX};
  Double_t maximum[4] = {-DBL_MAX,-DBL_MAX,-DBL_MAX,-DBL_MAX};
  Long64_t finiteCounts[4] = {0,0,0,0};
  Double_t vzWindowSumX = 0., vzWindowSumY = 0.;
  Long64_t vzWindowCount = 0;
  std::cout << "Stored tree entries=" << total << "; inspecting first " << requested << std::endl;
  for (Long64_t entry = 0; entry < requested; ++entry) {
    if (tree->GetEntry(entry) <= 0) {
      std::cerr << "ERROR: input read failed at " << entry << std::endl;
      file->Close(); delete file; gSystem->Exit(2); return;
    }
    Double_t x = leaves[0]->GetValue(), y = leaves[1]->GetValue(), z = leaves[2]->GetValue();
    Double_t vpd = leaves[3]->GetValue();
    Double_t ref = leaves[4]->GetValue() + leaves[5]->GetValue();
    Double_t dx = x - centerX, dy = y - centerY;
    Double_t vr = TMath::Sqrt(dx * dx + dy * dy);
    Double_t coordinates[4] = {x,y,z,vpd};
    for (j = 0; j < 4; ++j) if (TMath::Finite(coordinates[j])) {
      ++finiteCounts[j]; sums[j] += coordinates[j];
      if (coordinates[j] < minimum[j]) minimum[j] = coordinates[j];
      if (coordinates[j] > maximum[j]) maximum[j] = coordinates[j];
    }
    Bool_t finitePv = TMath::Finite(x) && TMath::Finite(y) && TMath::Finite(z) && TMath::Finite(vr);
    Bool_t failLowZ = finitePv && z < minVz;
    Bool_t failHighZ = finitePv && z > maxVz;
    Bool_t failVr = finitePv && vr > maxVr;
    Bool_t failRef = ref < minRef || ref > maxRef;
    Bool_t hasVpd = TMath::Finite(vpd) && TMath::Abs(vpd) < maxAbsVpd;
    Bool_t failVpd = hasVpd && TMath::Abs(z - vpd) > maxVzDiff;
    if (vpd == 0.) ++zeroVpd;
    if (hasVpd) ++validVpd; else ++invalidVpd;
    if (failVpd) ++vpdRejectBeforeOther;
    if (!finitePv) ++failures[0];
    if (failLowZ) ++failures[1];
    if (failHighZ) ++failures[2];
    if (failVr) ++failures[3];
    if (failRef) ++failures[4];
    if (failVpd) ++failures[5];
    if (entry < 8)
      std::cout << "event " << entry << " run=" << leaves[8]->GetValue()
                << " PV=(" << x << "," << y << "," << z << ") Vr=" << vr
                << " VPD=" << vpd << " ref=" << ref << " fxt=" << leaves[6]->GetValue()
                << " nTOF=" << leaves[7]->GetValue() << std::endl;
    if (!finitePv) continue;
    ++sequential[0];
    if (failLowZ || failHighZ) continue;
    ++sequential[1];
    ++vzWindowCount; vzWindowSumX += x; vzWindowSumY += y;
    if (failVr) continue;
    ++sequential[2];
    if (failRef) continue;
    ++sequential[3];
    if (failVpd) continue;
    ++sequential[4];
    ++sequential[5]; // maxNTr disabled, as explicitly checked above.
  }
  const char* coordinateNames[4] = {"PVx", "PVy", "PVz", "VzVPD"};
  for (j = 0; j < 4; ++j) if (finiteCounts[j])
    std::cout << coordinateNames[j] << ": finite=" << finiteCounts[j]
              << " range=[" << minimum[j] << "," << maximum[j] << "] mean="
              << sums[j] / Double_t(finiteCounts[j]) << std::endl;
  if (vzWindowCount)
    std::cout << "Within configured Vz window: " << vzWindowCount << " events; mean(x,y)=("
              << vzWindowSumX / Double_t(vzWindowCount) << ","
              << vzWindowSumY / Double_t(vzWindowCount) << ")" << std::endl;
  std::cout << "Independent failures (overlap): nonfinitePV=" << failures[0]
            << " lowVz=" << failures[1] << " highVz=" << failures[2]
            << " Vr=" << failures[3] << " refMult=" << failures[4] << " VPDdifference=" << failures[5]
            << "\nVPD values: exactly zero=" << zeroVpd << " valid-by-current-rule=" << validVpd
            << " invalid-or-outside-validity=" << invalidVpd << std::endl;
  const char* stepNames[6] = {"finite PV", "Vz window", "Vr", "refMult", "VPD difference", "maxNTr (disabled)"};
  for (j = 0; j < 6; ++j)
    std::cout << "Sequential survivors after " << stepNames[j] << ": " << sequential[j] << std::endl;
  std::cout << "This diagnostic does NOT apply bad-run, pileup, corrected-centrality, or KF cuts. "
               "Compare the last event-only count with the completed run's reconstructed-event count; "
               "no configuration has been changed." << std::endl;
  file->Close(); delete file;
  gSystem->Exit(requested > 0 ? 0 : 2);
}

