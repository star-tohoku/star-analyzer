#ifndef FEMTO_LAMBDA_RUN_SUPPORT_H
#define FEMTO_LAMBDA_RUN_SUPPORT_H
// I/O and StChain lifecycle only. Reconstruction and pair physics live in Makers.
#include "TChain.h"
#include "TFile.h"
#include "TStopwatch.h"
#include "TString.h"
#include "TSystem.h"
#include "StChain.h"
#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StMaker/StFemtoMaker/StFemtoMaker.h"
#include "ConfigManager.h"
#include "YamlParser.h"
#include "cuts/FemtoConfig.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <map>
#include <string>
#include <vector>

namespace FemtoLambdaRunSupport {
inline TString Absolute(const char* path) {
  const TString value(path ? path : "");
  if (value.IsNull() || value.BeginsWith("/") || value.Contains("://")) return value;
  return TString(gSystem->WorkingDirectory()) + "/" + value;
}

struct Run {
  TString input, output, config, job;
  Long64_t requested;
  std::vector<std::string> inputFiles;
  Run(const char* in, const char* out, const char* id, Long64_t n, const char* conf)
      : input(Absolute(in)), output(Absolute(out)), config(Absolute(conf)),
        job(id ? id : ""), requested(n) {}

  Bool_t Prepare(const char* expectedSpecies) {
    if (input.IsNull() || output.IsNull() || config.IsNull() || job.IsNull()
        || requested == 0 || requested < -1) {
      std::cerr << "ERROR: explicit input, output, jobid and mainconf arguments are required; "
                   "nEvents must be positive or -1. No environment/default fallback." << std::endl;
      return kFALSE;
    }
    std::cout << "[mainconf] (argument) " << config << "\n[input] (argument) " << input
              << "\n[output] (argument) " << output << "\n[jobid] (argument) " << job << std::endl;
    if (!gSystem->AccessPathName(output)) {
      std::cerr << "ERROR: refusing to overwrite existing output " << output << std::endl;
      return kFALSE;
    }
    std::map<std::string, std::string> keys;
    if (!YamlParser::ParseFile(config, keys)) return kFALSE;
    const char* required[] = {"analysis", "event", "centrality", "kf", "nuclearid", "maker",
                              "mixing", "femtoHist", "nuclearHist"};
    for (unsigned i = 0; i < sizeof(required)/sizeof(required[0]); ++i) {
      if (keys.find(required[i]) == keys.end() || keys[required[i]].empty()) {
        std::cerr << "ERROR: new KF Femto mainconf requires '" << required[i]
                  << "'. Legacy per-species Helix mainconfs are incompatible; "
                     "use main_auau<energy>_anaFemtoLambda_<species>_KFParticle_highpurity.yaml."
                  << std::endl;
        return kFALSE;
      }
    }
    if (!ConfigManager::GetInstance().LoadConfig(config)) return kFALSE;
    const FemtoConfig& femto = ConfigManager::GetInstance().GetFemtoConfig();
    if (femto.species.size() != 2 || femto.species.find("lambda") == femto.species.end()
        || femto.species.find(expectedSpecies) == femto.species.end()) {
      std::cerr << "ERROR: entry requires exactly lambda and " << expectedSpecies << std::endl;
      return kFALSE;
    }
    if (input.EndsWith(".root")) inputFiles.push_back(input.Data());
    else {
      std::ifstream stream(input.Data());
      if (!stream) { std::cerr << "ERROR: unreadable input list " << input << std::endl; return kFALSE; }
      std::string line;
      while (std::getline(stream, line)) {
        const std::string::size_type comment = line.find('#');
        if (comment != std::string::npos) line.erase(comment);
        std::istringstream fields(line);
        std::string path, extra;
        if (!(fields >> path)) continue;
        if ((fields >> extra) || path.find(".picoDst.root") == std::string::npos) {
          std::cerr << "ERROR: malformed input-list entry: " << line << std::endl; return kFALSE;
        }
        inputFiles.push_back(path);
      }
    }
    if (inputFiles.empty()) { std::cerr << "ERROR: empty input list" << std::endl; return kFALSE; }
    // Check every listed file, including files beyond an early event limit. Missing
    // trees/branches must not become an apparently complete job on a reduced list.
    Long64_t available = 0;
    for (unsigned i = 0; i < inputFiles.size(); ++i) {
      TFile* file = TFile::Open(inputFiles[i].c_str(), "READ");
      if (!file || file->IsZombie()) {
        std::cerr << "ERROR: unreadable PicoDst " << inputFiles[i] << std::endl;
        delete file; return kFALSE;
      }
      TTree* tree = dynamic_cast<TTree*>(file->Get("PicoDst"));
      if (!CheckTree(tree, inputFiles[i].c_str())) { delete file; return kFALSE; }
      const Long64_t entries = tree->GetEntries();
      std::cout << "[input-file] " << i << " entries=" << entries << " " << inputFiles[i] << std::endl;
      if (entries <= 0) { std::cerr << "ERROR: empty PicoDst" << std::endl; delete file; return kFALSE; }
      available += entries;
      delete file;
    }
    if (requested > available) {
      std::cerr << "ERROR: requested=" << requested << " exceeds available=" << available
                << "; refusing a silently truncated run" << std::endl;
      return kFALSE;
    }
    if (requested == -1) requested = available;
    // ROOT5 DirName returns a reusable buffer; copy before recursive mkdir.
    const TString outputDirectory = gSystem->DirName(output.Data());
    if (gSystem->mkdir(outputDirectory.Data(), kTRUE) != 0
        && gSystem->AccessPathName(outputDirectory.Data())) {
      std::cerr << "ERROR: cannot create output directory " << outputDirectory << std::endl;
      return kFALSE;
    }
    return kTRUE;
  }

  static Bool_t CheckTree(TTree* tree, const char* source) {
    if (!tree) { std::cerr << "ERROR: missing PicoDst tree in " << source << std::endl; return kFALSE; }
    // TOF QA remains part of the nuclear histogram contract even when final m2 PID is off.
    const char* required[] = {"Event", "Track", "TrackCovMatrix", "BTofPidTraits"};
    for (unsigned b = 0; b < sizeof(required)/sizeof(required[0]); ++b) {
      const TString dotted = TString(required[b]) + ".";
      if (!tree->GetBranch(required[b]) && !tree->GetBranch(dotted)) {
        std::cerr << "ERROR: missing " << required[b] << " in " << source << std::endl;
        return kFALSE;
      }
    }
    return kTRUE;
  }

  static void EnableBranches(StPicoDstMaker* pico) {
    pico->SetStatus("*", 0);
    pico->SetStatus("Event", 1);
    pico->SetStatus("Track", 1);
    pico->SetStatus("TrackCovMatrix", 1);
    pico->SetStatus("BTofHit", 1);
    pico->SetStatus("BTofPidTraits", 1);
  }

  Int_t Execute(StChain* chain, StPicoDstMaker* pico, StFemtoMaker* maker) {
    TStopwatch timer;
    timer.Start();
    if (chain->Init() != kStOK) {
      std::cerr << "ERROR: StChain Init failed" << std::endl;
      delete chain; return 2;
    }
    TChain* source = pico->chain();
    Int_t status = 0;
    if (!source || source->GetEntries() < requested || !source->GetListOfFiles()
        || source->GetListOfFiles()->GetEntries() != static_cast<Int_t>(inputFiles.size())) {
      std::cerr << "ERROR: Pico reader skipped files/entries" << std::endl; status = 2;
    }
    Long64_t completed = 0;
    Int_t previousTree = -1;
    for (Long64_t i = 0; !status && i < requested; ++i) {
      if (source->LoadTree(i) < 0 || !source->GetTree()) { status = 3; break; }
      if (source->GetTreeNumber() != previousTree) {
        previousTree = source->GetTreeNumber();
        if (!CheckTree(source->GetTree(), source->GetCurrentFile()->GetName())) { status = 3; break; }
      }
      if (i % 1000 == 0) std::cout << "Working on event " << i << std::endl;
      chain->Clear();
      const Int_t result = chain->Make(i);
      if (result != kStOK) {
        std::cerr << "ERROR: chain Make=" << result << " entry=" << i << std::endl;
        status = 3; break;
      }
      ++completed;
    }
    if (completed != requested || maker->GetReconstructedEvents() <= 0) status = status ? status : 3;
    maker->SetProcessingSucceeded(status == 0);
    if (chain->Finish() != kStOK) status = 4;
    timer.Stop();
    std::cout << "FemtoLambda job=" << job << " requested=" << requested << " completed=" << completed
              << " reconstructed=" << maker->GetReconstructedEvents() << " status=" << status
              << " RealTime=" << timer.RealTime() << " CpuTime=" << timer.CpuTime() << std::endl;
    delete chain;
    return status;
  }
};
}  // namespace FemtoLambdaRunSupport
#endif
