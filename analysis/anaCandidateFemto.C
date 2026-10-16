// Compiled StChain entry for StCandidateFemtoMaker.
// PicoDst input (a .picoDst.root file or a list of them): StPicoDstMaker -> StCandidateFemtoMaker. The
// reconstruction parts are created here and handed to the Maker (it owns them): the event selection and every
// CandidateSource the species of the candidateFemto YAML need. The Maker itself does not depend on KFParticle.
// With writeCandidateStore: true this is the reconstruction stage (the output is a candidate store).
// Candidate-store input (a store .root file or a list of them, every other .root entry): the pair stage, no
// PicoDst and no parts (StCandidateFemtoMaker::SetCandidateStoreInput). A list must not mix the two kinds.
#include "TChain.h"
#include "TStopwatch.h"
#include "TString.h"
#include "TSystem.h"
#include "StChain.h"
#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StMaker/StCandidateFemtoMaker/StCandidateFemtoMaker.h"
#include "StMaker/kfparticle/KfParticleCandidateSource.h"
#include "StMaker/kfparticle/OfficialEventSelection.h"
#include "StMaker/common/TrackPidCandidateSource.h"
#include "ConfigManager.h"
#include "cuts/KfParticleCutConfig.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

// Pair stage: the Maker reads the candidate stores (no PicoDst, no event selector, no reconstruction part).
Int_t anaCandidateFemtoFromStore(const Char_t* inputFile, const Char_t* outputFile, const Char_t* jobid,
                                 Long64_t nEventsMax, const TString& config, TStopwatch& timer) {
  // The stored reconstruction-stage histograms of every store are added to the output, so the pair stage reads whole stores:
  // a partial read would give histograms of more events than the pairs. For a short test give fewer stores.
  if (nEventsMax > 0) {
    std::cerr << "ERROR: pair stage (candidate store input) reads whole stores: nEvents must be -1 (give fewer store files "
                 "for a short test)" << std::endl;
    return 1;
  }
  gSystem->mkdir(gSystem->DirName(outputFile), kTRUE);
  StChain* chain = new StChain();
  StCandidateFemtoMaker* femto = new StCandidateFemtoMaker("candidateFemto", 0, outputFile);
  femto->SetMainConfigPath(config.Data());
  femto->SetCandidateStoreInput(inputFile);
  const Int_t init = chain->Init();
  if (init != kStOK) {
    std::cerr << "ERROR: chain Init failed: " << init << std::endl;
    delete chain;
    return 1;
  }
  Int_t status = 0;
  Long64_t completed = 0;
  for (Long64_t i = 0; nEventsMax < 0 || i < nEventsMax; ++i) {
    if (i % 10000 == 0) std::cout << "Working on stored event " << i << std::endl;
    chain->Clear();
    const Int_t result = chain->Make(i);
    if (femto->CandidateStoreInputFinished()) break;
    if (result != kStOK) {
      std::cerr << "ERROR: chain Make code " << result << " at stored event " << i << std::endl;
      status = 3;
      break;
    }
    ++completed;
  }
  if (!status && femto->GetAcceptedEvents() <= 0) {
    std::cerr << "ERROR: no event in the candidate stores" << std::endl;
    status = 2;
  }
  if (chain->Finish() != kStOK) status = 4;
  timer.Stop();
  std::cout << "CandidateFemto job=" << jobid << " (pair stage from candidate stores) requested="
            << (nEventsMax < 0 ? -1 : nEventsMax) << " completed=" << completed
            << " accepted=" << femto->GetAcceptedEvents() << " sePairs=" << femto->GetSameEventPairs()
            << " mePairs=" << femto->GetMixedEventPairs() << " status=" << status
            << " RealTime=" << timer.RealTime() << " CpuTime=" << timer.CpuTime() << std::endl;
  delete chain;
  return status;
}

Int_t anaCandidateFemto(const Char_t* inputFile, const Char_t* outputFile, const Char_t* jobid,
                        Long64_t nEventsMax, const Char_t* configPath) {
  TStopwatch timer;
  timer.Start();
  if (!inputFile || !*inputFile || !outputFile || !*outputFile || nEventsMax == 0 || nEventsMax < -1) {
    std::cerr << "ERROR: input/output required; nEvents must be positive or -1" << std::endl;
    return 1;
  }
  // The mainconf arrives as an argument (AGENT_RULES: no environment-only values, no default).
  if (!configPath || !*configPath) {
    std::cerr << "ERROR: anaCandidateFemto requires the mainconf path as an argument" << std::endl;
    return 1;
  }
  TString config(configPath);
  if (config(0) != '/') config = TString(gSystem->WorkingDirectory()) + "/" + config;
  std::cout << "[anaCandidateFemto] mainconf (argument): " << config << std::endl;
  if (!ConfigManager::GetInstance().LoadConfig(config.Data())) return 1;
  if (!gSystem->AccessPathName(outputFile)) {
    std::cerr << "ERROR: output already exists; choose a new file: " << outputFile << std::endl;
    return 1;
  }

  Long64_t expectedFiles = 0, storeFiles = 0;
  const TString inputName(inputFile);
  if (!inputName.Contains(".list")) {
    if (inputName.EndsWith(".root") && !inputName.EndsWith(".picoDst.root")) storeFiles = 1;
  } else {
    std::ifstream list(inputFile);
    if (!list) { std::cerr << "ERROR: cannot open input list " << inputFile << std::endl; return 2; }
    std::string line;
    while (std::getline(list, line)) {
      const std::string::size_type comment = line.find('#');
      if (comment != std::string::npos) line.erase(comment);
      std::istringstream fields(line);
      std::string file;
      if (!(fields >> file)) continue;
      const TString f(file.c_str());
      if (f.EndsWith(".picoDst.root")) ++expectedFiles;
      else if (f.EndsWith(".root")) ++storeFiles;
      else {
        std::cerr << "ERROR: unsupported input-list entry: " << file << std::endl;
        return 2;
      }
    }
    if (!expectedFiles && !storeFiles) { std::cerr << "ERROR: input list is empty" << std::endl; return 2; }
    if (expectedFiles && storeFiles) {
      std::cerr << "ERROR: the input list mixes PicoDst files and candidate stores" << std::endl;
      return 2;
    }
  }
  if (storeFiles) return anaCandidateFemtoFromStore(inputFile, outputFile, jobid, nEventsMax, config, timer);
  gSystem->mkdir(gSystem->DirName(outputFile), kTRUE);
  StChain* chain = new StChain();
  StPicoDstMaker* pico = new StPicoDstMaker(StPicoDstMaker::IoRead, inputFile, "picoDst");
  pico->SetStatus("*", 0);
  pico->SetStatus("Event", 1);
  pico->SetStatus("Track", 1);
  pico->SetStatus("TrackCovMatrix", 1);
  pico->SetStatus("BTofHit", 1);
  pico->SetStatus("BTofPidTraits", 1);
  StCandidateFemtoMaker* femto = new StCandidateFemtoMaker("candidateFemto", pico, outputFile);
  femto->SetMainConfigPath(config.Data());
  femto->SetEventSelector(createOfficialEventSelection());
  femto->AddCandidateSource(createKfParticleCandidateSource());
  femto->AddCandidateSource(createTrackPidCandidateSource());
  const Int_t init = chain->Init();
  if (init != kStOK) {
    std::cerr << "ERROR: chain Init failed: " << init << std::endl;
    delete chain;
    return 1;
  }
  TChain* input = pico->chain();
  const Long64_t total = input ? input->GetEntries() : 0;
  Long64_t requested = nEventsMax > 0 ? nEventsMax : total;
  if (requested > total) requested = total;
  Int_t status = 0;
  if (expectedFiles && (!input || input->GetListOfFiles()->GetEntries() != expectedFiles)) {
    std::cerr << "ERROR: Pico reader skipped input files; expected=" << expectedFiles
              << " attached=" << (input ? input->GetListOfFiles()->GetEntries() : 0) << std::endl;
    status = 2;
  }
  if (total <= 0) {
    std::cerr << "ERROR: no input events" << std::endl;
    status = 2;
  }
  // BTofPidTraits is always required: TrackPidCandidateSource needs TOF whatever the kf useTof setting.
  const Int_t nRequired = 4;
  Long64_t completed = 0;
  for (Long64_t i = 0; !status && i < requested; ++i) {
    if (input->LoadTree(i) < 0 || !input->GetTree()) {
      std::cerr << "ERROR: unreadable input tree at entry " << i << std::endl;
      status = 2;
      break;
    }
    const char* required[] = {"Event", "Track", "TrackCovMatrix", "BTofPidTraits"};
    for (Int_t b = 0; b < nRequired; ++b) {
      const TString dotted = TString(required[b]) + ".";
      if (!input->GetTree()->GetBranch(required[b]) && !input->GetTree()->GetBranch(dotted)) {
        std::cerr << "ERROR: input lacks required branch " << required[b] << " at entry " << i << std::endl;
        status = 2;
        break;
      }
    }
    if (status) break;
    if (i % 1000 == 0) std::cout << "Working on event " << i << std::endl;
    chain->Clear();
    const Int_t result = chain->Make(i);
    if (result != kStOK) {
      std::cerr << "ERROR: chain Make code " << result << " at event " << i << std::endl;
      status = 3;
      break;
    }
    ++completed;
  }
  if (!status && femto->GetAcceptedEvents() <= 0) {
    std::cerr << "ERROR: no event passed the event selection" << std::endl;
    status = 2;
  }
  if (chain->Finish() != kStOK) status = 4;
  timer.Stop();
  std::cout << "CandidateFemto job=" << jobid << " requested=" << requested << " completed=" << completed
            << " accepted=" << femto->GetAcceptedEvents() << " sePairs=" << femto->GetSameEventPairs()
            << " mePairs=" << femto->GetMixedEventPairs() << " status=" << status
            << " RealTime=" << timer.RealTime() << " CpuTime=" << timer.CpuTime() << std::endl;
  delete chain;
  return status;
}
