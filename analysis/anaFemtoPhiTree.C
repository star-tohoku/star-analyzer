// anaFemtoPhiTree.C - StChain loop for StFemtoPhiTreeMaker
// Usage: from project root via ./script/run_anaFemtoPhiTree.sh

#include "TROOT.h"
#include "TInterpreter.h"
#include "TSystem.h"
#include "TStopwatch.h"
#include "TString.h"
#include "TChain.h"
#include "StChain.h"
#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StMaker/StFemtoPhiTreeMaker/StFemtoPhiTreeMaker.h"
#include "ConfigManager.h"
#include <iostream>
#include <fstream>

StChain* chain = 0;
StFemtoPhiTreeMaker* treeMaker = 0;

void anaFemtoPhiTree(const Char_t* inputFile = "tmp/pico/bench1.picoDst.root",
                     const Char_t* outputFile = "rootfile/auau3p85fxt_anaFemtoPhiTree_dEdxOnly/temp.root",
                     const Char_t* jobid = "0",
                     Long64_t nEventsMax = -1,
                     const Char_t* configPath = 0,
                     Long64_t nSkip = 0)
{
  TStopwatch timer;
  timer.Start();

  Long64_t nEvents = (nEventsMax > 0) ? nEventsMax : 10000000;

  // Relative paths resolve against the process's real working directory, not against $PWD.
  // $PWD is a shell variable, and a SUMS job is a csh script that cd's into its runtime bundle;
  // csh does not update $PWD on cd, so on the farm it named a directory the job was not in.
  // See results/pilot-farm-20260915.md.
  TString cwd = gSystem->WorkingDirectory();
  const char* pwd = cwd.Data();

  // Choosing the configuration. The rule is that a configuration named on the command line is
  // used or the job fails: it must never be quietly replaced by another one. The 2026-09-15 farm
  // pilot produced 241 trees with the wrong selection precisely because the fallback chain below
  // used to run even when a path had been given.
  TString mainConfigPath;
  TString configSource;
  if (configPath && strlen(configPath) > 0) {
    mainConfigPath = configPath;
    if (mainConfigPath(0) != '/') mainConfigPath = cwd + "/" + mainConfigPath;
    configSource = "argument";
    if (gSystem->AccessPathName(mainConfigPath.Data())) {
      std::cerr << "ERROR: the mainconf given on the command line does not exist:\n"
                << "       " << mainConfigPath.Data() << "\n"
                << "       (working directory " << cwd.Data() << ")\n"
                << "       Refusing to fall back to a different configuration." << std::endl;
      return;
    }
  } else {
    const char* env_conf = gSystem->Getenv("STAR_ANA_MAINCONF");
    if (env_conf && strlen(env_conf) > 0) {
      mainConfigPath = env_conf;
      if (mainConfigPath(0) != '/') mainConfigPath = cwd + "/" + mainConfigPath;
      configSource = "STAR_ANA_MAINCONF";
    } else {
      TString fallbackPath = cwd + "/.current_mainconf";
      std::ifstream infile(fallbackPath.Data());
      std::string line;
      if (infile.is_open() && std::getline(infile, line) && !line.empty()) {
        mainConfigPath = line.c_str();
        if (mainConfigPath(0) != '/') mainConfigPath = cwd + "/" + mainConfigPath;
        configSource = ".current_mainconf";
      } else {
        std::cerr << "ERROR: no mainconf. Pass one as the 5th argument, set STAR_ANA_MAINCONF, or "
                  << "write one into .current_mainconf.\n"
                  << "       There is no built-in default: a job with no configuration must fail, "
                  << "not analyse with someone else's." << std::endl;
        return;
      }
    }
    if (gSystem->AccessPathName(mainConfigPath.Data())) {
      std::cerr << "ERROR: mainconf from " << configSource.Data() << " does not exist:\n"
                << "       " << mainConfigPath.Data() << std::endl;
      return;
    }
  }
  // Printed unconditionally: when a batch job analyses with the wrong configuration, this line is
  // what makes it visible in the log rather than only in the output months later.
  std::cout << "[anaFemtoPhiTree] mainconf (" << configSource.Data() << "): "
            << mainConfigPath.Data() << std::endl;

  if (!ConfigManager::GetInstance().LoadConfig(mainConfigPath.Data())) {
    std::cerr << "ERROR: Failed to load config: " << mainConfigPath.Data() << std::endl;
    return;
  }

  chain = new StChain();
  StPicoDstMaker* picoMaker = new StPicoDstMaker(StPicoDstMaker::IoRead, inputFile, "picoDst");
  // Only what StFemtoPhiTreeMaker actually reads. The tree maker uses StPicoEvent, StPicoTrack and
  // StPicoBTofPidTraits; CentralityHelper and StRefMultCorr take refMult and nBTOFMatch from the
  // event row, and nothing in this chain touches BbcHit / EpdHit / MtdHit / BTowHit /
  // ETofPidTraits / BTofHit (checked by grep over StMaker/, include/ and StRoot/StRefMultCorr).
  // At 2.28e9 events the I/O this saves is the dominant cost of production.
  picoMaker->SetStatus("*", 0);
  picoMaker->SetStatus("Event", 1);
  picoMaker->SetStatus("Track", 1);
  picoMaker->SetStatus("BTofPidTraits", 1);

  treeMaker = new StFemtoPhiTreeMaker("femtoPhiTree", picoMaker, outputFile);
  // Batch jobs reach this macro directly from the joblist, with no wrapper script and therefore
  // no STAR_ANA_JOBID in the environment; the jobid argument is the only copy that exists there.
  treeMaker->SetJobId(jobid);
  // Same reason: the tree-specific keys are re-parsed from the mainconf by the maker itself, and
  // on the farm there is no STAR_ANA_MAINCONF to find it by.
  treeMaker->SetMainconfPath(mainConfigPath.Data());
  chain->AddMaker(picoMaker);
  chain->AddMaker(treeMaker);

  if (chain->Init() == kStErr) {
    std::cerr << "ERROR: chain->Init() returned kStErr" << std::endl;
    delete chain;
    chain = 0;
    treeMaker = 0;
    return;
  }

  Long64_t totalEntries = picoMaker->chain() ? picoMaker->chain()->GetEntries() : 0;
  std::cout << "Total entries = " << totalEntries << std::endl;

  if (totalEntries <= 0) {
    std::cerr << "ERROR: no entries found. Check inputFile." << std::endl;
    chain->Finish();
    delete chain;
    chain = 0;
    treeMaker = 0;
    return;
  }

  if (nSkip < 0) nSkip = 0;
  Long64_t nLoop = nEvents + nSkip;
  if (nLoop > totalEntries) nLoop = totalEntries;
  if (nSkip >= nLoop) {
    std::cerr << "ERROR: nSkip >= available entries" << std::endl;
    nLoop = 0;
  }

  for (Long64_t i = 0; i < nLoop; i++) {
    if (i % 1000 == 0) std::cout << "Working on event " << i << std::endl;
    chain->Clear();
    Int_t iret = chain->Make(i);
    if (iret) {
      std::cerr << "Bad return code: " << iret << " at event " << i << std::endl;
      break;
    }
  }

  std::cout << "******************************************" << std::endl;
  std::cout << "Work done... chain->Finish()" << std::endl;
  std::cout << "******************************************" << std::endl;
  chain->Finish();

  timer.Stop();
  std::cout << "Processed events: " << nEvents << std::endl;
  std::cout << "RealTime: " << timer.RealTime() << " CpuTime: " << timer.CpuTime() << std::endl;

  delete chain;
  chain = 0;
  treeMaker = 0;
}
