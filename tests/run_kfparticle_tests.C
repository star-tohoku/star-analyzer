// Run compiled KF tests through ROOT 5's normal STAR library-loading path.
// Installed SL24y StarRoot has an unrelated, undefined THelix3d assignment
// symbol. Do not fake that method or relax executable link checks: run these
// shared test modules exactly as the production compiled Maker is loaded.
#include "TROOT.h"
#include "TSystem.h"
#include "TString.h"
#include <iostream>

// ROOT 5 CINT requires function-pointer typedefs at file scope.
typedef int (*StarAnalyzerKfpTestEntry)(const char*);

void run_kfparticle_tests(const char* testName = "full-chain",
                          const char* argument = "") {
  // A crash must fail the test process, not return to CINT and exit with zero.
  gSystem->ResetSignal(kSigSegmentationViolation, kTRUE);
  gSystem->ResetSignal(kSigBus, kTRUE);
  gSystem->ResetSignal(kSigIllegalInstruction, kTRUE);
  gSystem->ResetSignal(kSigFloatingException, kTRUE);
  const TString test(testName ? testName : "");
  const Bool_t pico = test == "pico-adapter";
  if (test != "full-chain" && !pico) {
    std::cerr << "ERROR: expected full-chain or pico-adapter test" << std::endl;
    gSystem->Exit(1); return;
  }
  if (gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C") < 0) {
    std::cerr << "ERROR: cannot load STAR loadSharedLibraries macro" << std::endl;
    gSystem->Exit(1); return;
  }
  Int_t error = 0;
  gROOT->ProcessLine("loadSharedLibraries();", &error);
  if (error) { gSystem->Exit(1); return; }
  // Old global KFParticle must be present before loading the isolated new core.
  if (gSystem->Load("StarRoot") < 0) { gSystem->Exit(1); return; }
  const TString directory = TString(gSystem->WorkingDirectory()) + "/lib/";
  if (pico) {
    const char* starLibraries[] = {"StBichsel", "StPicoEvent", "StPicoDstMaker"};
    for (UInt_t i = 0; i < sizeof(starLibraries)/sizeof(starLibraries[0]); ++i)
      if (gSystem->Load(starLibraries[i]) < 0) {
        std::cerr << "ERROR: loading " << starLibraries[i] << std::endl;
        gSystem->Exit(1); return;
      }
    if (gSystem->Load(directory + "libStarAnaConfig.so") < 0) {
      gSystem->Exit(1); return;
    }
  }
  if (gSystem->Load(directory + "libKFParticle.so") < 0) {
    gSystem->Exit(1); return;
  }
  if (pico && gSystem->Load(directory + "libStKfParticleCommon.so") < 0) {
    gSystem->Exit(1); return;
  }
  const TString module = directory + (pico ? "test_kfparticle_pico_adapter.so"
                                          : "test_kfparticle_full_chain.so");
  if (gSystem->Load(module) < 0) {
    std::cerr << "ERROR: loading test module " << module << std::endl;
    gSystem->Exit(1); return;
  }
  const char* symbol = pico ? "star_analyzer_kfp_pico_adapter_test"
                            : "star_analyzer_kfp_full_chain_test";
  StarAnalyzerKfpTestEntry entry =
      (StarAnalyzerKfpTestEntry)gSystem->DynFindSymbol(module.Data(), symbol);
  if (!entry) {
    std::cerr << "ERROR: missing compiled test entry " << symbol << std::endl;
    gSystem->Exit(1); return;
  }
  // Keep re-entrant CINT model/table initialization out of a declaration.
  int result = 1;
  result = entry(argument ? argument : "");
  std::cout << "[KF ROOT-native test] " << test << " status=" << result << std::endl;
  gSystem->Exit(result);
}
