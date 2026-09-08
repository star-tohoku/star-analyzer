// CINT only orchestrates ACLiC. Function pointers and test execution stay in
// compiled C++, as in the production two-macro pattern.
#include "TROOT.h"
#include "TSystem.h"
#include "TString.h"
#include "TUUID.h"
#include <iostream>

TString StarAnalyzerKfpTestQuote(const char* value) {
  TString quoted(value ? value : "");
  quoted.ReplaceAll("\\", "\\\\");
  quoted.ReplaceAll("\"", "\\\"");
  quoted.ReplaceAll("\n", "\\n");
  quoted.ReplaceAll("\r", "\\r");
  return TString("\"") + quoted + "\"";
}

void bootstrap_kfparticle_tests(const char* testName = "full-chain",
                                const char* argument = "") {
  gSystem->ResetSignal(kSigSegmentationViolation, kTRUE);
  gSystem->ResetSignal(kSigBus, kTRUE);
  gSystem->ResetSignal(kSigIllegalInstruction, kTRUE);
  gSystem->ResetSignal(kSigFloatingException, kTRUE);
  // Independent make targets may run concurrently: never share an ACLiC cache.
  TUUID identifier;
  TString buildDirectory = TString(gSystem->TempDirectory()) +
      "/star_analyzer_kfp_tests_" + identifier.AsString();
  if (gSystem->mkdir(buildDirectory.Data(), kTRUE) != 0) {
    std::cerr << "ERROR: creating KF test build directory " << buildDirectory << std::endl;
    gSystem->Exit(1); return;
  }
  TString source = TString(gSystem->WorkingDirectory()) + "/tests/run_kfparticle_tests.C";
  // No keep option: ACLiC may clean its generated temporary library at exit.
  if (!gSystem->CompileMacro(source.Data(), "", "", buildDirectory.Data())) {
    std::cerr << "ERROR: compiling KF test runner" << std::endl;
    gSystem->Exit(1); return;
  }
  TString call = TString("run_kfparticle_tests(") +
      StarAnalyzerKfpTestQuote(testName) + "," +
      StarAnalyzerKfpTestQuote(argument) + ")";
  Int_t error = 0;
  // Do not invoke re-entrant CINT from a const declaration initializer.
  gROOT->ProcessLine(call.Data(), &error);
  if (error) { gSystem->Exit(1); return; }
  // The compiled runner normally exits explicitly with the test's result.
  std::cerr << "ERROR: compiled KF runner returned without an explicit test status" << std::endl;
  gSystem->Exit(1);
}
