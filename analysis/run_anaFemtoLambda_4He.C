// ROOT5 interpreted loader; no vendor KF/SIMD headers are exposed to CINT.
#include "TROOT.h"
#include "TSystem.h"
#include "TString.h"
#include "TInterpreter.h"
#include <iostream>

TString FemtoLambda4HeQuote(const char* value) {
  TString quoted(value ? value : "");
  quoted.ReplaceAll("\\", "\\\\");
  quoted.ReplaceAll("\"", "\\\"");
  quoted.ReplaceAll("\n", "\\n");
  quoted.ReplaceAll("\r", "\\r");
  return TString("\"") + quoted + "\"";
}

void run_anaFemtoLambda_4He(const Char_t* inputFile, const Char_t* outputFile,
    const Char_t* jobid, Long64_t nEventsMax, const Char_t* configPath) {
  gSystem->ResetSignal(kSigSegmentationViolation, kTRUE);
  gSystem->ResetSignal(kSigBus, kTRUE);
  gSystem->ResetSignal(kSigIllegalInstruction, kTRUE);
  gSystem->ResetSignal(kSigFloatingException, kTRUE);
  if (!inputFile || !*inputFile || !outputFile || !*outputFile
      || !jobid || !*jobid || !configPath || !*configPath) {
    std::cerr << "ERROR: input/output/jobid/mainconf are mandatory arguments" << std::endl;
    gSystem->Exit(1); return;
  }
  TString cwd = gSystem->WorkingDirectory();
  if (gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C") < 0) {
    gSystem->Exit(1); return;
  }
  Int_t error = 0;
  gROOT->ProcessLine("loadSharedLibraries();", &error);
  if (error) { gSystem->Exit(1); return; }
  const char* starLibraries[] = {"StarRoot", "StBichsel", "StPicoEvent", "StPicoDstMaker"};
  for (UInt_t i = 0; i < sizeof(starLibraries)/sizeof(starLibraries[0]); ++i) {
    if (gSystem->Load(starLibraries[i]) < 0) { gSystem->Exit(1); return; }
  }
  const char* libraries[] = {"libStarAnaConfig.so", "libStRefMultCorr.so", "libKFParticle.so",
      "libStKfParticleCommon.so", "libStCommon.so", "libStFemtoMaker.so"};
  for (UInt_t i = 0; i < sizeof(libraries)/sizeof(libraries[0]); ++i) {
    TString path = cwd + "/lib/" + libraries[i];
    if (gSystem->Load(path) < 0) {
      std::cerr << "ERROR: loading " << path << std::endl;
      gSystem->Exit(1); return;
    }
  }
  gInterpreter->AddIncludePath(cwd.Data());
  gInterpreter->AddIncludePath((cwd + "/include").Data());
  gInterpreter->AddIncludePath((cwd + "/StMaker/common").Data());
  gInterpreter->AddIncludePath("$STAR/StRoot");
  gSystem->AddLinkedLibs(TString::Format(
      "-L%s/lib -lStFemtoMaker -lStKfParticleCommon -lKFParticle -lStCommon "
      "-lStarAnaConfig -lStRefMultCorr -Wl,-rpath,%s/lib", cwd.Data(), cwd.Data()));
  TString buildDir = TString::Format("%s/tmp/femto-lambda-aclic/%s-%d",
      cwd.Data(), gSystem->HostName(), gSystem->GetPid());
  gSystem->mkdir(buildDir, kTRUE);
  gSystem->SetBuildDir(buildDir, kTRUE);
  if (!gSystem->CompileMacro((cwd + "/analysis/anaFemtoLambda_4He.C").Data(), "kf")) {
    std::cerr << "ERROR: ACLiC compilation failed" << std::endl;
    gSystem->Exit(1); return;
  }
  TString call = TString::Format("anaFemtoLambda_4He(%s,%s,%s,%lld,%s)",
      FemtoLambda4HeQuote(inputFile).Data(), FemtoLambda4HeQuote(outputFile).Data(),
      FemtoLambda4HeQuote(jobid).Data(), nEventsMax, FemtoLambda4HeQuote(configPath).Data());
  Long_t result = 1;
  result = gROOT->ProcessLine(call.Data(), &error);
  if (error || result) gSystem->Exit(error ? 1 : static_cast<Int_t>(result));
}
