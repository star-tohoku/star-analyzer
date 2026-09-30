#include "TROOT.h"
#include "TSystem.h"
#include "TString.h"
#include "TInterpreter.h"
// Invoke with tests/run_femto_lambda_root5.sh after the batch-matched build.
void run_femto_lambda_preflight(const char* mainconf, const char* newFixtureDirectory) {
  if (gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C") < 0) { gSystem->Exit(1); return; }
  Int_t error = 0;
  gROOT->ProcessLine("loadSharedLibraries();", &error);
  if (error) { gSystem->Exit(1); return; }
  const char* star[] = {"StarRoot", "StBichsel", "StPicoEvent", "StPicoDstMaker"};
  for (UInt_t i = 0; i < sizeof(star)/sizeof(star[0]); ++i)
    if (gSystem->Load(star[i]) < 0) { gSystem->Exit(1); return; }
  const TString cwd = gSystem->WorkingDirectory();
  const char* libs[] = {"libStarAnaConfig.so", "libStRefMultCorr.so", "libStCommon.so", "libStFemtoMaker.so"};
  for (UInt_t i = 0; i < sizeof(libs)/sizeof(libs[0]); ++i)
    if (gSystem->Load(cwd + "/lib/" + libs[i]) < 0) { gSystem->Exit(1); return; }
  gInterpreter->AddIncludePath(cwd.Data());
  gInterpreter->AddIncludePath((cwd + "/include").Data());
  gInterpreter->AddIncludePath((cwd + "/StMaker/common").Data());
  gInterpreter->AddIncludePath("$STAR/StRoot");
  gSystem->AddLinkedLibs(TString::Format("-L%s/lib -lStFemtoMaker -lStCommon -lStarAnaConfig -lStRefMultCorr -Wl,-rpath,%s/lib", cwd.Data(), cwd.Data()));
  const TString build = TString::Format("%s/tmp/femto-lambda-preflight/%s-%d", cwd.Data(), gSystem->HostName(), gSystem->GetPid());
  gSystem->mkdir(build, kTRUE); gSystem->SetBuildDir(build, kTRUE);
  if (!gSystem->CompileMacro((cwd + "/tests/femto_lambda_preflight.C").Data(), "kf")) { gSystem->Exit(1); return; }
  TString config(mainconf), dir(newFixtureDirectory);
  config.ReplaceAll("\\", "\\\\"); config.ReplaceAll("\"", "\\\"");
  dir.ReplaceAll("\\", "\\\\"); dir.ReplaceAll("\"", "\\\"");
  const TString call = TString::Format("femto_lambda_preflight(\"%s\",\"%s\")", config.Data(), dir.Data());
  Long_t result = 1;
  result = gROOT->ProcessLine(call.Data(), &error);
  if (error || result) gSystem->Exit(error ? 1 : static_cast<Int_t>(result));
}
