#include "TROOT.h"
#include "TSystem.h"
#include "TString.h"
#include "TInterpreter.h"
#include <iostream>
TString FemtoLambdaQaQuote(const char* value) {
  TString quoted(value ? value : "");
  quoted.ReplaceAll("\\", "\\\\"); quoted.ReplaceAll("\"", "\\\"");
  quoted.ReplaceAll("\n", "\\n"); quoted.ReplaceAll("\r", "\\r");
  return TString("\"") + quoted + "\"";
}
void run_checkHistAnaFemtoLambda(const char* input, const char* mainconf, const char* outputPdf) {
  const TString cwd = gSystem->WorkingDirectory();
  if (gSystem->Load(cwd + "/lib/libStarAnaConfig.so") < 0) { gSystem->Exit(1); return; }
  gInterpreter->AddIncludePath((cwd + "/include").Data());
  gSystem->AddLinkedLibs(TString::Format("-L%s/lib -lStarAnaConfig -Wl,-rpath,%s/lib", cwd.Data(), cwd.Data()));
  const TString build = TString::Format("%s/tmp/femto-lambda-qa/%s-%d", cwd.Data(), gSystem->HostName(), gSystem->GetPid());
  gSystem->mkdir(build, kTRUE); gSystem->SetBuildDir(build, kTRUE);
  if (!gSystem->CompileMacro((cwd + "/common/macro/checkHistAnaFemtoLambda.C").Data(), "kf")) { gSystem->Exit(1); return; }
  TString call = TString::Format("checkHistAnaFemtoLambda(%s,%s,%s)",
      FemtoLambdaQaQuote(input).Data(), FemtoLambdaQaQuote(mainconf).Data(), FemtoLambdaQaQuote(outputPdf).Data());
  Int_t error = 0; Long_t result = 1;
  result = gROOT->ProcessLine(call.Data(), &error);
  if (error || result) gSystem->Exit(error ? 1 : static_cast<Int_t>(result));
}
