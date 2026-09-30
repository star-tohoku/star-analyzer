// Run independently of KF libraries to verify optional provider isolation.
#include "TROOT.h"
#include "TSystem.h"
#include "TString.h"
#include "TUUID.h"
#include "TInterpreter.h"
#include <iostream>

void femto_lambda_maker_bootstrap(const char* mainconf) {
  gSystem->ResetSignal(kSigSegmentationViolation,kTRUE);
  gSystem->ResetSignal(kSigBus,kTRUE);
  gSystem->ResetSignal(kSigIllegalInstruction,kTRUE);
  gSystem->ResetSignal(kSigFloatingException,kTRUE);
  if(!mainconf || !*mainconf) { gSystem->Exit(2); return; }
  if(gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C")<0) {
    gSystem->Exit(2); return;
  }
  Int_t error=0;
  gROOT->ProcessLine("loadSharedLibraries();",&error);
  if(error) {gSystem->Exit(2);return;}
  const char* star[]={"StarRoot","StBichsel","StPicoEvent","StPicoDstMaker"};
  for(unsigned i=0;i<sizeof(star)/sizeof(star[0]);++i)
    if(gSystem->Load(star[i])<0) {gSystem->Exit(2);return;}
  TString cwd=gSystem->WorkingDirectory();
  const char* libraries[]={"libStarAnaConfig.so","libStRefMultCorr.so","libStCommon.so",
    "libStFemtoMaker.so","libStLambdaMaker.so","libStNuclearIdMaker.so","libStLambdaNuclearMixMaker.so"};
  for(unsigned i=0;i<sizeof(libraries)/sizeof(libraries[0]);++i) {
    if(gSystem->Load(cwd+"/lib/"+libraries[i])<0) {gSystem->Exit(2);return;}
    std::cout << "CORE-ONLY LOAD " << libraries[i] << std::endl;
  }
  gInterpreter->AddIncludePath(cwd);
  gInterpreter->AddIncludePath(cwd+"/include");
  gInterpreter->AddIncludePath(cwd+"/StMaker/common");
  gInterpreter->AddIncludePath("$STAR/StRoot");
  gSystem->AddLinkedLibs(TString::Format("-L%s/lib -lStFemtoMaker -lStCommon -lStarAnaConfig -lStRefMultCorr -Wl,-rpath,%s/lib",cwd.Data(),cwd.Data()));
  TUUID id;
  TString build=TString(gSystem->TempDirectory())+"/star_femto_maker_aclic_"+id.AsString();
  if(gSystem->mkdir(build,kTRUE)!=0) {gSystem->Exit(2);return;}
  if(!gSystem->CompileMacro((cwd+"/tests/femto_lambda_maker_checks.cxx").Data(),"k","",build)) {
    std::cerr << "ERROR: cannot compile Maker checks" << std::endl;gSystem->Exit(2);return;
  }
  TString conf(mainconf);
  conf.ReplaceAll("\\","\\\\");conf.ReplaceAll("\"","\\\"");
  TString call=TString::Format("femto_lambda_maker_checks(\"%s\")",conf.Data());
  Long_t result=1;
  result=gROOT->ProcessLine(call,&error);
  gSystem->Exit(error?2:static_cast<Int_t>(result));
}
