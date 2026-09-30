// Compile/load the unchanged Phi/Kaon entry points with no optional KF module.
// This is a loader/ACLiC regression check, not an event-yield comparison.
#include "TROOT.h"
#include "TSystem.h"
#include "TString.h"
#include "TUUID.h"
#include "TInterpreter.h"
#include <iostream>

void femto_lambda_maker_core_load(const char* entry) {
  gSystem->ResetSignal(kSigSegmentationViolation,kTRUE);
  gSystem->ResetSignal(kSigBus,kTRUE);
  gSystem->ResetSignal(kSigIllegalInstruction,kTRUE);
  gSystem->ResetSignal(kSigFloatingException,kTRUE);
  TString name(entry ? entry : "");
  if(name!="anaFemtoPhiDeuteron" && name!="anaFemtoKaon") { gSystem->Exit(2); return; }
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
  const char* libraries[]={"libStarAnaConfig.so","libStRefMultCorr.so","libStCommon.so","libStFemtoMaker.so"};
  for(unsigned i=0;i<sizeof(libraries)/sizeof(libraries[0]);++i)
    if(gSystem->Load(cwd+"/lib/"+libraries[i])<0) {gSystem->Exit(2);return;}
  gInterpreter->AddIncludePath(cwd);
  gInterpreter->AddIncludePath(cwd+"/include");
  gInterpreter->AddIncludePath(cwd+"/StMaker/common");
  gInterpreter->AddIncludePath("$STAR/StRoot");
  gSystem->AddLinkedLibs(TString::Format("-L%s/lib -lStFemtoMaker -lStCommon -lStarAnaConfig -lStRefMultCorr -Wl,-rpath,%s/lib",cwd.Data(),cwd.Data()));
  TUUID id;
  TString build=TString(gSystem->TempDirectory())+"/star_femto_core_load_"+id.AsString();
  if(gSystem->mkdir(build,kTRUE)!=0) {gSystem->Exit(2);return;}
  if(!gSystem->CompileMacro((cwd+"/analysis/"+name+".C").Data(),"k","",build)) {
    std::cerr << "ERROR: compiling unchanged entry " << name << std::endl;gSystem->Exit(2);return;
  }
  if(gROOT->LoadMacro((cwd+"/analysis/run_"+name+".C").Data())<0) {
    std::cerr << "ERROR: loading unchanged runner " << name << std::endl;gSystem->Exit(2);return;
  }
  TString libs=gSystem->GetLibraries();
  if(libs.Contains("libKFParticle.so") || libs.Contains("libStKfParticleCommon.so")) {
    std::cerr << "ERROR: optional KF library was loaded" << std::endl;gSystem->Exit(2);return;
  }
  std::cout << "PASS unchanged " << name << " ACLiC + runner load; no optional KF dependencies" << std::endl;
  gSystem->Exit(0);
}
