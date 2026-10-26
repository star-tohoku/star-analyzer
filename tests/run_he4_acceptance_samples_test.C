#include "TROOT.h"
#include "TSystem.h"
#include "TInterpreter.h"
#include "TString.h"
TString ActualNuclearQuote(const char* text){TString v(text?text:"");v.ReplaceAll("\\","\\\\");v.ReplaceAll("\"","\\\"");return "\""+v+"\"";}
void run_he4_acceptance_samples_test(const char* source,const char* output){
  if(!source||!*source||!output||!*output){gSystem->Exit(1);return;}
  TString cwd=gSystem->WorkingDirectory();if(gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C")<0){gSystem->Exit(1);return;}Int_t error=0;gROOT->ProcessLine("loadSharedLibraries();",&error);if(error){gSystem->Exit(1);return;}
  const char* star[]={"StarRoot","StBichsel","StPicoEvent","StPicoDstMaker"};for(int i=0;i<4;++i)if(gSystem->Load(star[i])<0){gSystem->Exit(1);return;}
  const char* local[]={"libStarAnaConfig.so","libStRefMultCorr.so","libStCommon.so"};for(int i=0;i<3;++i)if(gSystem->Load(cwd+"/lib/"+local[i])<0){gSystem->Exit(1);return;}
  gInterpreter->AddIncludePath(cwd);gInterpreter->AddIncludePath(cwd+"/include");gInterpreter->AddIncludePath(cwd+"/StMaker/common");gInterpreter->AddIncludePath(cwd+"/src/third_party/yaml-cpp/include");gInterpreter->AddIncludePath("$STAR/StRoot");gSystem->AddLinkedLibs(TString::Format("-L%s/lib -lStCommon -lStarAnaConfig -lStRefMultCorr -Wl,-rpath,%s/lib",cwd.Data(),cwd.Data()));
  TString build=TString::Format("%s/tmp/he4-acceptance-samples-%d",cwd.Data(),gSystem->GetPid());gSystem->mkdir(build,kTRUE);gSystem->SetBuildDir(build,kTRUE);if(!gSystem->CompileMacro("tests/test_he4_acceptance_samples.C","kf")){gSystem->Exit(1);return;}
  TString call=TString::Format("test_he4_acceptance_samples(%s,%s)",ActualNuclearQuote(source).Data(),ActualNuclearQuote(output).Data());Long_t result=gROOT->ProcessLine(call,&error);if(error||result)gSystem->Exit(1);
}
