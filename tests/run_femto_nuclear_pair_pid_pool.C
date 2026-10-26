#include "TROOT.h"
#include "TSystem.h"
#include "TInterpreter.h"
#include "TString.h"
TString PairPidPoolQuote(const char* value){TString s(value?value:"");s.ReplaceAll("\\","\\\\");s.ReplaceAll("\"","\\\"");return "\""+s+"\"";}
void run_femto_nuclear_pair_pid_pool(const char* on,const char* off,const char* output){
  gInterpreter->AddIncludePath(TString(gSystem->WorkingDirectory())+"/src/third_party/yaml-cpp/include");
  if(!on||!*on||!off||!*off||!output||!*output){gSystem->Exit(2);return;}
  if(gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C")<0){gSystem->Exit(2);return;}
  Int_t error=0;gROOT->ProcessLine("loadSharedLibraries();",&error);if(error){gSystem->Exit(2);return;}
  const char* star[]={"StarRoot","StBichsel","StPicoEvent","StPicoDstMaker"};for(unsigned i=0;i<4;++i)if(gSystem->Load(star[i])<0){gSystem->Exit(2);return;}
  const TString cwd=gSystem->WorkingDirectory();const char* libs[]={"libStarAnaConfig.so","libStRefMultCorr.so","libStCommon.so","libStFemtoMaker.so"};for(unsigned i=0;i<4;++i)if(gSystem->Load(cwd+"/lib/"+libs[i])<0){gSystem->Exit(2);return;}
  gInterpreter->AddIncludePath(cwd);gInterpreter->AddIncludePath(cwd+"/include");gInterpreter->AddIncludePath(cwd+"/StMaker/common");gInterpreter->AddIncludePath("$STAR/StRoot");gSystem->AddLinkedLibs(TString::Format("-L%s/lib -lStFemtoMaker -lStCommon -lStarAnaConfig -lStRefMultCorr -Wl,-rpath,%s/lib",cwd.Data(),cwd.Data()));
  const TString build=TString::Format("%s/tmp/pair-pid-pool-%d",cwd.Data(),gSystem->GetPid());gSystem->mkdir(build,kTRUE);gSystem->SetBuildDir(build,kTRUE);if(!gSystem->CompileMacro("tests/femto_nuclear_pair_pid_pool.C","kf")){gSystem->Exit(2);return;}
  const TString call=TString::Format("femto_nuclear_pair_pid_pool(%s,%s,%s)",PairPidPoolQuote(on).Data(),PairPidPoolQuote(off).Data(),PairPidPoolQuote(output).Data());const Long_t result=gROOT->ProcessLine(call,&error);gSystem->Exit(error?2:int(result));
}
