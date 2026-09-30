#include "TROOT.h"
#include "TSystem.h"
#include "TInterpreter.h"
#include "TString.h"
#include <iostream>
#include <cstring>

TString LambdaLegacyTestQuote(const char* input) {
 TString q(input?input:"");q.ReplaceAll("\\","\\\\");q.ReplaceAll("\"","\\\"");return "\""+q+"\"";
}
void run_femto_lambda_legacy_test(const char* mainconf,const char* input="",int events=1000) {
 TString cwd=gSystem->WorkingDirectory();
 if(gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C")<0){gSystem->Exit(1);return;}
 Int_t error=0;gROOT->ProcessLine("loadSharedLibraries();",&error);
 if(error){gSystem->Exit(1);return;}
 const char* star[]={"StarRoot","StBichsel","StPicoEvent","StPicoDstMaker"};
 for(int i=0;i<4;++i)if(gSystem->Load(star[i])<0){gSystem->Exit(1);return;}
 const char* libs[]={"libStarAnaConfig.so","libStRefMultCorr.so","libStCommon.so","libStNuclearIdMaker.so"};
 for(int i=0;i<4;++i)if(gSystem->Load(cwd+"/lib/"+libs[i])<0){gSystem->Exit(1);return;}
 gInterpreter->AddIncludePath(cwd);gInterpreter->AddIncludePath(cwd+"/include");
 gInterpreter->AddIncludePath(cwd+"/StMaker/common");gInterpreter->AddIncludePath("$STAR/StRoot");
 gSystem->AddLinkedLibs(TString::Format("-L%s/lib -lStNuclearIdMaker -lStCommon -lStarAnaConfig -lStRefMultCorr -Wl,-rpath,%s/lib",cwd.Data(),cwd.Data()));
 TString build=TString::Format("%s/tmp/femto-legacy-test-%d",cwd.Data(),gSystem->GetPid());
 gSystem->mkdir(build,kTRUE);gSystem->SetBuildDir(build,kTRUE);
 Bool_t isData=kFALSE;
 if(input && strlen(input)>0) isData=kTRUE;
 std::cout << "[Legacy test] input=" << input << " data=" << isData << std::endl;
 TString path("tests/femto_lambda_legacy.C");
 if(isData) path="tests/compare_femto_lambda_nuclear.C";
 if(!gSystem->CompileMacro(path,"kf")){gSystem->Exit(1);return;}
 TString call;
 if(isData) call=TString::Format("compare_femto_lambda_nuclear(%s,%s,%d)",LambdaLegacyTestQuote(input).Data(),LambdaLegacyTestQuote(mainconf).Data(),events);
 else call=TString::Format("femto_lambda_legacy(%s)",LambdaLegacyTestQuote(mainconf).Data());
 Long_t result=gROOT->ProcessLine(call,&error);
 if(error||result)gSystem->Exit(1);
}

