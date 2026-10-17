// run_purityCorrectedCorrelation.C - load libStarAnaConfig, compile
// common/macro/candidate_femto/purityCorrectedCorrelation.C+ in a per-process build directory, and call it.
// Usage: root4star -b -q 'run_purityCorrectedCorrelation.C("file.root","mainconf.yaml","k0s_xi_minus","out/prefix",2,"4,5,6,7,8","")' (all arguments required)
#include "TROOT.h"
#include "TSystem.h"
#include "TString.h"
#include "TInterpreter.h"
#include <iostream>

void run_purityCorrectedCorrelation(const Char_t* rootFile, const Char_t* mainconfPath, const Char_t* channel,
                                    const Char_t* outPrefix, Int_t kstarRebin, const Char_t* cent9List,
                                    const Char_t* options) {
  const TString pwd = gSystem->WorkingDirectory();
  // libEG: TDatabasePDG (particle masses)
  if (gSystem->Load("libEG") < 0) {
    std::cerr << "ERROR: failed to load libEG" << std::endl;
    gSystem->Exit(1); return;
  }
  if (gSystem->Load(pwd + "/lib/libStarAnaConfig.so") < 0) {
    std::cerr << "ERROR: failed to load libStarAnaConfig.so" << std::endl;
    gSystem->Exit(1); return;
  }
  gInterpreter->AddIncludePath(pwd.Data());
  gInterpreter->AddIncludePath((pwd + "/include").Data());
  gInterpreter->AddIncludePath((pwd + "/StMaker/StCandidateFemtoMaker").Data());
  gInterpreter->AddIncludePath((pwd + "/common/macro/candidate_femto").Data());
  gSystem->AddLinkedLibs(TString::Format("-L%s/lib -lStarAnaConfig -Wl,-rpath,%s/lib", pwd.Data(), pwd.Data()));
  const TString buildDir = TString::Format("%s/tmp/aclic-purity/%s-%d", pwd.Data(), gSystem->HostName(), gSystem->GetPid());
  gSystem->mkdir(buildDir, kTRUE);
  gSystem->SetBuildDir(buildDir, kTRUE);
  if (!gSystem->CompileMacro((pwd + "/common/macro/candidate_femto/purityCorrectedCorrelation.C").Data(), "kf")) {
    std::cerr << "ERROR: compiling purityCorrectedCorrelation.C" << std::endl;
    gSystem->Exit(1); return;
  }
  Int_t error = 0;
  const Long_t result = gROOT->ProcessLine(TString::Format("purityCorrectedCorrelation(\"%s\",\"%s\",\"%s\",\"%s\",%d,\"%s\",\"%s\")", rootFile,
                                     mainconfPath, channel, outPrefix, kstarRebin, cent9List, options), &error);
  gSystem->Exec(TString::Format("rm -rf %s", buildDir.Data()));
  if (error || result != 0) {
    std::cerr << "ERROR: purityCorrectedCorrelation failed (status " << result << ")" << std::endl;
    gSystem->Exit(1);
  }
}
