// ROOT5 runner for anaCandidateFemto.C (StCandidateFemtoMaker + reconstruction parts).
#include "TROOT.h"
#include "TSystem.h"
#include "TString.h"
#include "TInterpreter.h"
#include <iostream>

TString CandidateFemtoQuote(const char* value) {
  TString quoted(value ? value : "");
  quoted.ReplaceAll("\\", "\\\\");
  quoted.ReplaceAll("\"", "\\\"");
  quoted.ReplaceAll("\n", "\\n");
  quoted.ReplaceAll("\r", "\\r");
  return TString("\"") + quoted + "\"";
}

// Every argument is required (no default job id or mainconf; docs/ai/AGENT_RULES.md "Batch safety").
void run_anaCandidateFemto(const Char_t* inputFile, const Char_t* outputFile, const Char_t* jobid,
                           Long64_t nEventsMax, const Char_t* configPath) {
  gSystem->ResetSignal(kSigSegmentationViolation, kTRUE);
  gSystem->ResetSignal(kSigBus, kTRUE);
  gSystem->ResetSignal(kSigIllegalInstruction, kTRUE);
  gSystem->ResetSignal(kSigFloatingException, kTRUE);
  const TString pwd = gSystem->WorkingDirectory();
  if (gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C") < 0) {
    gSystem->Exit(1); return;
  }
  Int_t loadError = 0;
  gROOT->ProcessLine("loadSharedLibraries();", &loadError);
  if (loadError) { gSystem->Exit(1); return; }
  // BDT scoring inside StCandidateFemtoMaker needs no extra library (CandidateTreeEnsemble.h, plain C++).
  const char* starLibraries[] = {"StarRoot", "StBichsel", "StPicoEvent", "StPicoDstMaker"};
  for (UInt_t i = 0; i < sizeof(starLibraries) / sizeof(starLibraries[0]); ++i)
    if (gSystem->Load(starLibraries[i]) < 0) {
      std::cerr << "ERROR: loading " << starLibraries[i] << std::endl;
      gSystem->Exit(1); return;
    }
  // KFParticle libraries are needed only by the parts created in anaCandidateFemto.C
  // (KfParticleCandidateSource, OfficialEventSelection), not by the Maker.
  const char* libraries[] = {"libStarAnaConfig.so", "libStRefMultCorr.so", "libKFParticle.so",
      "libStKfParticleCommon.so", "libStCommon.so", "libStCandidateFemtoMaker.so"};
  for (UInt_t i = 0; i < sizeof(libraries) / sizeof(libraries[0]); ++i) {
    const TString path = pwd + "/lib/" + libraries[i];
    if (gSystem->Load(path) < 0) {
      std::cerr << "ERROR: loading " << path << std::endl;
      gSystem->Exit(1); return;
    }
  }
  gInterpreter->AddIncludePath(pwd.Data());
  gInterpreter->AddIncludePath((pwd + "/include").Data());
  gInterpreter->AddIncludePath((pwd + "/StMaker/common").Data());
  gInterpreter->AddIncludePath((pwd + "/StMaker/kfparticle").Data());
  gInterpreter->AddIncludePath((pwd + "/StMaker").Data());
  gInterpreter->AddIncludePath("$STAR/StRoot");
  gSystem->AddLinkedLibs(TString::Format(
      "-L%s/lib -lStCandidateFemtoMaker -lStKfParticleCommon -lKFParticle -lStCommon -lStarAnaConfig "
      "-lStRefMultCorr -Wl,-rpath,%s/lib", pwd.Data(), pwd.Data()));
  const TString buildDir = TString::Format("%s/tmp/kf-aclic/%s-%d", pwd.Data(), gSystem->HostName(), gSystem->GetPid());
  gSystem->mkdir(buildDir, kTRUE);
  gSystem->SetBuildDir(buildDir, kTRUE);
  if (!gSystem->CompileMacro((pwd + "/analysis/anaCandidateFemto.C").Data(), "kf")) {
    std::cerr << "ERROR: compiling anaCandidateFemto.C" << std::endl;
    gSystem->Exit(1); return;
  }
  const TString call = TString::Format("anaCandidateFemto(%s,%s,%s,%lld,%s)",
      CandidateFemtoQuote(inputFile).Data(), CandidateFemtoQuote(outputFile).Data(),
      CandidateFemtoQuote(jobid).Data(), nEventsMax, CandidateFemtoQuote(configPath).Data());
  Int_t error = 0;
  Long_t result = gROOT->ProcessLine(call.Data(), &error);
  // This process's ACLiC build directory is not reused: remove it, or tmp/kf-aclic/ fills up.
  gSystem->Exec(TString::Format("rm -rf %s", buildDir.Data()));
  if (error || result) gSystem->Exit(error ? 1 : static_cast<Int_t>(result));
}
