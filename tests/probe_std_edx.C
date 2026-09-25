// Diagnostic only: mode 0 is STAR + StBichsel, without project KF libraries.
// Mode 1 additionally loads the project KF core; mode 2 loads the source row
// header first; mode 3 tests a simple interpreted spline-row field assignment.
// Mode 4 is STAR-only without any metadata lookup; mode 5 also loads the KF core.
// WARNING: mode 6 intentionally reproduces a segmentation fault (expected exit 139).
// Modes 6/7 isolate const declaration initialization versus separate assignment
// around the interpreter call which recursively interprets the spline table.
#include "TROOT.h"
#include "TSystem.h"
#include "TString.h"
#include "TClass.h"
#include "TDataMember.h"
#include <iostream>

void probe_std_edx(Int_t mode = 0) {
  gSystem->ResetSignal(kSigSegmentationViolation, kTRUE);
  gSystem->ResetSignal(kSigBus, kTRUE);
  gSystem->ResetSignal(kSigIllegalInstruction, kTRUE);
  gSystem->ResetSignal(kSigFloatingException, kTRUE);
  std::cout << "[dEdx probe] mode=" << mode << " ROOT=" << gROOT->GetVersion()
            << " STAR=" << gSystem->Getenv("STAR") << std::endl;
  if (gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C") < 0) {
    gSystem->Exit(1); return;
  }
  Int_t error = 0;
  gROOT->ProcessLine("loadSharedLibraries();", &error);
  if (error || gSystem->Load("StBichsel") < 0) { gSystem->Exit(1); return; }
  if (mode == 1 || mode == 5) {
    if (gSystem->Load(TString(gSystem->WorkingDirectory()) + "/lib/libKFParticle.so") < 0) {
      gSystem->Exit(1); return;
    }
  }
  if (mode < 4) {
  TClass* rowClass = gROOT->GetClass("spline3_st");
  TClass* tableClass = gROOT->GetClass("St_spline3");
  std::cout << "[dEdx probe] rowClass=" << rowClass << " tableClass=" << tableClass << std::endl;
  if (rowClass) {
    TDataMember* member = rowClass->GetDataMember("nknots");
    if (member) std::cout << "[dEdx probe] nknots type=" << member->GetFullTypeName()
                          << " property=" << member->Property()
                          << " offset=" << member->GetOffset() << std::endl;
  }
  }
  if (mode == 2) {
    gROOT->ProcessLine(".L $STAR/StRoot/StBichsel/spline3.h", &error);
    std::cout << "[dEdx probe] source row header load status=" << error << std::endl;
    if (error) { gSystem->Exit(1); return; }
  }
  if (mode == 3) {
    gROOT->ProcessLine("spline3_st star_kf_probe_row;", &error);
    std::cout << "[dEdx probe] row construction status=" << error << std::endl;
    if (error) { gSystem->Exit(1); return; }
    gROOT->ProcessLine("star_kf_probe_row.nknots = 14;", &error);
    std::cout << "[dEdx probe] row assignment status=" << error << std::endl;
    if (error) { gSystem->Exit(1); return; }
  }
  std::cout << "[dEdx probe] calling StdEdxModel::instance()" << std::endl;
  Long_t model = 0;
  if (mode == 6) {
    const Long_t nestedResult = gROOT->ProcessLine("(Long_t)StdEdxModel::instance();", &error);
    model = nestedResult;
  } else if (mode == 7) {
    Long_t nestedResult = 0;
    nestedResult = gROOT->ProcessLine("(Long_t)StdEdxModel::instance();", &error);
    model = nestedResult;
  } else {
    model = gROOT->ProcessLine("(Long_t)StdEdxModel::instance();", &error);
  }
  std::cout << "[dEdx probe] model=" << model << " status=" << error << std::endl;
  if (error || !model) { gSystem->Exit(1); return; }
  gROOT->ProcessLine("std::cout << StdEdxPull::EvalPred(4., 1, 1) << std::endl;", &error);
  std::cout << "[dEdx probe] prediction status=" << error << std::endl;
  gSystem->Exit(error ? 1 : 0);
}
