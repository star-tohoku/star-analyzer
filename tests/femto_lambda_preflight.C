// Deliberately simple ROOT branch fixtures: tests I/O guards, not KF reconstruction.
#include "analysis/FemtoLambdaRunSupport.h"
#include "TFile.h"
#include "TTree.h"
#include <fstream>
#include <stdexcept>

namespace FemtoLambdaPreflight {
void Require(Bool_t ok, const char* message) {
  if (!ok) throw std::runtime_error(message);
}
void Input(const TString& name, Bool_t covariance) {
  TFile file(name, "CREATE");
  Require(!file.IsZombie(), "fixture output exists or is not writable");
  TTree tree("PicoDst", "I/O-only preflight fixture");
  Int_t value = 0;
  tree.Branch("Event", &value, "Event/I");
  tree.Branch("Track", &value, "Track/I");
  tree.Branch("BTofPidTraits", &value, "BTofPidTraits/I");
  if (covariance) tree.Branch("TrackCovMatrix", &value, "TrackCovMatrix/I");
  tree.Fill(); tree.Write();
}
}

Int_t femto_lambda_preflight(const char* mainconf, const char* directory) {
  try {
    using namespace FemtoLambdaPreflight;
    Require(mainconf && *mainconf && directory && *directory, "explicit arguments required");
    const TString dir = FemtoLambdaRunSupport::Absolute(directory);
    Require(gSystem->AccessPathName(dir), "fixture directory already exists; use a fresh path");
    Require(gSystem->mkdir(dir, kTRUE) == 0, "cannot create fixture directory");
    const TString valid = dir + "/valid.picoDst.root";
    const TString valid2 = dir + "/valid2.picoDst.root";
    const TString bad = dir + "/missing_cov.picoDst.root";
    Input(valid, kTRUE); Input(valid2, kTRUE); Input(bad, kFALSE);
    const TString list = dir + "/second_missing.list", output = dir + "/analysis.root";
    { std::ofstream f(list.Data()); f << valid << "\n" << bad << "\n"; }
    FemtoLambdaRunSupport::Run missingSecond(list, output, "toy", 2, mainconf);
    Require(!missingSecond.Prepare("deuteron"), "second-file covariance omission was accepted");
    FemtoLambdaRunSupport::Run missingEvenOutsideLimit(list, output, "toy", 1, mainconf);
    Require(!missingEvenOutsideLimit.Prepare("deuteron"), "unread input file escaped whole-list preflight");
    FemtoLambdaRunSupport::Run tooMany(valid, output, "toy", 2, mainconf);
    Require(!tooMany.Prepare("deuteron"), "short input was silently clamped");
    FemtoLambdaRunSupport::Run noConfig(valid, output, "toy", 1, "");
    Require(!noConfig.Prepare("deuteron"), "missing mainconf was accepted");
    FemtoLambdaRunSupport::Run existing(valid, valid2, "toy", 1, mainconf);
    Require(!existing.Prepare("deuteron"), "existing output was accepted");
    FemtoLambdaRunSupport::Run wrongSpecies(valid, output, "toy", 1, mainconf);
    Require(!wrongSpecies.Prepare("he4"), "wrong entry species was accepted");
    FemtoLambdaRunSupport::Run good(valid, output, "toy", 1, mainconf);
    Require(good.Prepare("deuteron"), "valid preflight fixture unexpectedly rejected");
    Require(gSystem->AccessPathName(output), "preflight unexpectedly wrote analysis output");
    std::cout << "PASS: 6 failing guards + valid I/O fixture (no physics processing)" << std::endl;
    return 0;
  } catch (const std::exception& e) {
    std::cerr << "FAIL: " << e.what() << std::endl; return 1;
  }
}
