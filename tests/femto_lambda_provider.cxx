// ROOT5-linked deterministic provider/selector tests. Called by the integration
// test runner after loading the same libraries as the Lambda analysis.
#include "FemtoLambdaKfProvider.h"
#include "KfLambdaSelector.h"
#include "ConfigManager.h"
#include "cuts/KfParticleCutConfig.h"
#include "StPicoEvent/StPicoEvent.h"

#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>

namespace {
void Check(bool condition, const char* message) {
  if (!condition) throw std::runtime_error(message);
}
// Independent frozen pre-refactor predicate. Do not make this call the shared
// helper: its purpose is boundary/disabled-value compatibility verification.
Bool_t OriginalFinalCuts(const KfLambdaCandidate& c, const KfParticleCutConfig& k) {
  if (k.selectionProfile == "lambda_imp5" &&
      (!c.imp5PathValid || !std::isfinite(c.protonHelixPathLength) ||
       !std::isfinite(c.pionHelixPathLength) ||
       std::fabs(c.protonHelixPathLength) > k.imp5MaxPathLength ||
       std::fabs(c.pionHelixPathLength) > k.imp5MaxPathLength)) return kFALSE;
  if (c.mass < k.minMass || c.mass > k.maxMass) return kFALSE;
  if (k.maxMassError >= 0. && c.massError > k.maxMassError) return kFALSE;
  if (k.maxChi2Ndf >= 0. && c.chi2Ndf > k.maxChi2Ndf) return kFALSE;
  if (k.maxTopoChi2Ndf >= 0. && c.topoChi2Ndf > k.maxTopoChi2Ndf) return kFALSE;
  if (k.maxDaughterDistance >= 0. && c.daughterDistance > k.maxDaughterDistance) return kFALSE;
  if (k.maxDistanceToPv >= 0. && c.distanceToPv > k.maxDistanceToPv) return kFALSE;
  if (k.minDecayLength >= 0. && c.decayLength < k.minDecayLength) return kFALSE;
  if (k.minDecayLengthSignificance >= 0. && c.decayLengthSignificance < k.minDecayLengthSignificance) return kFALSE;
  if (k.minVertexLineSignificance >= 0. && c.vertexLineLengthSignificance < k.minVertexLineSignificance) return kFALSE;
  if (k.minCosPointing > -1. && c.cosPointing < k.minCosPointing) return kFALSE;
  return kTRUE;
}

KfLambdaCandidate Example() {
  KfLambdaCandidate c;
  c.mass = 1.125f; c.massError = .25f; c.chi2Ndf = .25f; c.topoChi2Ndf = .25f;
  c.daughterDistance = .25f; c.distanceToPv = .25f;
  c.decayLength = 1.f; c.decayLengthSignificance = 1.f; c.vertexLineLengthSignificance = 1.f;
  c.cosPointing = 1.f; c.imp5PathValid = true;
  c.protonHelixPathLength = 10.; c.pionHelixPathLength = -10.;
  c.protonId = 10003; c.pionId = 42; c.protonIndex = 0; c.pionIndex = 7; c.pdg = 3122;
  c.px = .75f; c.py = -.5f; c.pz = 1.5f;
  c.x = 1.f; c.y = 2.f; c.z = 3.f;
  return c;
}

void SelectorTests() {
  KfParticleCutConfig& k = KfParticleCutConfig::GetInstance();
  k.SetDefaults();
  k.selectionProfile = "lambda_imp5";
  k.imp5MaxPathLength = 10.;
  k.minMass = 1.; k.maxMass = 1.25;
  KfLambdaCandidate c = Example();
  Check(PassKfLambdaCandidateCuts(c, k), "baseline candidate rejected");
  struct Bound {
    Float_t KfLambdaCandidate::* value;
    Double_t KfParticleCutConfig::* cut;
    bool minimum;
  };
  const Bound bounds[] = {
    {&KfLambdaCandidate::massError, &KfParticleCutConfig::maxMassError, false},
    {&KfLambdaCandidate::chi2Ndf, &KfParticleCutConfig::maxChi2Ndf, false},
    {&KfLambdaCandidate::topoChi2Ndf, &KfParticleCutConfig::maxTopoChi2Ndf, false},
    {&KfLambdaCandidate::daughterDistance, &KfParticleCutConfig::maxDaughterDistance, false},
    {&KfLambdaCandidate::distanceToPv, &KfParticleCutConfig::maxDistanceToPv, false},
    {&KfLambdaCandidate::decayLength, &KfParticleCutConfig::minDecayLength, true},
    {&KfLambdaCandidate::decayLengthSignificance, &KfParticleCutConfig::minDecayLengthSignificance, true},
    {&KfLambdaCandidate::vertexLineLengthSignificance, &KfParticleCutConfig::minVertexLineSignificance, true},
    {&KfLambdaCandidate::cosPointing, &KfParticleCutConfig::minCosPointing, true}
  };
  for (size_t i = 0; i < sizeof(bounds)/sizeof(bounds[0]); ++i) {
    const Bound& b = bounds[i];
    const double saved = k.*(b.cut);
    k.*(b.cut) = .5;
    c.*(b.value) = .5f;
    Check(PassKfLambdaCandidateCuts(c,k), "inclusive final-cut boundary rejected");
    c.*(b.value) = b.minimum ? .25f : .75f;
    Check(!PassKfLambdaCandidateCuts(c,k), "outside final-cut boundary accepted");
    k.*(b.cut) = -1.;
    Check(PassKfLambdaCandidateCuts(c,k), "disabled final cut rejected candidate");
    k.*(b.cut) = saved;
    c = Example();
  }
  c.mass = 1.f; Check(PassKfLambdaCandidateCuts(c,k), "lower mass boundary rejected");
  c.mass = 1.25f; Check(PassKfLambdaCandidateCuts(c,k), "upper mass boundary rejected");
  c.mass = 1.5f; Check(!PassKfLambdaCandidateCuts(c,k), "out-of-mass accepted");
  c = Example(); c.imp5PathValid = false;
  Check(!PassKfLambdaCandidateCuts(c,k), "invalid Imp5 path accepted");
  c = Example(); c.pionHelixPathLength = std::numeric_limits<double>::quiet_NaN();
  Check(!PassKfLambdaCandidateCuts(c,k), "NaN path accepted");
  c = Example(); c.protonHelixPathLength = 10.0001;
  Check(!PassKfLambdaCandidateCuts(c,k), "above path bound accepted");
  // Cross-products exercise interaction and profile disabling, not just each
  // isolated predicate. Non-finite fits are adapter-level exclusions.
  for (int i = 0; i < 4096; ++i) {
    c = Example();
    k.selectionProfile = i%2 ? "lambda_imp5" : "kf_reference";
    c.imp5PathValid = i%3 != 0;
    c.mass = .9f + (i%7)*.0625f;
    c.protonHelixPathLength = (i%5)*5.;
    c.pionHelixPathLength = -(i%4)*5.;
    for (size_t b = 0; b < sizeof(bounds)/sizeof(bounds[0]); ++b) {
      k.*(bounds[b].cut) = (i >> (b%7))%3 == 0 ? -1. : .5;
      c.*(bounds[b].value) = ((i >> (b%9))%5)*.25f;
    }
    Check(PassKfLambdaCandidateCuts(c,k) == OriginalFinalCuts(c,k),
          "shared selector differs from original frozen predicate");
  }
}

void MappingTests() {
  KfLambdaCandidate raw = Example();
  FemtoLambdaCandidate output;
  std::string error;
  Check(ConvertKfLambdaForFemto(raw, 123, 1.115683, output, error), "conversion failed");
  Check(output.candidate.reso.dau1Index == 0 && output.candidate.reso.dau2Index == 7 &&
        output.candidate.reso.dau1EventIndex == 123 &&
        output.candidate.reso.dau2EventIndex == 123, "ID substituted for array index");
  Check(output.protonId == 10003 && output.pionId == 42 &&
        output.mass == raw.mass && output.candidate.reso.invMass == raw.mass &&
        output.candidate.px == raw.px && output.candidate.py == raw.py &&
        output.candidate.pz == raw.pz, "raw diagnostics/momentum changed");
  const float energy = output.candidate.energy;
  raw.mass = 1.20f;
  Check(ConvertKfLambdaForFemto(raw,123,1.115683,output,error) &&
        output.candidate.energy == energy && output.candidate.reso.invMass == raw.mass,
        "measured mass leaked into fixed-pair-mass energy");
  Check(std::fabs(output.candidate.P4().M() - 1.115683) < 1.e-6,
        "fixed pair mass not preserved to candidate precision");
  FemtoCandidate track;
  track.eventIndex = 123; track.trk.trackIndex = 7;
  Check(FemtoCandidatesShareTrack(output.candidate,track), "same-event shared daughter not detected");
  track.eventIndex = 124;
  Check(!FemtoCandidatesShareTrack(output.candidate,track), "cross-event local index treated as same track");
  raw.pionIndex = raw.protonIndex;
  Check(!ConvertKfLambdaForFemto(raw,123,1.115683,output,error), "identical daughter indices accepted");
  raw = Example(); raw.px = std::numeric_limits<float>::infinity();
  Check(!ConvertKfLambdaForFemto(raw,123,1.115683,output,error), "nonfinite momentum accepted");
}

void LifecycleTests(const char* mainconf) {
  std::auto_ptr<FemtoLambdaProvider> provider(createFemtoLambdaKfProvider());
  Check(!provider->Process(0,0) && !provider->LastError().empty(),
        "uninitialized provider processed event");
  std::ostringstream errors;
  Check(!provider->Init("",errors) && !provider->LastError().empty(),
        "missing mainconf accepted");
  StPicoEvent event;
  Check(!provider->AcceptEvent(event,0), "failed provider accepts events");
  if (mainconf && mainconf[0]) {
    Check(ConfigManager::GetInstance().LoadConfig(mainconf), "mainconf load failed");
    Check(provider->Init(mainconf,errors), "valid provider initialization failed");
    Check(!provider->Process(0,0) && !provider->LastError().empty(),
          "initialized provider accepts missing PicoDst");
    provider->Clear();
    Check(provider->Candidates().empty() && provider->Stats().selectedCandidates == 0 &&
          provider->LastError().empty(), "Clear leaves stale per-event state");
    Check(!provider->Init("",errors) && !provider->AcceptEvent(event,0),
          "failed reload retained active provider");
  }
}
}

extern "C" int run_femto_lambda_provider_tests(const char* mainconf) {
  try {
    SelectorTests();
    MappingTests();
    LifecycleTests(mainconf);
    std::cout << "PASS: Femto Lambda selector (4096 closure cases), mapping, fixed mass, lifecycle" << std::endl;
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "FAIL: Femto Lambda provider: " << error.what() << std::endl;
    return 1;
  }
}
