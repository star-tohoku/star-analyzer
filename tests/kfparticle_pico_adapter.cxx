// Integration fixture, not a generated PicoDst file or a physics-validation sample.
// Uses the actual SL24y Pico classes, ROOT TClonesArrays and full KF adapter.
#include "StPicoKFParticleInterface.h"
#include "KfParticleHelper.h"
#include "kfparticle_event_selection.h"
#include "cuts/KfParticleCutConfig.h"
#include "KFPTrack.h"
#include "StEvent/StDcaGeometry.h"
#include "StBichsel/StdEdxPull.h"
#include "StPicoEvent/StPicoArrays.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoTrackCovMatrix.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"
#include "TClonesArray.h"
#include "TSystem.h"
#include "TVector3.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iterator>
#include <unistd.h>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>

namespace kfp = star_analyzer_kfp;
namespace {
const double kFixtureMass = 1.122; // Deliberately different from the PDG Lambda mass.
const double kProtonMass = 0.938272;
const double kPionMass = 0.139570;

void Require(bool condition, const std::string& message) {
  if (!condition) throw std::runtime_error(message);
}

void Near(double actual, double expected, double tolerance, const char* quantity) {
  if (!std::isfinite(actual) || std::fabs(actual - expected) > tolerance) {
    std::ostringstream error;
    error << quantity << ": " << actual << " != " << expected << " +/- " << tolerance;
    throw std::runtime_error(error.str());
  }
}


class TemporaryConfig {
public:
  explicit TemporaryConfig(const std::string& text) {
    char pattern[] = "/tmp/star_analyzer_kf_config_XXXXXX";
    const int descriptor = mkstemp(pattern);
    Require(descriptor >= 0, "cannot create temporary KF configuration fixture");
    path = pattern;
    close(descriptor);
    std::ofstream output(path.c_str());
    output << text;
    output.close();
    if (!output) {
      unlink(path.c_str());
      throw std::runtime_error("cannot write temporary KF configuration fixture");
    }
  }
  ~TemporaryConfig() { unlink(path.c_str()); }
  std::string path;
private:
  TemporaryConfig(const TemporaryConfig&);
  TemporaryConfig& operator=(const TemporaryConfig&);
};

std::string ReplaceScalar(const std::string& yaml, const std::string& key,
                          const std::string& value) {
  const std::string marker = key + ":";
  size_t begin = yaml.find(marker);
  Require(begin != std::string::npos, "required test configuration key missing: " + key);
  const size_t end = yaml.find('\n', begin);
  std::string result = yaml;
  result.replace(begin, (end == std::string::npos ? yaml.size() : end) - begin,
                 marker + " " + value);
  return result;
}

// Preserve fixture YAML comments and existing keys; new optional profile keys
// may be absent from a historical baseline configuration.
std::string SetScalar(const std::string& yaml, const std::string& key,
                      const std::string& value) {
  const std::string marker = key + ":";
  size_t begin = yaml.find(marker);
  while (begin != std::string::npos && begin != 0 && yaml[begin - 1] != '\n')
    begin = yaml.find(marker, begin + marker.size());
  if (begin == std::string::npos) return yaml + "\n" + marker + " " + value + "\n";
  const size_t end = yaml.find('\n', begin);
  std::string result(yaml);
  result.replace(begin, (end == std::string::npos ? yaml.size() : end) - begin,
                 marker + " " + value);
  return result;
}

// Test input uses flat, one-line scalars/sequences, like the baseline KF YAML.
std::string RemoveScalar(const std::string& yaml, const std::string& key) {
  std::istringstream input(yaml);
  std::ostringstream output;
  std::string line;
  const std::string marker = key + ":";
  while (std::getline(input, line))
    if (line.compare(0, marker.size(), marker) != 0) output << line << '\n';
  return output.str();
}

void TestConfiguration(const char* cutsPath, KfParticleCutConfig& cuts) {
  std::ifstream input(cutsPath);
  Require(static_cast<bool>(input), "cannot read KF cut YAML for configuration tests");
  const std::string yaml((std::istreambuf_iterator<char>(input)), std::istreambuf_iterator<char>());
  {
    TemporaryConfig invalid(yaml + "\nunknownKfTestKey: 1\n");
    Require(!cuts.LoadFromFile(invalid.path.c_str()), "unknown KF configuration key accepted");
  }
  {
    TemporaryConfig invalid(yaml + "\npidProfile: pico_nsigma\n");
    Require(!cuts.LoadFromFile(invalid.path.c_str()), "duplicate KF configuration key accepted");
  }
  {
    TemporaryConfig invalid(ReplaceScalar(yaml, "schemaVersion", "1"));
    Require(!cuts.LoadFromFile(invalid.path.c_str()), "obsolete scalar schema accepted");
  }
  {
    TemporaryConfig invalid(ReplaceScalar(yaml, "finderLCut", "not_a_number"));
    Require(!cuts.LoadFromFile(invalid.path.c_str()), "nonnumeric Finder cut accepted");
  }
  {
    TemporaryConfig invalid(yaml + "\nminDCAProton: 0.7\n");
    Require(!cuts.LoadFromFile(invalid.path.c_str()), "obsolete scalar DCA key accepted");
  }
  Require(cuts.LoadFromFile(cutsPath), "cannot reload valid KF configuration");
  cuts.useTof = kTRUE;
  cuts.tofMean.resize(2);
  std::ostringstream errors;
  Require(!cuts.Validate(errors), "short TOF coefficient vector accepted");
  {
    StPicoKFParticleInterface invalid(cuts);
    Require(!invalid.ProcessEvent(0) && !invalid.LastError().empty(),
            "adapter failed to report invalid configuration safely");
  }
  Require(cuts.LoadFromFile(cutsPath), "cannot restore valid configuration after negative tests");
  Require(cuts.selectionProfile == "kf_reference" && cuts.imp5MinDCAProton < 0. &&
          cuts.imp5MinDCAPion < 0. && cuts.imp5MaxPathLength < 0.,
          "baseline YAML did not select the backwards-compatible disabled Imp5 defaults");
  std::string imp5 = SetScalar(yaml, "selectionProfile", "lambda_imp5");
  imp5 = SetScalar(imp5, "pidProfile", "pico_nsigma");
  imp5 = SetScalar(imp5, "useTof", "false");
  imp5 = SetScalar(imp5, "strictTofPid", "false");
  imp5 = SetScalar(imp5, "cleanKaonsWithTof", "false");
  imp5 = SetScalar(imp5, "useHftTracksOnly", "false");
  imp5 = SetScalar(imp5, "imp5MinDCAProton", "0");
  imp5 = SetScalar(imp5, "imp5MinDCAPion", "0");
  imp5 = SetScalar(imp5, "imp5MaxPathLength", "1000");
  const char* inactiveImp5Keys[] = {"tofCalibrationProfile", "minNHitsDedx", "minPt", "maxPt",
      "minEta", "maxEta", "minDedxError", "maxDedxError", "nSigmaKaon",
      "tofNSigma", "tofPMax", "tofKaonPMin", "tofKaonPMax", "tofMean", "tofSigma"};
  std::string denseImp5Dump;
  {
    TemporaryConfig valid(imp5);
    Require(cuts.LoadFromFile(valid.path.c_str()), "explicit lambda_imp5 YAML was rejected");
    Require(cuts.selectionProfile == "lambda_imp5" && cuts.pidProfile == "pico_nsigma" &&
            cuts.imp5MinDCAProton == 0. && cuts.imp5MinDCAPion == 0. &&
            cuts.imp5MaxPathLength == 1000., "Imp5 YAML keys did not reach the effective config");
    std::ostringstream dump; cuts.Dump(dump);
    Require(dump.str().find("lambda_imp5") != std::string::npos &&
            dump.str().find("imp5MaxPathLength") != std::string::npos,
            "Imp5 effective config dump omitted the opt-in profile or path guard");
    denseImp5Dump = dump.str();
    for (size_t i = 0; i < sizeof(inactiveImp5Keys) / sizeof(inactiveImp5Keys[0]); ++i)
      Require(denseImp5Dump.find(std::string(inactiveImp5Keys[i]) + ":") == std::string::npos,
              std::string("Imp5 dump advertises an unused threshold: ") + inactiveImp5Keys[i]);
  }
  std::string sparseImp5(imp5);
  for (size_t i = 0; i < sizeof(inactiveImp5Keys) / sizeof(inactiveImp5Keys[0]); ++i)
    sparseImp5 = RemoveScalar(sparseImp5, inactiveImp5Keys[i]);
  {
    TemporaryConfig valid(sparseImp5);
    Require(cuts.LoadFromFile(valid.path.c_str()), "sparse Imp5 YAML requires inactive quality/kaon/TOF fields");
    Require(cuts.tofMean.empty() && cuts.tofSigma.empty(),
            "sparse Imp5 reload inherited unused TOF coefficients from the previous YAML");
    std::ostringstream sparseDump; cuts.Dump(sparseDump);
    Require(sparseDump.str() == denseImp5Dump,
            "removing only inactive Imp5 fields changed the effective active configuration");
    const char* retainedKeys[] = {"selectionProfile", "imp5MinDCAProton", "imp5MinDCAPion",
        "imp5MaxPathLength", "minNHitsFit", "minNHitsRatio", "useHftTracksOnly", "nSigmaPion",
        "nSigmaProton", "useTof", "strictTofPid", "cleanKaonsWithTof", "finderLCut",
        "maxMassError", "maxChi2Ndf", "maxTopoChi2Ndf", "minDecayLength",
        "minDecayLengthSignificance", "minVertexLineSignificance", "minCosPointing"};
    for (size_t i = 0; i < sizeof(retainedKeys) / sizeof(retainedKeys[0]); ++i)
      Require(sparseDump.str().find(std::string(retainedKeys[i]) + ":") != std::string::npos,
              std::string("sparse Imp5 dump lost an active switch/cut or explicit disabled final cut: ") + retainedKeys[i]);
  }
  {
    TemporaryConfig invalid(RemoveScalar(sparseImp5, "useTof"));
    Require(!cuts.LoadFromFile(invalid.path.c_str()),
            "sparse Imp5 silently lost required explicit TOF disabling");
  }
  const char* invalidKeys[] = {"selectionProfile", "pidProfile", "useTof", "strictTofPid",
      "cleanKaonsWithTof", "useHftTracksOnly", "imp5MinDCAProton", "imp5MinDCAPion",
      "imp5MaxPathLength", "imp5MaxPathLength", "imp5MinDCAProton", "selectionProfile"};
  const char* invalidValues[] = {"unknown", "dedx_pull", "true", "true", "true", "true",
      "-1", "-1", "0", ".inf", ".nan", "kf_reference"};
  for (size_t i = 0; i < sizeof(invalidKeys) / sizeof(invalidKeys[0]); ++i) {
    TemporaryConfig invalid(SetScalar(imp5, invalidKeys[i], invalidValues[i]));
    Require(!cuts.LoadFromFile(invalid.path.c_str()),
            std::string("incompatible Imp5 configuration accepted: ") + invalidKeys[i] + "=" + invalidValues[i]);
  }
  Require(cuts.LoadFromFile(cutsPath), "cannot restore baseline after Imp5 config tests");
  Require(cuts.selectionProfile == "kf_reference" && cuts.imp5MinDCAProton < 0. &&
          cuts.imp5MinDCAPion < 0. && cuts.imp5MaxPathLength < 0.,
          "loading a baseline YAML retained an earlier Imp5 selection");
  std::ostringstream referenceDump; cuts.Dump(referenceDump);
  Require(referenceDump.str().find("minPt:") != std::string::npos &&
          referenceDump.str().find("nSigmaKaon:") != std::string::npos &&
          referenceDump.str().find("tofMean:") != std::string::npos,
          "Imp5 suppression leaked into the active reference-profile dump");
  Require(referenceDump.str().find("imp5MinDCAProton:") == std::string::npos &&
          referenceDump.str().find("imp5MinDCAPion:") == std::string::npos &&
          referenceDump.str().find("imp5MaxPathLength:") == std::string::npos,
          "reference dump advertises disabled Imp5-only thresholds");
  if (!cuts.cleanKaonsWithTof)
    Require(referenceDump.str().find("tofKaonPMin:") == std::string::npos &&
            referenceDump.str().find("tofKaonPMax:") == std::string::npos,
            "reference dump advertises an inactive kaon-cleaning range");
  cuts.useTof = kFALSE;
  cuts.cleanKaonsWithTof = kFALSE;
  std::ostringstream noTofDump; cuts.Dump(noTofDump);
  const char* inactiveTofKeys[] = {"tofCalibrationProfile", "tofNSigma", "tofPMax", "tofMean",
                                  "tofSigma", "tofKaonPMin", "tofKaonPMax"};
  for (size_t i = 0; i < sizeof(inactiveTofKeys) / sizeof(inactiveTofKeys[0]); ++i)
    Require(noTofDump.str().find(std::string(inactiveTofKeys[i]) + ":") == std::string::npos,
            std::string("TOF-disabled reference dump advertises an unused field: ") + inactiveTofKeys[i]);
  Require(noTofDump.str().find("useTof:") != std::string::npos &&
          noTofDump.str().find("minPt:") != std::string::npos,
          "TOF-disabled reference dump lost an explicit switch or active track cut");
  // Existing reference kaon-cleaning logic uses this range even without a TOF
  // measurement; its active range must not be omitted just because useTof=false.
  cuts.cleanKaonsWithTof = kTRUE;
  std::ostringstream cleanKaonDump; cuts.Dump(cleanKaonDump);
  Require(cleanKaonDump.str().find("tofKaonPMin:") != std::string::npos &&
          cleanKaonDump.str().find("tofKaonPMax:") != std::string::npos &&
          cleanKaonDump.str().find("tofMean:") == std::string::npos,
          "kaon-cleaning range was not retained independently of TOF-polynomial usage");
  Require(cuts.LoadFromFile(cutsPath), "cannot restore reference after conditional-dump tests");
  std::cout << "PASS strict config: unknown/duplicate/obsolete/schema/type/TOF dimensions, "
               "explicit/sparse Imp5 compatibility, inactive-field dump omission and baseline restoration" << std::endl;
}

void ConfigureFixture(KfParticleCutConfig& cuts) {
  // Explicit test-only settings; production YAML is neither changed nor saved.
  cuts.selectionProfile = "kf_reference";
  cuts.imp5MinDCAProton = -1.;
  cuts.imp5MinDCAPion = -1.;
  cuts.imp5MaxPathLength = -1.;
  cuts.pidProfile = "pico_nsigma";
  cuts.useTof = kFALSE;
  cuts.strictTofPid = kFALSE;
  cuts.cleanKaonsWithTof = kFALSE;
  cuts.useHftTracksOnly = kFALSE;
  cuts.minNHitsFit = 15;
  cuts.minNHitsDedx = 5;
  cuts.minNHitsRatio = 0.52;
  cuts.minPt = 0.15;
  cuts.maxPt = 10.;
  cuts.minEta = -1.;
  cuts.maxEta = 1.;
  cuts.minDedxError = 0.04;
  cuts.maxDedxError = 0.12;
  cuts.nSigmaPion = 3.;
  cuts.nSigmaKaon = 2.;
  cuts.nSigmaProton = 3.;
  cuts.interfaceChiPrimaryCut = 3.;
  cuts.primaryProbCut = 0.0001;
  cuts.finderMaxDaughterDistance = 1.5;
  cuts.finderLCut = 1.;
  cuts.finderChiPrimary2D = 3.;
  cuts.finderChi2Ndf2D = 10.;
  cuts.finderLdL2D = 3.;
  cuts.reconstructAntiLambda = kTRUE;
  cuts.rejectBadCovariance = kTRUE;
  cuts.maxPositionVariance = 100.;
  cuts.maxMomentumVariance = 1.;
  Require(cuts.Validate(std::cerr), "synthetic test configuration is invalid");
}

struct Fixture {
  Fixture() {
    for (int i = 0; i < StPicoArrays::NAllPicoArrays; ++i) arrays[i] = 0;
    arrays[StPicoArrays::Event] = new TClonesArray("StPicoEvent", 1);
    arrays[StPicoArrays::Track] = new TClonesArray("StPicoTrack", 8);
    arrays[StPicoArrays::TrackCovMatrix] = new TClonesArray("StPicoTrackCovMatrix", 8);
    arrays[StPicoArrays::BTofPidTraits] = new TClonesArray("StPicoBTofPidTraits", 8);
    StPicoDst::set(arrays);
    Reset();
  }
  ~Fixture() {
    StPicoDst::unset();
    for (int i = 0; i < StPicoArrays::NAllPicoArrays; ++i) delete arrays[i];
  }
  void Reset() {
    for (int i = 0; i < StPicoArrays::NAllPicoArrays; ++i)
      if (arrays[i]) arrays[i]->Clear("C");
    StPicoEvent* event = new ((*arrays[StPicoArrays::Event])[0]) StPicoEvent;
    event->setPrimaryVertexPosition(0.f, 0.f, 0.f);
    event->setPrimaryVertexPositionError(0.01f, 0.01f, 0.015f);
    event->setBField(0.); // Straight-line fixture in a valid homogeneous zero field.
  }
  void AddTrack(int index, int id, int charge, const TVector3& momentum,
                const TVector3& decay, bool isProton) {
    StPicoTrack* track = new ((*arrays[StPicoArrays::Track])[index]) StPicoTrack;
    track->setId(id);
    track->setNHitsFit(35 * charge);
    track->setNHitsMax(40);
    track->setNHitsDedx(25);
    track->setDedx(2.e-6f); // Setter receives GeV/cm; getter returns keV/cm.
    track->setDedxError(0.08f);
    track->setNSigmaPion(isProton ? 10.f : 0.125f);
    track->setNSigmaKaon(10.f);
    track->setNSigmaProton(isProton ? -0.25f : 10.f);
    track->setBTofPidTraitsIndex(-1);
    track->setGlobalMomentum(momentum);

    // StDcaGeometry origin is the closest point to the beam axis. For B=0,
    // project the displaced daughter line analytically to that point.
    const double psi = momentum.Phi();
    const double flight = -(decay.X() * momentum.X() + decay.Y() * momentum.Y()) /
                           momentum.Perp2();
    const TVector3 origin = decay + flight * momentum;
    track->setOrigin(origin);
    Float_t parameters[6] = {
      static_cast<Float_t>(decay.Y() * std::cos(psi) - decay.X() * std::sin(psi)),
      static_cast<Float_t>(origin.Z()), static_cast<Float_t>(psi),
      static_cast<Float_t>(-charge / momentum.Pt()),
      static_cast<Float_t>(momentum.Z() / momentum.Pt()), 0.f};
    Float_t sigmas[5] = {0.03f, 0.04f, 0.002f, 0.002f, 0.003f};
    Float_t correlations[10] = {0.01f, -0.02f, 0.015f, -0.01f, 0.025f,
                               -0.015f, 0.012f, -0.018f, 0.022f, -0.008f};
    StPicoTrackCovMatrix* covariance =
        new ((*arrays[StPicoArrays::TrackCovMatrix])[index]) StPicoTrackCovMatrix;
    covariance->setParams(parameters);
    covariance->setSigmas(sigmas);
    covariance->setCorrelations(correlations);
  }
  void Lambda(int sign = 1, int protonId = 700, int pionId = 19) {
    Reset();
    const double m2 = kFixtureMass * kFixtureMass;
    const double sum = kProtonMass + kPionMass;
    const double difference = kProtonMass - kPionMass;
    const double restP = std::sqrt((m2 - sum * sum) * (m2 - difference * difference)) /
                         (2. * kFixtureMass);
    const double protonEnergy = std::sqrt(kProtonMass * kProtonMass + restP * restP);
    const double pionEnergy = std::sqrt(kPionMass * kPionMass + restP * restP);
    // Boost a physical two-body decay to parent px=1 GeV/c. Small daughter-z
    // residual yields positive, nonzero fit chi2 (upstream explicitly requires it).
    const TVector3 protonP(protonEnergy / kFixtureMass, restP, 0.);
    const TVector3 pionP(pionEnergy / kFixtureMass, -restP, 0.);
    AddTrack(0, pionId, -sign, pionP, TVector3(5., 0., 0.003), false);
    AddTrack(1, protonId, sign, protonP, TVector3(5., 0., 0.), true);
  }
  StPicoTrack* Track(int index) { return dst.track(index); }
  StPicoTrackCovMatrix* Covariance(int index) { return dst.trackCovMatrix(index); }
  TClonesArray* arrays[StPicoArrays::NAllPicoArrays];
  StPicoDst dst;
};

const KfLambdaCandidate& FindLambda(const StPicoKFParticleInterface& interface,
                                  int sign, int protonId, int pionId) {
  const std::vector<KfLambdaCandidate>& candidates = interface.Candidates();
  for (size_t i = 0; i < candidates.size(); ++i) {
    const KfLambdaCandidate& c = candidates[i];
    if (c.pdg == sign * 3122 && c.protonId == protonId && c.pionId == pionId) return c;
  }
  const KfParticleEventStats& s = interface.Stats();
  std::ostringstream error;
  error << "expected Lambda fixture missing: quality=" << s.qualityTracks
        << " covariance=" << s.covarianceTracks << " pid=" << s.pidHypotheses
        << " primary=" << s.primaryHypotheses << " finder=" << s.finderParticles
        << " lambda=" << s.lambdaParticles << " invalid=" << s.invalidCandidates;
  throw std::runtime_error(error.str());
}

void TestCovariance(Fixture& fixture, KfParticleCutConfig& cuts) {
  KfParticleHelper helper(cuts);
  for (int index = 0; index < 2; ++index) {
    const StPicoTrackCovMatrix* c = fixture.Covariance(index);
    Float_t covariance[15] = {0.f};
    const int diagonal[5] = {0, 2, 5, 9, 14};
    const int offDiagonal[10] = {1, 3, 4, 6, 7, 8, 10, 11, 12, 13};
    const int rows[10] = {1, 2, 2, 3, 3, 3, 4, 4, 4, 4};
    const int cols[10] = {0, 0, 1, 0, 1, 2, 0, 1, 2, 3};
    for (int i = 0; i < 5; ++i) covariance[diagonal[i]] = c->sigmas()[i] * c->sigmas()[i];
    for (int i = 0; i < 10; ++i)
      covariance[offDiagonal[i]] = c->correlations()[i] * c->sigmas()[rows[i]] * c->sigmas()[cols[i]];
    StDcaGeometry reference;
    reference.set(c->params(), covariance);
    Double_t parameters[6], cartesian[21];
    reference.GetXYZ(parameters, cartesian);
    kfp::KFPTrack converted;
    Require(helper.BuildTrack(&fixture.dst, index, converted), "BuildTrack rejected a valid helix");
    for (int i = 0; i < 6; ++i)
      Near(converted.GetParameter(i), parameters[i], 2.e-6, "Cartesian parameter");
    for (int i = 0; i < 21; ++i)
      Near(converted.GetCovariance(i), cartesian[i], 2.e-9, "Cartesian covariance");
    Require(converted.GetID() == fixture.Track(index)->id(), "original track ID not retained");
    Require(converted.Charge() == fixture.Track(index)->charge(), "track charge changed");
    Require(converted.GetNDF() == 1, "reference GetTrack NDF changed");
    Near(converted.GetChi2(), 0., 0., "reference GetTrack chi2");
  }
  kfp::KFPTrack dummy;
  Require(!helper.BuildTrack(&fixture.dst, -1, dummy), "negative track index accepted");
  Require(!helper.BuildTrack(&fixture.dst, 2, dummy), "out-of-range track index accepted");
  std::cout << "PASS covariance: 6 Cartesian parameters + 21 covariance entries, both signs" << std::endl;
}

void TestEvents(Fixture& fixture, KfParticleCutConfig& cuts) {
  StPicoKFParticleInterface interface(cuts);
  Require(!interface.ProcessEvent(0), "null PicoDst accepted");
  Require(interface.ProcessEvent(&fixture.dst), interface.LastError());
  const KfLambdaCandidate& lambda = FindLambda(interface, 1, 700, 19);
  Require(lambda.protonIndex == 1 && lambda.pionIndex == 0, "sorted/noncontiguous daughter mapping failed");
  Near(lambda.mass, kFixtureMass, 0.001, "unconstrained Lambda mass");
  Require(std::fabs(lambda.mass - 1.115683) > 0.003, "parent mass appears artificially constrained");
  Near(lambda.protonPidPull, -0.25, 1.e-5, "stored proton PID");
  Near(lambda.pionPidPull, 0.125, 1.e-5, "stored pion PID");
  Require(lambda.topoChi2Ndf < 3.f && lambda.massError > 0.f,
          "invalid topology or mass uncertainty");
  std::cout << "PASS Lambda: mass=" << lambda.mass << " topoChi2Ndf=" << lambda.topoChi2Ndf << std::endl;

  fixture.Lambda(-1, 77, 1205);
  Require(interface.ProcessEvent(&fixture.dst), interface.LastError());
  const KfLambdaCandidate& anti = FindLambda(interface, -1, 77, 1205);
  Near(anti.mass, kFixtureMass, 0.001, "anti-Lambda mass");
  Require(interface.Stats().rawTracks == 2, "event counters accumulated instead of reset");
  fixture.Reset();
  Require(interface.ProcessEvent(&fixture.dst), interface.LastError());
  Require(interface.Candidates().empty() && interface.Stats().finderParticles == 0,
          "previous event candidates survived an empty event");

  fixture.Lambda();
  fixture.Track(1)->setNSigmaPion(0.25f);
  Require(interface.ProcessEvent(&fixture.dst), interface.LastError());
  Require(interface.Stats().pidHypotheses == 3 && interface.Stats().pionHypotheses == 2,
          "multiple PID hypotheses were collapsed");
  FindLambda(interface, 1, 700, 19);
  fixture.Track(1)->setId(19);
  Require(!interface.ProcessEvent(&fixture.dst), "duplicate original IDs accepted");
  Require(interface.Stats().invalidTrackIds == 1 && interface.Candidates().empty(),
          "duplicate-ID failure retained candidates");

  fixture.Lambda();
  TClonesArray* saved = fixture.arrays[StPicoArrays::TrackCovMatrix];
  fixture.arrays[StPicoArrays::TrackCovMatrix] = 0;
  const bool missingAccepted = interface.ProcessEvent(&fixture.dst);
  fixture.arrays[StPicoArrays::TrackCovMatrix] = saved;
  Require(!missingAccepted, "missing covariance branch accepted");
  saved->Clear("C");
  Require(!interface.ProcessEvent(&fixture.dst), "mismatched covariance count accepted");
  fixture.Lambda();
  Float_t invalidSigmas[5] = {0.03f, 0.04f, 0.002f, 0.002f,
                             std::numeric_limits<Float_t>::quiet_NaN()};
  fixture.Covariance(0)->setSigmas(invalidSigmas);
  Require(interface.ProcessEvent(&fixture.dst), interface.LastError());
  Require(interface.Stats().invalidCovariance == 1 && interface.Candidates().empty(),
          "nonfinite covariance was not rejected");
  std::cout << "PASS anti-Lambda, multi-PID, event reset, duplicate IDs and covariance errors" << std::endl;
}


void ConfigureImp5Fixture(KfParticleCutConfig& cuts) {
  ConfigureFixture(cuts);
  cuts.selectionProfile = "lambda_imp5";
  cuts.imp5MinDCAProton = 0.;
  cuts.imp5MinDCAPion = 0.;
  cuts.imp5MaxPathLength = 1000.;
  Require(cuts.Validate(std::cerr), "synthetic Imp5 configuration is invalid");
}

void TestImp5(Fixture& fixture, KfParticleCutConfig& cuts) {
  ConfigureImp5Fixture(cuts);
  fixture.Lambda();
  fixture.Track(0)->setNSigmaPion(3.f);
  fixture.Track(1)->setNSigmaProton(-3.f);
  const double pionDca = fixture.Track(0)->gDCA(0.f, 0.f, 0.f);
  const double protonDca = fixture.Track(1)->gDCA(0.f, 0.f, 0.f);
  Require(pionDca > protonDca + 0.01 && protonDca > 0., "unexpected fixture DCA geometry");
  cuts.imp5MinDCAPion = pionDca;
  cuts.imp5MinDCAProton = protonDca;
  {
    StPicoKFParticleInterface interface(cuts);
    Require(interface.ProcessEvent(&fixture.dst), interface.LastError());
    const KfLambdaCandidate& candidate = FindLambda(interface, 1, 700, 19);
    Require(interface.Stats().pidHypotheses == 2, "Imp5 inclusive PID/DCA boundaries rejected daughters");
    Near(candidate.pionPidPull, 3., 1.e-6, "inclusive positive stored pion nSigma");
    Near(candidate.protonPidPull, -3., 1.e-6, "inclusive negative stored proton nSigma");
    Near(candidate.pionDcaToPv, pionDca, 1.e-6, "original Pico pion DCA");
    Near(candidate.protonDcaToPv, protonDca, 1.e-6, "original Pico proton DCA");
    Require(candidate.imp5PathValid && std::isfinite(candidate.protonHelixPathLength) &&
            std::isfinite(candidate.pionHelixPathLength), "original Imp5 helix path QA is invalid");
    const double expectedProtonPath = (TVector3(5., 0., 0.) - fixture.Track(1)->origin()).Mag();
    const double expectedPionPath = (TVector3(5., 0., 0.003) - fixture.Track(0)->origin()).Mag();
    Near(candidate.protonHelixPathLength, expectedProtonPath, 1.e-5, "original straight proton helix path");
    Near(candidate.pionHelixPathLength, expectedPionPath, 1.e-5, "original straight pion helix path");
    Near(candidate.mass, kFixtureMass, 0.001, "Imp5 retains unconstrained KF mass");
  }
  fixture.Track(0)->setNSigmaPion(3.1f);
  {
    StPicoKFParticleInterface interface(cuts);
    Require(interface.ProcessEvent(&fixture.dst), interface.LastError());
    Require(interface.Stats().pionHypotheses == 0 && interface.Stats().protonHypotheses == 1 &&
            interface.Candidates().empty(), "Imp5 outside-nSigma pion was accepted");
  }

  ConfigureImp5Fixture(cuts);
  fixture.Lambda();
  // Both species hypotheses are independent. The proton track also passes pion
  // nSigma, but a proton-specific DCA cut may remove only its proton hypothesis.
  fixture.Track(1)->setNSigmaPion(0.f);
  cuts.imp5MinDCAProton = protonDca + 0.01;
  {
    StPicoKFParticleInterface interface(cuts);
    Require(interface.ProcessEvent(&fixture.dst), interface.LastError());
    Require(interface.Stats().protonHypotheses == 0 && interface.Stats().pionHypotheses == 2,
            "proton DCA cut vetoed an otherwise valid independent pion hypothesis");
  }
  cuts.imp5MinDCAProton = 0.;
  cuts.imp5MinDCAPion = protonDca + 0.01;
  {
    StPicoKFParticleInterface interface(cuts);
    Require(interface.ProcessEvent(&fixture.dst), interface.LastError());
    Require(interface.Stats().protonHypotheses == 1 && interface.Stats().pionHypotheses == 1,
            "pion DCA cut vetoed the valid proton or failed to remove its pion hypothesis");
    FindLambda(interface, 1, 700, 19);
  }

  ConfigureImp5Fixture(cuts);
  fixture.Lambda();
  // Production Imp5 does not apply these unrelated quality cuts. Keep only the
  // original hit count/ratio selection and the finite, nonzero momentum safety.
  cuts.minNHitsFit = 35;
  cuts.minNHitsDedx = 99;
  cuts.minPt = 5.; cuts.maxPt = 6.;
  cuts.minEta = 2.; cuts.maxEta = 3.;
  cuts.minDedxError = 0.1; cuts.maxDedxError = 0.2;
  for (int index = 0; index < 2; ++index) {
    fixture.Track(index)->setNHitsMax(0);
    fixture.Track(index)->setNHitsDedx(0);
    fixture.Track(index)->setDedxError(std::numeric_limits<float>::quiet_NaN());
    // Stored nSigma is quantized Short_t/1000, so no safe NaN fixture exists.
    // An extreme unrelated kaon pull must not veto a valid pion/proton.
    fixture.Track(index)->setNSigmaKaon(-32.f);
    fixture.Track(index)->setBTofPidTraitsIndex(123);
  }
  {
    StPicoKFParticleInterface interface(cuts);
    TClonesArray* savedTof = fixture.arrays[StPicoArrays::BTofPidTraits];
    fixture.arrays[StPicoArrays::BTofPidTraits] = 0;
    const bool processed = interface.ProcessEvent(&fixture.dst);
    fixture.arrays[StPicoArrays::BTofPidTraits] = savedTof;
    Require(processed, interface.LastError());
    const KfLambdaCandidate& candidate = FindLambda(interface, 1, 700, 19);
    Require(interface.Stats().qualityTracks == 2 && interface.Stats().pidHypotheses == 2 &&
            interface.Stats().tofTracks == 0 && !candidate.protonHasTof && !candidate.pionHasTof,
            "Imp5 applied unrelated quality/kaon/TOF selection or rejected zero hit-ratio denominator");
  }
  fixture.Track(1)->setNHitsFit(34);
  {
    StPicoKFParticleInterface interface(cuts);
    Require(interface.ProcessEvent(&fixture.dst), interface.LastError());
    Require(interface.Stats().qualityTracks == 1, "Imp5 minNHitsFit was bypassed with the other quality cuts");
  }
  fixture.Track(1)->setNHitsFit(35);
  fixture.Track(0)->setGlobalMomentum(TVector3(0., 0., 0.));
  {
    StPicoKFParticleInterface interface(cuts);
    Require(interface.ProcessEvent(&fixture.dst), interface.LastError());
    Require(interface.Stats().qualityTracks == 1, "Imp5 zero-momentum safety was bypassed");
  }
  fixture.Track(0)->setGlobalMomentum(TVector3(std::numeric_limits<float>::quiet_NaN(), 0., 0.));
  {
    StPicoKFParticleInterface interface(cuts);
    Require(interface.ProcessEvent(&fixture.dst), interface.LastError());
    Require(interface.Stats().qualityTracks == 1, "Imp5 nonfinite-momentum safety was bypassed");
  }

  ConfigureImp5Fixture(cuts);
  fixture.Lambda();
  cuts.minNHitsRatio = 0.875; // Exactly 35 / 40.
  fixture.Track(0)->setNSigmaKaon(0.f);
  fixture.Track(1)->setNSigmaKaon(0.f);
  {
    StPicoKFParticleInterface interface(cuts);
    Require(interface.ProcessEvent(&fixture.dst), interface.LastError());
    Require(interface.Stats().qualityTracks == 2 && interface.Stats().pidHypotheses == 2 &&
            interface.Stats().protonHypotheses == 1 && interface.Stats().pionHypotheses == 1,
            "Imp5 ratio boundary changed or kaon hypotheses leaked into the Lambda-only profile");
    FindLambda(interface, 1, 700, 19);
  }
  fixture.Track(1)->setNHitsFit(34);
  {
    StPicoKFParticleInterface interface(cuts);
    Require(interface.ProcessEvent(&fixture.dst), interface.LastError());
    Require(interface.Stats().qualityTracks == 1, "positive-denominator hit-ratio rejection was ignored");
  }

  ConfigureImp5Fixture(cuts);
  fixture.Lambda(-1, 77, 1205);
  {
    StPicoKFParticleInterface interface(cuts);
    Require(interface.ProcessEvent(&fixture.dst), interface.LastError());
    Near(FindLambda(interface, -1, 77, 1205).mass, kFixtureMass, 0.001,
         "signed Imp5 anti-Lambda mass");
  }
  cuts.reconstructAntiLambda = kFALSE;
  {
    StPicoKFParticleInterface interface(cuts);
    Require(interface.ProcessEvent(&fixture.dst), interface.LastError());
    Require(interface.Candidates().empty(), "Imp5 ignored disabled anti-Lambda reconstruction");
  }

  // Restore and actively check the historical behavior before the existing
  // reference-profile TOF and dEdx tests, not merely the profile string.
  ConfigureFixture(cuts);
  fixture.Lambda();
  fixture.Track(0)->setNHitsMax(0);
  {
    StPicoKFParticleInterface interface(cuts);
    Require(interface.ProcessEvent(&fixture.dst), interface.LastError());
    Require(interface.Stats().qualityTracks == 1, "Imp5 zero-denominator bypass leaked into kf_reference");
  }
  fixture.Lambda();
  fixture.Track(0)->setNSigmaPion(3.f);
  fixture.Track(1)->setNSigmaProton(-3.f);
  {
    StPicoKFParticleInterface interface(cuts);
    Require(interface.ProcessEvent(&fixture.dst), interface.LastError());
    Require(interface.Stats().unknownHypotheses == 2 && interface.Candidates().empty(),
            "Imp5 inclusive nSigma boundary leaked into the strict reference PID");
  }
  fixture.Lambda();
  std::cout << "PASS opt-in Imp5: independent inclusive PID/DCA, Lambda-only hypotheses, quality/TOF "
               "scope, hit boundaries, signed reconstruction, original helix QA and default restoration" << std::endl;
}

double TofMean(const KfParticleCutConfig& cuts, int row, double p) {
  const double x = std::min(p, cuts.tofPMax);
  double result = 0.;
  for (int term = 0; term < 5; ++term) result += cuts.tofMean[row * 5 + term] * std::pow(x, term);
  return result;
}

void TestTof(Fixture& fixture, KfParticleCutConfig& cuts) {
  fixture.Lambda();
  cuts.useTof = kTRUE;
  Require(cuts.Validate(std::cerr), "TOF coefficients missing from provided KF YAML");
  StPicoKFParticleInterface interface(cuts);
  KfParticleHelper helper(cuts);
  kfp::KFPTrack pion;
  Require(helper.BuildTrack(&fixture.dst, 0, pion), "cannot construct TOF pion fixture");
  StPicoBTofPidTraits* tof =
      new ((*fixture.arrays[StPicoArrays::BTofPidTraits])[0]) StPicoBTofPidTraits;
  tof->setTrackIndex(0);
  const double massSquared = TofMean(cuts, 3, pion.GetP()); // pi-minus coefficient row
  tof->setBeta(pion.GetP() / std::sqrt(pion.GetP() * pion.GetP() + massSquared));
  fixture.Track(0)->setBTofPidTraitsIndex(0);
  Require(interface.ProcessEvent(&fixture.dst), interface.LastError());
  const KfLambdaCandidate& candidate = FindLambda(interface, 1, 700, 19);
  Require(interface.Stats().tofTracks == 1 && candidate.pionHasTof && !candidate.protonHasTof,
          "TOF index zero was discarded or attached to the wrong daughter");
  const double beta = tof->btofBeta(); // Account for actual Pico UShort_t beta quantization.
  Near(candidate.pionTofM2, pion.GetP() * pion.GetP() * (1. / (beta * beta) - 1.),
       1.e-6, "TOF mass-squared");
  tof->setBeta(0.f);
  Require(interface.ProcessEvent(&fixture.dst), interface.LastError());
  Require(!FindLambda(interface, 1, 700, 19).pionHasTof,
          "invalid beta became a valid TOF measurement");
  fixture.Track(0)->setBTofPidTraitsIndex(4);
  Require(!interface.ProcessEvent(&fixture.dst), "out-of-range referenced TOF index accepted");
  std::cout << "PASS TOF index zero, quantized beta, missing beta and invalid reference" << std::endl;
}

void TestDedx(Fixture& fixture, KfParticleCutConfig& cuts) {
  fixture.Lambda();
  cuts.useTof = kFALSE;
  cuts.pidProfile = "dedx_pull";
  for (int index = 0; index < 2; ++index) {
    StPicoTrack* track = fixture.Track(index);
    const Float_t mass = index == 0 ? kPionMass : kProtonMass;
    const Float_t momentum = track->gMom().Mag();
    const Float_t betaGamma = momentum / mass;
    track->setDedx(StdEdxPull::EvalPred(betaGamma, 1, 1));
    // Saved nSigma intentionally disagrees: only the selected profile may be used.
    track->setNSigmaPion(10.f);
    track->setNSigmaProton(10.f);
  }
  StPicoKFParticleInterface interface(cuts);
  Require(interface.ProcessEvent(&fixture.dst), interface.LastError());
  const KfLambdaCandidate& candidate = FindLambda(interface, 1, 700, 19);
  for (int index = 0; index < 2; ++index) {
    const StPicoTrack* track = fixture.Track(index);
    const Float_t mass = index == 0 ? kPionMass : kProtonMass;
    const Float_t momentum = track->gMom().Mag();
    const Float_t betaGamma = momentum / mass;
    const Float_t measured = 1.e-6 * track->dEdx();
    const Float_t expected = StdEdxPull::Eval(measured, track->dEdxError(), betaGamma, 1, 1);
    Near(index == 0 ? candidate.pionPidPull : candidate.protonPidPull,
         expected, 1.e-5, "SL24y dEdxPull public-API equivalent");
  }
  fixture.Track(0)->setDedx(-1.f);
  Require(interface.ProcessEvent(&fixture.dst), interface.LastError());
  Require(interface.Stats().invalidPid == 1 && interface.Candidates().empty(),
          "invalid dEdx silently fell back to saved nSigma");
  std::cout << "PASS dEdxPull numerical wrapper and no saved-nSigma fallback" << std::endl;
}
}

extern "C" int star_analyzer_kfp_pico_adapter_test(const char* cutsPath) {
  if (!cutsPath || !*cutsPath) return 2;
  try {
    // The executable links STAR libraries, including the old global StarRoot
    // KF types. The adapter must still call only namespaced vendored KF symbols.
    Require(gSystem->Load("StPicoEvent") >= 0, "cannot load SL24y Pico dictionaries");
    TestKfEventSelection();
    KfParticleCutConfig& cuts = KfParticleCutConfig::GetInstance();
    Require(cuts.LoadFromFile(cutsPath), "cannot load the provided KF configuration");
    TestConfiguration(cutsPath, cuts);
    ConfigureFixture(cuts);
    Fixture fixture;
    fixture.Lambda();
    TestCovariance(fixture, cuts);
    TestEvents(fixture, cuts);
    TestImp5(fixture, cuts);
    TestTof(fixture, cuts);
    TestDedx(fixture, cuts);
    std::cout << "PASS full Pico adapter synthetic integration (not real-data validation)" << std::endl;
  } catch (const std::exception& error) {
    std::cerr << "FAIL Pico adapter: " << error.what() << std::endl;
    return 1;
  }
  return 0;
}

#ifndef STAR_ANALYZER_KFP_ROOT_TEST
int main(int argc, char** argv) {
  if (argc != 2) {
    std::cerr << "usage: " << argv[0] << " <KF cut YAML>" << std::endl;
    return 2;
  }
  return star_analyzer_kfp_pico_adapter_test(argv[1]);
}
#endif




