// Synthetic regression test for the full, namespaced KFParticle engine.
// These exact fixtures/tolerances are test inputs, not analysis cut defaults.
// This does not validate PicoDst conversion, PID calibration, or mass spectra.
#include "KFParticle.h"
#include "KFParticleSIMD.h"
#include "KFParticleTopoReconstructor.h"
#include "KFPTrack.h"
#include "KFPVertex.h"
#include "KFVertex.h"

#include "TSystem.h"

#include <cmath>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace kfp = star_analyzer_kfp;

namespace {

const double kProtonMass = 0.9382720813;
const double kPionMass = 0.13957039;
const double kLambdaMass = 1.115683;
// Intentionally away from the PDG mass: a hidden parent mass constraint must fail.
const double kFixtureMass = kLambdaMass + 0.005;

void Require(bool condition, const std::string& message)
{
  if (!condition) throw std::runtime_error(message);
}

std::string EventLabel(int pdg, float field)
{
  std::ostringstream text;
  text << "PDG=" << pdg << ", Bz=" << field << ": ";
  return text.str();
}

kfp::KFPTrack Track(int id, int charge, float px, float py, float z)
{
  kfp::KFPTrack track;
  float parameters[6] = {5.f, 0.f, z, px, py, 0.f};
  float covariance[21] = {};
  // Independent 0.02 cm position and 0.001 GeV/c momentum uncertainties.
  covariance[0] = covariance[2] = covariance[5] = 0.0004f;
  covariance[9] = covariance[14] = covariance[20] = 0.000001f;
  track.SetParameters(parameters);
  track.SetCovarianceMatrix(covariance);
  track.SetCharge(charge);
  track.SetID(id);
  track.SetChi2(1.f);
  track.SetNDF(1);
  return track;
}

void AddHypothesis(const kfp::KFPTrack& track, int pdg,
                   std::vector<kfp::KFParticle>& particles,
                   std::vector<int>& pdgs)
{
  kfp::KFParticle particle(track, pdg);
  // KFParticle(KFPTrack, PDG) does not transfer the input track ID.
  particle.SetId(track.GetID());
  particle.AddDaughterId(track.GetID());
  particles.push_back(particle);
  pdgs.push_back(pdg);
}

void Fixture(int sign, int protonId, int pionId,
             std::vector<kfp::KFParticle>& particles, std::vector<int>& pdgs)
{
  const double mass2 = kFixtureMass * kFixtureMass;
  const double sum = kProtonMass + kPionMass;
  const double difference = kProtonMass - kPionMass;
  const double q = std::sqrt((mass2 - sum * sum) *
                            (mass2 - difference * difference)) / (2 * kFixtureMass);
  const double protonEnergy = std::sqrt(kProtonMass * kProtonMass + q * q);
  const double pionEnergy = std::sqrt(kPionMass * kPionMass + q * q);
  // Two-body decay boosted along x to parent momentum 1 GeV/c.
  const double gammaBeta = 1.0 / kFixtureMass;
  const kfp::KFPTrack proton = Track(protonId, sign,
      float(gammaBeta * protonEnergy), float(q), 0.003f);
  const kfp::KFPTrack pion = Track(pionId, -sign,
      float(gammaBeta * pionEnergy), float(-q), -0.003f);

  // Non-contiguous IDs, deliberately unsorted input, and a second PID
  // hypothesis of the SAME proton track exercise sorting/identity separation.
  AddHypothesis(pion, -sign * 211, particles, pdgs);
  AddHypothesis(proton, sign * 321, particles, pdgs);
  AddHypothesis(proton, sign * 2212, particles, pdgs);
}

kfp::KFVertex PrimaryVertex()
{
  kfp::KFPVertex input;
  input.SetXYZ(0.f, 0.f, 0.f);
  input.SetCovarianceMatrix(0.0001f, 0.f, 0.0001f, 0.f, 0.f, 0.0001f);
  input.SetChi2(1.f);
  input.SetNDF(17);
  input.SetNContributors(10);
  return kfp::KFVertex(input);
}

void Configure(kfp::KFParticleTopoReconstructor& topo, bool includeAnti = true)
{
  kfp::KFParticleFinder* finder = topo.GetKFParticleFinder();
  finder->AddDecayToReconstructionList(3122);
  if (includeAnti) finder->AddDecayToReconstructionList(-3122);
  finder->SetMaxDistanceBetweenParticlesCut(1.5f);
  finder->SetLCut(1.f);
  finder->SetChiPrimaryCut2D(3.f);
  finder->SetChi2Cut2D(10.f);
  finder->SetLdLCut2D(3.f);
}

void CheckField(float expected)
{
  kfp::KFParticle scalar;
  float position[3] = {}, field[3] = {};
  scalar.GetFieldValue(position, field);
  Require(field[0] == 0.f && field[1] == 0.f && field[2] == expected,
          "Scalar magnetic field not updated");

  kfp::KFParticleSIMD simd;
  kfp::float32_v simdPosition[3] = {0.f, 0.f, 0.f};
  kfp::float32_v simdField[3];
  simd.GetFieldValue(simdPosition, simdField);
  for (int lane = 0; lane < kfp::SimdLen; ++lane)
    Require(simdField[2][lane] == expected, "SIMD magnetic field not updated");
}

void Reconstruct(kfp::KFParticleTopoReconstructor& topo,
                 std::vector<kfp::KFParticle>& input, std::vector<int>& pdgs,
                 float field, bool addPv = true)
{
  topo.Clear();
  Require(topo.GetParticles().empty() && topo.NPrimaryVertices() == 0,
          "Clear retained particles or PVs");
  topo.SetField(field);
  CheckField(field);
  topo.Init(input, &pdgs);
  if (addPv) {
    const std::vector<int> primaryTrackIndices;
    topo.AddPV(PrimaryVertex(), primaryTrackIndices);
    topo.FillPVIndices();
    Require(topo.NPrimaryVertices() == 1, "Input PV not registered exactly once");
  }
  topo.SortTracks();
  topo.ReconstructParticles();
}

int CountLambda(const kfp::KFParticleTopoReconstructor& topo)
{
  int count = 0;
  const std::vector<kfp::KFParticle>& particles = topo.GetParticles();
  for (size_t i = 0; i < particles.size(); ++i)
    if (std::abs(particles[i].GetPDG()) == 3122) ++count;
  return count;
}

float CheckLambda(const kfp::KFParticleTopoReconstructor& topo,
                  int expectedPdg, int protonId, int pionId, float field)
{
  const std::string context = EventLabel(expectedPdg, field);
  const std::vector<kfp::KFParticle>& particles = topo.GetParticles();
  Require(CountLambda(topo) == 1, context + "expected exactly one selected Lambda");
  for (size_t i = 0; i < particles.size(); ++i) {
    const kfp::KFParticle& raw = particles[i];
    if (std::abs(raw.GetPDG()) != 3122) continue;
    Require(raw.GetPDG() == expectedPdg, context + "wrong parent PDG sign");
    Require(raw.NDaughters() == 2, context + "wrong daughter count");
    int foundProtonId = -1, foundPionId = -1;
    for (int daughter = 0; daughter < raw.NDaughters(); ++daughter) {
      const int particleIndex = raw.DaughterIds()[daughter];
      Require(particleIndex >= 0 && size_t(particleIndex) < particles.size(),
              context + "composite daughter is not a particle-array index");
      const kfp::KFParticle& leaf = particles[particleIndex];
      Require(leaf.NDaughters() == 1, context + "track particle must have one input ID");
      const int trackId = leaf.DaughterIds()[0];
      const int sign = expectedPdg > 0 ? 1 : -1;
      if (leaf.GetPDG() == sign * 2212) foundProtonId = trackId;
      else if (leaf.GetPDG() == -sign * 211) foundPionId = trackId;
      else Require(false, context + "unexpected daughter PDG hypothesis");
    }
    Require(foundProtonId == protonId && foundPionId == pionId && protonId != pionId,
            context + "sorting/multiple PID corrupted daughter track identity");

    float mass = 0.f, massError = 0.f;
    Require(raw.GetMass(mass, massError) == 0 && std::isfinite(mass) &&
            std::isfinite(massError) && massError > 0.f,
            context + "invalid unconstrained mass/error");
    Require(std::fabs(mass - kFixtureMass) < 0.002,
            context + "raw mass does not match synthetic two-body kinematics");
    Require(std::fabs(mass - kLambdaMass) > 0.001,
            context + "raw parent mass appears constrained to the PDG value");
    Require(raw.NDF() > 0 && std::isfinite(raw.Chi2()), context + "invalid raw fit quality");

    float parameters[8], covariance[36];
    for (int j = 0; j < 8; ++j) parameters[j] = raw.GetParameter(j);
    for (int j = 0; j < 36; ++j) covariance[j] = raw.GetCovariance(j);

    kfp::KFParticle withPv = raw;
    withPv.SetProductionVertex(topo.GetPrimVertex());
    float length = 0.f, lengthError = 0.f;
    Require(withPv.NDF() > 0 && std::isfinite(withPv.Chi2()) &&
            withPv.GetDecayLength(length, lengthError) == 0 &&
            std::isfinite(length) && std::isfinite(lengthError) &&
            length > 0.f && lengthError > 0.f,
            context + "invalid PV-constrained copy");
    kfp::KFParticle withMass = raw;
    withMass.SetNonlinearMassConstraint(float(kLambdaMass));
    float constrainedMass = 0.f, constrainedError = 0.f;
    const int constrainedStatus = withMass.GetMass(constrainedMass, constrainedError);
    Require(std::isfinite(constrainedMass) &&
            std::fabs(constrainedMass - kLambdaMass) < 0.0001 &&
            std::fabs(constrainedMass - mass) > 0.001 &&
            withMass.NDF() == raw.NDF() + 1 && std::isfinite(withMass.Chi2()),
            context + "explicit mass constraint did not work on the copy");
    for (int j = 0; j < 8; ++j)
      Require(std::isfinite(withMass.GetParameter(j)), context + "non-finite constrained state");
    for (int j = 0; j < 36; ++j)
      Require(std::isfinite(withMass.GetCovariance(j)), context + "non-finite constrained covariance");

    // An exact mass constraint removes the mass-variance direction. GetMass
    // still returns the correct mass but reports status=1/error=1000 if its
    // float evaluation of g^T C g rounds just below zero. Do not demand a
    // positive mass error after an EXACT constraint (unlike the raw fit above).
    // Independently require zero variance within float roundoff, relative to
    // the actual covariance-term scale, so a genuinely invalid covariance fails.
    double gradient[4] = {-double(withMass.GetParameter(3)),
                          -double(withMass.GetParameter(4)),
                          -double(withMass.GetParameter(5)),
                           double(withMass.GetParameter(6))};
    double massVarianceNumerator = 0., covarianceTermScale = 0.;
    for (int row = 0; row < 4; ++row) {
      for (int column = 0; column <= row; ++column) {
        const int index = (row + 3) * (row + 4) / 2 + column + 3;
        const double term = (row == column ? 1. : 2.) * gradient[row] *
                            gradient[column] * withMass.GetCovariance(index);
        massVarianceNumerator += term;
        covarianceTermScale += std::fabs(term);
      }
    }
    const double varianceRoundoff = 64. * std::numeric_limits<float>::epsilon() *
                                    covarianceTermScale;
    Require(std::isfinite(massVarianceNumerator) && covarianceTermScale > 0. &&
            std::fabs(massVarianceNumerator) <= varianceRoundoff,
            context + "exact mass constraint did not remove the mass variance");
    Require((constrainedStatus == 0 && std::isfinite(constrainedError) &&
             constrainedError >= 0.f &&
             constrainedError <= std::sqrt(varianceRoundoff) / constrainedMass) ||
            (constrainedStatus == 1 && constrainedError == 1000.f),
            context + "unexpected GetMass status/error for exact constrained copy");
    std::cout << "PASS mass-constrained copy " << context << "mass=" << constrainedMass
              << ", GetMass status=" << constrainedStatus
              << ", gCg=" << massVarianceNumerator
              << ", float roundoff bound=" << varianceRoundoff << '\n';

    for (int j = 0; j < 8; ++j)
      Require(parameters[j] == raw.GetParameter(j), context + "copy constraints changed raw parameters");
    for (int j = 0; j < 36; ++j)
      Require(covariance[j] == raw.GetCovariance(j), context + "copy constraints changed raw covariance");
    std::cout << "PASS " << context << "mass=" << mass << " +/- " << massError
              << ", track IDs=" << foundProtonId << ',' << foundPionId
              << ", decay length=" << length << " +/- " << lengthError << '\n';
    return mass;
  }
  throw std::runtime_error(context + "no Lambda found");
}

} // namespace

static int RunStarAnalyzerFullChainTest(const char* starRootArgument)
{
  try {
    const char* starLibraryDirectory = gSystem->Getenv("STAR_LIB");
    const bool hasArgument = starRootArgument && *starRootArgument;
    std::string starRoot = hasArgument ? starRootArgument : "libStarRoot";
    if (!hasArgument && starLibraryDirectory && *starLibraryDirectory)
      starRoot = std::string(starLibraryDirectory) + "/StarRoot.so";
    Require(gSystem->Load(starRoot.c_str()) >= 0,
            "Cannot load legacy StarRoot: " + starRoot);
    std::cout << "Loaded legacy StarRoot before full-chain fixture execution: " << starRoot << '\n';

    kfp::KFParticleTopoReconstructor topo;
    Configure(topo);
    std::vector<kfp::KFParticle> input;
    std::vector<int> pdgs;
    Fixture(1, 101, 7, input, pdgs);
    Reconstruct(topo, input, pdgs, 5.f);
    const float firstMass = CheckLambda(topo, 3122, 101, 7, 5.f);

    input.clear(); pdgs.clear();
    Fixture(-1, 503, 43, input, pdgs);
    Reconstruct(topo, input, pdgs, -5.f);
    CheckLambda(topo, -3122, 503, 43, -5.f);

    kfp::KFParticleTopoReconstructor positiveOnly;
    Configure(positiveOnly, false);
    Reconstruct(positiveOnly, input, pdgs, -5.f);
    Require(CountLambda(positiveOnly) == 0, "Signed decay list unexpectedly reconstructed anti-Lambda");

    input.clear(); pdgs.clear();
    Fixture(1, 909, 909, input, pdgs);
    Reconstruct(topo, input, pdgs, 5.f);
    Require(CountLambda(topo) == 0, "Candidate reused one physical input track twice");

    input.clear(); pdgs.clear();
    Reconstruct(topo, input, pdgs, 0.f);
    Require(topo.GetParticles().empty(), "Empty event retained previous candidates");

    Fixture(1, 101, 7, input, pdgs);
    Reconstruct(topo, input, pdgs, 5.f, false);
    Require(topo.GetParticles().empty(), "Event without PV retained/reconstructed candidates");
    Reconstruct(topo, input, pdgs, 5.f);
    const float repeatedMass = CheckLambda(topo, 3122, 101, 7, 5.f);
    Require(std::fabs(firstMass - repeatedMass) < 0.000001,
            "Repeating an event after reset changed its mass");
    std::cout << "PASS full Topo/Finder synthetic regression: signs, sorting, multiple PID, "
                 "track IDs, raw/copy constraints, duplicate rejection, empty/no-PV reset, field changes.\n"
                 "PicoDst/PID calibration and real-event physics are NOT validated by this test.\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "FAIL full Topo/Finder synthetic regression: " << error.what() << '\n';
    return 1;
  }
}

// ROOT5 loads this entry through TSystem::DynFindSymbol. Keeping the test in a
// compiled module exercises the same lazy-loaded STAR environment as analysis.
extern "C" int star_analyzer_kfp_full_chain_test(const char* starRootArgument)
{
  return RunStarAnalyzerFullChainTest(starRootArgument);
}

#ifndef STAR_ANALYZER_KFP_ROOT_TEST
int main(int argc, char** argv)
{
  return RunStarAnalyzerFullChainTest(argc > 1 ? argv[1] : 0);
}
#endif
