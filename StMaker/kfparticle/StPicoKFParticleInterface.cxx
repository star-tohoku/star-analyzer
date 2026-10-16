#include "StPicoKFParticleInterface.h"

#include "KFParticle.h"
#include "KFParticleFinder.h"
#include "KFParticleTopoReconstructor.h"
#include "KFParticleSnapshot.h"
#include "KFPTrack.h"
#include "KFPVertex.h"
#include "KFVertex.h"
#include "cuts/KfParticleCutConfig.h"
#include "StBichsel/StdEdxPull.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "TVector3.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <sstream>

namespace kfp = star_analyzer_kfp;

namespace {
const int kPdg[3] = {211, 321, 2212};
// PDG hypotheses used by the reference StPicoTrack::dEdxPull calls, not cuts.
const Float_t kMass[3] = {0.139570f, 0.493677f, 0.938272f};

bool FiniteVector(const TVector3& v) {
  return std::isfinite(v.X()) && std::isfinite(v.Y()) && std::isfinite(v.Z());
}

bool FiniteParticle(const kfp::KFParticle& p) {
  for (int i = 0; i < 8; ++i)
    if (!std::isfinite(p.Parameters()[i])) return false;
  for (int i = 0; i < 36; ++i)
    if (!std::isfinite(p.CovarianceMatrix()[i])) return false;
  return std::isfinite(p.GetChi2());
}

double Polynomial(const std::vector<double>& coefficients, int row, double p) {
  double result = coefficients[row * 5 + 4];
  for (int term = 3; term >= 0; --term)
    result = result * p + coefficients[row * 5 + term];
  return result;
}

struct TrackPid {
  TrackPid() : index(-1), hasTof(false), m2(0.), dcaToPv(0.) {
    for (int i = 0; i < 3; ++i) tpc[i] = tof[i] = 0.;
  }
  int index;
  bool hasTof;
  double m2, dcaToPv, tpc[3], tof[3];
};
}

KfParticleEventStats::KfParticleEventStats()
    : rawTracks(0), qualityTracks(0), covarianceTracks(0), pidTracks(0),
      pidHypotheses(0), protonHypotheses(0), pionHypotheses(0),
      unknownHypotheses(0), primaryHypotheses(0), finderParticles(0),
      lambdaParticles(0), validCandidates(0), xiParticles(0),
      validXiCandidates(0), invalidXiCandidates(0), xiDaughterRecovered(0),
      k0Particles(0), validK0Candidates(0), invalidK0Candidates(0),
      invalidCandidates(0), invalidCovariance(0), invalidPid(0),
      invalidTrackIds(0), tofTracks(0) {}

class StPicoKFParticleInterface::Impl {
public:
  explicit Impl(const KfParticleCutConfig& configuration)
      : cuts(configuration), helper(configuration) {
    std::ostringstream errors;
    if (!cuts.Validate(errors)) {
      configurationError = errors.str();
      return;
    }
    // For NDF=2 the inverse chi-square survival probability is exactly
    // -2 ln(prob). SetChi2PrimaryCut also sets Finder ChiPrimary2D, hence the
    // explicit Finder value MUST be applied after this probability setting.
    topo.SetChi2PrimaryCut(-2. * std::log(cuts.primaryProbCut));
    topo.SetLambdaTopoChi2NdfCut(cuts.topoChi2NdfCut);
    kfp::KFParticleFinder* finder = topo.GetKFParticleFinder();
    finder->SetMaxDistanceBetweenParticlesCut(cuts.finderMaxDaughterDistance);
    finder->SetLCut(cuts.finderLCut);
    finder->SetChiPrimaryCut2D(cuts.finderChiPrimary2D);
    finder->SetChi2Cut2D(cuts.finderChi2Ndf2D);
    finder->SetLdLCut2D(cuts.finderLdL2D);
    finder->SetApplyLambdaGeometryCuts(cuts.applyLambdaGeometryCuts);
    finder->SetLdLCutXiOmega(cuts.finderLdLXiOmega);
    finder->SetChi2TopoCutXiOmega(cuts.finderChi2TopoXiOmega);
    finder->SetChi2CutXiOmega(cuts.finderChi2XiOmega);
    // Whitelist: empty list reconstructs everything; any entry restricts.
    // Ξ requires Λ in the Finder list. K0S is PDG 310 only (self-conjugate).
    if (cuts.reconstructLambda) finder->AddDecayToReconstructionList(3122);
    if (cuts.reconstructAntiLambda || cuts.reconstructAntiXi)
      finder->AddDecayToReconstructionList(-3122);
    if (cuts.reconstructXi) finder->AddDecayToReconstructionList(3312);
    if (cuts.reconstructAntiXi) finder->AddDecayToReconstructionList(-3312);
    if (cuts.reconstructK0short) finder->AddDecayToReconstructionList(310);
  }

  void Clear() {
    stats = KfParticleEventStats();
    candidates.clear();
    xiCandidates.clear();
    k0Candidates.clear();
    error.clear();
    trackPid.clear();
    topo.Clear();
    rawTrackPt.clear();
    rawTrackEta.clear();
  }

  bool Fail(const std::string& reason) {
    candidates.clear();
    xiCandidates.clear();
    k0Candidates.clear();
    error = reason;
    return false;
  }

  bool Quality(const StPicoTrack& track) const {
    if (cuts.selectionProfile == "lambda_imp5") {
      // Match StLambdaMaker, not the generic (unused there) track YAML.
      if (!track.charge() || track.nHitsFit() < cuts.minNHitsFit) return false;
      if (track.nHitsMax() > 0 &&
          double(track.nHitsFit()) / double(track.nHitsMax()) < cuts.minNHitsRatio)
        return false;
      // Numerical safety for KF conversion, not an extra pT/eta acceptance.
      const TVector3 p = track.gMom();
      return FiniteVector(p) && p.Mag2() > 0.;
    }
    if (!track.charge() || track.nHitsFit() < cuts.minNHitsFit ||
        track.nHitsDedx() < cuts.minNHitsDedx) return false;
    if (cuts.minNHitsRatio > 0. &&
        (track.nHitsMax() <= 0 ||
         track.nHitsFit() < cuts.minNHitsRatio * track.nHitsMax())) return false;
    const double dedxError = track.dEdxError();
    if (!std::isfinite(dedxError) || dedxError < cuts.minDedxError ||
        dedxError > cuts.maxDedxError) return false;
    const TVector3 p = track.gMom();
    if (!FiniteVector(p) || !std::isfinite(p.Eta()) || p.Mag2() <= 0.) return false;
    return p.Pt() >= cuts.minPt && p.Pt() <= cuts.maxPt &&
           p.Eta() >= cuts.minEta && p.Eta() <= cuts.maxEta;
  }

  bool Pid(StPicoDst* dst, const StPicoTrack& track, double momentum,
           TrackPid& result, std::vector<int>& hypotheses) {
    if (cuts.selectionProfile == "lambda_imp5") {
      // Independent stored nSigma tests, inclusive at the legacy boundary.
      // No kaon/PID competition or TOF lookup, even with a stored TOF index.
      result.tpc[0] = track.nSigmaPion();
      result.tpc[2] = track.nSigmaProton();
      const TVector3 pv = dst->event()->primaryVertex();
      result.dcaToPv = track.gDCA(pv.X(), pv.Y(), pv.Z());
      if (!std::isfinite(result.dcaToPv)) return false;
      for (int species = 0; species < 3; species += 2) {
        if (!cuts.reconstructAntiLambda &&
            ((species == 0 && track.charge() > 0) ||
             (species == 2 && track.charge() < 0))) continue;
        const double limit = species == 0 ? cuts.nSigmaPion : cuts.nSigmaProton;
        const double minimumDca = species == 0 ? cuts.imp5MinDCAPion : cuts.imp5MinDCAProton;
        if (!std::isfinite(result.tpc[species]) ||
            std::fabs(result.tpc[species]) > limit || result.dcaToPv < minimumDca) continue;
        hypotheses.push_back(kPdg[species] * (track.charge() > 0 ? 1 : -1));
      }
      if (hypotheses.empty()) hypotheses.push_back(-1);
      return true;
    }
    if (cuts.pidProfile == "dedx_pull") {
      if (!std::isfinite(track.dEdx()) || track.dEdx() <= 0. ||
          !std::isfinite(track.dEdxError()) || track.dEdxError() <= 0.) return false;
      // Identical units, fit=1 and |charge|=1 to the SL24y __TFG__ guarded
      // StPicoTrack::dEdxPull, called through public APIs only.
      const Float_t picoMomentum = track.gMom().Mag();
      const Float_t measured = 1.e-6 * track.dEdx();
      const Float_t resolution = track.dEdxError();
      for (int species = 0; species < 3; ++species) {
        const Float_t betaGamma = picoMomentum / kMass[species];
        result.tpc[species] = static_cast<Float_t>(
            StdEdxPull::Eval(measured, resolution, betaGamma, 1, 1));
      }
    } else if (cuts.pidProfile == "pico_nsigma") {
      result.tpc[0] = track.nSigmaPion();
      result.tpc[1] = track.nSigmaKaon();
      result.tpc[2] = track.nSigmaProton();
    } else {
      error = "Unknown KF PID profile: " + cuts.pidProfile;
      return false;
    }
    for (int species = 0; species < 3; ++species)
      if (!std::isfinite(result.tpc[species]) || result.tpc[species] == -999.)
        return false; // -999 is the SL24y failed-pull sentinel, not a measurement.

    const int tofIndex = track.bTofPidTraitsIndex();
    if (cuts.useTof && tofIndex >= 0) {
      if (static_cast<unsigned int>(tofIndex) >= dst->numberOfBTofPidTraits()) {
        error = "BTofPidTraits index is out of range; check the enabled branch and input consistency";
        return false;
      }
      const StPicoBTofPidTraits* tof = dst->btofPidTraits(tofIndex);
      if (!tof) {
        error = "BTofPidTraits contains a null entry referenced by a track";
        return false;
      }
      const double beta = tof->btofBeta();
      // A missing/invalid beta is not a measurement. In soft mode it does not
      // veto an otherwise valid TPC hypothesis; valid beta can exceed one.
      if (std::isfinite(beta) && beta > 0. && beta * beta > 1.e-6) {
        result.m2 = momentum * momentum * (1. / (beta * beta) - 1.);
        result.hasTof = std::isfinite(result.m2);
      }
    }

    bool tofAccepted[3] = {false, false, false};
    if (result.hasTof) {
      ++stats.tofTracks;
      const double x = std::min(momentum, cuts.tofPMax);
      int bestSpecies = 0;
      for (int species = 0; species < 3; ++species) {
        const int row = species + (track.charge() < 0 ? 3 : 0);
        const double mean = Polynomial(cuts.tofMean, row, x);
        const double sigma = Polynomial(cuts.tofSigma, row, x);
        if (!std::isfinite(mean) || !std::isfinite(sigma) || sigma <= 0.) {
          error = "TOF calibration gives an invalid width; check the selected profile and coefficients";
          return false;
        }
        result.tof[species] = (result.m2 - mean) / sigma;
        if (!std::isfinite(result.tof[species])) return false;
        tofAccepted[species] = std::fabs(result.tof[species]) < cuts.tofNSigma;
        if (std::fabs(result.tof[species]) < std::fabs(result.tof[bestSpecies]))
          bestSpecies = species;
      }
      if (cuts.strictTofPid)
        for (int species = 0; species < 3; ++species)
          if (species != bestSpecies) tofAccepted[species] = false;
    }
    const double tpcCut[3] = {cuts.nSigmaPion, cuts.nSigmaKaon, cuts.nSigmaProton};
    const bool requireKaonTof = cuts.cleanKaonsWithTof &&
        momentum > cuts.tofKaonPMin && momentum < cuts.tofKaonPMax;
    for (int species = 0; species < 3; ++species) {
      if (std::fabs(result.tpc[species]) >= tpcCut[species]) continue;
      if (species == 1 && requireKaonTof && !tofAccepted[1]) continue;
      if (result.hasTof && !tofAccepted[species]) continue;
      hypotheses.push_back(kPdg[species] * (track.charge() > 0 ? 1 : -1));
    }
    // Preserve the reference Interface's unknown-PID input hypothesis. It is
    // not counted as accepted p/pi PID and does not silently become a proton.
    if (hypotheses.empty()) hypotheses.push_back(-1);
    return true;
  }

  bool Candidate(const kfp::KFParticle& raw,
                 const std::vector<kfp::KFParticle>& particles,
                 const kfp::KFVertex& pv, const TVector3& picoPv,
                 KfLambdaCandidate& out, StPicoDst* dst) {
    if (raw.NDaughters() != 2 || !FiniteParticle(raw) || raw.GetNDF() <= 0 ||
        raw.GetChi2() < 0.f) return false;
    const kfp::KFParticle* proton = 0;
    const kfp::KFParticle* pion = 0;
    const int sign = raw.GetPDG() > 0 ? 1 : -1;
    for (int daughter = 0; daughter < 2; ++daughter) {
      const int particleIndex = raw.DaughterIds()[daughter];
      if (particleIndex < 0 || static_cast<size_t>(particleIndex) >= particles.size())
        return false;
      const kfp::KFParticle& particle = particles[particleIndex];
      if (particle.NDaughters() != 1) return false;
      if (particle.GetPDG() == 2212 * sign && particle.GetQ() == sign)
        proton = &particle;
      else if (particle.GetPDG() == -211 * sign && particle.GetQ() == -sign)
        pion = &particle;
      else return false;
    }
    if (!proton || !pion) return false;
    // A single-daughter wrapper particle's DaughterIds() can be empty (same condition
    // documented for the Xi bachelor path, see ResolveXiDaughters()); guard identically
    // rather than indexing an empty vector.
    if (proton->DaughterIds().empty() || pion->DaughterIds().empty()) return false;
    out.protonId = proton->DaughterIds()[0];
    out.pionId = pion->DaughterIds()[0];
    const std::map<int, TrackPid>::const_iterator protonInput = trackPid.find(out.protonId);
    const std::map<int, TrackPid>::const_iterator pionInput = trackPid.find(out.pionId);
    if (out.protonId == out.pionId || protonInput == trackPid.end() ||
        pionInput == trackPid.end()) return false;
    const TrackPid& protonPid = protonInput->second;
    const TrackPid& pionPid = pionInput->second;
    out.protonIndex = protonPid.index;
    out.pionIndex = pionPid.index;
    if (cuts.selectionProfile == "lambda_imp5") {
      out.protonDcaToPv = protonPid.dcaToPv;
      out.pionDcaToPv = pionPid.dcaToPv;
      helper.FillImp5PathLengths(dst, out);
    }
    if (out.protonIndex == out.pionIndex) return false;
    out.protonPidPull = protonPid.tpc[2];
    out.pionPidPull = pionPid.tpc[0];
    out.protonHasTof = protonPid.hasTof;
    out.pionHasTof = pionPid.hasTof;
    out.protonTofPull = protonPid.tof[2];
    out.pionTofPull = pionPid.tof[0];
    out.protonTofM2 = protonPid.m2;
    out.pionTofM2 = pionPid.m2;
    out.pdg = raw.GetPDG();
    out.x = raw.GetX(); out.y = raw.GetY(); out.z = raw.GetZ();
    out.px = raw.GetPx(); out.py = raw.GetPy(); out.pz = raw.GetPz();
    out.chi2Ndf = raw.GetChi2() / raw.GetNDF();
    if (raw.GetMass(out.mass, out.massError) != 0 || out.massError <= 0.f) return false;
    out.daughterDistance = proton->GetDistanceFromParticle(*pion);
    out.distanceToPv = raw.GetDistanceFromVertex(pv);
    out.protonChi2Prim = proton->GetDeviationFromVertex(pv);
    out.pionChi2Prim = pion->GetDeviationFromVertex(pv);
    out.protonDistanceToPv = proton->GetDistanceFromVertex(pv);
    out.pionDistanceToPv = pion->GetDistanceFromVertex(pv);
    raw.GetDistanceToVertexLine(pv, out.vertexLineLength, out.vertexLineLengthError);
    if (out.vertexLineLengthError <= 0.f) return false;
    out.vertexLineLengthSignificance = out.vertexLineLength / out.vertexLineLengthError;
    // Never constrain the parent mass used for the spectrum. Constrain only a
    // distinct copy to the input production vertex for topology/lifetime QA.
    kfp::KFParticle topological = raw;
    topological.SetProductionVertex(pv);
    if (!FiniteParticle(topological) || topological.GetNDF() <= 0 ||
        topological.GetChi2() < 0.f) return false;
    out.topoChi2Ndf = topological.GetChi2() / topological.GetNDF();
    if (topological.GetDecayLength(out.decayLength, out.decayLengthError) != 0 ||
        out.decayLengthError <= 0.f) return false;
    out.decayLengthSignificance = out.decayLength / out.decayLengthError;
    const TVector3 flight(out.x - picoPv.X(), out.y - picoPv.Y(), out.z - picoPv.Z());
    const TVector3 momentum(out.px, out.py, out.pz);
    const double denominator = flight.Mag() * momentum.Mag();
    if (!std::isfinite(denominator) || denominator <= 0.) return false;
    out.cosPointing = std::max(-1., std::min(1., flight.Dot(momentum) / denominator));
    const float values[] = {out.mass, out.massError, out.chi2Ndf, out.topoChi2Ndf,
      out.daughterDistance, out.distanceToPv, out.vertexLineLength,
      out.vertexLineLengthError, out.vertexLineLengthSignificance, out.decayLength,
      out.decayLengthError, out.decayLengthSignificance, out.cosPointing,
      out.protonChi2Prim, out.pionChi2Prim, out.protonDistanceToPv, out.pionDistanceToPv};
    for (size_t i = 0; i < sizeof(values) / sizeof(values[0]); ++i)
      if (!std::isfinite(values[i])) return false;
    return true;
  }

  // K0S → π+ π-. Reuses KfLambdaCandidate: proton* = π+, pion* = π-.
  bool K0shortCandidate(const kfp::KFParticle& raw,
                        const std::vector<kfp::KFParticle>& particles,
                        const kfp::KFVertex& pv, const TVector3& picoPv,
                        KfLambdaCandidate& out) {
    if (raw.NDaughters() != 2 || !FiniteParticle(raw) || raw.GetNDF() <= 0 ||
        raw.GetChi2() < 0.f || raw.GetPDG() != 310) return false;
    const kfp::KFParticle* pionPlus = 0;
    const kfp::KFParticle* pionMinus = 0;
    for (int daughter = 0; daughter < 2; ++daughter) {
      const int particleIndex = raw.DaughterIds()[daughter];
      if (particleIndex < 0 || static_cast<size_t>(particleIndex) >= particles.size())
        return false;
      const kfp::KFParticle& particle = particles[particleIndex];
      if (particle.NDaughters() != 1) return false;
      if (particle.GetPDG() == 211 && particle.GetQ() > 0)
        pionPlus = &particle;
      else if (particle.GetPDG() == -211 && particle.GetQ() < 0)
        pionMinus = &particle;
      else return false;
    }
    if (!pionPlus || !pionMinus) return false;
    // See the identical guard/comment in Candidate() above.
    if (pionPlus->DaughterIds().empty() || pionMinus->DaughterIds().empty()) return false;
    out.protonId = pionPlus->DaughterIds()[0];
    out.pionId = pionMinus->DaughterIds()[0];
    const std::map<int, TrackPid>::const_iterator plusInput = trackPid.find(out.protonId);
    const std::map<int, TrackPid>::const_iterator minusInput = trackPid.find(out.pionId);
    if (out.protonId == out.pionId || plusInput == trackPid.end() ||
        minusInput == trackPid.end()) return false;
    const TrackPid& plusPid = plusInput->second;
    const TrackPid& minusPid = minusInput->second;
    out.protonIndex = plusPid.index;
    out.pionIndex = minusPid.index;
    if (out.protonIndex == out.pionIndex) return false;
    out.protonPidPull = plusPid.tpc[0];
    out.pionPidPull = minusPid.tpc[0];
    out.protonHasTof = plusPid.hasTof;
    out.pionHasTof = minusPid.hasTof;
    out.protonTofPull = plusPid.tof[0];
    out.pionTofPull = minusPid.tof[0];
    out.protonTofM2 = plusPid.m2;
    out.pionTofM2 = minusPid.m2;
    out.pdg = raw.GetPDG();
    out.x = raw.GetX(); out.y = raw.GetY(); out.z = raw.GetZ();
    out.px = raw.GetPx(); out.py = raw.GetPy(); out.pz = raw.GetPz();
    out.chi2Ndf = raw.GetChi2() / raw.GetNDF();
    if (raw.GetMass(out.mass, out.massError) != 0 || out.massError <= 0.f) return false;
    out.daughterDistance = pionPlus->GetDistanceFromParticle(*pionMinus);
    out.distanceToPv = raw.GetDistanceFromVertex(pv);
    out.protonChi2Prim = pionPlus->GetDeviationFromVertex(pv);
    out.pionChi2Prim = pionMinus->GetDeviationFromVertex(pv);
    out.protonDistanceToPv = pionPlus->GetDistanceFromVertex(pv);
    out.pionDistanceToPv = pionMinus->GetDistanceFromVertex(pv);
    raw.GetDistanceToVertexLine(pv, out.vertexLineLength, out.vertexLineLengthError);
    if (out.vertexLineLengthError <= 0.f) return false;
    out.vertexLineLengthSignificance = out.vertexLineLength / out.vertexLineLengthError;
    kfp::KFParticle topological = raw;
    topological.SetProductionVertex(pv);
    if (!FiniteParticle(topological) || topological.GetNDF() <= 0 ||
        topological.GetChi2() < 0.f) return false;
    out.topoChi2Ndf = topological.GetChi2() / topological.GetNDF();
    if (topological.GetDecayLength(out.decayLength, out.decayLengthError) != 0 ||
        out.decayLengthError <= 0.f) return false;
    out.decayLengthSignificance = out.decayLength / out.decayLengthError;
    const TVector3 flight(out.x - picoPv.X(), out.y - picoPv.Y(), out.z - picoPv.Z());
    const TVector3 momentum(out.px, out.py, out.pz);
    const double denominator = flight.Mag() * momentum.Mag();
    if (!std::isfinite(denominator) || denominator <= 0.) return false;
    out.cosPointing = std::max(-1., std::min(1., flight.Dot(momentum) / denominator));
    const float values[] = {out.mass, out.massError, out.chi2Ndf, out.topoChi2Ndf,
      out.daughterDistance, out.distanceToPv, out.vertexLineLength,
      out.vertexLineLengthError, out.vertexLineLengthSignificance, out.decayLength,
      out.decayLengthError, out.decayLengthSignificance, out.cosPointing,
      out.protonChi2Prim, out.pionChi2Prim, out.protonDistanceToPv, out.pionDistanceToPv};
    for (size_t i = 0; i < sizeof(values) / sizeof(values[0]); ++i)
      if (!std::isfinite(values[i])) return false;
    return true;
  }

  bool IsBachelorPion(const kfp::KFParticle& part, int sign) const {
    return part.NDaughters() == 1 && part.GetPDG() == -211 * sign &&
           part.GetQ() == -sign;
  }

  bool IsLambdaParticle(const kfp::KFParticle& part, int sign) const {
    return part.NDaughters() == 2 && part.GetPDG() == 3122 * sign;
  }

  // Prefer Finder convention: both DaughterIds are particles[] indices (as for Λ).
  // Fallback: one slot is Λ index, the other is bachelor track id.
  bool ResolveXiDaughters(const kfp::KFParticle& raw,
                          const std::vector<kfp::KFParticle>& particles, int sign,
                          const kfp::KFParticle*& lambda,
                          const kfp::KFParticle*& bachelorPart,
                          int& bachelorTrackId) const {
    lambda = 0;
    bachelorPart = 0;
    bachelorTrackId = -1;
    if (raw.NDaughters() != 2) return false;
    const int id0 = raw.DaughterIds()[0];
    const int id1 = raw.DaughterIds()[1];
    const bool id0Ok = id0 >= 0 && static_cast<size_t>(id0) < particles.size();
    const bool id1Ok = id1 >= 0 && static_cast<size_t>(id1) < particles.size();

    // Primary: both DaughterIds are particles[] indices (Construct uses Daughter.Id()).
    if (id0Ok && id1Ok) {
      const kfp::KFParticle& p0 = particles[id0];
      const kfp::KFParticle& p1 = particles[id1];
      // Accept Λ by PDG (±3122) or by 2-daughter topology if PDG was remapped.
      const bool lam0 = (std::abs(p0.GetPDG()) == 3122 && p0.NDaughters() == 2) ||
                        (p0.NDaughters() == 2 && p1.NDaughters() == 1 &&
                         std::abs(p0.GetPDG()) != 211 && std::abs(p0.GetPDG()) != 321);
      const bool lam1 = (std::abs(p1.GetPDG()) == 3122 && p1.NDaughters() == 2) ||
                        (p1.NDaughters() == 2 && p0.NDaughters() == 1 &&
                         std::abs(p1.GetPDG()) != 211 && std::abs(p1.GetPDG()) != 321);
      const bool bach0 = p0.NDaughters() == 1;
      const bool bach1 = p1.NDaughters() == 1;
      if (lam0 && bach1) {
        lambda = &p0;
        bachelorPart = &p1;
      } else if (lam1 && bach0) {
        lambda = &p1;
        bachelorPart = &p0;
      } else if (bach0 && bach1) {
        // Both look like tracks: one id may already be a track id — try PDG.
        if (IsBachelorPion(p0, sign) && IsLambdaParticle(p1, sign)) {
          bachelorPart = &p0;
          lambda = &p1;
        } else if (IsBachelorPion(p1, sign) && IsLambdaParticle(p0, sign)) {
          bachelorPart = &p1;
          lambda = &p0;
        }
      }
      if (lambda && bachelorPart) {
        bachelorTrackId = bachelorPart->DaughterIds().empty()
                              ? -1
                              : bachelorPart->DaughterIds()[0];
        // If bachelor Id itself is the track (NDaughters==1 but DaughterIds empty),
        // fall back to particle Id as track key used in some snapshots.
        if (bachelorTrackId < 0) bachelorTrackId = bachelorPart->Id();
        return bachelorTrackId >= 0;
      }
    }

    // Legacy fallback: Λ via PDG at a particles[] index; other id = track id.
    int lambdaIndex = -1;
    int trackId = -1;
    if (id1Ok && (IsLambdaParticle(particles[id1], sign) ||
                  (particles[id1].NDaughters() == 2 &&
                   std::abs(particles[id1].GetPDG()) == 3122))) {
      lambdaIndex = id1;
      trackId = id0;
    } else if (id0Ok && (IsLambdaParticle(particles[id0], sign) ||
                         (particles[id0].NDaughters() == 2 &&
                          std::abs(particles[id0].GetPDG()) == 3122))) {
      lambdaIndex = id0;
      trackId = id1;
    }
    if (lambdaIndex < 0 || trackId < 0) return false;
    lambda = &particles[lambdaIndex];
    for (size_t ip = 0; ip < particles.size(); ++ip) {
      const kfp::KFParticle& part = particles[ip];
      if (part.NDaughters() == 1 && !part.DaughterIds().empty() &&
          part.DaughterIds()[0] == trackId) {
        bachelorPart = &part;
        bachelorTrackId = trackId;
        return true;
      }
      // Also match particle Id == trackId (Id set to Particles.size at creation).
      if (part.NDaughters() == 1 && part.Id() == trackId) {
        bachelorPart = &part;
        bachelorTrackId =
            part.DaughterIds().empty() ? trackId : part.DaughterIds()[0];
        return true;
      }
    }
    // Last resort: bachelor KFParticle may not be in the vector; DCA(Λ,bach)
    // needs a particle. Cannot recover without it.
    return false;
  }

  bool XiCandidate(const kfp::KFParticle& raw,
                   const std::vector<kfp::KFParticle>& particles,
                   const kfp::KFVertex& pv, const TVector3& picoPv,
                   KfXiCandidate& out, StPicoDst* dst, bool& daughtersRecovered) {
    // Parent kinematics always; optional daughter recovery fills DCA / Λ L.
    daughtersRecovered = false;
    if (!FiniteParticle(raw) || raw.GetNDF() <= 0 || raw.GetChi2() < 0.f) return false;
    const int sign = raw.GetPDG() > 0 ? 1 : -1;
    out = KfXiCandidate();
    out.pdg = raw.GetPDG();
    out.lambdaPdg = 3122 * sign;
    out.x = raw.GetX(); out.y = raw.GetY(); out.z = raw.GetZ();
    out.px = raw.GetPx(); out.py = raw.GetPy(); out.pz = raw.GetPz();
    out.chi2Ndf = raw.GetChi2() / raw.GetNDF();
    if (raw.GetMass(out.mass, out.massError) != 0 || out.massError <= 0.f) return false;
    out.distanceToPv = raw.GetDistanceFromVertex(pv);
    raw.GetDistanceToVertexLine(pv, out.vertexLineLength, out.vertexLineLengthError);
    if (out.vertexLineLengthError <= 0.f) return false;
    out.vertexLineLengthSignificance = out.vertexLineLength / out.vertexLineLengthError;
    kfp::KFParticle topological = raw;
    topological.SetProductionVertex(pv);
    if (!FiniteParticle(topological) || topological.GetNDF() <= 0 ||
        topological.GetChi2() < 0.f) return false;
    out.topoChi2Ndf = topological.GetChi2() / topological.GetNDF();
    if (topological.GetDecayLength(out.decayLength, out.decayLengthError) != 0 ||
        out.decayLengthError <= 0.f) return false;
    out.decayLengthSignificance = out.decayLength / out.decayLengthError;
    const TVector3 flight(out.x - picoPv.X(), out.y - picoPv.Y(), out.z - picoPv.Z());
    const TVector3 momentum(out.px, out.py, out.pz);
    const double denominator = flight.Mag() * momentum.Mag();
    if (!std::isfinite(denominator) || denominator <= 0.) return false;
    out.cosPointing = std::max(-1., std::min(1., flight.Dot(momentum) / denominator));
    out.daughterDistance = -1.f;
    out.lambdaMass = 0.f;
    out.lambdaMassError = -1.f;
    out.lambdaDecayLength = -1.f;
    out.lambdaDecayLengthSignificance = -1.f;
    out.lambdaDaughterDistance = -1.f;

    if (dst) {
      const kfp::KFParticle* lambda = 0;
      const kfp::KFParticle* bachelorPart = 0;
      int bachelorTrackId = -1;
      if (ResolveXiDaughters(raw, particles, sign, lambda, bachelorPart, bachelorTrackId) &&
          lambda && bachelorPart && FiniteParticle(*bachelorPart)) {
        float unused = 0.f;
        if (lambda->GetMass(out.lambdaMass, unused) == 0) {
          out.lambdaMassError = unused;
        }
        kfp::KFParticle lambdaTopo = *lambda;
        lambdaTopo.SetProductionVertex(pv);
        float lambdaLErr = 0.f;
        if (FiniteParticle(lambdaTopo) && lambdaTopo.GetNDF() > 0 &&
            lambdaTopo.GetChi2() >= 0.f &&
            lambdaTopo.GetDecayLength(out.lambdaDecayLength, lambdaLErr) == 0 &&
            lambdaLErr > 0.f) {
          out.lambdaDecayLengthSignificance = out.lambdaDecayLength / lambdaLErr;
        } else {
          out.lambdaDecayLength = -1.f;
          out.lambdaDecayLengthSignificance = -1.f;
        }
        out.daughterDistance = lambda->GetDistanceFromParticle(*bachelorPart);
        if (!std::isfinite(out.daughterDistance)) out.daughterDistance = -1.f;

        out.bachelorChi2Prim = bachelorPart->GetDeviationFromVertex(pv);
        out.bachelorDistanceToPv = bachelorPart->GetDistanceFromVertex(pv);
        if (!std::isfinite(out.bachelorChi2Prim)) out.bachelorChi2Prim = -1.f;
        if (!std::isfinite(out.bachelorDistanceToPv)) out.bachelorDistanceToPv = -1.f;

        const kfp::KFParticle* proton = 0;
        const kfp::KFParticle* pion = 0;
        if (lambda->NDaughters() == 2) {
          for (int daughter = 0; daughter < 2; ++daughter) {
            const int particleIndex = lambda->DaughterIds()[daughter];
            if (particleIndex < 0 || static_cast<size_t>(particleIndex) >= particles.size())
              continue;
            const kfp::KFParticle& particle = particles[particleIndex];
            if (particle.NDaughters() != 1) continue;
            if (particle.GetPDG() == 2212 * sign && particle.GetQ() == sign)
              proton = &particle;
            else if (particle.GetPDG() == -211 * sign && particle.GetQ() == -sign)
              pion = &particle;
          }
        }
        if (proton && pion) {
          out.lambdaDaughterDistance = proton->GetDistanceFromParticle(*pion);
          if (!std::isfinite(out.lambdaDaughterDistance)) out.lambdaDaughterDistance = -1.f;
          out.protonChi2Prim = proton->GetDeviationFromVertex(pv);
          out.pionChi2Prim = pion->GetDeviationFromVertex(pv);
          out.protonDistanceToPv = proton->GetDistanceFromVertex(pv);
          out.pionDistanceToPv = pion->GetDistanceFromVertex(pv);
          if (!std::isfinite(out.protonChi2Prim)) out.protonChi2Prim = -1.f;
          if (!std::isfinite(out.pionChi2Prim)) out.pionChi2Prim = -1.f;
          if (!std::isfinite(out.protonDistanceToPv)) out.protonDistanceToPv = -1.f;
          if (!std::isfinite(out.pionDistanceToPv)) out.pionDistanceToPv = -1.f;
          const std::map<int, TrackPid>::const_iterator bachelorInput =
              trackPid.find(bachelorTrackId);
          if (bachelorInput != trackPid.end() &&
              !proton->DaughterIds().empty() && !pion->DaughterIds().empty()) {
            out.protonId = proton->DaughterIds()[0];
            out.pionId = pion->DaughterIds()[0];
            out.bachelorId = bachelorTrackId;
            const std::map<int, TrackPid>::const_iterator protonInput =
                trackPid.find(out.protonId);
            const std::map<int, TrackPid>::const_iterator pionInput =
                trackPid.find(out.pionId);
            if (protonInput != trackPid.end() && pionInput != trackPid.end()) {
              out.protonIndex = protonInput->second.index;
              out.pionIndex = pionInput->second.index;
              out.bachelorIndex = bachelorInput->second.index;
              out.protonPidPull = protonInput->second.tpc[2];
              out.pionPidPull = pionInput->second.tpc[0];
              out.bachelorPidPull = bachelorInput->second.tpc[0];
            }
          }
        }
        daughtersRecovered =
            out.daughterDistance >= 0.f && out.lambdaDaughterDistance >= 0.f;
      }
    }

    const float values[] = {out.mass, out.massError, out.chi2Ndf, out.topoChi2Ndf,
      out.distanceToPv, out.vertexLineLength, out.vertexLineLengthError,
      out.vertexLineLengthSignificance, out.decayLength, out.decayLengthError,
      out.decayLengthSignificance, out.cosPointing};
    for (size_t i = 0; i < sizeof(values) / sizeof(values[0]); ++i)
      if (!std::isfinite(values[i])) return false;
    return true;
  }

  bool ProcessEvent(StPicoDst* dst) {
    Clear();
    if (!configurationError.empty()) return Fail(configurationError);
    std::ostringstream errors;
    if (!cuts.Validate(errors)) return Fail(errors.str());
    if (!dst || !dst->picoArray(StPicoArrays::Event) || !dst->event())
      return Fail("PicoDst or Event is missing");
    if (!dst->picoArray(StPicoArrays::Track) ||
        !dst->picoArray(StPicoArrays::TrackCovMatrix))
      return Fail("Track and TrackCovMatrix branch arrays are required");
    if (cuts.useTof && !dst->picoArray(StPicoArrays::BTofPidTraits))
      return Fail("The enabled TOF profile requires the BTofPidTraits branch array");
    stats.rawTracks = dst->numberOfTracks();
    if (dst->numberOfTrackCovMatrices() != stats.rawTracks)
      return Fail("TrackCovMatrix count differs from Track count; the covariance branch is required");
    if (cuts.useTof && (cuts.tofMean.size() != 30 || cuts.tofSigma.size() != 30))
      return Fail("KF TOF calibration requires 30 mean and 30 width coefficients");
    const StPicoEvent* event = dst->event();
    const TVector3 position = event->primaryVertex();
    const TVector3 uncertainty = event->primaryVertexError();
    if (!FiniteVector(position) || !FiniteVector(uncertainty) ||
        uncertainty.X() <= 0. || uncertainty.Y() <= 0. || uncertainty.Z() <= 0. ||
        !std::isfinite(event->bField()))
      return Fail("Input primary vertex, its errors, or the magnetic field is invalid");
    // Set the current event's field BEFORE track-to-PV deviations are computed.
    // It is in PicoDst/STAR units (kG), with no additional Tesla conversion.
    topo.SetField(event->bField());
    kfp::KFPVertex vertex;
    vertex.SetXYZ(position.X(), position.Y(), position.Z());
    vertex.SetCovarianceMatrix(uncertainty.X() * uncertainty.X(), 0.f,
        uncertainty.Y() * uncertainty.Y(), 0.f, 0.f, uncertainty.Z() * uncertainty.Z());
    const kfp::KFVertex pv(vertex);
    std::vector<kfp::KFParticle> inputs;
    std::vector<int> pdgs, hftHits, primaryTracks;
    inputs.reserve(stats.rawTracks * 3);
    pdgs.reserve(stats.rawTracks * 3);
    hftHits.reserve(stats.rawTracks * 3);
    std::map<int, int> originalIds;
    for (UInt_t index = 0; index < stats.rawTracks; ++index) {
      const StPicoTrack* track = dst->track(index);
      if (!track) return Fail("Track branch contains a null entry");
      // IDs need not be contiguous or equal to the PicoDst array index.
      if (!originalIds.insert(std::make_pair(track->id(), index)).second) {
        ++stats.invalidTrackIds;
        return Fail("Duplicate StPicoTrack::id(): cannot unambiguously recover daughters");
      }
      // Diagnostic only: raw kinematics before the Quality() cut (minPt/maxPt/
      // minEta/maxEta), see the Impl members' comment above.
      {
        const TVector3 rawP = track->gMom();
        if (FiniteVector(rawP) && rawP.Mag2() > 0. && std::isfinite(rawP.Eta())) {
          rawTrackPt.push_back(static_cast<Float_t>(rawP.Pt()));
          rawTrackEta.push_back(static_cast<Float_t>(rawP.Eta()));
        }
      }
      if (!Quality(*track)) continue;
      const int nHft = static_cast<int>(track->hasPxl1Hit()) +
          static_cast<int>(track->hasPxl2Hit()) + static_cast<int>(track->hasIstHit());
      if (cuts.useHftTracksOnly && nHft < 3) continue;
      ++stats.qualityTracks;
      kfp::KFPTrack input;
      if (!helper.BuildTrack(dst, index, input)) { ++stats.invalidCovariance; continue; }
      ++stats.covarianceTracks;
      TrackPid pid;
      pid.index = index;
      std::vector<int> hypotheses;
      if (!Pid(dst, *track, input.GetP(), pid, hypotheses)) {
        ++stats.invalidPid;
        if (!error.empty()) return Fail(error);
        continue;
      }
      trackPid.insert(std::make_pair(track->id(), pid));
      if (hypotheses[0] != -1) ++stats.pidTracks;
      for (size_t h = 0; h < hypotheses.size(); ++h) {
        const int pdg = hypotheses[h];
        kfp::KFParticle particle(input, pdg);
        const float chiPrimary = particle.GetDeviationFromVertex(pv);
        if (!FiniteParticle(particle) || !std::isfinite(chiPrimary) || chiPrimary < 0.f) {
          ++stats.invalidPid;
          continue;
        }
        // As in reference AddTrackToParticleList: do NOT square this API's
        // return value. The opt-in Imp5 profile separately selects original
        // daughter gDCA in PID; the default reference has no such prefilter.
        if (chiPrimary < cuts.interfaceChiPrimaryCut) {
          primaryTracks.push_back(inputs.size());
          ++stats.primaryHypotheses;
        }
        particle.SetId(track->id());
        inputs.push_back(particle);
        pdgs.push_back(pdg);
        hftHits.push_back(nHft);
        ++stats.pidHypotheses;
        if (std::abs(pdg) == 2212) ++stats.protonHypotheses;
        if (std::abs(pdg) == 211) ++stats.pionHypotheses;
        if (pdg == -1) ++stats.unknownHypotheses;
      }
    }
    // Same STAR Interface sequence, with the external Pico PV. No PV refit,
    // custom p x pi loop, extra competition pass, or scalar fallback is performed.
    // Upstream ReconstructParticles itself runs SelectParticleCandidates:
    // Lambda must pass the configured PV-copy chi2/NDF cut (default 3).
    topo.CleanPV();
    topo.Init(inputs, &pdgs, &hftHits);
    topo.GetKFParticleFinder()->Init(topo.NPrimaryVertices());
    topo.FillPVIndices();
    topo.AddPV(pv, primaryTracks);
    topo.FillPVIndices();
    topo.SortTracks();
    topo.ReconstructParticles();
    const std::vector<kfp::KFParticle>& particles = topo.GetParticles();
    stats.finderParticles = particles.size();
    for (size_t i = 0; i < particles.size(); ++i) {
      const int pdg = particles[i].GetPDG();
      if (std::abs(pdg) == 3122) {
        ++stats.lambdaParticles;
        KfLambdaCandidate candidate;
        if (!Candidate(particles[i], particles, pv, position, candidate, dst)) {
          ++stats.invalidCandidates;
          continue;
        }
        candidates.push_back(candidate);
      } else if (std::abs(pdg) == 3312) {
        ++stats.xiParticles;
        KfXiCandidate candidate;
        bool daughtersRecovered = false;
        if (!XiCandidate(particles[i], particles, pv, position, candidate, dst,
                         daughtersRecovered)) {
          ++stats.invalidXiCandidates;
          continue;
        }
        if (daughtersRecovered) ++stats.xiDaughterRecovered;
        xiCandidates.push_back(candidate);
      } else if (pdg == 310) {
        ++stats.k0Particles;
        KfLambdaCandidate candidate;
        if (!K0shortCandidate(particles[i], particles, pv, position, candidate)) {
          ++stats.invalidK0Candidates;
          continue;
        }
        k0Candidates.push_back(candidate);
      }
    }
    stats.validCandidates = candidates.size();
    stats.validXiCandidates = xiCandidates.size();
    stats.validK0Candidates = k0Candidates.size();
    return true;
  }

  const KfParticleCutConfig& cuts;
  KfParticleHelper helper;
  kfp::KFParticleTopoReconstructor topo;
  KfParticleEventStats stats;
  std::vector<KfLambdaCandidate> candidates;
  std::vector<KfXiCandidate> xiCandidates;
  std::vector<KfLambdaCandidate> k0Candidates;
  std::map<int, TrackPid> trackPid;
  std::string error;
  std::string configurationError;
  // Diagnostic only: raw StPicoTrack gPt/gEta filled before the Quality()
  // (minPt/maxPt/minEta/maxEta) cut so Makers can histogram what Stage A rejects.
  // No selection change; filled regardless of Quality() outcome.
  std::vector<Float_t> rawTrackPt;
  std::vector<Float_t> rawTrackEta;
};

StPicoKFParticleInterface::StPicoKFParticleInterface(const KfParticleCutConfig& cuts)
    : mImpl(new Impl(cuts)) {}
StPicoKFParticleInterface::~StPicoKFParticleInterface() { delete mImpl; }
Bool_t StPicoKFParticleInterface::ProcessEvent(StPicoDst* dst) { return mImpl->ProcessEvent(dst); }
void StPicoKFParticleInterface::Clear() { mImpl->Clear(); }
const std::vector<KfLambdaCandidate>& StPicoKFParticleInterface::Candidates() const {
  return mImpl->candidates;
}
const std::vector<KfXiCandidate>& StPicoKFParticleInterface::XiCandidates() const {
  return mImpl->xiCandidates;
}
const std::vector<KfLambdaCandidate>& StPicoKFParticleInterface::K0shortCandidates() const {
  return mImpl->k0Candidates;
}
const KfParticleEventStats& StPicoKFParticleInterface::Stats() const { return mImpl->stats; }
const std::string& StPicoKFParticleInterface::LastError() const { return mImpl->error; }
const std::vector<Float_t>& StPicoKFParticleInterface::RawTrackPt() const { return mImpl->rawTrackPt; }
const std::vector<Float_t>& StPicoKFParticleInterface::RawTrackEta() const { return mImpl->rawTrackEta; }

const char* StPicoKFParticleInterface::BackendDescription() {
  return "finder_topo; core=.DEV2 sha256=" STAR_ANALYZER_KFP_SOURCE_ID
      "; snapshot=" STAR_ANALYZER_KFP_SNAPSHOT_DATE
      "; namespace=star_analyzer_kfp; __ROOT__,KFParticleStandalone,HomogeneousField,SSE4.1"
      "; Pico source=xwu2/KFTree_lambda/StKFParticleInterface"
      "; local=lambda-study-controls-v1; Lambda PV chi2/NDF and Finder geometry use effective config"
      "; Lambda, optional Xi (+/-3312), optional K0S (310); parent mass unconstrained";
}
