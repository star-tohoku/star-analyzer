#include "KfCandidateCuts.h"
#include "KfParticleHelper.h"
#include "cuts/KfParticleCutConfig.h"

#include <algorithm>
#include <cmath>

Bool_t PassKfXiCandidateCuts(const KfXiCandidate& c, const KfParticleCutConfig& k) {
  if (c.mass < k.minMass || c.mass > k.maxMass) return kFALSE;
  if (k.maxMassError >= 0. && c.massError > k.maxMassError) return kFALSE;
  if (k.maxChi2Ndf >= 0. && c.chi2Ndf > k.maxChi2Ndf) return kFALSE;
  if (k.maxTopoChi2Ndf >= 0. && c.topoChi2Ndf > k.maxTopoChi2Ndf) return kFALSE;
  if (k.maxDaughterDistance >= 0. && c.daughterDistance >= 0. &&
      c.daughterDistance > k.maxDaughterDistance) return kFALSE;
  // DCA between the Lambda daughters (p, pi) [cm]; a missing value fails when the cut is enabled.
  if (k.maxLambdaDaughterDistance >= 0.) {
    if (c.lambdaDaughterDistance < 0.) return kFALSE;
    if (c.lambdaDaughterDistance > k.maxLambdaDaughterDistance) return kFALSE;
  }
  if (k.maxDistanceToPv >= 0. && c.distanceToPv > k.maxDistanceToPv) return kFALSE;
  if (k.minDecayLength >= 0. && c.decayLength < k.minDecayLength) return kFALSE;
  if (k.minDecayLengthSignificance >= 0. && c.decayLengthSignificance < k.minDecayLengthSignificance)
    return kFALSE;
  if (k.minVertexLineSignificance >= 0. && c.vertexLineLengthSignificance < k.minVertexLineSignificance)
    return kFALSE;
  if (k.minCosPointing > -1. && c.cosPointing < k.minCosPointing) return kFALSE;
  // requireDecayLengthOrder: decay length order L_Lambda >= L_Xi.
  // Missing Lambda length fails closed when enabled.
  if (k.requireDecayLengthOrder) {
    if (c.lambdaDecayLength < 0.f || c.decayLength < 0.f) return kFALSE;
    if (c.lambdaDecayLength < c.decayLength) return kFALSE;
  }
  // Fake-Lambda / fake-Xi veto (same concept as K0short's daughter cuts).
  if (k.minLambdaDaughterChi2Prim >= 0.) {
    if (c.protonChi2Prim < 0. || c.pionChi2Prim < 0.) return kFALSE;
    if (c.protonChi2Prim < k.minLambdaDaughterChi2Prim ||
        c.pionChi2Prim < k.minLambdaDaughterChi2Prim) return kFALSE;
  }
  if (k.minLambdaDaughterDistanceToPv >= 0.) {
    if (c.protonDistanceToPv < 0. || c.pionDistanceToPv < 0.) return kFALSE;
    if (c.protonDistanceToPv < k.minLambdaDaughterDistanceToPv ||
        c.pionDistanceToPv < k.minLambdaDaughterDistanceToPv) return kFALSE;
  }
  if (k.minBachelorChi2Prim >= 0.) {
    if (c.bachelorChi2Prim < 0.) return kFALSE;
    if (c.bachelorChi2Prim < k.minBachelorChi2Prim) return kFALSE;
  }
  if (k.minBachelorDistanceToPv >= 0.) {
    if (c.bachelorDistanceToPv < 0.) return kFALSE;
    if (c.bachelorDistanceToPv < k.minBachelorDistanceToPv) return kFALSE;
  }
  return kTRUE;
}

Bool_t PassKfK0CandidateCuts(const KfLambdaCandidate& c, const KfParticleCutConfig& k) {
  if (c.mass < k.k0MinMass || c.mass > k.k0MaxMass) return kFALSE;
  if (k.k0MaxMassError >= 0. && c.massError > k.k0MaxMassError) return kFALSE;
  if (k.k0MaxChi2Ndf >= 0. && c.chi2Ndf > k.k0MaxChi2Ndf) return kFALSE;
  if (k.k0MaxTopoChi2Ndf >= 0. && c.topoChi2Ndf > k.k0MaxTopoChi2Ndf) return kFALSE;
  if (k.k0MaxDaughterDistance >= 0. && c.daughterDistance > k.k0MaxDaughterDistance) return kFALSE;
  if (k.k0MaxDistanceToPv >= 0. && c.distanceToPv > k.k0MaxDistanceToPv) return kFALSE;
  if (k.k0MinDecayLength >= 0. && c.decayLength < k.k0MinDecayLength) return kFALSE;
  if (k.k0MinDecayLengthSignificance >= 0. &&
      c.decayLengthSignificance < k.k0MinDecayLengthSignificance) return kFALSE;
  if (k.k0MinVertexLineSignificance >= 0. &&
      c.vertexLineLengthSignificance < k.k0MinVertexLineSignificance) return kFALSE;
  if (k.k0MinCosPointing > -1. && c.cosPointing < k.k0MinCosPointing) return kFALSE;
  if (k.k0MinDaughterChi2Prim >= 0. &&
      (c.protonChi2Prim < k.k0MinDaughterChi2Prim || c.pionChi2Prim < k.k0MinDaughterChi2Prim))
    return kFALSE;
  if (k.k0MinDaughterDistanceToPv >= 0. &&
      (c.protonDistanceToPv < k.k0MinDaughterDistanceToPv ||
       c.pionDistanceToPv < k.k0MinDaughterDistanceToPv))
    return kFALSE;
  // K0 reuses KfLambdaCandidate: proton* = pi+, pion* = pi-.
  if (k.k0MaxPionPlusPidPull >= 0. && std::fabs(c.protonPidPull) > k.k0MaxPionPlusPidPull) return kFALSE;
  if (k.k0MaxPionMinusPidPull >= 0. && std::fabs(c.pionPidPull) > k.k0MaxPionMinusPidPull) return kFALSE;
  return kTRUE;
}

Double_t KfK0LambdaVetoDeviation(Double_t pxPlus, Double_t pyPlus, Double_t pzPlus,
                                 Double_t pxMinus, Double_t pyMinus, Double_t pzMinus,
                                 Double_t nSigmaProtonPlus, Double_t nSigmaProtonMinus,
                                 const KfParticleCutConfig& k) {
  const Double_t kMassProton = 0.938272, kMassPion = 0.139570, kMassLambda = 1.115683;
  const Double_t pp2 = pxPlus * pxPlus + pyPlus * pyPlus + pzPlus * pzPlus;
  const Double_t pm2 = pxMinus * pxMinus + pyMinus * pyMinus + pzMinus * pzMinus;
  const Double_t sx = pxPlus + pxMinus, sy = pyPlus + pyMinus, sz = pzPlus + pzMinus;
  const Double_t s2 = sx * sx + sy * sy + sz * sz;
  Double_t best = 1e9;
  if (std::fabs(nSigmaProtonPlus) < k.k0LambdaVetoMaxNSigmaProton) {
    const Double_t e = std::sqrt(pp2 + kMassProton * kMassProton) + std::sqrt(pm2 + kMassPion * kMassPion);
    best = std::min(best, std::fabs(std::sqrt(std::max(0., e * e - s2)) - kMassLambda));
  }
  if (k.k0LambdaVetoAntiLambda && std::fabs(nSigmaProtonMinus) < k.k0LambdaVetoMaxNSigmaProton) {
    const Double_t e = std::sqrt(pp2 + kMassPion * kMassPion) + std::sqrt(pm2 + kMassProton * kMassProton);
    best = std::min(best, std::fabs(std::sqrt(std::max(0., e * e - s2)) - kMassLambda));
  }
  return best;
}

Bool_t PassKfK0LambdaVeto(Double_t deviation, const KfParticleCutConfig& k) {
  return !k.HasK0LambdaVeto() || deviation >= k.k0LambdaVetoMassWindow;
}
