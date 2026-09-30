#ifndef STAR_ANALYZER_KF_LAMBDA_SELECTOR_H
#define STAR_ANALYZER_KF_LAMBDA_SELECTOR_H

#include "KfParticleHelper.h"
#include "cuts/KfParticleCutConfig.h"
#include <cmath>

// Exact final predicate formerly in StLambdaKFParticleMaker. The full Pico
// adapter already rejects non-finite fits before these additional cuts. Keep
// inequalities/disabled thresholds unchanged for standalone/Femto closure.
inline Bool_t PassKfLambdaCandidateCuts(const KfLambdaCandidate& c,
                                       const KfParticleCutConfig& k) {
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

#endif
