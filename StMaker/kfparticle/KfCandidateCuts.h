#ifndef KF_CANDIDATE_CUTS_H
#define KF_CANDIDATE_CUTS_H

#include "Rtypes.h"

struct KfXiCandidate;
struct KfLambdaCandidate;
class KfParticleCutConfig;

// Stage B candidate selections for Xi- and K0short (rectangular cuts of the kf YAML), applied by the
// KFParticle reconstruction part (KfParticleCandidateSource).
Bool_t PassKfXiCandidateCuts(const KfXiCandidate& c, const KfParticleCutConfig& k);
Bool_t PassKfK0CandidateCuts(const KfLambdaCandidate& c, const KfParticleCutConfig& k);

// K0S Lambda veto: inputs are pi+/pi- momenta at the K0S decay vertex [GeV/c] and
// |nSigmaProton|. Returns the smallest |M(p pi) - M_Lambda| [GeV] over hypotheses
// enabled by the config (pi+ as p; pi- as pbar only when k0LambdaVetoAntiLambda),
// skipping a daughter whose |nSigmaProton| is not below k0LambdaVetoMaxNSigmaProton;
// returns a large value if no hypothesis applies.
Double_t KfK0LambdaVetoDeviation(Double_t pxPlus, Double_t pyPlus, Double_t pzPlus,
                                 Double_t pxMinus, Double_t pyMinus, Double_t pzMinus,
                                 Double_t nSigmaProtonPlus, Double_t nSigmaProtonMinus,
                                 const KfParticleCutConfig& k);
// kTRUE if the veto is disabled or the deviation is outside the veto window.
Bool_t PassKfK0LambdaVeto(Double_t deviation, const KfParticleCutConfig& k);

#endif
