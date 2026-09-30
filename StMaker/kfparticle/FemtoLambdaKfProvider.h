#ifndef FEMTO_LAMBDA_KF_PROVIDER_H
#define FEMTO_LAMBDA_KF_PROVIDER_H

#include "FemtoLambdaProvider.h"
struct KfLambdaCandidate;

// Only Lambda runners link this factory from libStKfParticleCommon. The
// general StFemtoMaker depends on the neutral abstract API, not this factory.
FemtoLambdaProvider* createFemtoLambdaKfProvider();

// Pure value conversion used by the provider and deterministic closure tests.
// No selection or mass constraint is introduced here.
bool ConvertKfLambdaForFemto(const KfLambdaCandidate& source, int eventIndex,
                            double pairMass, FemtoLambdaCandidate& output,
                            std::string& error);
#endif
