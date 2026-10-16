#ifndef OFFICIAL_EVENT_SELECTION_H
#define OFFICIAL_EVENT_SELECTION_H

#include "CandidateSource.h"

// Official Run21 3 GeV event selection as a FemtoEventSelector part (see include/CandidateSource.h):
// bad-run list, trigger ids, vertex/event cuts (KfEventSelection, event YAML), positive
// multiplicities, per-run-batch pile-up band and centrality (cent9 from fxtMult). Settings are the
// official* keys of the kf YAML (KfParticleCutConfig, officialEventSelection: true). It does not use
// KFParticle itself and lives in libStKfParticleCommon only because KfEventSelection does.
// Created by the analysis macro and handed to StCandidateFemtoMaker, which owns it.
FemtoEventSelector* createOfficialEventSelection();

#endif
