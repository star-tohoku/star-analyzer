#ifndef KF_PARTICLE_CANDIDATE_SOURCE_H
#define KF_PARTICLE_CANDIDATE_SOURCE_H

#include "CandidateSource.h"

// KFParticle reconstruction part for StCandidateFemtoMaker (see include/CandidateSource.h).
// Runs StPicoKFParticleInterface once per event (one Finder pass for every particle enabled
// in the kf YAML) and returns the Stage-B selected candidates (KfCandidateCuts) per particle
// key. Particle keys:
//   xi_minus_kf  Xi- -> Lambda pi-   (kf: reconstructXi)       daughters proton, pion, bachelor
//   xi_plus_kf   Xi+ (anti)          (kf: reconstructAntiXi)   daughters proton, pion, bachelor (the antiproton,
//                                                              pi+ and pi+ of the Xi+; same names as Xi-)
//   k0s_kf       K0S -> pi+ pi-      (kf: reconstructK0short)  daughters pionPlus, pionMinus
// To add a particle: enable it in the interface/kf YAML, add a key here (SpeciesKeys,
// DaughterNames, VariableNames, ProcessEvent).
// Created by the analysis macro and handed to the Maker, which owns it.
CandidateSource* createKfParticleCandidateSource();

#endif
