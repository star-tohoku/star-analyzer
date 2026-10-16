#ifndef TRACK_PID_CANDIDATE_SOURCE_H
#define TRACK_PID_CANDIDATE_SOURCE_H

#include "CandidateSource.h"

// Identified-track part for StCandidateFemtoMaker (see include/CandidateSource.h): charged
// tracks selected by charge, kinematics, DCA, nHitsFit, TPC nSigma and (optionally) TOF. Every
// particle key and its cuts come from the YAML referenced by the mainconf key "trackPid":
//   particleKeys: kaon_plus_tof
//   kaon_plus_tof_charge: 1                 # +1 or -1
//   kaon_plus_tof_minPt: 0.1                # global pT >= minPt [GeV/c]
//   kaon_plus_tof_maxP: 1.2                 # global |p| < maxP [GeV/c]
//   kaon_plus_tof_maxDca: 3.0               # global DCA to the PV < maxDca [cm]
//   kaon_plus_tof_minNHitsFit: 15
//   kaon_plus_tof_nSigmaHypothesis: kaon    # pion | kaon | proton | electron
//   kaon_plus_tof_maxAbsNSigma: 3.0         # |nSigma| < maxAbsNSigma
//   kaon_plus_tof_requireTof: true          # TOF match with beta > 0
// massValue = TOF m^2 = p^2 (1/beta^2 - 1) from the global momentum (-999 without TOF). The
// momentum used for k* is the primary momentum for primary tracks, the global one otherwise.
CandidateSource* createTrackPidCandidateSource();

#endif
