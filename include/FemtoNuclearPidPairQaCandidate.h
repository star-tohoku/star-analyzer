#ifndef FEMTO_NUCLEAR_PID_PAIR_QA_CANDIDATE_H
#define FEMTO_NUCLEAR_PID_PAIR_QA_CANDIDATE_H
#include "FemtoCandidate.h"

// Independent diagnostic payload; never admits an event to the physical pool.
struct FemtoNuclearPidPairQaCandidate {
  FemtoCandidate candidate;
  // Original Double-precision PID result before storing Float kinematics/pulls.
  bool selected;
  FemtoNuclearPidPairQaCandidate() : selected(false) {}
};
typedef std::map<std::string, std::vector<FemtoNuclearPidPairQaCandidate> >
    FemtoNuclearPidPairQaStore;
#endif
