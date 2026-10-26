#ifndef FEMTO_NUCLEAR_PAIR_PID_SPARSE_H
#define FEMTO_NUCLEAR_PAIR_PID_SPARSE_H

#include "THnSparse.h"
#include <iosfwd>
#include <string>

class TDirectory;

// Optional, additive pair QA. No candidate selection or mixing is performed here.
class FemtoNuclearPairPidSparse {
 public:
  FemtoNuclearPairPidSparse();
  ~FemtoNuclearPairPidSparse();
  bool Init(const char* yamlPath, std::ostream& diagnostics);
  // representation: he4=0, he3=1, he3as4=2. massClass is exclusive:
  // SBNeg=0, signalLow=1, signalHigh=2, SBPos=3.
  void Fill(double kstar, double rigidity, double pull, int cent9,
            int representation, bool mixed, int massClass, bool prePid);
  THnSparseD* Histogram() const { return mHistogram; }
  const std::string& Metadata() const { return mMetadata; }
  bool Write(TDirectory* directory) const;

 private:
  FemtoNuclearPairPidSparse(const FemtoNuclearPairPidSparse&);
  FemtoNuclearPairPidSparse& operator=(const FemtoNuclearPairPidSparse&);
  THnSparseD* mHistogram;
  std::string mMetadata;
  unsigned long long mAccepted;
  unsigned long long mRejectedNonfinite;
  unsigned long long mRejectedCategory;
};

#endif
