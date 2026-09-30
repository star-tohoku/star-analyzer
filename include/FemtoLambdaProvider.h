#ifndef FEMTO_LAMBDA_PROVIDER_H
#define FEMTO_LAMBDA_PROVIDER_H

#include "FemtoCandidate.h"
#include <iosfwd>
#include <string>
#include <vector>

class StPicoDst;
class StPicoEvent;

// Value-only bridge: this header has no KFParticle/SIMD or helper-library
// dependency. The Maker owns its injected provider; mixed pools own copies of
// FemtoCandidate, never a provider, Pico pointer, or KF object.
struct FemtoLambdaCandidate {
  FemtoCandidate candidate;
  float x, y, z, mass, massError, chi2Ndf, topoChi2Ndf;
  float daughterDistance, distanceToPv;
  float decayLength, decayLengthError, decayLengthSignificance;
  float vertexLineLength, vertexLineLengthError, vertexLineLengthSignificance;
  float cosPointing, protonPidPull, pionPidPull;
  float protonTofPull, pionTofPull, protonTofM2, pionTofM2;
  double protonDcaToPv, pionDcaToPv;
  double protonHelixPathLength, pionHelixPathLength;
  bool imp5PathValid, protonHasTof, pionHasTof;
  int protonId, pionId, protonIndex, pionIndex, pdg;

  FemtoLambdaCandidate()
      : x(0), y(0), z(0), mass(0), massError(0), chi2Ndf(0), topoChi2Ndf(0),
        daughterDistance(0), distanceToPv(0), decayLength(0), decayLengthError(0),
        decayLengthSignificance(0), vertexLineLength(0), vertexLineLengthError(0),
        vertexLineLengthSignificance(0), cosPointing(0), protonPidPull(0),
        pionPidPull(0), protonTofPull(0), pionTofPull(0), protonTofM2(0), pionTofM2(0),
        protonDcaToPv(0), pionDcaToPv(0), protonHelixPathLength(0), pionHelixPathLength(0),
        imp5PathValid(false), protonHasTof(false), pionHasTof(false),
        protonId(-1), pionId(-1), protonIndex(-1), pionIndex(-1), pdg(0) {}
};

struct FemtoLambdaEventStats {
  unsigned rawTracks, qualityTracks, covarianceTracks, pidTracks;
  unsigned pidHypotheses, protonHypotheses, pionHypotheses, unknownHypotheses;
  unsigned primaryHypotheses, finderParticles, lambdaParticles, validCandidates;
  unsigned invalidCandidates, invalidCovariance, invalidPid, invalidTrackIds;
  unsigned tofTracks, selectedCandidates;
  FemtoLambdaEventStats()
      : rawTracks(0), qualityTracks(0), covarianceTracks(0), pidTracks(0),
        pidHypotheses(0), protonHypotheses(0), pionHypotheses(0), unknownHypotheses(0),
        primaryHypotheses(0), finderParticles(0), lambdaParticles(0), validCandidates(0),
        invalidCandidates(0), invalidCovariance(0), invalidPid(0), invalidTrackIds(0),
        tofTracks(0), selectedCandidates(0) {}
};

class FemtoLambdaProvider {
public:
  virtual ~FemtoLambdaProvider() {}
  virtual bool Init(const char* mainconf, std::ostream& errors) = 0;
  // Event acceptance and centrality are applied once by the Maker. Process
  // reconstructs the full PicoDst; it neither reapplies event cuts nor centrality.
  virtual bool AcceptEvent(const StPicoEvent& event, int nTracks) const = 0;
  virtual double ComputeVr(double x, double y) const = 0;
  virtual bool Process(StPicoDst* picoDst, int eventIndex) = 0;
  virtual void Clear() = 0;
  virtual const std::vector<FemtoLambdaCandidate>& Candidates() const = 0;
  virtual const FemtoLambdaEventStats& Stats() const = 0;
  virtual const std::string& LastError() const = 0;
};

#endif
