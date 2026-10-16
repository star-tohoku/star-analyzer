#ifndef CANDIDATE_SOURCE_H
#define CANDIDATE_SOURCE_H

#include "Rtypes.h"
#include <string>
#include <vector>

class HistManager;
class StPicoDst;

// Reconstruction-tool-neutral interfaces between particle reconstruction parts ("sources")
// and the femto pairing Maker (StCandidateFemtoMaker). Nothing here depends on KFParticle or
// on any particular reconstruction method: a source turns one PicoDst event into lists of
// candidates per species key, and the Maker only sees these plain structs. The candidate store
// (StMaker/StCandidateFemtoMaker/CandidateStore.h) replays stored candidates through the same
// interface, so the pair stage is identical with and without the PicoDst.

// Kinematics of one final-state track as the PicoDst stores it (global momentum at the DCA to the
// PV, charge). The Maker fills it from the PicoDst; the candidate store keeps it so that the pair
// stage can run without the PicoDst. A missing or unusable track is all zeros.
struct FemtoDaughterTrack {
  FemtoDaughterTrack() : gx(0.f), gy(0.f), gz(0.f), charge(0) {}
  Float_t gx, gy, gz;
  Int_t charge;
};

// One reconstructed particle candidate (composite or track).
struct FemtoSourceCandidate {
  FemtoSourceCandidate()
      : px(0.f), py(0.f), pz(0.f), massValue(0.f), x(0.f), y(0.f), z(0.f), decayLength(-1.f) {}
  // Momentum used for k* [GeV/c].
  Float_t px, py, pz;
  // Mass-like variable of this species: invariant mass for composites, TOF m^2 for tracks.
  Float_t massValue;
  // Decay vertex [cm] (composites; the primary vertex for tracks) and decay length [cm].
  Float_t x, y, z;
  Float_t decayLength;
  // PicoDst track indices of the final-state tracks (a track candidate lists itself).
  // Order follows CandidateSource::DaughterNames().
  std::vector<Int_t> daughters;
  // Values of CandidateSource::VariableNames() (BDT inputs and training-tree branches); exactly
  // one value per name (the Maker stops otherwise). BDT scores are computed by the Maker, not here.
  std::vector<Float_t> variables;
  // One entry per daughter (daughters order). Reconstruction parts leave it empty (the Maker reads
  // the PicoDst); the candidate store fills it.
  std::vector<FemtoDaughterTrack> daughterTracks;
};

// A particle reconstruction part. The Maker owns the sources it is given, calls Init once,
// then ProcessEvent / Clear once per accepted event, and Finish at the end of the job.
class CandidateSource {
public:
  virtual ~CandidateSource() {}
  // Short label for logs, e.g. "kfparticle", "track_pid".
  virtual const char* Name() const = 0;
  // mainConfigPath: the mainconf of this job; hist: the Maker's HistManager (QA histograms
  // missing from the hist YAML are skipped). Returns kFALSE on any configuration error.
  virtual Bool_t Init(const char* mainConfigPath, HistManager* hist) = 0;
  // Species keys this source can produce (e.g. "xi_minus", "k0s").
  virtual std::vector<std::string> SpeciesKeys() const = 0;
  virtual Bool_t ProvidesSpecies(const std::string& key) const = 0;
  // Final-state track names of a species, in FemtoSourceCandidate::daughters order.
  virtual std::vector<std::string> DaughterNames(const std::string& key) const = 0;
  // Names of FemtoSourceCandidate::variables of a species.
  virtual std::vector<std::string> VariableNames(const std::string& key) const = 0;
  // Reconstruct the current event. cent9 is the event-selection centrality bin. kFALSE = error.
  virtual Bool_t ProcessEvent(StPicoDst* picoDst, Int_t cent9) = 0;
  virtual const std::vector<FemtoSourceCandidate>& Candidates(const std::string& key) const = 0;
  virtual void Clear() = 0;
  // Called before the Maker writes its histograms (fill counters into histograms, print summary).
  virtual void Finish() {}
  virtual const std::string& LastError() const = 0;
};

// Event-level selection part (bad runs, triggers, vertex, pile-up, centrality).
class FemtoEventSelector {
public:
  virtual ~FemtoEventSelector() {}
  virtual const char* Name() const = 0;
  virtual Bool_t Init(const char* mainConfigPath, HistManager* hist) = 0;
  // kTRUE if the event is accepted; cent9 is set (0 = most peripheral bin, 8 = most central).
  virtual Bool_t Accept(StPicoDst* picoDst, Int_t& cent9) = 0;
  virtual void Finish() {}
  // Human-readable dump of the effective configuration (written to the output file).
  virtual std::string Describe() const = 0;
};

#endif
