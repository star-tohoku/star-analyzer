#ifndef ST_CANDIDATE_FEMTO_MAKER_H
#define ST_CANDIDATE_FEMTO_MAKER_H

#include "StMaker.h"
#include "TString.h"

class CandidateSource;
class FemtoEventSelector;
class StPicoDstMaker;
class StCandidateFemtoMakerImpl;

class StCandidateFemtoMaker;
StCandidateFemtoMaker* createStCandidateFemtoMaker(const char* name, StPicoDstMaker* picoMaker,
                                                   const char* outName);

// Femtoscopic pairing of candidates from reconstruction parts. The Maker does not know how a
// candidate was reconstructed (KFParticle, track PID, ...): the analysis macro hands it one
// FemtoEventSelector and any number of CandidateSource parts (include/CandidateSource.h), and
// the YAML of the mainconf key "candidateFemto" (CandidateFemtoConfig.h) declares species and
// channels. Per event: event selection -> sources -> BDT scores (mainconf key "candidateBdt",
// optional) -> same-event pairs -> mixed-event pairs against the pool -> store the event.
// Output: THnSparseF hPairSE_<channel> / hPairME_<channel> (k*, cent9, mass classes, masses,
// BDT scores) plus QA histograms from the hist YAML (mainconf key "<MakerName>Hist").
// Windows, sidebands, BDT thresholds and purity corrections are chosen in the checkHist step.
//
// Two stages with one PicoDst pass (CandidateStore.h): with writeCandidateStore: true the Maker is the
// reconstruction stage (event selection and parts on the PicoDst, every candidate of every part written with the
// event records into the job output; no pairs). SetCandidateStoreInput makes it the pair stage: no PicoDst,
// no event selector, no parts; the stored parts are replayed through the same CandidateSource interface, so the
// output equals that of one PicoDst pass over the same files (reconstruction-stage histograms are added from
// the stores). Without either, one PicoDst pass does everything.
class StCandidateFemtoMaker : public StMaker {
public:
  StCandidateFemtoMaker(const char* name, StPicoDstMaker* picoMaker, const char* outName);
  virtual ~StCandidateFemtoMaker();
  virtual Int_t Init();
  virtual Int_t Make();
  virtual void Clear(Option_t* option = "");
  virtual Int_t Finish();

  void SetMainConfigPath(const char* path) { mMainConfigPath = path ? path : ""; }
  // Ownership passes to the Maker.
  void SetEventSelector(FemtoEventSelector* selector);
  void AddCandidateSource(CandidateSource* source);
  // Pair stage: text file listing candidate stores (reconstruction-stage outputs), read in order. picoMaker must
  // be 0 and no event selector or source may be given.
  void SetCandidateStoreInput(const char* storeList) { mStoreList = storeList ? storeList : ""; }
  Bool_t CandidateStoreInputFinished() const;
  Long64_t GetAcceptedEvents() const;
  Long64_t GetSameEventPairs() const;
  Long64_t GetMixedEventPairs() const;

private:
  StPicoDstMaker* mPicoDstMaker;
  TString mOutName;
  TString mMainConfigPath;
  TString mStoreList;
  StCandidateFemtoMakerImpl* mImpl;

  StCandidateFemtoMaker(const StCandidateFemtoMaker&);
  StCandidateFemtoMaker& operator=(const StCandidateFemtoMaker&);
};

#endif
