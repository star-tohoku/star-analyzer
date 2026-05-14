#ifndef StPhiMaker_h
#define StPhiMaker_h

#include "StMaker.h"
#include "StarClassLibrary/StPhysicalHelixD.hh"

#include <deque>
#include <map>
#include <vector>

class StPicoDst;
class StPicoDstMaker;
class StPicoEvent;
class StPicoTrack;
class StPicoBTofPidTraits;
class TString;
class HistManager;
class TVector3;

class StPhiMaker;

// Factory for CINT (new StPhiMaker fails with "declared but not defined")
StPhiMaker* createStPhiMaker(const char* name, StPicoDstMaker* picoMaker, const char* outName);

// Extern "C" wrapper for CINT symbol resolution
extern "C" void* createStPhiMakerC(const char* name, void* picoMaker, const char* outName);

class StPhiMaker : public StMaker {
 public:
  StPhiMaker(const char* name, StPicoDstMaker* picoMaker, const char* outName);
  virtual ~StPhiMaker();

  virtual Int_t Init();
  virtual Int_t Make();
  virtual void Clear(Option_t* opt = "");
  virtual Int_t Finish();

  void WriteHistograms();

 private:
  StPicoDstMaker* mPicoDstMaker;
  StPicoDst* mPicoDst;
  TString mOutName;
  Int_t mEventCounter;
  HistManager* m_histManager;

  // Track structure for KK pair reconstruction
  struct Track_t {
    Float_t pT, eta, phi;
    Short_t charge;
    Float_t nSigmaKaon;
    Float_t DCA;
    Short_t nHitsFit, nHitsMax, nHitsDedx;
    Float_t chi2;
    Bool_t tofMatch;
    Float_t mass2;
    Float_t deltaOneOverBeta;
    Float_t originX, originY, originZ;
    Float_t momentumX, momentumY, momentumZ;
    Float_t BField;
  };

  struct PhiMixingEvent {
    std::vector<Track_t> kaonsPlus;
    std::vector<Track_t> kaonsMinus;
  };

  std::map<Int_t, std::deque<PhiMixingEvent> > m_phiMixingPool;

  Bool_t PassTofKaonPid(const Track_t& trk) const;
  void FillTofInfo(Track_t& track, StPicoTrack* trk, const TVector3& pMom, Int_t btofIndex);
  Int_t GetMixingVzBin(Float_t vz) const;
  void FillMixedEventPairs(const std::vector<Track_t>& kaonsPlus, const std::vector<Track_t>& kaonsMinus, Float_t vz);
  void StoreEventForMixing(const std::vector<Track_t>& kaonsPlus, const std::vector<Track_t>& kaonsMinus, Float_t vz);
  void FillPhiPairHistograms(Double_t invMass, const TVector3& phiMom, Double_t openingAngle, Double_t pairRapidity,
                             Int_t nKaonPlus, Int_t nKaonMinus, Bool_t applySignalCuts);
  void FinalizeBackgroundSubtractedHistogram();

  // Helper methods
  Bool_t PassEventCuts(Float_t vz, Float_t vr, Int_t refMult, Float_t vzVpd);
  Bool_t PassTrackCuts(StPicoTrack* trk, TVector3& pVtx);
  Bool_t PassKaonCuts(StPicoTrack* trk, TVector3& pVtx);
  Bool_t PassTrackCuts(const Track_t& trk);
  Bool_t IsKaon(const Track_t& trk);
  void BuildTrack(Track_t& track, StPicoTrack* pico, StPicoEvent* event, TVector3& pVtx);
  StPhysicalHelixD BuildHelix(const Track_t& trk);
  Double_t CalculateDCA(const Track_t& trk1, const Track_t& trk2, TVector3& dcaPos1, TVector3& dcaPos2);
  Bool_t ReconstructPhi(const Track_t& kPlus, const Track_t& kMinus, Double_t& invMass, TVector3& phiMom, TVector3& dcaPosPlus, TVector3& dcaPosMinus);
  Double_t CalculateInvariantMass(const Track_t& trk1, const Track_t& trk2, Double_t mass1, Double_t mass2);
  Double_t CalculateOpeningAngle(const Track_t& trk1, const Track_t& trk2);
  Double_t CalculatePairRapidity(Double_t invMass, const TVector3& phiMom);
};

#endif
