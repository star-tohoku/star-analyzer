#ifndef StFemtoMaker_h
#define StFemtoMaker_h

#include "StMaker.h"
#include "StPhiKKReconstruction.h"
#include "FemtoCandidate.h"
#include "cuts/FemtoConfig.h"
#include "StarClassLibrary/StPhysicalHelixD.hh"

#include <deque>
#include <map>
#include <string>
#include <vector>

class StPicoDst;
class StPicoDstMaker;
class StPicoEvent;
class StPicoTrack;
class TString;
class HistManager;
class CentralityHelper;
class TVector3;

class StFemtoMaker;

StFemtoMaker* createStFemtoMaker(const char* name, StPicoDstMaker* picoMaker, const char* outName);
extern "C" void* createStFemtoMakerC(const char* name, void* picoMaker, const char* outName);

class StFemtoMaker : public StMaker {
 public:
  StFemtoMaker(const char* name, StPicoDstMaker* picoMaker, const char* outName);
  virtual ~StFemtoMaker();

  virtual Int_t Init();
  virtual Int_t Make();
  virtual void Clear(Option_t* opt = "");
  virtual Int_t Finish();

  void WriteHistograms();

 private:
  struct TrackState {
    Float_t pT, eta, phi;
    Short_t charge;
    Float_t nSigmaKaon;
    Float_t nSigmaProton;
    Float_t DCA;
    Short_t nHitsFit, nHitsMax, nHitsDedx;
    Float_t chi2;
    Bool_t tofMatch;
    Float_t mass2;
    Float_t deltaOneOverBeta;
    Float_t originX, originY, originZ;
    Float_t momentumX, momentumY, momentumZ;
    Float_t BField;
    Int_t trackIndex;
  };

  struct FemtoMixingEvent {
    FemtoCandidateStore candidates;
  };

  StPicoDstMaker* mPicoDstMaker;
  StPicoDst* mPicoDst;
  TString mOutName;
  Int_t mEventCounter;
  HistManager* m_histManager;
  CentralityHelper* m_centrality;
  Int_t m_cent9;
  Int_t m_cent16;
  Double_t m_refMultCorr;
  Double_t m_centWeight;
  Double_t m_centralityPercent;
  Double_t m_psi2;

  FemtoCandidateStore m_eventCandidates;
  std::map<Int_t, std::deque<FemtoMixingEvent> > m_mixingPool;

  Bool_t PassEventCuts(Float_t vz, Float_t vr, Int_t refMult, Float_t vzVpd);
  Bool_t PassTrackCuts(StPicoTrack* trk, TVector3& pVtx);
  Bool_t PassKaonCuts(StPicoTrack* trk, TVector3& pVtx);
  Bool_t PassProtonCuts(StPicoTrack* trk, TVector3& pVtx);

  static PhiKkTrackState ToPhiKkTrack(const TrackState& trk);
  void BuildTrackState(TrackState& track, StPicoTrack* pico, StPicoEvent* event, TVector3& pVtx, Int_t index);
  void FillTofInfo(TrackState& track, StPicoTrack* trk, const TVector3& pMom, Int_t btofIndex);
  Bool_t PassTofKaonPid(const TrackState& trk) const;
  Bool_t PassTofProtonPid(const TrackState& trk) const;
  Bool_t IsKaon(const TrackState& trk);
  Bool_t IsProton(const TrackState& trk);

  TVector3 TrackMomentum(const TrackState& trk) const;
  StPhysicalHelixD BuildHelix(const TrackState& trk);
  Double_t CalculateDCA(const TrackState& trk1, const TrackState& trk2, TVector3& dcaPos1, TVector3& dcaPos2);
  Bool_t ReconstructPhi(const TrackState& kPlus, const TrackState& kMinus, Double_t& invMass, TVector3& phiMom,
                        TVector3& dcaPosPlus, TVector3& dcaPosMinus);
  Double_t CalculateOpeningAngle(const TrackState& trk1, const TrackState& trk2);
  Double_t CalculatePairRapidity(Double_t invMass, const TVector3& phiMom);
  Double_t ApplyRapidityFrame(Double_t yLab) const;
  Bool_t PassPairTofCut(const TrackState& kPlus, const TrackState& kMinus) const;

  void BuildTrackPidCandidates(const std::string& speciesKey, const std::string& particleKey,
                               const std::vector<TrackState>& tracks, Int_t eventIndex);
  void BuildResonanceCandidates(const std::string& speciesKey, const std::string& particleKey,
                                const std::vector<TrackState>& kaonsPlus, const std::vector<TrackState>& kaonsMinus,
                                Int_t eventIndex);
  FemtoCandidate MakeProtonCandidate(const TrackState& trk, Int_t eventIndex, const std::string& speciesKey) const;
  FemtoCandidate MakePhiCandidate(const TrackState& kPlus, const TrackState& kMinus, Double_t invMass,
                                  const TVector3& phiMom, Double_t openingAngle, Double_t pairRapidity,
                                  Double_t dcaKK, Int_t eventIndex, const std::string& speciesKey) const;

  Double_t ComputeKStar(const TLorentzVector& pA, const TLorentzVector& pB) const;
  TLorentzVector ProtonP4(const TVector3& p) const;
  TLorentzVector CandidateP4(const FemtoCandidate& cand) const;
  Bool_t TracksOverlap(const FemtoCandidate& phiCand, const FemtoCandidate& trkCand) const;

  Int_t GetMixingBin(Float_t vz, Int_t cent9, Double_t psi2) const;
  void FillSameEventPairs(const FemtoConfig::ChannelDef& ch);
  void FillMixedEventPairs(const FemtoConfig::ChannelDef& ch, Float_t vz, Int_t cent9, Double_t psi2);
  void StoreEventForMixing(Float_t vz, Int_t cent9, Double_t psi2);
  void FillCandidateQA();
  void FillCentralityEventQA(Int_t cent9, Int_t rawMult, Double_t refMultCorr, Int_t nTracks, Int_t nBTOFMatch,
                             Int_t nKaonPlus, Int_t nKaonMinus, Int_t nPhiCandidates, Int_t nProtons);
  void FinalizeCorrelationFunctions();
  std::string HistName(const std::string& prefix, const std::string& channelName) const;
};

#endif
