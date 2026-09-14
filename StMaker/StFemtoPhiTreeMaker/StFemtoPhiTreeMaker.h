#ifndef StFemtoPhiTreeMaker_h
#define StFemtoPhiTreeMaker_h

#include "StMaker.h"
#include "StPhiKKReconstruction.h"
#include "FemtoPhiTreeSchema.h"
#include <string>
#include <vector>

class StPicoDst;
class StPicoDstMaker;
class StPicoEvent;
class StPicoTrack;
class TFile;
class TTree;
class TString;
class CentralityHelper;
class TVector3;

class StFemtoPhiTreeMaker;

StFemtoPhiTreeMaker* createStFemtoPhiTreeMaker(const char* name, StPicoDstMaker* picoMaker,
                                               const char* outName);
extern "C" void* createStFemtoPhiTreeMakerC(const char* name, void* picoMaker, const char* outName);

class StFemtoPhiTreeMaker : public StMaker {
 public:
  StFemtoPhiTreeMaker(const char* name, StPicoDstMaker* picoMaker, const char* outName);
  virtual ~StFemtoPhiTreeMaker();

  virtual Int_t Init();
  virtual Int_t Make();
  virtual void Clear(Option_t* opt = "");
  virtual Int_t Finish();

 private:
  struct TrackState {
    Float_t pT, eta, phi;
    Short_t charge;
    Float_t nSigmaKaon;
    Float_t nSigmaPion;
    Float_t nSigmaProton;
    Float_t nSigmaDeuteron;
    Float_t dEdx;
    Float_t DCA;
    Short_t nHitsFit, nHitsMax, nHitsDedx;
    Float_t chi2;
    Char_t tofMatch;
    Float_t mass2;
    Float_t deltaOneOverBeta;
    Float_t tofBeta;
    Float_t originX, originY, originZ;
    Float_t momentumX, momentumY, momentumZ;
    Float_t BField;
    Int_t trackIndex;
    UInt_t selFlags;
  };

  StPicoDstMaker* mPicoDstMaker;
  StPicoDst* mPicoDst;
  TString mOutName;
  TFile* mOutFile;
  TTree* mEventTree;
  TTree* mTrackTree;
  TTree* mPairTree;
  CentralityHelper* m_centrality;

  femto_phi_tree::EventRow mEvt;
  femto_phi_tree::TrackRow mTrk;
  femto_phi_tree::PhiPairRow mPair;
  femto_phi_tree::EventRowV2 mEvt2;
  femto_phi_tree::TrackRowV2 mTrk2;
  femto_phi_tree::PackStats mPackStats;

  Bool_t mWriteEvent;
  Bool_t mWriteTrack;
  Bool_t mWritePhiPair;
  Bool_t mStoreProtons;
  Bool_t mStoreDeuterons;
  Bool_t mStoreKaons;
  Int_t mCompressLevel;
  Int_t mAutoFlush;
  UInt_t mSchemaVersion;
  UInt_t mSubjobId;
  std::vector<UInt_t> mTriggerIds;          // configured triggers, from YAML
  std::vector<std::string> mSourceFiles;    // sourceFileIndex -> PicoDst path
  std::string mPidCorrectionState;
  Long64_t mSkipRemaining;
  std::string mMakerYamlPath;
  std::string mMainconfPath;

  Double_t mEnvMinNHitsRatio;
  Int_t mEnvMinNHitsFit;
  Int_t mEnvMinNHitsDedx;
  Double_t mEnvMaxDca;
  Double_t mEnvMinPt;
  Double_t mEnvMaxPt;
  Double_t mEnvMinEta;
  Double_t mEnvMaxEta;
  Double_t mEnvMaxChi2;
  Double_t mEnvMaxAbsNSigmaKaon;
  Double_t mEnvMaxAbsNSigmaDeuteron;
  Double_t mEnvMaxAbsNSigmaProton;
  Double_t mEnvMaxDcaKaon;
  Double_t mEnvDeuteronMaxDca;
  Double_t mEnvDeuteronMinPMom;
  Double_t mEnvDeuteronMaxPMom;
  Double_t mEnvProtonMaxDca;
  Double_t mEnvProtonMinPt;
  Int_t mEnvMinNHitsDedxNuclear;

  Long64_t mNInput;
  Long64_t mNNoPico;
  Long64_t mNBadRun;
  Long64_t mNFailEventCuts;
  Long64_t mNPileup;
  Long64_t mNFailCent;
  Long64_t mNFailMaxNTr;
  Long64_t mNAccepted;
  Long64_t mNKp;
  Long64_t mNKm;
  Long64_t mNDeuteron;
  Long64_t mNProton;
  Long64_t mNPhiPair;

  Bool_t LoadTreeConfig();
  void BookTrees();
  void BookTreesV1();
  void BookTreesV2();
  UShort_t SourceFileIndex(const TString& path);
  void FillTriggerInfo(StPicoEvent* event);
  void WriteSourceFileTable();
  void WriteMetadata();
  void WriteCounterHistogram();
  Bool_t PassEventCuts(Float_t vz, Float_t vr, Int_t refMult, Float_t vzVpd) const;
  Bool_t PassEnvTrackCuts(StPicoTrack* trk, TVector3& pVtx) const;
  Bool_t PassNominalTrackCuts(StPicoTrack* trk, TVector3& pVtx) const;
  void BuildTrackState(TrackState& track, StPicoTrack* pico, StPicoEvent* event, TVector3& pVtx,
                       Int_t index);
  void FillTofInfo(TrackState& track, StPicoTrack* trk, const TVector3& pMom, Int_t btofIndex);
  static PhiKkTrackState ToPhiKkTrack(const TrackState& trk);
  Bool_t PassTofKaonPid(const TrackState& trk) const;
  Bool_t PassPhiDaughterTofPid(const TrackState& trk) const;
  Bool_t PassTofProtonPid(const TrackState& trk) const;
  Bool_t PassNominalKaonCuts(StPicoTrack* trk, TVector3& pVtx) const;
  Bool_t PassFemtoDeuteronCuts(const TrackState& trk) const;
  Bool_t PassFemtoProtonCuts(const TrackState& trk) const;
  Double_t ApplyRapidityFrame(Double_t yLab) const;
  Double_t DeuteronRapidityCm(const TrackState& trk) const;
  Double_t ProtonRapidityCm(const TrackState& trk) const;
  TVector3 TrackMomentum(const TrackState& trk) const;
  void FillTrackTree(const TrackState& src, UChar_t species, ULong64_t eventUID);
  void FillTrackTreeV2(const TrackState& src, UChar_t species, ULong64_t eventUID);
  void FillOptionalPhiPairs(const std::vector<TrackState>& kp, const std::vector<TrackState>& km,
                            ULong64_t eventUID);
  void MixBins(Float_t vz, Int_t cent9, Double_t psi2, Int_t& vzBin, Int_t& centBin, Int_t& epBin,
               Int_t& mixBin) const;
  TString CurrentPicoPath() const;
  Long64_t CurrentEntry() const;
};

#endif
