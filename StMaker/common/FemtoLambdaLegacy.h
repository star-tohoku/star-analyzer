#ifndef FEMTO_LAMBDA_LEGACY_H
#define FEMTO_LAMBDA_LEGACY_H

#include "FemtoCandidate.h"
#include "FemtoNuclearPidPairQaCandidate.h"
#include "FemtoLambdaPairGeometry.h"
#include "TVector3.h"
#include <iosfwd>
#include <map>
#include <string>
#include <vector>

class HistManager;
class StPicoDst;
class StPicoTrack;
class TDirectory;
struct FemtoLambdaCandidate;
class FemtoNuclearDedxCalibration;
class FemtoNuclearPairPidSparse;

// Lambda-only compatibility support. It owns its histograms; no KF/SIMD linkage.
// StFemtoMaker owns event acceptance, centrality, pair calculation and mixing.
class FemtoLambdaLegacy {
 public:
  FemtoLambdaLegacy();
  ~FemtoLambdaLegacy();
  bool Init(const char* mainconf, std::ostream& diagnostics, bool producerQaOnly=false,
            bool collectPrePidControls=true);
  HistManager* RootHistograms() const { return mRoot; }
  HistManager* NuclearHistograms() const { return mNuclear; }
  // Optional selected/pre-TPC-PID pair QA. Null for historical mainconfs.
  HistManager* NuclearPidPairQaHistograms() const { return mNuclearPidPairQa; }
  FemtoNuclearPairPidSparse* NuclearPidPairSparse() const { return mNuclearPidPairSparse; }
  // This flag controls the pre-PID bank, not the selected-pair histogram view.
  bool NuclearPidPairQaEnabled() const { return mCollectNuclearPidPairQa; }
  const std::string& SpeciesKey() const { return mSpecies; }
  bool IsPairNucleus(int type) const;
  void FillNucleusAcceptance(const FemtoCandidate&, int type);
  double NuclearDedxPull(int type, double rigidity, double dedx) const;
  const std::string& NuclearDedxCalibrationMetadata() const;
  void CollectNuclei(StPicoDst*, int eventIndex, int cent9, FemtoCandidateStore&,
                     FemtoNuclearPidPairQaStore* pidQa=0);
  void FillLambda(const FemtoLambdaCandidate&, int cent9, double refmultcorr,
                  const TVector3& primaryVertex);
  // Signal-candidate acceptance, not weighted by the number of nuclear pairs.
  bool FillLambdaAcceptance(const FemtoLambdaCandidate&, StPicoDst*);
  struct TrackGeometry {
    int trackIndex, charge;
    TVector3 globalMomentum, origin;
    TrackGeometry() : trackIndex(-1), charge(0) {}
  };
  // Owned current-event geometry; original Pico array indices remain the keys.
  void SetReplayTrackGeometry(const std::vector<TrackGeometry>&, double fieldKilogauss);
  void ClearReplayTrackGeometry();
  bool FillLambdaAcceptanceFromReplay(const FemtoLambdaCandidate&);
  bool PairGeometryEnabled() const { return mPairGeometryConfigured&&!mProducerQaOnly; }
  void BuildPairGeometry(int eventIndex,const FemtoCandidateStore&,StPicoDst*,
                        femto_lambda_pair_geometry::EventGeometry&) const;
  // Diagnostic EP label of the event currently being processed, not the
  // Lambda event in both mixed-event directions. Negative means unavailable.
  void SetPairEventPlaneBin(int bin) { mPairEventPlaneBin=bin; }
  int PairEventPlaneBins() const { return mPairEventPlaneNBins; }
  bool FillPair(const FemtoCandidate& lambda, const FemtoCandidate& nucleus,
                double kstar, double qlab, int cent9, bool mixed,
                StPicoDst* currentDst = 0,
                const femto_lambda_pair_geometry::EventGeometry* lambdaGeometry=0,
                const femto_lambda_pair_geometry::EventGeometry* nucleusGeometry=0,
                bool* rejectedByMerging=0);
  bool FillPrePidNuclearPairQa(const FemtoCandidate& lambda, const FemtoCandidate& nucleus,
                              double kstar, int cent9, bool mixed);
  bool Write(TDirectory* destination);

  // 0=outside windows, 1=signal, 2=left sideband, 3=right sideband.
  int MassRegion(double mass) const;
  // 0=not signal, 1=low half (mass<mean), 2=high half (mass>=mean).
  int SignalHalf(double mass) const;
  int SelectNuclearType(double rigidity, const double pulls[4],
                        bool validTof, double mass2) const;
  // TOF uses raw rigidity: mass2 is (m/q)^2, before the Z=2 P4 correction.
  bool PassHe4HighRigidityTof(double rigidity, bool validTof, double mass2) const;
  bool PassNuclearHighRigidityTof(int type, double rigidity,
                                bool validTof, double mass2) const;
  bool PassNuclearPidQaDomain(int type, double rigidity, bool validTof, double mass2) const;
  TLorentzVector NuclearP4(const TVector3& rawRigidity, int type) const;
  // First forward cylinder intersection of the STAR global-track helix.
  static bool MergingTrackPositionAtRadius(const StPicoTrack*, double fieldKilogauss,
                                          double radiusMeters, double tpcHalfLengthCm,
                                          double radiusToleranceCm, TVector3& position);
  static bool MergingTrackPositionAtRadius(const TrackGeometry&, double fieldKilogauss,
                                          double radiusMeters, double tpcHalfLengthCm,
                                          double radiusToleranceCm, TVector3& position);

 private:
  FemtoLambdaLegacy(const FemtoLambdaLegacy&);
  FemtoLambdaLegacy& operator=(const FemtoLambdaLegacy&);
  FemtoCandidate MapHe3AsHe4(const FemtoCandidate&) const;
  bool FillPairHistograms(const FemtoCandidate&, const FemtoCandidate&, double kstar,
                         double qlab, int cent9, bool mixed, StPicoDst*,
                         const std::string& species, const std::string& legacy, int massIndex,
                         int selectedNuclearType);
  void FillNuclearPidPairQa(const FemtoCandidate& lambda, const FemtoCandidate& nucleus,
                           double kstar, int cent9, bool mixed,
                           const std::string& legacy, int selectedNuclearType, bool prePid);
  void FillMerging(const FemtoCandidate&, const FemtoCandidate&, StPicoDst*,
                   const std::string& legacy, int massIndex);
  const TrackGeometry* ReplayTrack(int originalIndex) const;
  bool AcceptMergingPair(const FemtoCandidate&,const FemtoCandidate&,bool mixed,
                        const femto_lambda_pair_geometry::EventGeometry*,
                        const femto_lambda_pair_geometry::EventGeometry*);
  bool FillLambdaAcceptanceKinematics(const FemtoLambdaCandidate&, const TVector3&,
                                     const TVector3&);
  int He4LabRapidityRegion(double yLab) const;
  void FillHe4RapiditySplit(const FemtoCandidate&, double kstar, int cent9, bool mixed,
                           const std::string& massSuffix);
  void FillHe4RapidityBinning(const FemtoCandidate&, double kstar, int cent9, bool mixed,
                             const std::string& massSuffix);
  void FillHe4PairRapidity(const TLorentzVector&, double kstar, int cent9, bool mixed,
                           const std::string& massSuffix);
  void FillPairEventPlane(const std::string& species, double kstar, int cent9, bool mixed,
                          const std::string& massSuffix);
  bool ValidateHistograms(std::ostream&);
  bool He4TofQaEnabled() const;
  double He4TofThreshold() const;
  HistManager* mRoot;
  HistManager* mNuclear;
  HistManager* mNuclearPidPairQa;
  FemtoNuclearPairPidSparse* mNuclearPidPairSparse;
  FemtoNuclearDedxCalibration* mDedxCalibration;
  bool mProducerQaOnly, mCollectNuclearPidPairQa, mReplayTrackGeometryEnabled;
  bool mPairGeometryConfigured, mMergingPairCutEnabled, mMergingPairQaEnabled;
  double mMergingPairPhiMax, mMergingPairEtaMax;
  double mReplayFieldKilogauss;
  std::vector<TrackGeometry> mReplayTracks;
  std::map<int, std::vector<TrackGeometry>::size_type> mReplayTrackSlots;
  std::vector<std::string> mRootKeys, mNuclearKeys;
  std::vector<std::string> mNuclearPidPairQaKeys;
  std::string mSpecies, mLegacySpecies;
  int mSpeciesIndex;
  bool mHe3ControlEnabled;
  std::vector<int> mPairSpeciesIndices;
  double mMean, mWindow, mOuterFactor, mFieldTesla, mRadiusMeters;
  bool mMergingIncludePion, mHe4RapiditySplitEnabled, mHe4RapidityBinningEnabled;
  bool mPairRapiditySplitEnabled;
  bool mPairEventPlaneSplitEnabled;
  double mMergingTpcHalfLengthCm, mMergingRadiusToleranceCm, mHe4RapidityThreshold;
  int mHe4RapidityNBins, mPairRapidityNBins, mPairEventPlaneBin, mPairEventPlaneNBins;
  double mHe4RapidityMin, mHe4RapidityMax, mPairRapidityMin, mPairRapidityMax;
  double mDaughterProtonMass, mDaughterPionMass, mLambdaPairMass;
  int mMinNHitsDedx;
  double mMinPt, mNuclearMass[4], mTofMean[3], mTofSigma[3];
  double mNSigmaFill, mNSigmaExclude, mMaxNSigma, mM2SigmaCut;
  double mMinPTofQa, mMinM2, mMaxM2, mMaxRigidity;
  bool mM2Selection;
  bool mHe4HighRigidityTof;
  double mHe4TofThreshold, mHe4TofMass2Min, mHe4TofMass2Max;
  // Optional CF-domain gates. Absent blocks retain the historical selection.
  bool mSpeciesMinRigidityEnabled, mHighRigidityTof, mExtendHe4BeyondMaxRigidity;
  bool mMinRigidityExclusive[4];
  double mMinRigidity[4], mHighTofThreshold, mHighTofMass2Min[4], mHighTofMass2Max[4];
  struct RigidityVetoRange { double minimum, maximum; };
  bool mSpeciesRigidityVetoConfigured;
  std::vector<RigidityVetoRange> mRigidityVeto[4];
  std::string mPidDomainMetadata;
};

#endif
