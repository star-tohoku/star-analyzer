#ifndef KF_PARTICLE_CUT_CONFIG_H
#define KF_PARTICLE_CUT_CONFIG_H

#include "Rtypes.h"
#include <iosfwd>
#include <string>
#include <vector>

// Plain configuration only: this header must never depend on the KF backend.
class KfParticleCutConfig {
public:
  static KfParticleCutConfig& GetInstance();
  Bool_t LoadFromFile(const Char_t* filename);
  void SetDefaults();
  Bool_t Validate(std::ostream& errors) const;
  void Dump(std::ostream& output) const;

  Int_t schemaVersion;
  std::string backend;
  // Explicit, backward-compatible non-KF selection profile.
  std::string selectionProfile;
  Double_t imp5MinDCAProton;
  Double_t imp5MinDCAPion;
  Double_t imp5MaxPathLength;
  std::string pidProfile;
  std::string tofCalibrationProfile;
  Int_t minNHitsFit;
  Int_t minNHitsDedx;
  Double_t minNHitsRatio;
  Double_t minPt;
  Double_t maxPt;
  Double_t minEta;
  Double_t maxEta;
  Double_t minDedxError;
  Double_t maxDedxError;
  Bool_t useHftTracksOnly;
  Double_t nSigmaProton;
  Double_t nSigmaPion;
  Double_t nSigmaKaon;
  Bool_t useTof;
  Bool_t strictTofPid;
  Bool_t cleanKaonsWithTof;
  Double_t tofNSigma;
  Double_t tofPMax;
  Double_t tofKaonPMin;
  Double_t tofKaonPMax;
  std::vector<double> tofMean;
  std::vector<double> tofSigma;
  Bool_t rejectBadCovariance;
  Double_t maxPositionVariance;
  Double_t maxMomentumVariance;
  Double_t interfaceChiPrimaryCut;
  Double_t primaryProbCut;
  Double_t finderMaxDaughterDistance;
  Double_t finderLCut;
  Double_t finderChiPrimary2D;
  Double_t finderChi2Ndf2D;
  Double_t finderLdL2D;
  Double_t finderLdLXiOmega;
  Double_t finderChi2TopoXiOmega;
  Double_t finderChi2XiOmega;
  Double_t minMass;
  Double_t maxMass;
  Double_t maxMassError;
  Double_t maxChi2Ndf;
  Double_t maxTopoChi2Ndf;
  Double_t maxDaughterDistance;
  Double_t maxLambdaDaughterDistance;
  Double_t maxDistanceToPv;
  Double_t minDecayLength;
  Double_t minDecayLengthSignificance;
  Double_t minVertexLineSignificance;
  Double_t minCosPointing;
  // Require L_Lambda >= L_Xi when both lengths are available (decay-order cut).
  Bool_t requireDecayLengthOrder;
  // Fake-Lambda / fake-Xi veto (Xi's own daughter-level chi2prim/DCA-to-PV cuts, same
  // concept as K0short's k0MinDaughter* keys). Applied to the Lambda's own daughters
  // (proton, pion) and separately to the bachelor pion. Defaults are disabled (<0).
  Double_t minLambdaDaughterChi2Prim;
  Double_t minLambdaDaughterDistanceToPv;
  Double_t minBachelorChi2Prim;
  Double_t minBachelorDistanceToPv;
  // Optional K0S Stage B (when set, used instead of the shared mass/topo keys for PDG 310).
  // Defaults are disabled (<0); K0-only jobs may keep using the shared keys.
  Double_t k0MinMass;
  Double_t k0MaxMass;
  Double_t k0MaxMassError;
  Double_t k0MaxChi2Ndf;
  Double_t k0MaxTopoChi2Ndf;
  Double_t k0MaxDaughterDistance;
  Double_t k0MaxDistanceToPv;
  Double_t k0MinDecayLength;
  Double_t k0MinDecayLengthSignificance;
  Double_t k0MinVertexLineSignificance;
  Double_t k0MinCosPointing;
  // Daughter-level PSN0892-style cuts, applied to the smaller of (pi+, pi-).
  // Defaults are disabled (<0).
  Double_t k0MinDaughterChi2Prim;
  Double_t k0MinDaughterDistanceToPv;
  // Per-daughter |PID pull| upper limits for K0S (pi+, pi- separately); disabled (<0) by default.
  Double_t k0MaxPionPlusPidPull;
  Double_t k0MaxPionMinusPidPull;
  // K0S Lambda veto: reject if re-assigning a daughter the proton mass gives
  // |M(p pi) - M_Lambda| < k0LambdaVetoMassWindow [GeV]. Momenta at the K0S decay
  // vertex; only daughters with |nSigmaProton| < k0LambdaVetoMaxNSigmaProton.
  // pi+ -> p always; pi- -> pbar only if k0LambdaVetoAntiLambda. Disabled (<0).
  Double_t k0LambdaVetoMassWindow;
  Double_t k0LambdaVetoMaxNSigmaProton;
  Bool_t k0LambdaVetoAntiLambda;
  Bool_t reconstructLambda;
  Bool_t reconstructAntiLambda;
  // Appended study controls; omitted YAML keys retain the established cuts.
  Double_t topoChi2NdfCut;
  Bool_t applyLambdaGeometryCuts;

  Bool_t reconstructXi;
  Bool_t reconstructAntiXi;
  Bool_t reconstructK0short;

  // Official Run21 3 GeV event selection, applied by OfficialEventSelection (StMaker/kfparticle/) instead of
  // CentralityHelper when officialEventSelection is true (default false = unchanged):
  //  - officialBadRunFile: one run number per line (empty = none)
  //  - officialTriggerIds: event must have at least one of these trigger ids (empty = none)
  //  - run batches: officialBatchRunEdges (ascending run numbers); batch = number of edges <= runId
  //    (empty = a single batch). Pileup and centrality parameters are given per batch.
  //  - pileup: keep low(n) < fxtMult < high(n), n = nBTOFMatch; officialPileupLow/High are lists of
  //    polynomial coefficients (any order, c0 first), one list per batch (empty = no pileup cut)
  //  - centrality: cent9 = largest k with fxtMult > officialCentralityEdges[batch][k] (9 lower
  //    bounds, cent9 0 = 70-80% ... 8 = 0-5%); fxtMult <= edges[0] gives cent9 = -1 (not accepted)
  //  - officialRequirePositiveMult: reject nBTOFMatch <= 0 or fxtMult <= 0
  // Events that pass bad run, trigger and the vertex/event cuts are reconstructed and written to
  // the candidate trees (with flags); only events that also pass pileup and centrality are
  // "accepted" (histograms, selected candidates, femto pairing).
  Bool_t officialEventSelection;
  std::string officialBadRunFile;
  std::vector<int> officialTriggerIds;
  std::vector<int> officialBatchRunEdges;
  std::vector<std::vector<double> > officialPileupLow;
  std::vector<std::vector<double> > officialPileupHigh;
  std::vector<std::vector<double> > officialCentralityEdges;
  Bool_t officialRequirePositiveMult;
  Int_t OfficialBatch(Int_t runId) const {
    Int_t b = 0;
    for (size_t i = 0; i < officialBatchRunEdges.size(); ++i) if (runId >= officialBatchRunEdges[i]) b = static_cast<Int_t>(i) + 1;
    return b;
  }

  Bool_t HasK0LambdaVeto() const { return k0LambdaVetoMassWindow > 0.; }

  Bool_t HasK0StageBOverride() const {
    return k0MinMass >= 0. && k0MaxMass > k0MinMass;
  }

private:
  KfParticleCutConfig();
  ~KfParticleCutConfig();
  KfParticleCutConfig(const KfParticleCutConfig&);
  KfParticleCutConfig& operator=(const KfParticleCutConfig&);
};

#endif
