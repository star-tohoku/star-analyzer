#include "KfParticleCandidateSource.h"

#include "ConfigManager.h"
#include "HistManager.h"
#include "KfCandidateCuts.h"
#include "KfParticleHelper.h"
#include "StPicoKFParticleInterface.h"
#include "cuts/KfParticleCutConfig.h"

#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoPhysicalHelix.h"
#include "StPicoEvent/SystemOfUnits.h"

#include "TH1.h"
#include "TMath.h"
#include "TString.h"
#include "TVector3.h"

#include <algorithm>
#include <iostream>
#include <map>
#include <sstream>

namespace {
const char* kXiMinusKey = "xi_minus_kf";
const char* kXiPlusKey = "xi_plus_kf";
const char* kK0sKey = "k0s_kf";

Double_t CandidateRapidity(const TVector3& momentum, Double_t mass) {
  const Double_t energy = TMath::Sqrt(momentum.Mag2() + mass * mass);
  const Double_t denominator = energy - momentum.Z();
  return denominator > 0. ? 0.5 * TMath::Log((energy + momentum.Z()) / denominator) : 0.;
}

void FillStageLabels(TH1* stages, const char* const labels[14]) {
  if (!stages) return;
  for (Int_t i = 0; i < 14; ++i) stages->GetXaxis()->SetBinLabel(i + 1, labels[i]);
}

class KfParticleCandidateSource : public CandidateSource {
public:
  KfParticleCandidateSource()
      : mHist(0), mInterface(0), mPicoDst(0), mEvents(0), mXiSelected(0), mK0Selected(0) {
    for (Int_t i = 0; i < 14; ++i) mXiStages[i] = mK0Stages[i] = 0;
  }
  virtual ~KfParticleCandidateSource() { delete mInterface; }

  virtual const char* Name() const { return "kfparticle"; }

  virtual Bool_t Init(const char* /*mainConfigPath*/, HistManager* hist) {
    mHist = hist;
    const KfParticleCutConfig& cuts = ConfigManager::GetInstance().GetKfParticleCuts();
    if (!cuts.Validate(std::cerr)) { mError = "invalid kf configuration"; return kFALSE; }
    if (!cuts.reconstructXi && !cuts.reconstructAntiXi && !cuts.reconstructK0short) {
      mError = "kf YAML enables neither Xi nor K0short reconstruction";
      return kFALSE;
    }
    if (cuts.reconstructK0short && !cuts.HasK0StageBOverride()) {
      mError = "K0short reconstruction requires the k0MinMass/k0MaxMass Stage B override";
      return kFALSE;
    }
    if (cuts.reconstructXi) mKeys.push_back(kXiMinusKey);
    if (cuts.reconstructAntiXi) mKeys.push_back(kXiPlusKey);
    if (cuts.reconstructK0short) mKeys.push_back(kK0sKey);
    for (size_t i = 0; i < mKeys.size(); ++i) mCandidates[mKeys[i]];
    mInterface = new StPicoKFParticleInterface(cuts);
    cuts.Dump(std::cout);
    std::cout << StPicoKFParticleInterface::BackendDescription() << std::endl;
    const char* xiLabels[] = {"read events", "reconstructed events", "raw tracks", "quality tracks",
        "covariance tracks", "PID tracks", "PID hypotheses", "primary hypotheses",
        "Topo particle slots", "Topo Xis", "invalid Xi", "valid raw Xis",
        "selected Xis", "selected anti-Xis"};
    const char* k0Labels[] = {"read events", "reconstructed events", "raw tracks", "quality tracks",
        "covariance tracks", "PID tracks", "PID hypotheses", "primary hypotheses",
        "Topo particle slots", "Topo K0S", "invalid K0S", "valid raw K0S",
        "selected K0S", "unused"};
    if (mHist) {
      FillStageLabels(mHist->Get("hKfStages"), xiLabels);
      FillStageLabels(mHist->Get("hKfK0Stages"), k0Labels);
    }
    std::cout << "[KfParticleCandidateSource] particle keys:";
    for (size_t i = 0; i < mKeys.size(); ++i) std::cout << " " << mKeys[i];
    std::cout << std::endl;
    return kTRUE;
  }

  virtual std::vector<std::string> SpeciesKeys() const { return mKeys; }

  virtual Bool_t ProvidesSpecies(const std::string& key) const {
    return std::find(mKeys.begin(), mKeys.end(), key) != mKeys.end();
  }

  virtual std::vector<std::string> DaughterNames(const std::string& key) const {
    std::vector<std::string> n;
    if (key == kXiMinusKey || key == kXiPlusKey) {
      n.push_back("proton"); n.push_back("pion"); n.push_back("bachelor");
    } else if (key == kK0sKey) {
      n.push_back("pionPlus"); n.push_back("pionMinus");
    }
    return n;
  }

  virtual std::vector<std::string> VariableNames(const std::string& key) const {
    std::vector<std::string> n;
    if (key == kXiMinusKey || key == kXiPlusKey) {
      const char* v[] = {"massError", "chi2Ndf", "topoChi2Ndf", "daughterDistance", "distanceToPv",
          "decayLength", "decayLengthSignificance", "vertexLineLengthSignificance", "cosPointing",
          "lambdaMass", "lambdaDecayLength", "lambdaDecayLengthSignificance", "lambdaDaughterDistance",
          "protonChi2Prim", "pionChi2Prim", "bachelorChi2Prim",
          "protonDistanceToPv", "pionDistanceToPv", "bachelorDistanceToPv",
          "protonPidPull", "pionPidPull", "bachelorPidPull",
          "protonNHitsFit", "pionNHitsFit", "bachelorNHitsFit"};
      n.assign(v, v + sizeof(v) / sizeof(v[0]));
    } else if (key == kK0sKey) {
      const char* v[] = {"massError", "chi2Ndf", "topoChi2Ndf", "daughterDistance", "distanceToPv",
          "decayLength", "decayLengthSignificance", "vertexLineLengthSignificance", "cosPointing",
          "pionPlusChi2Prim", "pionMinusChi2Prim", "pionPlusDistanceToPv", "pionMinusDistanceToPv",
          "pionPlusPidPull", "pionMinusPidPull", "pionPlusNHitsFit", "pionMinusNHitsFit",
          "lambdaVetoDeviation"};
      n.assign(v, v + sizeof(v) / sizeof(v[0]));
    }
    return n;
  }

  virtual Bool_t ProcessEvent(StPicoDst* picoDst, Int_t cent9) {
    mPicoDst = picoDst;
    if (!mInterface || !picoDst || !picoDst->event()) { mError = "not initialised or no event"; return kFALSE; }
    ++mEvents;
    ++mXiStages[0];
    ++mK0Stages[0];
    if (!mInterface->ProcessEvent(picoDst)) {
      std::ostringstream s;
      s << "run=" << picoDst->event()->runId() << " event=" << picoDst->event()->eventId() << ": "
        << mInterface->LastError();
      mError = s.str();
      return kFALSE;
    }
    const KfParticleCutConfig& k = ConfigManager::GetInstance().GetKfParticleCuts();
    const KfParticleEventStats& stats = mInterface->Stats();
    const Long64_t xiInc[] = {0, 1, stats.rawTracks, stats.qualityTracks,
        stats.covarianceTracks, stats.pidTracks, stats.pidHypotheses, stats.primaryHypotheses,
        stats.finderParticles, stats.xiParticles, stats.invalidXiCandidates, stats.validXiCandidates, 0, 0};
    const Long64_t k0Inc[] = {0, 1, stats.rawTracks, stats.qualityTracks,
        stats.covarianceTracks, stats.pidTracks, stats.pidHypotheses, stats.primaryHypotheses,
        stats.finderParticles, stats.k0Particles, stats.invalidK0Candidates, stats.validK0Candidates, 0, 0};
    for (Int_t i = 0; i < 14; ++i) { mXiStages[i] += xiInc[i]; mK0Stages[i] += k0Inc[i]; }
    if (mHist) {
      const std::vector<Float_t>& rawPt = mInterface->RawTrackPt();
      const std::vector<Float_t>& rawEta = mInterface->RawTrackEta();
      for (size_t i = 0; i < rawPt.size(); ++i) mHist->Fill("hRawTrackPt", rawPt[i]);
      for (size_t i = 0; i < rawEta.size(); ++i) mHist->Fill("hRawTrackEta", rawEta[i]);
      const std::vector<KfLambdaCandidate>& lambdas = mInterface->Candidates();
      for (size_t i = 0; i < lambdas.size(); ++i)
        mHist->Fill(lambdas[i].pdg > 0 ? "hKfLambdaMassRaw" : "hKfAntiLambdaMassRaw", lambdas[i].mass);
    }

    const std::vector<KfXiCandidate>& xis = mInterface->XiCandidates();
    for (size_t i = 0; i < xis.size(); ++i) {
      const KfXiCandidate& c = xis[i];
      if (mHist) {
        // Xi- and Xi+ never share a histogram (hXiMinus_* / hXiPlus_*)
        const char* q = c.pdg > 0 ? "hXiMinus_" : "hXiPlus_";
        mHist->Fill(c.pdg > 0 ? "hKfXiMassRaw" : "hKfAntiXiMassRaw", c.mass);
        mHist->Fill(Form("%sTopoChi2NdfRaw", q), c.topoChi2Ndf);
        if (c.lambdaMass > 0.f) mHist->Fill(Form("%sLambdaMassRaw", q), c.lambdaMass);
      }
      if (!PassKfXiCandidateCuts(c, k)) continue;
      const char* key = c.pdg > 0 ? kXiMinusKey : kXiPlusKey;
      std::map<std::string, std::vector<FemtoSourceCandidate> >::iterator out = mCandidates.find(key);
      if (out == mCandidates.end()) continue;  // charge not requested in the kf YAML
      ++mXiSelected;
      ++mXiStages[c.pdg > 0 ? 12 : 13];
      FillXiQa(c, cent9);
      FemtoSourceCandidate f;
      f.px = c.px; f.py = c.py; f.pz = c.pz;
      f.massValue = c.mass;
      f.x = c.x; f.y = c.y; f.z = c.z;
      f.decayLength = c.decayLength;
      f.daughters.push_back(c.protonIndex);
      f.daughters.push_back(c.pionIndex);
      f.daughters.push_back(c.bachelorIndex);
      const Float_t v[] = {c.massError, c.chi2Ndf, c.topoChi2Ndf, c.daughterDistance, c.distanceToPv,
          c.decayLength, c.decayLengthSignificance, c.vertexLineLengthSignificance, c.cosPointing,
          c.lambdaMass, c.lambdaDecayLength, c.lambdaDecayLengthSignificance, c.lambdaDaughterDistance,
          c.protonChi2Prim, c.pionChi2Prim, c.bachelorChi2Prim,
          c.protonDistanceToPv, c.pionDistanceToPv, c.bachelorDistanceToPv,
          c.protonPidPull, c.pionPidPull, c.bachelorPidPull,
          NHitsFit(c.protonIndex), NHitsFit(c.pionIndex), NHitsFit(c.bachelorIndex)};
      f.variables.assign(v, v + sizeof(v) / sizeof(v[0]));
      out->second.push_back(f);
    }

    const std::vector<KfLambdaCandidate>& k0s = mInterface->K0shortCandidates();
    std::map<std::string, std::vector<FemtoSourceCandidate> >::iterator k0Out = mCandidates.find(kK0sKey);
    for (size_t i = 0; i < k0s.size() && k0Out != mCandidates.end(); ++i) {
      const KfLambdaCandidate& c = k0s[i];  // K0S layout: proton* = pi+, pion* = pi-
      if (mHist) {
        mHist->Fill("hKfK0shortMassRaw", c.mass);
        mHist->Fill("hK0short_TopoChi2NdfRaw", c.topoChi2Ndf);
      }
      const Bool_t passTopo = PassKfK0CandidateCuts(c, k);
      Double_t dauP[6];
      Double_t nSigmaProton[2];
      const Int_t idx[2] = {c.protonIndex, c.pionIndex};
      for (Int_t d = 0; d < 2; ++d) {
        const StPicoTrack* t = Track(idx[d]);
        // Momentum at the K0 decay vertex (helix propagated to its point of closest approach
        // to the KF vertex), so the daughter 4-vectors reproduce the K0 mass.
        TVector3 g(0., 0., 0.);
        if (t) {
          const Float_t bField = picoDst->event()->bField();
          const StPicoPhysicalHelix helix = t->helix(bField);
          const Double_t sPath = helix.pathLength(TVector3(c.x, c.y, c.z));
          g = helix.momentumAt(sPath, bField * units::kilogauss);
        }
        dauP[3 * d] = g.X(); dauP[3 * d + 1] = g.Y(); dauP[3 * d + 2] = g.Z();
        nSigmaProton[d] = t ? t->nSigmaProton() : -999.;
      }
      const Double_t vetoDev = KfK0LambdaVetoDeviation(dauP[0], dauP[1], dauP[2], dauP[3], dauP[4], dauP[5],
                                                       nSigmaProton[0], nSigmaProton[1], k);
      if (passTopo && k.HasK0LambdaVeto() && mHist)
        mHist->Fill("hK0short_InvMass_vs_LambdaVetoDeviation", std::min(vetoDev, 0.1999), c.mass);
      if (!passTopo || !PassKfK0LambdaVeto(vetoDev, k)) continue;
      ++mK0Selected;
      ++mK0Stages[12];
      FillK0Qa(c, cent9);
      FemtoSourceCandidate f;
      f.px = c.px; f.py = c.py; f.pz = c.pz;
      f.massValue = c.mass;
      f.x = c.x; f.y = c.y; f.z = c.z;
      f.decayLength = c.decayLength;
      f.daughters.push_back(c.protonIndex);
      f.daughters.push_back(c.pionIndex);
      const Float_t v[] = {c.massError, c.chi2Ndf, c.topoChi2Ndf, c.daughterDistance, c.distanceToPv,
          c.decayLength, c.decayLengthSignificance, c.vertexLineLengthSignificance, c.cosPointing,
          c.protonChi2Prim, c.pionChi2Prim, c.protonDistanceToPv, c.pionDistanceToPv,
          c.protonPidPull, c.pionPidPull, NHitsFit(c.protonIndex), NHitsFit(c.pionIndex),
          static_cast<Float_t>(std::min(vetoDev, 1e3))};
      f.variables.assign(v, v + sizeof(v) / sizeof(v[0]));
      k0Out->second.push_back(f);
    }
    return kTRUE;
  }

  virtual const std::vector<FemtoSourceCandidate>& Candidates(const std::string& key) const {
    std::map<std::string, std::vector<FemtoSourceCandidate> >::const_iterator it = mCandidates.find(key);
    return it != mCandidates.end() ? it->second : mEmpty;
  }

  virtual void Clear() {
    if (mInterface) mInterface->Clear();
    for (std::map<std::string, std::vector<FemtoSourceCandidate> >::iterator it = mCandidates.begin();
         it != mCandidates.end(); ++it)
      it->second.clear();
  }

  virtual void Finish() {
    if (mHist) {
      TH1* xiStages = mHist->Get("hKfStages");
      if (xiStages) for (Int_t i = 0; i < 14; ++i) xiStages->SetBinContent(i + 1, mXiStages[i]);
      TH1* k0Stages = mHist->Get("hKfK0Stages");
      if (k0Stages) for (Int_t i = 0; i < 14; ++i) k0Stages->SetBinContent(i + 1, mK0Stages[i]);
    }
    std::cout << "[KfParticleCandidateSource] events=" << mEvents << " selectedXi=" << mXiSelected
              << " selectedK0=" << mK0Selected << std::endl;
  }

  virtual const std::string& LastError() const { return mError; }

private:
  const StPicoTrack* Track(Int_t index) const {
    if (!mPicoDst || index < 0 || index >= static_cast<Int_t>(mPicoDst->numberOfTracks())) return 0;
    return mPicoDst->track(index);
  }
  Float_t NHitsFit(Int_t index) const {
    const StPicoTrack* t = Track(index);
    return t ? static_cast<Float_t>(t->nHitsFit()) : -1.f;
  }

  // Selected Xi candidates, one set of histograms per charge (hXiMinus_* for pdg > 0, hXiPlus_* otherwise).
  void FillXiQa(const KfXiCandidate& c, Int_t cent9) {
    if (!mHist) return;
    const TVector3 momentum(c.px, c.py, c.pz);
    const TString q = c.pdg > 0 ? "hXiMinus_" : "hXiPlus_";
    mHist->Fill(q + "InvMass", c.mass);
    mHist->Fill(c.pdg > 0 ? "hKfXiMassSelected" : "hKfAntiXiMassSelected", c.mass);
    if (c.lambdaMass > 0.f) mHist->Fill(q + "LambdaMassSelected", c.lambdaMass);
    mHist->Fill(q + "Pt", momentum.Pt());
    mHist->Fill(q + "Eta", momentum.PseudoRapidity());
    mHist->Fill(q + "DCA12", c.daughterDistance);
    mHist->Fill(q + "DCAV0", c.distanceToPv);
    mHist->Fill(q + "CosPointing", c.cosPointing);
    mHist->Fill(q + "InvMass_vs_Pt", momentum.Pt(), c.mass);
    mHist->Fill(q + "InvMass_vs_DecayLength", c.decayLength, c.mass);
    mHist->Fill(q + "InvMass_vs_Y", CandidateRapidity(momentum, c.mass), c.mass);
    if (cent9 >= 0) mHist->Fill(q + "InvMass_vs_Cent9", cent9, c.mass);
    mHist->Fill(q + "Chi2Ndf", c.chi2Ndf);
    mHist->Fill(q + "TopoChi2Ndf", c.topoChi2Ndf);
    mHist->Fill(q + "DecayLengthSignificance", c.decayLengthSignificance);
    mHist->Fill("hKfPdg", c.pdg);
  }

  void FillK0Qa(const KfLambdaCandidate& c, Int_t cent9) {
    if (!mHist) return;
    const TVector3 momentum(c.px, c.py, c.pz);
    mHist->Fill("hK0short_InvMass", c.mass);
    mHist->Fill("hKfK0shortMassSelected", c.mass);
    mHist->Fill("hK0short_Pt", momentum.Pt());
    mHist->Fill("hK0short_Eta", momentum.PseudoRapidity());
    mHist->Fill("hK0short_DCA12", c.daughterDistance);
    mHist->Fill("hK0short_DCAV0", c.distanceToPv);
    mHist->Fill("hK0short_CosPointing", c.cosPointing);
    mHist->Fill("hK0short_Chi2Ndf", c.chi2Ndf);
    mHist->Fill("hK0short_TopoChi2Ndf", c.topoChi2Ndf);
    mHist->Fill("hK0short_DecayLengthSignificance", c.decayLengthSignificance);
    mHist->Fill("hK0short_Pdg", c.pdg);
    mHist->Fill("hK0short_InvMass_vs_Pt", momentum.Pt(), c.mass);
    mHist->Fill("hK0short_InvMass_vs_DecayLength", c.decayLength, c.mass);
    mHist->Fill("hK0short_InvMass_vs_Y", CandidateRapidity(momentum, c.mass), c.mass);
    if (cent9 >= 0) mHist->Fill("hK0short_InvMass_vs_Cent9", cent9, c.mass);
    mHist->Fill("hK0short_InvMass_vs_Chi2Ndf", c.chi2Ndf, c.mass);
    mHist->Fill("hK0short_InvMass_vs_TopoChi2Ndf", c.topoChi2Ndf, c.mass);
    mHist->Fill("hK0short_InvMass_vs_DecayLengthSignificance", c.decayLengthSignificance, c.mass);
    mHist->Fill("hK0short_InvMass_vs_CosPointing", c.cosPointing, c.mass);
    mHist->Fill("hK0short_InvMass_vs_DaughterChi2Prim", std::min(c.protonChi2Prim, c.pionChi2Prim), c.mass);
    mHist->Fill("hK0short_InvMass_vs_DaughterDistanceToPv",
                std::min(c.protonDistanceToPv, c.pionDistanceToPv), c.mass);
  }

  HistManager* mHist;
  StPicoKFParticleInterface* mInterface;
  StPicoDst* mPicoDst;
  std::vector<std::string> mKeys;
  std::map<std::string, std::vector<FemtoSourceCandidate> > mCandidates;
  std::vector<FemtoSourceCandidate> mEmpty;
  std::string mError;
  Long64_t mEvents;
  Long64_t mXiSelected;
  Long64_t mK0Selected;
  Long64_t mXiStages[14];
  Long64_t mK0Stages[14];
};
}  // namespace

CandidateSource* createKfParticleCandidateSource() { return new KfParticleCandidateSource(); }
