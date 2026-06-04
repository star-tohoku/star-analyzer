#include "StFemtoMaker.h"
#include "ConfigManager.h"
#include "HistManager.h"
#include "kinematics.h"
#include "cuts/EventCutConfig.h"
#include "cuts/TrackCutConfig.h"
#include "cuts/PIDCutConfig.h"
#include "cuts/PhiCutConfig.h"
#include "cuts/MixingConfig.h"
#include "cuts/CentralityCutConfig.h"
#include "cuts/FemtoConfig.h"
#include "CentralityHelper.h"
#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"
#include "StarClassLibrary/StPhysicalHelixD.hh"
#include "StarClassLibrary/StThreeVectorF.hh"
#include "StarClassLibrary/SystemOfUnits.h"

#include "TFile.h"
#include "TH1.h"
#include "TString.h"
#include "TMath.h"
#include "TSystem.h"
#include "TVector3.h"
#include "TVector2.h"
#include "TLorentzVector.h"
#include "TRandom.h"
#include <iostream>
#include <sstream>
#include <vector>

namespace {
const Double_t kKaonMass = 0.493677;
const Double_t kProtonMass = 0.938272;
const Double_t kPhiMass = 1.019461;
}  // namespace

StFemtoMaker::StFemtoMaker(const char* name, StPicoDstMaker* picoMaker, const char* outName)
    : StMaker(name),
      mPicoDstMaker(picoMaker),
      mPicoDst(0),
      mOutName(outName),
      mEventCounter(0),
      m_histManager(0),
      m_centrality(0),
      m_cent9(-1),
      m_cent16(-1),
      m_refMultCorr(-1.0),
      m_centWeight(1.0),
      m_centralityPercent(-1.0),
      m_psi2(-1.0) {}

StFemtoMaker::~StFemtoMaker() {
  if (m_centrality) {
    delete m_centrality;
    m_centrality = 0;
  }
  if (m_histManager) {
    delete m_histManager;
    m_histManager = 0;
  }
}

StFemtoMaker* createStFemtoMaker(const char* name, StPicoDstMaker* picoMaker, const char* outName) {
  return new StFemtoMaker(name, picoMaker, outName);
}

extern "C" void* createStFemtoMakerC(const char* name, void* picoMaker, const char* outName) {
  return (void*)createStFemtoMaker(name, (StPicoDstMaker*)picoMaker, outName);
}

Int_t StFemtoMaker::Init() {
  std::string histPath = ConfigManager::GetInstance().GetHistConfigPath(GetName());
  if (histPath.empty()) {
    std::cerr << "[StFemtoMaker] GetHistConfigPath() returned empty; no histograms will be filled." << std::endl;
    m_histManager = 0;
    return kStOK;
  }
  m_histManager = new HistManager();
  if (!m_histManager->LoadFromFile(histPath.c_str())) {
    std::cerr << "[StFemtoMaker] Failed to load hist config from " << histPath << std::endl;
    delete m_histManager;
    m_histManager = 0;
    return kStOK;
  }

  m_centrality = new CentralityHelper();
  if (!m_centrality->Init(ConfigManager::GetInstance().GetCentralityCuts())) {
    std::cerr << "[StFemtoMaker] CentralityHelper init failed" << std::endl;
  }

  if (!ConfigManager::GetInstance().GetPhiCuts().FinalizeRapidityFrame(
          ConfigManager::GetInstance().GetCentralityCuts())) {
    std::cerr << "[StFemtoMaker] FinalizeRapidityFrame failed; fix maker/centrality YAML." << std::endl;
    return kStErr;
  }

  if (!ConfigManager::GetInstance().GetFemtoConfig().Validate()) {
    std::cerr << "[StFemtoMaker] FemtoConfig validation failed." << std::endl;
    return kStErr;
  }

  return kStOK;
}

void StFemtoMaker::Clear(Option_t* opt) {
  (void)opt;
  m_eventCandidates.clear();
}

Int_t StFemtoMaker::Make() {
  if (!mPicoDstMaker) return kStWarn;
  mPicoDst = mPicoDstMaker->picoDst();
  if (!mPicoDst) return kStWarn;

  StPicoEvent* event = mPicoDst->event();
  if (!event) return kStWarn;

  mEventCounter++;
  m_eventCandidates.clear();

  TVector3 pVtx = event->primaryVertex();
  Float_t vzVpd = event->vzVpd();
  Int_t refMult = event->refMult();
  EventCutConfig& evCuts = ConfigManager::GetInstance().GetEventCuts();
  Float_t vr = evCuts.ComputeVr(pVtx.X(), pVtx.Y());
  const Double_t vz = pVtx.Z();
  const Int_t runId = event->runId();
  const Int_t nBTOFMatch = event->nBTOFMatch();

  CentralityCutConfig& centCfg = ConfigManager::GetInstance().GetCentralityCuts();
  Int_t rawMult = refMult;
  if (centCfg.enabled) {
    TString mode(centCfg.mode.c_str());
    mode.ToLower();
    if (mode == "fxtmult") rawMult = event->fxtMult();
  }

  m_cent9 = -1;
  m_cent16 = -1;
  m_refMultCorr = -1.0;
  m_centWeight = 1.0;
  m_centralityPercent = -1.0;
  m_psi2 = -1.0;

  if (m_histManager) {
    m_histManager->Fill("hRefMultVsNTOFMatch", (Double_t)nBTOFMatch, (Double_t)rawMult);
    m_histManager->Fill("hVz", pVtx.Z());
    m_histManager->Fill("hVr", vr);
    m_histManager->Fill("hVxVy", pVtx.X(), pVtx.Y());
    m_histManager->Fill("hRefMult", refMult);
    m_histManager->Fill("hVzVsRun", (Double_t)event->runId(), pVtx.Z());
    m_histManager->Fill("hRefMultVsVz", pVtx.Z(), refMult);
    EventCutConfig& ev = ConfigManager::GetInstance().GetEventCuts();
    if (TMath::Abs(vzVpd) < ev.maxAbsVzVpd) {
      m_histManager->Fill("hVzDiff", pVtx.Z() - vzVpd);
    }
    std::vector<unsigned int> triggerIds = event->triggerIds();
    for (size_t i = 0; i < triggerIds.size(); i++) {
      m_histManager->Fill("hTriggerIds", triggerIds[i]);
    }
  }

  CentralityRejectReason centReason = kCentralityOk;
  if (m_centrality && m_centrality->IsEnabled()) {
    if (!m_centrality->CheckBadRun(runId, centReason)) return kStOK;
  }

  if (!PassEventCuts(pVtx.Z(), vr, refMult, vzVpd)) return kStOK;

  if (m_histManager) {
    m_histManager->Fill("hVz_After", pVtx.Z());
    m_histManager->Fill("hVr_After", vr);
    m_histManager->Fill("hVxVy_After", pVtx.X(), pVtx.Y());
    m_histManager->Fill("hRefMult_After", refMult);
    EventCutConfig& ev = ConfigManager::GetInstance().GetEventCuts();
    if (TMath::Abs(vzVpd) < ev.maxAbsVzVpd) {
      m_histManager->Fill("hVzDiff_After", pVtx.Z() - vzVpd);
    }
  }

  if (m_histManager && centCfg.fillCentralityQA) {
    m_histManager->Fill("hRawMult", (Double_t)rawMult);
  }

  if (m_centrality && m_centrality->IsEnabled()) {
    if (!m_centrality->CheckPileup(rawMult, nBTOFMatch, vz, centReason)) return kStOK;
    if (!m_centrality->ComputeBins(event, rawMult, vz, m_cent9, m_cent16, m_refMultCorr, m_centWeight, centReason)) {
      return kStOK;
    }

    if (m_histManager) {
      m_histManager->Fill("hCentralityRaw", (Double_t)m_cent9);
      m_histManager->Fill("hRefMultCorr", m_refMultCorr);
      m_histManager->Fill("hCentralityVsVz", vz, (Double_t)m_cent9);
      m_histManager->Fill("hRefMultWeight", m_centWeight);
      m_histManager->Fill("hRefMultVsNTOFMatchAfter", (Double_t)nBTOFMatch, (Double_t)rawMult);
    }

    if (!m_centrality->AcceptCentBin(m_cent9, centReason)) return kStOK;

    m_centralityPercent = CentralityHelper::Cent9ToPercentile(m_cent9);
    if (m_histManager) {
      const Double_t w = centCfg.useWeight ? m_centWeight : 1.0;
      TH1* hCent = m_histManager->Get("hCentrality");
      if (hCent) hCent->Fill((Double_t)m_cent9, w);
      TH1* hCent16 = m_histManager->Get("hCentrality16");
      if (hCent16) hCent16->Fill((Double_t)m_cent16, w);
    }
  }

  PhiCutConfig& phiCfg = ConfigManager::GetInstance().GetPhiCuts();
  const Int_t kMaxTracks = 4000;
  std::vector<TrackState> kaonsPlus;
  std::vector<TrackState> kaonsMinus;
  std::vector<TrackState> protons;
  kaonsPlus.reserve(500);
  kaonsMinus.reserve(500);
  protons.reserve(500);

  Double_t Qx = 0.0, Qy = 0.0;
  Int_t nTracks = mPicoDst->numberOfTracks();
  if (m_histManager) m_histManager->Fill("hNTracks", (Double_t)nTracks);
  if (phiCfg.maxNTr > 0 && nTracks > phiCfg.maxNTr) return kStOK;

  for (Int_t itrk = 0; itrk < nTracks; itrk++) {
    StPicoTrack* trk = mPicoDst->track(itrk);
    if (!trk) continue;
    TVector3 pMom = trk->pMom();
    if (m_histManager) {
      Float_t ptRaw = pMom.Perp();
      Float_t etaRaw = pMom.PseudoRapidity();
      m_histManager->Fill("hPt_Raw", ptRaw);
      m_histManager->Fill("hEta_Raw", etaRaw);
      m_histManager->Fill("hPhi_Raw", pMom.Phi());
      m_histManager->Fill("hNHitsFit_Raw", trk->nHitsFit());
      if (trk->nHitsMax() > 0) {
        m_histManager->Fill("hNHitsRatio_Raw", (Float_t)trk->nHitsFit() / (Float_t)trk->nHitsMax());
      }
      m_histManager->Fill("hNHitsDedx_Raw", trk->nHitsDedx());
      m_histManager->Fill("hChi2_Raw", trk->chi2());
      m_histManager->Fill("hDCA_Raw", trk->gDCA(pVtx).Mag());
    }
    if (!PassTrackCuts(trk, pVtx)) continue;

    Float_t pt = pMom.Perp();
    Float_t eta = pMom.PseudoRapidity();
    Float_t phi = pMom.Phi();
    Int_t btofIndex = trk->bTofPidTraitsIndex();

    if (m_histManager) {
      m_histManager->Fill("hPt", pt);
      m_histManager->Fill("hEta", eta);
      m_histManager->Fill("hPhi", phi);
      m_histManager->Fill("hNHitsFit", trk->nHitsFit());
      m_histManager->Fill("hNHitsRatio", (Float_t)trk->nHitsFit() / (Float_t)trk->nHitsMax());
      m_histManager->Fill("hNHitsDedx", trk->nHitsDedx());
      m_histManager->Fill("hDCA", trk->gDCA(pVtx).Mag());
      m_histManager->Fill("hCharge", trk->charge());
      m_histManager->Fill("hChi2", trk->chi2());
      m_histManager->Fill("hNSigmaKaon_Raw", trk->nSigmaKaon());
      m_histManager->Fill("hNSigmaProtonVsPt", pt, trk->nSigmaProton());
      if (trk->charge() > 0) {
        m_histManager->Fill("hNSigmaProtonVsPt_Pos", pt, trk->nSigmaProton());
      }
    }

    if (pt >= phiCfg.minPtEp && pt <= phiCfg.maxPtEp && TMath::Abs(eta) < phiCfg.maxEtaEp) {
      Qx += TMath::Cos(2.0 * phi);
      Qy += TMath::Sin(2.0 * phi);
    }

    if (PassKaonCuts(trk, pVtx)) {
      TrackState kTrack;
      BuildTrackState(kTrack, trk, event, pVtx, itrk);
      FillTofInfo(kTrack, trk, pMom, btofIndex);
      if (IsKaon(kTrack)) {
      if (m_histManager) {
        m_histManager->Fill("hK_Pt", kTrack.pT);
        m_histManager->Fill("hK_Eta", kTrack.eta);
        m_histManager->Fill("hK_NSigma", kTrack.nSigmaKaon);
      }
      if (kTrack.charge > 0 && (Int_t)kaonsPlus.size() < kMaxTracks) {
          kaonsPlus.push_back(kTrack);
        } else if (kTrack.charge < 0 && (Int_t)kaonsMinus.size() < kMaxTracks) {
          kaonsMinus.push_back(kTrack);
        }
      }
    }

    if (PassProtonCuts(trk, pVtx)) {
      TrackState pTrack;
      BuildTrackState(pTrack, trk, event, pVtx, itrk);
      FillTofInfo(pTrack, trk, pMom, btofIndex);
      if (IsProton(pTrack) && pTrack.charge > 0 && (Int_t)protons.size() < kMaxTracks) {
        if (m_histManager && pTrack.tofMatch) {
          m_histManager->Fill("hP_Mass2", pTrack.mass2);
        }
        protons.push_back(pTrack);
      }
    }
  }

  if (m_histManager) {
    m_histManager->Fill("hTofMatchMult", (Double_t)nBTOFMatch);
    m_histManager->Fill("hNKaonPlus", (Double_t)kaonsPlus.size());
    m_histManager->Fill("hNKaonMinus", (Double_t)kaonsMinus.size());
    m_histManager->Fill("hNKaonPlusVsNKaonMinus", (Double_t)kaonsPlus.size(), (Double_t)kaonsMinus.size());
  }

  TVector2 Q(Qx, Qy);
  if (Q.Mod() > 0) {
    m_psi2 = 0.5 * TMath::ATan2(Qy, Qx);
    if (m_psi2 < 0) m_psi2 += TMath::Pi();
  }

  FemtoConfig& femtoCfg = ConfigManager::GetInstance().GetFemtoConfig();
  for (std::map<std::string, FemtoConfig::SpeciesDef>::const_iterator it = femtoCfg.species.begin();
       it != femtoCfg.species.end(); ++it) {
    const FemtoConfig::SpeciesDef& sp = it->second;
    if (sp.builderType == "track") {
      BuildTrackPidCandidates(sp.key, sp.particleKey, protons, mEventCounter);
    } else if (sp.builderType == "resonance") {
      BuildResonanceCandidates(sp.key, sp.particleKey, kaonsPlus, kaonsMinus, mEventCounter);
    }
  }

  FillCandidateQA();

  for (size_t ic = 0; ic < femtoCfg.channels.size(); ic++) {
    const FemtoConfig::ChannelDef& ch = femtoCfg.channels[ic];
    if (!ch.enabled) continue;
    FillSameEventPairs(ch);
    if (ch.doMixing) FillMixedEventPairs(ch, vz, m_cent9, m_psi2);
  }

  StoreEventForMixing(vz, m_cent9, m_psi2);

  Int_t nPhiCandidates = 0;
  FemtoCandidateStore::const_iterator itPhi = m_eventCandidates.find("phi");
  if (itPhi != m_eventCandidates.end()) nPhiCandidates = (Int_t)itPhi->second.size();

  if (m_histManager && centCfg.fillCentralityQA && m_cent9 >= 0) {
    FillCentralityEventQA(m_cent9, rawMult, m_refMultCorr, nTracks, nBTOFMatch, (Int_t)kaonsPlus.size(),
                          (Int_t)kaonsMinus.size(), nPhiCandidates, (Int_t)protons.size());
  }

  if (m_histManager) {
    TVector2 Q(Qx, Qy);
    m_histManager->Fill("hQxQy", Qx, Qy);
    if (Q.Mod() > 0) {
      Double_t psi2 = 0.5 * TMath::ATan2(Qy, Qx);
      if (psi2 < 0) psi2 += TMath::Pi();
      m_histManager->Fill("hPsi2", psi2);
    }
    m_histManager->Fill("hN", 0);
  }
  return kStOK;
}

Int_t StFemtoMaker::Finish() {
  FinalizeCorrelationFunctions();
  if (mOutName != "") {
    TFile* fout = new TFile(mOutName.Data(), "RECREATE");
    fout->cd();
    WriteHistograms();
    if (m_histManager) m_histManager->ReleaseOwnership();
    fout->Close();
    delete fout;
  }
  std::cout << "StFemtoMaker::Finish() processed " << mEventCounter << " events" << std::endl;
  if (m_centrality && m_centrality->IsEnabled()) {
    m_centrality->Finish();
  }
  return kStOK;
}

void StFemtoMaker::WriteHistograms() {
  if (m_histManager) m_histManager->Write();
}

Bool_t StFemtoMaker::PassEventCuts(Float_t vz, Float_t vr, Int_t refMult, Float_t vzVpd) {
  EventCutConfig& ev = ConfigManager::GetInstance().GetEventCuts();
  if (vz < ev.minVz || vz > ev.maxVz) return kFALSE;
  if (vr > ev.maxVr) return kFALSE;
  if (refMult < ev.minRefMult) return kFALSE;
  if (refMult > ev.maxRefMult) return kFALSE;
  if (TMath::Abs(vz - vzVpd) > ev.maxVzDiff && TMath::Abs(vzVpd) < ev.maxAbsVzVpd) return kFALSE;
  return kTRUE;
}

Bool_t StFemtoMaker::PassTrackCuts(StPicoTrack* trk, TVector3& pVtx) {
  TrackCutConfig& tr = ConfigManager::GetInstance().GetTrackCuts();
  if (trk->nHitsFit() < tr.minNHitsFit) return kFALSE;
  if ((Float_t)trk->nHitsFit() / (Float_t)trk->nHitsMax() < tr.minNHitsRatio) return kFALSE;
  if (trk->nHitsDedx() < tr.minNHitsDedx) return kFALSE;
  TVector3 pMom = trk->pMom();
  if (pMom.Mag() < 1e-4) return kFALSE;
  Float_t pt = pMom.Perp();
  Float_t eta = pMom.PseudoRapidity();
  if (pt < tr.minPt || pt > tr.maxPt) return kFALSE;
  if (eta < tr.minEta || eta > tr.maxEta) return kFALSE;
  if (trk->gDCA(pVtx).Mag() > tr.maxDCA) return kFALSE;
  if (tr.requirePrimaryTrack && !trk->isPrimary()) return kFALSE;
  if (trk->chi2() > tr.maxChi2) return kFALSE;
  return kTRUE;
}

Bool_t StFemtoMaker::PassKaonCuts(StPicoTrack* trk, TVector3& pVtx) {
  if (!PassTrackCuts(trk, pVtx)) return kFALSE;
  PhiCutConfig& phi = ConfigManager::GetInstance().GetPhiCuts();
  if (trk->gDCA(pVtx).Mag() > phi.maxDCAKaon) return kFALSE;
  if (TMath::Abs(trk->nSigmaKaon()) > phi.nSigmaKaon) return kFALSE;
  return kTRUE;
}

Bool_t StFemtoMaker::PassProtonCuts(StPicoTrack* trk, TVector3& pVtx) {
  if (!PassTrackCuts(trk, pVtx)) return kFALSE;
  PIDCutConfig& pid = ConfigManager::GetInstance().GetPIDCuts();
  if (TMath::Abs(trk->nSigmaProton()) > pid.nSigmaProton) return kFALSE;
  return kTRUE;
}

void StFemtoMaker::BuildTrackState(TrackState& track, StPicoTrack* pico, StPicoEvent* event, TVector3& pVtx,
                                   Int_t index) {
  TVector3 pmom = pico->pMom();
  TVector3 org = pico->origin();
  track.originX = org.X();
  track.originY = org.Y();
  track.originZ = org.Z();
  track.momentumX = pmom.X();
  track.momentumY = pmom.Y();
  track.momentumZ = pmom.Z();
  track.BField = event->bField();
  track.pT = pmom.Perp();
  track.eta = pmom.PseudoRapidity();
  track.phi = pmom.Phi();
  track.charge = pico->charge();
  track.nHitsFit = pico->nHitsFit();
  track.nHitsMax = pico->nHitsMax();
  track.nHitsDedx = pico->nHitsDedx();
  track.DCA = pico->gDCA(pVtx).Mag();
  track.chi2 = pico->chi2();
  track.nSigmaKaon = pico->nSigmaKaon();
  track.nSigmaProton = pico->nSigmaProton();
  track.trackIndex = index;
  track.tofMatch = kFALSE;
  track.mass2 = -999.0f;
  track.deltaOneOverBeta = 999.0f;
}

void StFemtoMaker::FillTofInfo(TrackState& track, StPicoTrack* trk, const TVector3& pMom, Int_t btofIndex) {
  track.mass2 = -999.0;
  track.deltaOneOverBeta = 999.0;
  track.tofMatch = kFALSE;
  if (btofIndex < 0) return;
  StPicoBTofPidTraits* tof = mPicoDst->btofPidTraits(btofIndex);
  if (!tof) return;
  Double_t beta = tof->btofBeta();
  if (beta <= 1e-4) return;
  Double_t pMag = pMom.Mag();
  track.mass2 = pMom.Mag2() * (1.0 / (beta * beta) - 1.0);
  track.deltaOneOverBeta = DeltaOneOverBeta(beta, kKaonMass, pMag);
  track.tofMatch = kTRUE;
}

Bool_t StFemtoMaker::PassTofKaonPid(const TrackState& trk) const {
  const PIDCutConfig& pid = ConfigManager::GetInstance().GetPIDCuts();
  if (!pid.requireTOF) return kTRUE;
  TString fallbackMode(pid.tofFallbackMode.c_str());
  fallbackMode.ToLower();
  if (fallbackMode.IsNull()) fallbackMode = "acceptlowpt";
  if (fallbackMode == "acceptlowpt" && trk.pT <= pid.pTofFallbackMax) return kTRUE;
  if (!trk.tofMatch) {
    if (fallbackMode == "tpconly") return kTRUE;
    return kFALSE;
  }
  Bool_t pass = kTRUE;
  if (pid.tofUseMass2Cut) {
    pass = pass && (trk.mass2 >= pid.minMass2Kaon && trk.mass2 <= pid.maxMass2Kaon);
  }
  if (pid.tofUseDeltaInvBetaCut) {
    pass = pass && (TMath::Abs(trk.deltaOneOverBeta) <= pid.maxAbsDeltaOneOverBetaKaon);
  }
  return pass;
}

Bool_t StFemtoMaker::PassTofProtonPid(const TrackState& trk) const {
  const PIDCutConfig& pid = ConfigManager::GetInstance().GetPIDCuts();
  if (!pid.requireTOF) return kTRUE;
  if (!trk.tofMatch) {
    TString fallbackMode(pid.tofFallbackMode.c_str());
    fallbackMode.ToLower();
    if (fallbackMode == "tpconly") return kTRUE;
    return kFALSE;
  }
  if (pid.tofUseMass2Cut) {
    return (trk.mass2 >= pid.minMass2Proton && trk.mass2 <= pid.maxMass2Proton);
  }
  return kTRUE;
}

Bool_t StFemtoMaker::IsKaon(const TrackState& trk) { return PassTofKaonPid(trk); }

Bool_t StFemtoMaker::IsProton(const TrackState& trk) { return PassTofProtonPid(trk); }

TVector3 StFemtoMaker::TrackMomentum(const TrackState& trk) const {
  return TVector3(trk.momentumX, trk.momentumY, trk.momentumZ);
}

StPhysicalHelixD StFemtoMaker::BuildHelix(const TrackState& trk) {
  StThreeVectorF gmomSt(trk.momentumX, trk.momentumY, trk.momentumZ);
  StThreeVectorF orgSt(trk.originX, trk.originY, trk.originZ);
  return StPhysicalHelixD(gmomSt, orgSt, trk.BField * units::kilogauss, static_cast<float>(trk.charge));
}

Double_t StFemtoMaker::CalculateDCA(const TrackState& trk1, const TrackState& trk2, TVector3& dcaPos1,
                                    TVector3& dcaPos2) {
  StPhysicalHelixD helix1 = BuildHelix(trk1);
  StPhysicalHelixD helix2 = BuildHelix(trk2);
  std::pair<Double_t, Double_t> pathLengths = helix1.pathLengths(helix2);
  StThreeVectorF pos1 = helix1.at(pathLengths.first);
  StThreeVectorF pos2 = helix2.at(pathLengths.second);
  dcaPos1 = TVector3(pos1.x(), pos1.y(), pos1.z());
  dcaPos2 = TVector3(pos2.x(), pos2.y(), pos2.z());
  return (dcaPos1 - dcaPos2).Mag();
}

Bool_t StFemtoMaker::ReconstructPhi(const TrackState& kPlus, const TrackState& kMinus, Double_t& invMass,
                                  TVector3& phiMom, TVector3& dcaPosPlus, TVector3& dcaPosMinus) {
  PhiCutConfig& phi = ConfigManager::GetInstance().GetPhiCuts();
  Double_t dca = CalculateDCA(kPlus, kMinus, dcaPosPlus, dcaPosMinus);
  if (dca > phi.maxDCAKK) return kFALSE;
  const TVector3 pPlus = TrackMomentum(kPlus);
  const TVector3 pMinus = TrackMomentum(kMinus);
  Double_t EPlus = TMath::Sqrt(kKaonMass * kKaonMass + pPlus.Mag2());
  Double_t EMinus = TMath::Sqrt(kKaonMass * kKaonMass + pMinus.Mag2());
  phiMom = pPlus + pMinus;
  invMass = TMath::Sqrt((EPlus + EMinus) * (EPlus + EMinus) - phiMom.Mag2());
  return kTRUE;
}

Double_t StFemtoMaker::CalculateOpeningAngle(const TrackState& trk1, const TrackState& trk2) {
  TVector3 p1 = TrackMomentum(trk1);
  TVector3 p2 = TrackMomentum(trk2);
  Double_t mag1 = p1.Mag();
  Double_t mag2 = p2.Mag();
  if (mag1 < 1e-10 || mag2 < 1e-10) return TMath::Pi();
  Double_t cosTheta = p1.Dot(p2) / (mag1 * mag2);
  if (cosTheta > 1.0) cosTheta = 1.0;
  if (cosTheta < -1.0) cosTheta = -1.0;
  return TMath::ACos(cosTheta);
}

Double_t StFemtoMaker::CalculatePairRapidity(Double_t invMass, const TVector3& phiMom) {
  Double_t E = TMath::Sqrt(invMass * invMass + phiMom.Mag2());
  Double_t pz = phiMom.Z();
  if (E <= TMath::Abs(pz)) return 0.0;
  return 0.5 * TMath::Log((E + pz) / (E - pz));
}

Double_t StFemtoMaker::ApplyRapidityFrame(Double_t yLab) const {
  return ConfigManager::GetInstance().GetPhiCuts().ApplyAnalysisRapidity(yLab);
}

Bool_t StFemtoMaker::PassPairTofCut(const TrackState& kPlus, const TrackState& kMinus) const {
  const PIDCutConfig& pid = ConfigManager::GetInstance().GetPIDCuts();
  const Float_t pLow = pid.pMomKaonPID;
  const Float_t pKplus = TrackMomentum(kPlus).Mag();
  const Float_t pKminus = TrackMomentum(kMinus).Mag();
  auto inKaonMass2 = [&](Float_t mass2) {
    return (mass2 > pid.minMass2Kaon && mass2 < pid.maxMass2Kaon);
  };
  Bool_t passPlus = kFALSE;
  if (pKplus <= pLow) {
    passPlus = !kPlus.tofMatch || (kPlus.tofMatch && inKaonMass2(kPlus.mass2));
  } else {
    passPlus = kPlus.tofMatch && inKaonMass2(kPlus.mass2);
  }
  Bool_t passMinus = kFALSE;
  if (pKminus <= pLow) {
    passMinus = kMinus.tofMatch && inKaonMass2(kMinus.mass2);
  } else {
    passMinus = inKaonMass2(kMinus.mass2);
  }
  return passPlus && passMinus;
}

FemtoCandidate StFemtoMaker::MakeProtonCandidate(const TrackState& trk, Int_t eventIndex,
                                                 const std::string& speciesKey) const {
  FemtoCandidate cand;
  cand.eventIndex = eventIndex;
  cand.source = kFemtoCandTrack;
  cand.speciesKey = speciesKey;
  cand.charge = trk.charge;
  TVector3 p = TrackMomentum(trk);
  TLorentzVector p4 = ProtonP4(p);
  cand.SetP4(p4);
  cand.trk.trackIndex = trk.trackIndex;
  cand.trk.nSigmaProton = trk.nSigmaProton;
  cand.trk.nSigmaKaon = trk.nSigmaKaon;
  cand.trk.mass2 = trk.mass2;
  cand.trk.dca = trk.DCA;
  cand.trk.nHitsFit = trk.nHitsFit;
  return cand;
}

FemtoCandidate StFemtoMaker::MakePhiCandidate(const TrackState& kPlus, const TrackState& kMinus, Double_t invMass,
                                              const TVector3& phiMom, Double_t openingAngle, Double_t pairRapidity,
                                              Double_t dcaKK, Int_t eventIndex, const std::string& speciesKey) const {
  FemtoCandidate cand;
  cand.eventIndex = eventIndex;
  cand.source = kFemtoCandResonance;
  cand.speciesKey = speciesKey;
  cand.charge = 0;
  Double_t E = TMath::Sqrt(invMass * invMass + phiMom.Mag2());
  TLorentzVector p4(phiMom.X(), phiMom.Y(), phiMom.Z(), E);
  cand.SetP4(p4);
  cand.y = (Float_t)pairRapidity;
  cand.reso.invMass = (Float_t)invMass;
  cand.reso.dcaDaughters = (Float_t)dcaKK;
  cand.reso.dau1Index = kPlus.trackIndex;
  cand.reso.dau2Index = kMinus.trackIndex;
  (void)openingAngle;
  return cand;
}

void StFemtoMaker::BuildTrackPidCandidates(const std::string& speciesKey, const std::string& particleKey,
                                           const std::vector<TrackState>& tracks, Int_t eventIndex) {
  if (particleKey != "proton") {
    std::cerr << "[StFemtoMaker] TrackPidBuilder: unsupported particleKey '" << particleKey << "'" << std::endl;
    return;
  }
  std::vector<FemtoCandidate>& out = m_eventCandidates[speciesKey];
  for (size_t i = 0; i < tracks.size(); i++) {
    out.push_back(MakeProtonCandidate(tracks[i], eventIndex, speciesKey));
  }
}

void StFemtoMaker::BuildResonanceCandidates(const std::string& speciesKey, const std::string& particleKey,
                                            const std::vector<TrackState>& kaonsPlus,
                                            const std::vector<TrackState>& kaonsMinus, Int_t eventIndex) {
  if (particleKey != "phi") {
    std::cerr << "[StFemtoMaker] ResonanceBuilder: unsupported particleKey '" << particleKey << "'" << std::endl;
    return;
  }
  PhiCutConfig& phiCfg = ConfigManager::GetInstance().GetPhiCuts();
  TrackCutConfig& trkCfg = ConfigManager::GetInstance().GetTrackCuts();
  std::vector<FemtoCandidate>& out = m_eventCandidates[speciesKey];
  Int_t nPhiPairStage0 = 0;
  Int_t nPhiPairTofStrict = 0;

  for (size_t iPlus = 0; iPlus < kaonsPlus.size(); iPlus++) {
    for (size_t iMinus = 0; iMinus < kaonsMinus.size(); iMinus++) {
      TVector3 dcaPosPlus, dcaPosMinus;
      Double_t invMass;
      TVector3 phiMom;
      if (!ReconstructPhi(kaonsPlus[iPlus], kaonsMinus[iMinus], invMass, phiMom, dcaPosPlus, dcaPosMinus)) continue;

      Double_t openingAngle = CalculateOpeningAngle(kaonsPlus[iPlus], kaonsMinus[iMinus]);
      Double_t yLab = CalculatePairRapidity(invMass, phiMom);
      Double_t pairRapidity = ApplyRapidityFrame(yLab);
      Double_t dcaKK = CalculateDCA(kaonsPlus[iPlus], kaonsMinus[iMinus], dcaPosPlus, dcaPosMinus);

      const Bool_t passStage =
          (TMath::Abs(kaonsPlus[iPlus].DCA) <= phiCfg.maxDCAKaon) &&
          (TMath::Abs(kaonsMinus[iMinus].DCA) <= phiCfg.maxDCAKaon) &&
          (kaonsPlus[iPlus].pT >= trkCfg.minPt && kaonsPlus[iPlus].pT <= trkCfg.maxPt) &&
          (kaonsMinus[iMinus].pT >= trkCfg.minPt && kaonsMinus[iMinus].pT <= trkCfg.maxPt) &&
          (kaonsPlus[iPlus].eta >= trkCfg.minEta && kaonsPlus[iPlus].eta <= trkCfg.maxEta) &&
          (kaonsMinus[iMinus].eta >= trkCfg.minEta && kaonsMinus[iMinus].eta <= trkCfg.maxEta) &&
          (TMath::Abs(kaonsPlus[iPlus].nSigmaKaon) <= phiCfg.nSigmaKaon) &&
          (TMath::Abs(kaonsMinus[iMinus].nSigmaKaon) <= phiCfg.nSigmaKaon) &&
          (openingAngle >= phiCfg.minOpeningAngle && openingAngle <= phiCfg.maxOpeningAngle) &&
          (pairRapidity >= phiCfg.minPairRapidity && pairRapidity <= phiCfg.maxPairRapidity);
      if (!passStage) continue;
      if (m_histManager) {
        nPhiPairStage0++;
        m_histManager->Fill("hPhiPair_Mass_stage0", invMass);
        m_histManager->Fill("hPhiPair_Pt_stage0", phiMom.Pt());
        m_histManager->Fill("hPhiPair_Rapidity_stage0", pairRapidity);
        m_histManager->Fill("hPhiPair_OpeningAngle_stage0", openingAngle);
      }

      if (!PassPairTofCut(kaonsPlus[iPlus], kaonsMinus[iMinus])) continue;

      if (m_histManager) {
        nPhiPairTofStrict++;
        m_histManager->Fill("hPhiPair_Mass_tofStrict", invMass);
        m_histManager->Fill("hPhiPair_Pt_tofStrict", phiMom.Pt());
        m_histManager->Fill("hPhiPair_Rapidity_tofStrict", pairRapidity);
        m_histManager->Fill("hPhiPair_OpeningAngle_tofStrict", openingAngle);
        m_histManager->Fill("hPhiPair_MassVsPt_tofStrict", phiMom.Pt(), invMass);
        if (m_cent9 >= 0 && m_cent9 <= 8) {
          m_histManager->Fill("hPhiPair_MassVsCent_tofStrict", (Double_t)m_cent9, invMass);
        }
      }

      out.push_back(MakePhiCandidate(kaonsPlus[iPlus], kaonsMinus[iMinus], invMass, phiMom, openingAngle,
                                     pairRapidity, dcaKK, eventIndex, speciesKey));
    }
  }
  if (m_histManager) {
    m_histManager->Fill("hPhiPair_NPairs_stage0", (Double_t)nPhiPairStage0);
    m_histManager->Fill("hPhiPair_NPairs_tofStrict", (Double_t)nPhiPairTofStrict);
  }
}

void StFemtoMaker::FillCentralityEventQA(Int_t cent9, Int_t rawMult, Double_t refMultCorr, Int_t nTracks,
                                         Int_t nBTOFMatch, Int_t nKaonPlus, Int_t nKaonMinus, Int_t nPhiCandidates,
                                         Int_t nProtons) {
  if (!m_histManager || cent9 < 0) return;
  const Double_t centX = (Double_t)cent9;
  m_histManager->Fill("hRawMult_vs_Cent9", centX, (Double_t)rawMult);
  if (refMultCorr >= 0.0) {
    m_histManager->Fill("hRefMultCorr_vs_Cent9", centX, refMultCorr);
    m_histManager->Fill("hRawMult_vs_RefMultCorr", refMultCorr, (Double_t)rawMult);
  }
  m_histManager->Fill("hNTracks_vs_Cent9", centX, (Double_t)nTracks);
  m_histManager->Fill("hTofMatchMult_vs_Cent9", centX, (Double_t)nBTOFMatch);
  m_histManager->Fill("hNKaonPlus_vs_Cent9", centX, (Double_t)nKaonPlus);
  m_histManager->Fill("hNKaonMinus_vs_Cent9", centX, (Double_t)nKaonMinus);
  m_histManager->Fill("hNPhiPairs_vs_Cent9", centX, (Double_t)nPhiCandidates);
  m_histManager->Fill("hNProton_vs_Cent9", centX, (Double_t)nProtons);
}

TLorentzVector StFemtoMaker::ProtonP4(const TVector3& p) const {
  Double_t e = TMath::Sqrt(kProtonMass * kProtonMass + p.Mag2());
  return TLorentzVector(p.X(), p.Y(), p.Z(), e);
}

TLorentzVector StFemtoMaker::CandidateP4(const FemtoCandidate& cand) const {
  if (cand.source == kFemtoCandResonance) {
    return cand.P4();
  }
  TVector3 p(cand.px, cand.py, cand.pz);
  return ProtonP4(p);
}

Double_t StFemtoMaker::ComputeKStar(const TLorentzVector& pA, const TLorentzVector& pB) const {
  TLorentzVector q = pA - pB;
  TLorentzVector pair = pA + pB;
  q.Boost(-pair.BoostVector());
  return 0.5 * q.Vect().Mag();
}

Bool_t StFemtoMaker::TracksOverlap(const FemtoCandidate& phiCand, const FemtoCandidate& trkCand) const {
  if (phiCand.source != kFemtoCandResonance) return kFALSE;
  Int_t idx = trkCand.trk.trackIndex;
  return (idx == phiCand.reso.dau1Index || idx == phiCand.reso.dau2Index);
}

std::string StFemtoMaker::HistName(const std::string& prefix, const std::string& channelName) const {
  return prefix + "_" + channelName;
}

void StFemtoMaker::FillSameEventPairs(const FemtoConfig::ChannelDef& ch) {
  FemtoCandidateStore::const_iterator itA = m_eventCandidates.find(ch.partA);
  FemtoCandidateStore::const_iterator itB = m_eventCandidates.find(ch.partB);
  if (itA == m_eventCandidates.end() || itB == m_eventCandidates.end()) return;

  const std::vector<FemtoCandidate>& candsA = itA->second;
  const std::vector<FemtoCandidate>& candsB = itB->second;
  if (candsA.empty() || candsB.empty()) return;

  std::string hName = HistName("hKstarSE", ch.name);
  for (size_t i = 0; i < candsA.size(); i++) {
    const FemtoCandidate& a = candsA[i];
    if (a.source == kFemtoCandResonance) {
      if (a.reso.invMass < ch.signalMin || a.reso.invMass > ch.signalMax) continue;
    }
    for (size_t j = 0; j < candsB.size(); j++) {
      const FemtoCandidate& b = candsB[j];
      if (TracksOverlap(a, b)) continue;
      Double_t kstar = ComputeKStar(CandidateP4(a), CandidateP4(b));
      if (m_histManager) m_histManager->Fill(hName.c_str(), kstar);
    }
  }
}

Int_t StFemtoMaker::GetMixingBin(Float_t vz, Int_t cent9, Double_t psi2) const {
  const EventCutConfig& ev = ConfigManager::GetInstance().GetEventCuts();
  const MixingConfig& mix = ConfigManager::GetInstance().GetMixingConfig();
  Int_t vzBin = 0;
  if (mix.nVzBins > 0) {
    Double_t vzSpan = ev.maxVz - ev.minVz;
    if (vzSpan > 0) {
      vzBin = (Int_t)((vz - ev.minVz) / vzSpan * mix.nVzBins);
      if (vzBin < 0) vzBin = 0;
      if (vzBin >= mix.nVzBins) vzBin = mix.nVzBins - 1;
    }
  }
  Int_t centBin = 0;
  if (mix.nCentralityBins > 0 && cent9 >= 0) {
    centBin = cent9;
    if (centBin >= mix.nCentralityBins) centBin = mix.nCentralityBins - 1;
  }
  Int_t epBin = 0;
  if (mix.nEventPlaneBins > 1 && psi2 >= 0) {
    epBin = (Int_t)(psi2 / TMath::Pi() * mix.nEventPlaneBins);
    if (epBin < 0) epBin = 0;
    if (epBin >= mix.nEventPlaneBins) epBin = mix.nEventPlaneBins - 1;
  }
  return vzBin + mix.nVzBins * (centBin + mix.nCentralityBins * epBin);
}

void StFemtoMaker::FillMixedEventPairs(const FemtoConfig::ChannelDef& ch, Float_t vz, Int_t cent9, Double_t psi2) {
  FemtoCandidateStore::const_iterator itA = m_eventCandidates.find(ch.partA);
  FemtoCandidateStore::const_iterator itB = m_eventCandidates.find(ch.partB);
  if (itA == m_eventCandidates.end() || itB == m_eventCandidates.end()) return;
  const std::vector<FemtoCandidate>& candsA = itA->second;
  const std::vector<FemtoCandidate>& candsB = itB->second;
  if (candsA.empty() || candsB.empty()) return;

  Int_t mixBin = GetMixingBin(vz, cent9, psi2);
  std::map<Int_t, std::deque<FemtoMixingEvent> >::const_iterator poolIt = m_mixingPool.find(mixBin);
  if (poolIt == m_mixingPool.end() || poolIt->second.empty()) return;

  std::string hName = HistName("hKstarME", ch.name);
  Int_t nPairs = (Int_t)candsA.size() * (Int_t)candsB.size();
  if (nPairs > 500) nPairs = 500;
  if (nPairs < 1) nPairs = 1;

  for (Int_t ip = 0; ip < nPairs; ip++) {
    const FemtoMixingEvent& mixEvt = poolIt->second[(Int_t)gRandom->Uniform(0, poolIt->second.size())];
    FemtoCandidateStore::const_iterator mixB = mixEvt.candidates.find(ch.partB);
    if (mixB == mixEvt.candidates.end() || mixB->second.empty()) continue;

    const FemtoCandidate& a = candsA[(Int_t)gRandom->Uniform(0, candsA.size())];
    if (a.source == kFemtoCandResonance) {
      if (a.reso.invMass < ch.signalMin || a.reso.invMass > ch.signalMax) continue;
    }
    const FemtoCandidate& b = mixB->second[(Int_t)gRandom->Uniform(0, mixB->second.size())];
    Double_t kstar = ComputeKStar(CandidateP4(a), CandidateP4(b));
    if (m_histManager) m_histManager->Fill(hName.c_str(), kstar);
  }
}

void StFemtoMaker::StoreEventForMixing(Float_t vz, Int_t cent9, Double_t psi2) {
  if (m_eventCandidates.empty()) return;
  Int_t mixBin = GetMixingBin(vz, cent9, psi2);
  FemtoMixingEvent evt;
  evt.candidates = m_eventCandidates;
  std::deque<FemtoMixingEvent>& pool = m_mixingPool[mixBin];
  pool.push_back(evt);
  const MixingConfig& mix = ConfigManager::GetInstance().GetMixingConfig();
  while ((Int_t)pool.size() > mix.bufferSize) pool.pop_front();
}

void StFemtoMaker::FillCandidateQA() {
  if (!m_histManager) return;
  FemtoCandidateStore::const_iterator itP = m_eventCandidates.find("proton");
  if (itP != m_eventCandidates.end()) {
    for (size_t i = 0; i < itP->second.size(); i++) {
      const FemtoCandidate& c = itP->second[i];
      m_histManager->Fill("hP_Pt", c.pt);
      m_histManager->Fill("hP_Eta", c.eta);
      m_histManager->Fill("hP_Phi", c.phi);
      m_histManager->Fill("hP_NSigmaProton", c.trk.nSigmaProton);
      m_histManager->Fill("hP_DCA", c.trk.dca);
      if (c.trk.mass2 > -900) m_histManager->Fill("hP_Mass2", c.trk.mass2);
    }
    m_histManager->Fill("hP_NCand", (Double_t)itP->second.size());
  }
  FemtoCandidateStore::const_iterator itPhi = m_eventCandidates.find("phi");
  if (itPhi != m_eventCandidates.end()) {
    for (size_t i = 0; i < itPhi->second.size(); i++) {
      const FemtoCandidate& c = itPhi->second[i];
      m_histManager->Fill("hPhi_MKK", c.reso.invMass);
      m_histManager->Fill("hPhi_Pt", c.pt);
      m_histManager->Fill("hPhi_Rapidity", c.y);
    }
    m_histManager->Fill("hPhi_NCand", (Double_t)itPhi->second.size());
  }
}

void StFemtoMaker::FinalizeCorrelationFunctions() {
  if (!m_histManager) return;
  FemtoConfig& femtoCfg = ConfigManager::GetInstance().GetFemtoConfig();
  for (size_t ic = 0; ic < femtoCfg.channels.size(); ic++) {
    const FemtoConfig::ChannelDef& ch = femtoCfg.channels[ic];
    if (!ch.enabled) continue;
    std::string seName = HistName("hKstarSE", ch.name);
    std::string meName = HistName("hKstarME", ch.name);
    std::string cfName = HistName("hCF", ch.name);
    TH1* hSE = m_histManager->Get(seName.c_str());
    TH1* hME = m_histManager->Get(meName.c_str());
    TH1* hCF = m_histManager->Get(cfName.c_str());
    if (!hSE || !hME || !hCF) continue;

    Int_t binLo = hSE->FindBin(ch.normQMin + 1e-9);
    Int_t binHi = hSE->FindBin(ch.normQMax - 1e-9);
    Double_t seNorm = hSE->Integral(binLo, binHi);
    Double_t meNorm = hME->Integral(binLo, binHi);
    if (seNorm <= 0 || meNorm <= 0) continue;

    hCF->Reset();
    hCF->Add(hSE);
    TH1* hMEScaled = (TH1*)hME->Clone("_me_scaled");
    hMEScaled->Scale(seNorm / meNorm);
    hCF->Divide(hMEScaled);
    delete hMEScaled;
  }
}
