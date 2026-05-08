#include "StLambdaMaker.h"
#include "ConfigManager.h"
#include "HistManager.h"
#include "cuts/EventCutConfig.h"
#include "cuts/LambdaCutConfig.h"
#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StarClassLibrary/StPhysicalHelixD.hh"
#include "StarClassLibrary/StThreeVectorF.hh"
#include "StarClassLibrary/SystemOfUnits.h"

#include "TFile.h"
#include "TMath.h"
#include "TVector3.h"
#include "TLorentzVector.h"

#include <iostream>
#include <utility>

namespace {
  const Double_t kProtonMass = 0.938272;
  const Double_t kPionMass   = 0.139570;
}

//-----------------------------------------------------------------------------
StLambdaMaker::StLambdaMaker(const char* name, StPicoDstMaker* picoMaker, const char* outName)
  : StMaker(name),
    mPicoDstMaker(picoMaker),
    mPicoDst(0),
    mOutName(outName),
    mEventCounter(0),
    m_histManager(0) {}

//-----------------------------------------------------------------------------
StLambdaMaker::~StLambdaMaker() {
  if (m_histManager) {
    delete m_histManager;
    m_histManager = 0;
  }
}

//-----------------------------------------------------------------------------
StLambdaMaker* createStLambdaMaker(const char* name, StPicoDstMaker* picoMaker, const char* outName) {
  return new StLambdaMaker(name, picoMaker, outName);
}

extern "C" void* createStLambdaMakerC(const char* name, void* picoMaker, const char* outName) {
  return (void*)createStLambdaMaker(name, (StPicoDstMaker*)picoMaker, outName);
}

//-----------------------------------------------------------------------------
Int_t StLambdaMaker::Init() {
  std::string histPath = ConfigManager::GetInstance().GetHistConfigPath();
  if (histPath.empty()) {
    std::cerr << "[StLambdaMaker] GetHistConfigPath() returned empty; no histograms will be filled." << std::endl;
    m_histManager = 0;
    return kStOK;
  }
  m_histManager = new HistManager();
  if (!m_histManager->LoadFromFile(histPath.c_str())) {
    std::cerr << "[StLambdaMaker] Failed to load hist config from " << histPath << std::endl;
    delete m_histManager;
    m_histManager = 0;
    return kStOK;
  }
  return kStOK;
}

//-----------------------------------------------------------------------------
void StLambdaMaker::Clear(Option_t* opt) {}

//-----------------------------------------------------------------------------
Bool_t StLambdaMaker::PassEventCuts(Int_t nTracks) {
  EventCutConfig& ev = ConfigManager::GetInstance().GetEventCuts();
  if (ev.maxNTr > 0 && nTracks > ev.maxNTr) return kFALSE;
  return kTRUE;
}

//-----------------------------------------------------------------------------
Bool_t StLambdaMaker::PassProtonCuts(StPicoTrack* trk, const TVector3& pVtx) {
  if (!trk || trk->charge() <= 0) return kFALSE;
  LambdaCutConfig& lam = ConfigManager::GetInstance().GetLambdaCuts();
  if (TMath::Abs(trk->nSigmaProton()) > lam.nSigmaProton) return kFALSE;
  Double_t dca = trk->gDCA(pVtx.X(), pVtx.Y(), pVtx.Z());
  if (dca < lam.minDCAProton) return kFALSE;
  return kTRUE;
}

//-----------------------------------------------------------------------------
Bool_t StLambdaMaker::PassPionCuts(StPicoTrack* trk, const TVector3& pVtx) {
  if (!trk || trk->charge() >= 0) return kFALSE;
  LambdaCutConfig& lam = ConfigManager::GetInstance().GetLambdaCuts();
  if (TMath::Abs(trk->nSigmaPion()) > lam.nSigmaPion) return kFALSE;
  Double_t dca = trk->gDCA(pVtx.X(), pVtx.Y(), pVtx.Z());
  if (dca < lam.minDCAPion) return kFALSE;
  return kTRUE;
}

//-----------------------------------------------------------------------------
StPhysicalHelixD StLambdaMaker::MakeHelix(StPicoTrack* trk, Double_t bField) {
  StThreeVectorF p(trk->gMom().X(), trk->gMom().Y(), trk->gMom().Z());
  StThreeVectorF o(trk->origin().X(), trk->origin().Y(), trk->origin().Z());
  return StPhysicalHelixD(p, o, bField * units::kilogauss, (Float_t)trk->charge());
}

//-----------------------------------------------------------------------------
Bool_t StLambdaMaker::MakeLambdaHelix(StPicoTrack* p, StPicoTrack* pi, Double_t bField,
                                      TVector3& v0, TVector3& momP, TVector3& momPi, Double_t& dca12) {
  LambdaCutConfig& lam = ConfigManager::GetInstance().GetLambdaCuts();
  StPhysicalHelixD hp  = MakeHelix(p,  bField);
  StPhysicalHelixD hpi = MakeHelix(pi, bField);

  std::pair<Double_t, Double_t> s = hp.pathLengths(hpi);
  if (TMath::Abs(s.first) > lam.maxPathLength || TMath::Abs(s.second) > lam.maxPathLength)
    return kFALSE;

  StThreeVectorD dcaA = hp.at(s.first);
  StThreeVectorD dcaB = hpi.at(s.second);

  StThreeVectorD v0_( (dcaA.x() + dcaB.x()) * 0.5,
                      (dcaA.y() + dcaB.y()) * 0.5,
                      (dcaA.z() + dcaB.z()) * 0.5 );
  dca12 = (dcaA - dcaB).mag();

  if (dca12 < 0 || dca12 > lam.maxDaughterDCA) return kFALSE;

  StThreeVectorD pp  = hp.momentumAt(s.first,  bField * units::kilogauss);
  StThreeVectorD ppi = hpi.momentumAt(s.second, bField * units::kilogauss);

  momP.SetXYZ(pp.x(), pp.y(), pp.z());
  momPi.SetXYZ(ppi.x(), ppi.y(), ppi.z());
  v0.SetXYZ(v0_.x(), v0_.y(), v0_.z());

  return kTRUE;
}

//-----------------------------------------------------------------------------
Int_t StLambdaMaker::Make() {
  if (!mPicoDstMaker) return kStWarn;
  mPicoDst = mPicoDstMaker->picoDst();
  if (!mPicoDst) return kStWarn;

  StPicoEvent* event = mPicoDst->event();
  if (!event) return kStWarn;

  mEventCounter++;

  TVector3 pVtx = event->primaryVertex();
  Int_t nTr = mPicoDst->numberOfTracks();

  if (m_histManager) {
    m_histManager->Fill("hVz", pVtx.Z());
    m_histManager->Fill("hRefMult", event->refMult());
  }

  if (!PassEventCuts(nTr)) return kStOK;

  Double_t bField = event->bField();

  for (Int_t ip = 0; ip < nTr; ip++) {
    StPicoTrack* p = mPicoDst->track(ip);
    if (!p) continue;
    if (!PassProtonCuts(p, pVtx)) continue;

    for (Int_t ii = 0; ii < nTr; ii++) {
      if (ii == ip) continue;
      StPicoTrack* pi = mPicoDst->track(ii);
      if (!pi) continue;
      if (!PassPionCuts(pi, pVtx)) continue;

      TVector3 v0, momP, momPi;
      Double_t dca12 = 0;
      if (!MakeLambdaHelix(p, pi, bField, v0, momP, momPi, dca12)) continue;

      TVector3 pLam = momP + momPi;
      Double_t pLamMag = pLam.Mag();
      if (pLamMag < 1e-5) continue;

      TVector3 pLamUnit = pLam * (1.0 / pLamMag);
      TVector3 diff = pVtx - v0;
      Double_t dcaV0 = (diff.Cross(pLamUnit)).Mag();

      LambdaCutConfig& lam = ConfigManager::GetInstance().GetLambdaCuts();
      if (dcaV0 > lam.maxDCAV0) continue;

      TVector3 flight = v0 - pVtx;
      Double_t cosPoint = flight.Dot(pLam) / (flight.Mag() * pLam.Mag() + 1e-10);
      if (cosPoint < lam.minCosPointing) continue;

      TLorentzVector lp, lpi;
      lp.SetVectM(momP,  kProtonMass);
      lpi.SetVectM(momPi, kPionMass);
      Double_t invMass = (lp + lpi).M();

      if (m_histManager) {
        m_histManager->Fill("hLambda_InvMass", invMass);
        m_histManager->Fill("hLambda_Pt", pLam.Pt());
        m_histManager->Fill("hLambda_Eta", pLam.PseudoRapidity());
        m_histManager->Fill("hLambda_Phi", pLam.Phi());
        m_histManager->Fill("hDCA12", dca12);
        m_histManager->Fill("hDCAV0", dcaV0);
        m_histManager->Fill("hCosPointing", cosPoint);
        m_histManager->Fill("hNSigmaProton", p->nSigmaProton());
        m_histManager->Fill("hNSigmaPion", pi->nSigmaPion());
        m_histManager->Fill("hLambda_InvMass_vs_Pt", pLam.Pt(), invMass);
        m_histManager->Fill("hDCAV0_vs_InvMass", invMass, dcaV0);
        m_histManager->Fill("hCosPointing_vs_InvMass", invMass, cosPoint);
      }
    }
  }

  if (m_histManager) m_histManager->Fill("hN", 0);
  return kStOK;
}

//-----------------------------------------------------------------------------
Int_t StLambdaMaker::Finish() {
  // No HistManager => no histograms to write; skip creating empty ROOT output.
  if (mOutName != "" && m_histManager) {
    TFile* fout = new TFile(mOutName.Data(), "RECREATE");
    fout->cd();
    WriteHistograms();
    m_histManager->ReleaseOwnership();
    fout->Close();
    delete fout;
  }
  std::cout << "StLambdaMaker::Finish() processed " << mEventCounter << " events" << std::endl;
  return kStOK;
}

//-----------------------------------------------------------------------------
void StLambdaMaker::WriteHistograms() {
  if (m_histManager) m_histManager->Write();
}
