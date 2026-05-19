#include "StNuclearIdMaker.h"
#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"
#include "ConfigManager.h"
#include "cuts/NuclearIdCutConfig.h"
#include "HistManager.h"
#include "TH2D.h"
#include "TFile.h"
#include "TVector3.h"
#include "TMath.h"
#include <iostream>

// Factory functions
StNuclearIdMaker* createStNuclearIdMaker(const char* name, StPicoDstMaker* picoMaker, const char* outName) {
  return new StNuclearIdMaker(name, picoMaker, outName);
}
extern "C" void* createStNuclearIdMakerC(const char* name, void* picoMaker, const char* outName) {
  return (void*)createStNuclearIdMaker(name, (StPicoDstMaker*)picoMaker, outName);
}

// NuclearPID struct inside the implementation
struct NuclearPID {
  enum ParticleType { kDeuteron = 0, kTriton, kHe3, kHe4 };
  double fMeanParams[4][4];
  double fSigmaParams[4][4];
  double fRange[4][2];

  NuclearPID() {
    double mean[4][4] = {
      {10.2875, -36.6639, 67.63, -45.2483},
      {19.2177, -70.3198, 113.876, -67.0591},
      {28.8747, -68.8637, 104.455, -56.9342},
      {39.974, -103.365, 139.182, -67.1446}
    };
    double sigma[4][4] = {
      {0.771189, -2.63609, 2.8044, 0.0},
      {1.40996, -4.97147, 4.97608, 0.0},
      {1.636, -4.35588, 4.40594, 0.0},
      {2.24127, -7.08328, 8.07124, 0.0}
    };
    double range[4][2] = {
      {0.2, 6.0}, {0.2, 6.0}, {0.4, 5.0}, {0.4, 5.0}
    };
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        fMeanParams[i][j] = mean[i][j];
        fSigmaParams[i][j] = sigma[i][j];
      }
      fRange[i][0] = range[i][0];
      fRange[i][1] = range[i][1];
    }
  }

  bool IsValid(ParticleType type, double p) const {
    int t = (int)type;
    return (p >= fRange[t][0] && p <= fRange[t][1]);
  }

  double LogPoly(double p, const double* par) const {
    if (p <= 0) return 0;
    double lp = log10(p);
    return par[0] + par[1]*lp + par[2]*lp*lp + par[3]*lp*lp*lp;
  }

  double GetMean(ParticleType type, double p) const {
    if (!IsValid(type, p)) return 0;
    return LogPoly(p, fMeanParams[(int)type]);
  }

  double GetSigma(ParticleType type, double p) const {
    if (!IsValid(type, p)) return 0;
    return LogPoly(p, fSigmaParams[(int)type]);
  }

  double GetNSigma(ParticleType type, double p, double dedx) const {
    if (!IsValid(type, p)) return 999.0;
    double mean = GetMean(type, p);
    double sigma = GetSigma(type, p);
    if (sigma <= 0) return 999.0;
    return (dedx - mean) / sigma;
  }
};

StNuclearIdMaker::StNuclearIdMaker(const char* name, StPicoDstMaker* picoMaker, const char* outName)
    : StMaker(name), mPicoDstMaker(picoMaker), mPicoDst(nullptr), mOutName(outName), m_histManager(0) {}

StNuclearIdMaker::~StNuclearIdMaker() {
  if (m_histManager) {
    delete m_histManager;
    m_histManager = 0;
  }
}

Int_t StNuclearIdMaker::Init() {
  std::string histPath = ConfigManager::GetInstance().GetHistConfigPath();
  if (histPath.empty()) {
    std::cerr << "[StNuclearIdMaker] GetHistConfigPath() returned empty; no histograms will be filled." << std::endl;
    m_histManager = 0;
    return kStOK;
  }
  m_histManager = new HistManager();
  if (!m_histManager->LoadFromFile(histPath.c_str())) {
    std::cerr << "[StNuclearIdMaker] Failed to load hist config from " << histPath << std::endl;
    delete m_histManager;
    m_histManager = 0;
    return kStOK;
  }

  return kStOK;
}

void StNuclearIdMaker::Clear(Option_t* opt) {
  StMaker::Clear(opt);
}

Int_t StNuclearIdMaker::Make() {
  mPicoDst = mPicoDstMaker->picoDst();
  if (!mPicoDst) return kStWarn;
  StPicoEvent* event = mPicoDst->event();
  if (!event) return kStWarn;

  NuclearPID pid;
  NuclearIdCutConfig& cuts = ConfigManager::GetInstance().GetNuclearIdCuts();

  double d_mean   = 3.48096e+00;
  double d_sigma  = 1.41458e-01;
  double He3_mean = 1.92385e+00;
  double He3_sigma= 0.94e-01;
  double t_mean   = 7.76906e+00;
  double t_sigma  = 3.41755e-01;

  Int_t nTracks = mPicoDst->numberOfTracks();
  for (Int_t it = 0; it < nTracks; it++) {
    StPicoTrack* trk = mPicoDst->track(it);
    if (!trk) continue;

    if (trk->nHitsDedx() < 15) continue;
    if (trk->gPt() < 0.1) continue;

    TVector3 pMom = trk->pMom();
    Double_t p = pMom.Mag();

    Double_t dedx = trk->dEdx();
    if (dedx <= 0) continue;

    if (m_histManager) {
      m_histManager->Fill("hDedxP", p, dedx);

      if (fabs(trk->nSigmaElectron()) < cuts.nSigmaFill1D) m_histManager->Fill("hDedxP_e", p, dedx);
      if (fabs(trk->nSigmaPion()) < cuts.nSigmaFill1D) m_histManager->Fill("hDedxP_pi", p, dedx);
      if (fabs(trk->nSigmaKaon()) < cuts.nSigmaFill1D) m_histManager->Fill("hDedxP_K", p, dedx);
      if (fabs(trk->nSigmaProton()) < cuts.nSigmaFill1D) m_histManager->Fill("hDedxP_p", p, dedx);

      if (fabs(trk->nSigmaElectron()) > cuts.nSigmaExclude &&
          fabs(trk->nSigmaPion()) > cuts.nSigmaExclude &&
          fabs(trk->nSigmaKaon()) > cuts.nSigmaExclude &&
          fabs(trk->nSigmaProton()) > cuts.nSigmaExclude) {
        m_histManager->Fill("hDedxP_else", p, dedx);
      }
    }

    double nSigma_d   = pid.GetNSigma(NuclearPID::kDeuteron, p, dedx);
    double nSigma_t   = pid.GetNSigma(NuclearPID::kTriton,  p, dedx);
    double nSigma_3He = pid.GetNSigma(NuclearPID::kHe3,     p, dedx);
    double nSigma_4He = pid.GetNSigma(NuclearPID::kHe4,     p, dedx);

    if (m_histManager) {
      if (fabs(nSigma_d)   < cuts.maxNSigmaNuclear) m_histManager->Fill("hDedxP_d", p, dedx);
      if (fabs(nSigma_t)   < cuts.maxNSigmaNuclear) m_histManager->Fill("hDedxP_t", p, dedx);
      if (fabs(nSigma_3He) < cuts.maxNSigmaNuclear) m_histManager->Fill("hDedxP_3He", p, dedx);
      if (fabs(nSigma_4He) < cuts.maxNSigmaNuclear) m_histManager->Fill("hDedxP_4He", p, dedx);
    }

    const double M_d   = 1.87561;
    const double M_t   = 2.80892;
    const double M_3He = 2.80839;
    const double M_4He = 3.72742;
    double eta = pMom.PseudoRapidity();
    double pT  = pMom.Perp();
    double pz  = pMom.Z();

    if (fabs(nSigma_d) < cuts.maxNSigmaNuclear) {
      double E = TMath::Sqrt(M_d*M_d + p*p);
      double y = 0.5 * TMath::Log((E + pz) / (E - pz));
      if (m_histManager) {
        m_histManager->Fill("hPvsEta_d", p, eta);
        m_histManager->Fill("hPvsY_d", p, y);
        m_histManager->Fill("hPvsPt_d", p, pT);
      }
    }
    if (fabs(nSigma_t) < cuts.maxNSigmaNuclear) {
      double E = TMath::Sqrt(M_t*M_t + p*p);
      double y = 0.5 * TMath::Log((E + pz) / (E - pz));
      if (m_histManager) {
        m_histManager->Fill("hPvsEta_t", p, eta);
        m_histManager->Fill("hPvsY_t", p, y);
        m_histManager->Fill("hPvsPt_t", p, pT);
      }
    }
    if (fabs(nSigma_3He) < cuts.maxNSigmaNuclear) {
      double p_act = 2.0 * p;
      double pz_act = 2.0 * pz;
      double pT_act = 2.0 * pT;
      double E = TMath::Sqrt(M_3He*M_3He + p_act*p_act);
      double y = 0.5 * TMath::Log((E + pz_act) / (E - pz_act));
      if (m_histManager) {
        m_histManager->Fill("hPvsEta_3He", p_act, eta);
        m_histManager->Fill("hPvsY_3He", p_act, y);
        m_histManager->Fill("hPvsPt_3He", p_act, pT_act);
      }
    }
    if (fabs(nSigma_4He) < cuts.maxNSigmaNuclear) {
      double p_act = 2.0 * p;
      double pz_act = 2.0 * pz;
      double pT_act = 2.0 * pT;
      double E = TMath::Sqrt(M_4He*M_4He + p_act*p_act);
      double y = 0.5 * TMath::Log((E + pz_act) / (E - pz_act));
      if (m_histManager) {
        m_histManager->Fill("hPvsEta_4He", p_act, eta);
        m_histManager->Fill("hPvsY_4He", p_act, y);
        m_histManager->Fill("hPvsPt_4He", p_act, pT_act);
      }
    }

    if (!trk->isTofTrack()) {
      if (m_histManager) m_histManager->Fill("hDedxP_cut", p, dedx);
      continue;
    }

    int idx = trk->bTofPidTraitsIndex();
    if (idx < 0) {
      if (m_histManager) m_histManager->Fill("hDedxP_cut", p, dedx);
      continue;
    }
    StPicoBTofPidTraits* btof = mPicoDst->btofPidTraits(idx);
    if (!btof) {
      if (m_histManager) m_histManager->Fill("hDedxP_cut", p, dedx);
      continue;
    }
    double beta = btof->btofBeta();
    if (beta <= 0) {
      if (m_histManager) m_histManager->Fill("hDedxP_cut", p, dedx);
      continue;
    }

    double m2 = p*p*(1.0/(beta*beta) - 1.0);
    if (m_histManager) m_histManager->Fill("hM2P", p, m2);

    if (m_histManager) {
      if (m2 > d_mean   - cuts.m2SigmaCut * d_sigma   && m2 < d_mean   + cuts.m2SigmaCut * d_sigma)   m_histManager->Fill("hDedxP_d_m2", p, dedx);
      if (m2 > t_mean   - cuts.m2SigmaCut * t_sigma   && m2 < t_mean   + cuts.m2SigmaCut * t_sigma)   m_histManager->Fill("hDedxP_t_m2", p, dedx);
      if (m2 > He3_mean - cuts.m2SigmaCut * He3_sigma && m2 < He3_mean + cuts.m2SigmaCut * He3_sigma) m_histManager->Fill("hDedxP_3He_m2", p, dedx);
    }

    if (m2 > cuts.minM2_M2cut && m2 < cuts.maxM2_M2cut && p > cuts.minP_M2cut) continue;
    else if (m_histManager) m_histManager->Fill("hDedxP_cut", p, dedx);
  }

  return kStOK;
}

Int_t StNuclearIdMaker::Finish() {
  if (mOutName != "" && m_histManager) {
    TFile* fout = new TFile(mOutName.Data(), "RECREATE");
    fout->cd();
    WriteHistograms();
    m_histManager->ReleaseOwnership();
    fout->Close();
    delete fout;
  }
  std::cout << "StNuclearIdMaker::Finish() processed events" << std::endl;
  return kStOK;
}

void StNuclearIdMaker::WriteHistograms() {
  if (m_histManager) m_histManager->Write();
}
