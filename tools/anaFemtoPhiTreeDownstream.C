// anaFemtoPhiTreeDownstream.C
// Rebuild SE/ME phi-d diagnostics from a reduced tree. Does not modify StFemtoMaker.
// C++98 / ROOT 5 ACLiC.

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TMath.h"
#include "TString.h"
#include "TNamed.h"
#include "TSystem.h"
#include "TRandom3.h"
#include "ConfigManager.h"
#include "cuts/PhiCutConfig.h"
#include "cuts/MixingConfig.h"
#include "cuts/FemtoConfig.h"
#include "FemtoPhiTreeSchema.h"
#include "StPhiKKReconstruction.h"

#include <deque>
#include <iostream>
#include <map>
#include <set>
#include <vector>

namespace {
const Double_t kKaonMass = 0.493677;
const Double_t kDeuteronMass = 1.875612;

struct Cand {
  Int_t trackIndex;
  TLorentzVector p4;
  Float_t mKK;
  Int_t dau1;
  Int_t dau2;
};

struct MixEvent {
  ULong64_t eventUID;
  std::vector<Cand> phis;
  std::vector<Cand> deuterons;
};

struct MixRef {
  size_t ib;
  Int_t reverse;
  size_t i;
  size_t j;
};

Double_t KStar(const TLorentzVector& a, const TLorentzVector& b) {
  TLorentzVector q = a - b;
  TLorentzVector pair = a + b;
  q.Boost(-pair.BoostVector());
  return 0.5 * q.Vect().Mag();
}

PhiKkTrackState ToKk(const femto_phi_tree::TrackRow& t) {
  PhiKkTrackState s;
  TVector3 p(t.px, t.py, t.pz);
  s.pT = (Float_t)p.Pt();
  s.eta = (Float_t)p.Eta();
  s.phi = (Float_t)p.Phi();
  s.charge = t.charge;
  s.originX = t.originX;
  s.originY = t.originY;
  s.originZ = t.originZ;
  s.momentumX = t.px;
  s.momentumY = t.py;
  s.momentumZ = t.pz;
  s.BField = t.bField;
  s.tofMatch = t.tofMatch ? kTRUE : kFALSE;
  s.mass2 = t.mass2;
  s.deltaOneOverBeta = t.deltaOneOverBeta;
  return s;
}

Bool_t SharedTrack(const Cand& phi, const Cand& d) {
  return d.trackIndex >= 0 && (d.trackIndex == phi.dau1 || d.trackIndex == phi.dau2);
}

Bool_t PassKaonForPhi(const femto_phi_tree::TrackRow& t) {
  const UInt_t need = femto_phi_tree::kSelKaonCutsNom | femto_phi_tree::kSelLoosePid;
  return (t.selFlags & need) == need;
}

Bool_t PassDeuteronNominal(const femto_phi_tree::TrackRow& t) {
  const UInt_t need = femto_phi_tree::kSelNominalPid | femto_phi_tree::kSelNominalFemto;
  return (t.selFlags & need) == need;
}

Bool_t PassDeuteronVariation(const femto_phi_tree::TrackRow& t, Double_t nSigmaDeuteronMax,
                             Double_t dcaDeuteronMax) {
  if (t.speciesCode != femto_phi_tree::kSpeciesDeuteron) return kFALSE;
  if (!(t.selFlags & femto_phi_tree::kSelNominalPid)) return kFALSE;
  if (nSigmaDeuteronMax < 0 && dcaDeuteronMax < 0) return PassDeuteronNominal(t);
  const Bool_t tightenOnly =
      (nSigmaDeuteronMax < 0 || nSigmaDeuteronMax <= 2.0) && (dcaDeuteronMax < 0 || dcaDeuteronMax <= 1.0);
  if (tightenOnly) {
    if (!PassDeuteronNominal(t)) return kFALSE;
    if (nSigmaDeuteronMax > 0 && TMath::Abs(t.nSigmaDeuteron) > nSigmaDeuteronMax) return kFALSE;
    if (dcaDeuteronMax > 0 && t.dca >= dcaDeuteronMax) return kFALSE;
    return kTRUE;
  }
  if (nSigmaDeuteronMax > 0 && TMath::Abs(t.nSigmaDeuteron) > nSigmaDeuteronMax) return kFALSE;
  if (dcaDeuteronMax > 0 && t.dca >= dcaDeuteronMax) return kFALSE;
  return (t.selFlags & femto_phi_tree::kSelTrackQualityNom) ? kTRUE : kFALSE;
}
}  // namespace

void anaFemtoPhiTreeDownstream(const Char_t* treeFile, const Char_t* outFile, const Char_t* mainconf,
                               Int_t bufferSizeOverride = -1, const Char_t* mixingModeOverride = "bufferAll",
                               Double_t nSigmaDeuteronMax = -1.0, Double_t dcaDeuteronMax = -1.0,
                               Double_t signalMin = -1.0, Double_t signalMax = -1.0) {
  if (!ConfigManager::GetInstance().LoadConfig(mainconf)) {
    std::cerr << "ERROR: failed to load " << mainconf << std::endl;
    return;
  }
  PhiCutConfig& phiCfg = ConfigManager::GetInstance().GetPhiCuts();
  if (!phiCfg.FinalizeRapidityFrame(ConfigManager::GetInstance().GetCentralityCuts())) {
    std::cerr << "ERROR: FinalizeRapidityFrame failed" << std::endl;
    return;
  }
  MixingConfig& mix = ConfigManager::GetInstance().GetMixingConfig();
  const FemtoConfig& femtoCfg = ConfigManager::GetInstance().GetFemtoConfig();
  const Int_t bufferSize = (bufferSizeOverride > 0) ? bufferSizeOverride : mix.bufferSize;
  TString mode = mixingModeOverride ? mixingModeOverride : mix.mixingMode.c_str();
  Double_t sigMin = 1.012;
  Double_t sigMax = 1.026;
  const FemtoConfig::ChannelDef* ch = femtoCfg.FindChannel("phi_deuteron_signal");
  if (ch) {
    sigMin = ch->signalMin;
    sigMax = ch->signalMax;
  }
  if (signalMin > 0) sigMin = signalMin;
  if (signalMax > 0) sigMax = signalMax;

  TFile* fin = TFile::Open(treeFile, "READ");
  if (!fin || fin->IsZombie()) {
    std::cerr << "ERROR: cannot open " << treeFile << std::endl;
    return;
  }
  TTree* ev = (TTree*)fin->Get("FemtoEventTree");
  TTree* tr = (TTree*)fin->Get("FemtoTrackTree");
  if (!ev || !tr) {
    std::cerr << "ERROR: missing FemtoEventTree or FemtoTrackTree" << std::endl;
    return;
  }

  femto_phi_tree::EventRow e;
  femto_phi_tree::TrackRow t;
  ev->SetBranchAddress("schemaVersion", &e.schemaVersion);
  ev->SetBranchAddress("eventUID", &e.eventUID);
  ev->SetBranchAddress("runId", &e.runId);
  ev->SetBranchAddress("eventId", &e.eventId);
  ev->SetBranchAddress("cent9", &e.cent9);
  ev->SetBranchAddress("mixBin", &e.mixBin);
  ev->SetBranchAddress("vz", &e.vz);
  ev->SetBranchAddress("psi2", &e.psi2);

  tr->SetBranchAddress("eventUID", &t.eventUID);
  tr->SetBranchAddress("trackIndex", &t.trackIndex);
  tr->SetBranchAddress("speciesCode", &t.speciesCode);
  tr->SetBranchAddress("charge", &t.charge);
  tr->SetBranchAddress("px", &t.px);
  tr->SetBranchAddress("py", &t.py);
  tr->SetBranchAddress("pz", &t.pz);
  tr->SetBranchAddress("nSigmaDeuteron", &t.nSigmaDeuteron);
  tr->SetBranchAddress("tofMatch", &t.tofMatch);
  tr->SetBranchAddress("mass2", &t.mass2);
  tr->SetBranchAddress("deltaOneOverBeta", &t.deltaOneOverBeta);
  tr->SetBranchAddress("dca", &t.dca);
  tr->SetBranchAddress("originX", &t.originX);
  tr->SetBranchAddress("originY", &t.originY);
  tr->SetBranchAddress("originZ", &t.originZ);
  tr->SetBranchAddress("bField", &t.bField);
  tr->SetBranchAddress("selFlags", &t.selFlags);

  std::map<ULong64_t, std::vector<Long64_t> > tracksByEvent;
  const Long64_t nTr = tr->GetEntries();
  for (Long64_t i = 0; i < nTr; ++i) {
    tr->GetEntry(i);
    tracksByEvent[t.eventUID].push_back(i);
  }

  TFile* fout = new TFile(outFile, "RECREATE");
  TH1D* hMkk = new TH1D("hMkk", "M_{KK} tofStrict-like;M_{KK} (GeV/c^{2});counts", 80, 0.98, 1.06);
  TH1D* hKstarSE = new TH1D("hKstarSE_phi_deuteron_signal", "SE k*;k* (GeV/c);counts", 50, 0.0, 1.0);
  TH1D* hKstarME = new TH1D("hKstarME_phi_deuteron_signal", "ME k*;k* (GeV/c);counts", 50, 0.0, 1.0);
  TH1D* hNPhi = new TH1D("hNPhi", "nominal #phi / event", 21, -0.5, 20.5);
  TH1D* hND = new TH1D("hND", "nominal d / event", 21, -0.5, 20.5);
  TH1D* hRejectShared = new TH1D("hRejectShared", "shared-track rejects", 2, -0.5, 1.5);
  TH1D* hSameEventME = new TH1D("hSameEventME", "ME pairs from same eventUID (must be 0)", 2, -0.5, 1.5);
  TH1D* hFlow = new TH1D("hFlow", "counters", 10, 0.5, 10.5);
  hFlow->GetXaxis()->SetBinLabel(1, "events");
  hFlow->GetXaxis()->SetBinLabel(2, "nom_kp");
  hFlow->GetXaxis()->SetBinLabel(3, "nom_km");
  hFlow->GetXaxis()->SetBinLabel(4, "nom_phi");
  hFlow->GetXaxis()->SetBinLabel(5, "nom_d");
  hFlow->GetXaxis()->SetBinLabel(6, "se_pairs");
  hFlow->GetXaxis()->SetBinLabel(7, "me_pairs");
  hFlow->GetXaxis()->SetBinLabel(8, "shared_reject");
  hFlow->GetXaxis()->SetBinLabel(9, "env_d_not_nom");
  hFlow->GetXaxis()->SetBinLabel(10, "env_d_nsig_2to5");

  std::map<Int_t, std::deque<MixEvent> > pool;
  Long64_t nSE = 0, nME = 0, nShared = 0, nPhiTot = 0, nDTot = 0;
  Long64_t nNomKp = 0, nNomKm = 0, nSameEventME = 0, nEnvDNotNom = 0, nEnvDNsig2to5 = 0;
  const Long64_t nEv = ev->GetEntries();
  const Int_t maxMixed = mix.maxMixedPairsPerEvent;
  TRandom3 rng(1);

  for (Long64_t ie = 0; ie < nEv; ++ie) {
    ev->GetEntry(ie);
    hFlow->Fill(1);
    std::vector<Cand> kp, km, deuterons, phis;
    std::map<Int_t, femto_phi_tree::TrackRow> kmap;
    const std::vector<Long64_t>& idxs = tracksByEvent[e.eventUID];
    for (size_t k = 0; k < idxs.size(); ++k) {
      tr->GetEntry(idxs[k]);
      if (t.speciesCode == femto_phi_tree::kSpeciesKp || t.speciesCode == femto_phi_tree::kSpeciesKm) {
        if (PassKaonForPhi(t)) {
          kmap[t.trackIndex] = t;
          Cand c;
          c.trackIndex = t.trackIndex;
          c.mKK = 0;
          c.dau1 = c.dau2 = -1;
          const Double_t p2 = t.px * t.px + t.py * t.py + t.pz * t.pz;
          c.p4 = TLorentzVector(t.px, t.py, t.pz, TMath::Sqrt(kKaonMass * kKaonMass + p2));
          if (t.speciesCode == femto_phi_tree::kSpeciesKp) {
            kp.push_back(c);
            nNomKp++;
            hFlow->Fill(2);
          } else {
            km.push_back(c);
            nNomKm++;
            hFlow->Fill(3);
          }
        }
      } else if (t.speciesCode == femto_phi_tree::kSpeciesDeuteron) {
        if (!(t.selFlags & femto_phi_tree::kSelNominalPid) || !(t.selFlags & femto_phi_tree::kSelNominalFemto)) {
          nEnvDNotNom++;
          hFlow->Fill(9);
        }
        if (TMath::Abs(t.nSigmaDeuteron) > 2.0 && TMath::Abs(t.nSigmaDeuteron) <= 5.0) {
          nEnvDNsig2to5++;
          hFlow->Fill(10);
        }
        if (!PassDeuteronVariation(t, nSigmaDeuteronMax, dcaDeuteronMax)) continue;
        Cand c;
        c.trackIndex = t.trackIndex;
        c.mKK = 0;
        c.dau1 = c.dau2 = -1;
        const Double_t p2 = t.px * t.px + t.py * t.py + t.pz * t.pz;
        c.p4 = TLorentzVector(t.px, t.py, t.pz, TMath::Sqrt(kDeuteronMass * kDeuteronMass + p2));
        deuterons.push_back(c);
      }
    }

    for (size_t i = 0; i < kp.size(); ++i) {
      for (size_t j = 0; j < km.size(); ++j) {
        const femto_phi_tree::TrackRow& a = kmap[kp[i].trackIndex];
        const femto_phi_tree::TrackRow& b = kmap[km[j].trackIndex];
        Double_t invMass = 0.0;
        TVector3 phiMom, dcaP, dcaM;
        if (!StPhiKKReconstruction::ReconstructPhi(ToKk(a), ToKk(b), invMass, phiMom, dcaP, dcaM)) continue;
        if (!StPhiKKReconstruction::PassPairTofCut(ToKk(a), ToKk(b))) continue;
        Double_t opening = StPhiKKReconstruction::CalculateOpeningAngle(ToKk(a), ToKk(b));
        Double_t yLab = StPhiKKReconstruction::CalculatePairRapidity(invMass, phiMom);
        Double_t yPair = phiCfg.ApplyAnalysisRapidity(yLab);
        if (opening < phiCfg.minOpeningAngle || opening > phiCfg.maxOpeningAngle) continue;
        if (yPair < phiCfg.minPairRapidity || yPair > phiCfg.maxPairRapidity) continue;
        hMkk->Fill(invMass);
        Cand phi;
        phi.trackIndex = -1;
        phi.dau1 = a.trackIndex;
        phi.dau2 = b.trackIndex;
        phi.mKK = (Float_t)invMass;
        const Double_t ePhi = TMath::Sqrt(invMass * invMass + phiMom.Mag2());
        phi.p4 = TLorentzVector(phiMom.X(), phiMom.Y(), phiMom.Z(), ePhi);
        phis.push_back(phi);
        nPhiTot++;
        hFlow->Fill(4);
      }
    }
    nDTot += deuterons.size();
    hNPhi->Fill(phis.size());
    hND->Fill(deuterons.size());
    for (size_t i = 0; i < deuterons.size(); ++i) hFlow->Fill(5);

    for (size_t i = 0; i < phis.size(); ++i) {
      if (phis[i].mKK < sigMin || phis[i].mKK > sigMax) continue;
      for (size_t j = 0; j < deuterons.size(); ++j) {
        if (SharedTrack(phis[i], deuterons[j])) {
          nShared++;
          hRejectShared->Fill(1);
          hFlow->Fill(8);
          continue;
        }
        hKstarSE->Fill(KStar(phis[i].p4, deuterons[j].p4));
        nSE++;
        hFlow->Fill(6);
      }
    }

    std::deque<MixEvent>& binPool = pool[e.mixBin];
    std::vector<MixRef> refs;
    for (size_t ib = 0; ib < binPool.size(); ++ib) {
      const MixEvent& buf = binPool[ib];
      for (size_t i = 0; i < phis.size(); ++i) {
        if (phis[i].mKK < sigMin || phis[i].mKK > sigMax) continue;
        for (size_t j = 0; j < buf.deuterons.size(); ++j) {
          MixRef r;
          r.ib = ib;
          r.reverse = 0;
          r.i = i;
          r.j = j;
          refs.push_back(r);
        }
      }
      for (size_t i = 0; i < buf.phis.size(); ++i) {
        if (buf.phis[i].mKK < sigMin || buf.phis[i].mKK > sigMax) continue;
        for (size_t j = 0; j < deuterons.size(); ++j) {
          MixRef r;
          r.ib = ib;
          r.reverse = 1;
          r.i = i;
          r.j = j;
          refs.push_back(r);
        }
      }
    }

    std::vector<size_t> pick;
    const Bool_t isRandom = (mode == "randomSample");
    if (!isRandom || maxMixed <= 0 || (Int_t)refs.size() <= maxMixed) {
      for (size_t ir = 0; ir < refs.size(); ++ir) pick.push_back(ir);
    } else {
      rng.SetSeed((UInt_t)(e.eventUID & 0xffffffffULL) ^ (UInt_t)(e.eventUID >> 32) ^ 0x9e3779b9u);
      std::set<Long64_t> selected;
      const Long64_t pop = (Long64_t)refs.size();
      const Long64_t sample = maxMixed;
      for (Long64_t j = pop - sample; j < pop; ++j) {
        Long64_t candidate = (Long64_t)rng.Integer((UInt_t)(j + 1));
        if (!selected.insert(candidate).second) selected.insert(j);
      }
      for (std::set<Long64_t>::const_iterator it = selected.begin(); it != selected.end(); ++it) {
        pick.push_back((size_t)(*it));
      }
    }

    for (size_t ip = 0; ip < pick.size(); ++ip) {
      const MixRef& r = refs[pick[ip]];
      const MixEvent& buf = binPool[r.ib];
      if (buf.eventUID == e.eventUID) {
        nSameEventME++;
        hSameEventME->Fill(1);
        continue;
      }
      hSameEventME->Fill(0);
      Double_t ks = 0.0;
      if (r.reverse) {
        ks = KStar(buf.phis[r.i].p4, deuterons[r.j].p4);
      } else {
        ks = KStar(phis[r.i].p4, buf.deuterons[r.j].p4);
      }
      hKstarME->Fill(ks);
      nME++;
      hFlow->Fill(7);
    }

    MixEvent stored;
    stored.eventUID = e.eventUID;
    stored.phis = phis;
    stored.deuterons = deuterons;
    if (!stored.phis.empty() || !stored.deuterons.empty()) {
      binPool.push_back(stored);
      while ((Int_t)binPool.size() > bufferSize) binPool.pop_front();
    }
  }

  fout->cd();
  hMkk->Write();
  hKstarSE->Write();
  hKstarME->Write();
  hNPhi->Write();
  hND->Write();
  hRejectShared->Write();
  hSameEventME->Write();
  hFlow->Write();
  TNamed tf("treeFile", treeFile);
  tf.Write();
  TNamed mm("mixingMode", mode.Data());
  mm.Write();
  TNamed ns("nSE", TString::Format("%lld", (long long)nSE).Data());
  ns.Write();
  TNamed nm("nME", TString::Format("%lld", (long long)nME).Data());
  nm.Write();
  TNamed np("nPhi", TString::Format("%lld", (long long)nPhiTot).Data());
  np.Write();
  TNamed nd("nD", TString::Format("%lld", (long long)nDTot).Data());
  nd.Write();
  TNamed nSame("nSameEventME", TString::Format("%lld", (long long)nSameEventME).Data());
  nSame.Write();
  TNamed nEnv("nEnvDNotNom", TString::Format("%lld", (long long)nEnvDNotNom).Data());
  nEnv.Write();
  TNamed nEdge("nEnvDNsig2to5", TString::Format("%lld", (long long)nEnvDNsig2to5).Data());
  nEdge.Write();
  fout->Close();
  fin->Close();
  std::cout << "[downstream] events=" << nEv << " nomKp=" << nNomKp << " nomKm=" << nNomKm
            << " phi=" << nPhiTot << " d=" << nDTot << " SE=" << nSE << " ME=" << nME
            << " sharedReject=" << nShared << " sameEventME=" << nSameEventME
            << " sig=[" << sigMin << "," << sigMax << "] mode=" << mode.Data() << std::endl;
}
