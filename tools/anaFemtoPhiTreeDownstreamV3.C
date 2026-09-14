// anaFemtoPhiTreeDownstreamV3.C
// Rebuild SE/ME phi-d diagnostics from a schema 3 reduced tree.
// The schema 1 reader (anaFemtoPhiTreeDownstream.C) still handles the older trees; this file
// exists because schema 3 stores packed rows and drops fields that are recomputable:
//
//   px, py, pz        -> from pT / eta / phi
//   charge            -> sign of nHitsFit
//   tofMatch          -> kSelTofMatch bit
//   mass2             -> p^2 (1/beta^2 - 1)
//   deltaOneOverBeta  -> 1/beta - sqrt(m^2 + p^2)/p
//   mixBin            -> read from the tree for the nominal (the maker's own assignment) and
//                        recomputed with the same formula on request, so the binning can still
//                        be varied without re-reading PicoDst
//   originX/Y/Z, bField(track) -> NOT stored. The KK decay DCA therefore cannot be computed,
//                        so this reader refuses to run if maxDCAKK is set to an active value
//                        instead of silently feeding the cut a garbage number.
//
// C++98 / ROOT 5 ACLiC.
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TMath.h"
#include "TString.h"
#include "TNamed.h"
#include "TRandom3.h"
#include "ConfigManager.h"
#include "cuts/PhiCutConfig.h"
#include "cuts/MixingConfig.h"
#include "cuts/FemtoConfig.h"
#include "cuts/EventCutConfig.h"
#include "cuts/PIDCutConfig.h"
#include "cuts/NuclearIdCutConfig.h"
#include "FemtoPhiTreeSchema.h"

#include <deque>
#include <iostream>
#include <map>
#include <vector>

namespace {
const Double_t kKaonMass = 0.493677;
const Double_t kDeuteronMass = 1.875612;
// A physical KK decay DCA is a fraction of a cm; anything at or above this means the cut was
// deliberately parked in the "off" position.
const Double_t kDcaKKDisabledAbove = 50.0;

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

struct MixRef { size_t ib; Int_t reverse; size_t i; size_t j; };

Double_t KStar(const TLorentzVector& a, const TLorentzVector& b) {
  TLorentzVector q = a - b;
  TLorentzVector pair = a + b;
  q.Boost(-pair.BoostVector());
  return 0.5 * q.Vect().Mag();
}

TVector3 Momentum(const femto_phi_tree::TrackRowV2& t) {
  TVector3 p;
  p.SetPtEtaPhi(t.Pt(), t.Eta(), t.Phi());
  return p;
}

// The maker's StFemtoPhiTreeMaker::MixBins, reproduced from the same config values.
Int_t MixBinOf(Double_t vz, Int_t cent9, Double_t psi2) {
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

// Production phi-daughter PID, evaluated on a packed row with mass2 / deltaOneOverBeta
// recomputed rather than read.
//
// Use this only for a PID *variation*. For the nominal, read kSelNominalPid instead: the maker
// evaluated PassPhiDaughterTofPid on the unpacked track and stored the answer, so the flag is
// exact, while recomputing runs the packed momentum into hard cut edges. Measured on the 5-file
// pilot: two TPC-only K+ sit at p = 0.4999, and packing pT to 1/6000 GeV/c rounds them just past
// the p <= pMomKaonPID = 0.5 edge, which cost one phi candidate out of 249. The shift itself is
// harmless (1.7e-4 GeV/c); what is not harmless is letting it decide a boolean the maker already
// decided. A variation moves the edge by far more than the quantization, so recomputing is fine
// there.
Bool_t RecomputeDaughterPid(const femto_phi_tree::TrackRowV2& t) {
  const PIDCutConfig& pid = ConfigManager::GetInstance().GetPIDCuts();
  const Double_t p = t.P();
  const Bool_t tof = t.HasTof();
  const Short_t q = (Short_t)t.Charge();
  if (!tof) {
    if (q < 0 && pid.phiDaughterKaonMinusRequireTof) return kFALSE;
    return p <= pid.pMomKaonPID;
  }
  Bool_t pass = kTRUE;
  if (pid.tofUseMass2Cut) {
    const Double_t m2 = t.Mass2();
    pass = pass && (m2 >= pid.minMass2Kaon && m2 <= pid.maxMass2Kaon);
  }
  if (pid.tofUseDeltaInvBetaCut) {
    const Double_t d = t.DeltaOneOverBeta(kKaonMass);
    pass = pass && (TMath::Abs(d) <= pid.maxAbsDeltaOneOverBetaKaon);
  }
  return pass;
}

Bool_t PassDaughterPid(const femto_phi_tree::TrackRowV2& t, Bool_t recompute) {
  if (!recompute) return (t.selFlags & femto_phi_tree::kSelNominalPid) != 0;
  return RecomputeDaughterPid(t);
}

Bool_t PassKaonForPhi(const femto_phi_tree::TrackRowV2& t) {
  const UInt_t need = femto_phi_tree::kSelKaonCutsNom | femto_phi_tree::kSelLoosePid;
  return (t.selFlags & need) == need;
}

// StFemtoMaker.cxx:494 gates the nuclear-ID deuteron collection on
//   nHitsDedx >= nuclearId.minNHitsDedxNuclear
// before any PID or femtoscopy cut. That is an analysis cut, not a storage envelope: the
// tree maker's envMinNHitsDedxNuclear only decides which rows are written, and is set
// looser (10) than the analysis value (15) so the cut can be varied here. It therefore has
// to be re-applied downstream, or the tree yields deuterons the maker-direct chain rejects.
Bool_t PassNuclearDedxHits(const femto_phi_tree::TrackRowV2& t, Int_t minNHitsDedxOverride) {
  const Int_t need = (minNHitsDedxOverride > 0)
                         ? minNHitsDedxOverride
                         : (Int_t)ConfigManager::GetInstance().GetNuclearIdCuts().minNHitsDedxNuclear;
  return (Int_t)t.nHitsDedx >= need;
}

// StFemtoMaker.cxx:323 gates the whole track loop on PassTrackCuts (the `track` YAML block:
// nHitsFit, nHits ratio, nHitsDedx, pT, eta, global DCA, chi2) before any species branch runs.
// The tree maker records that same predicate as kSelTrackQualityNom instead of applying it, so
// the nominal deuteron selection has to require it here; without it the tree yields deuterons
// with chi2 or nHitsFit outside the nominal track cuts that the maker-direct chain never sees.
// Keeping it a flag rather than an envelope is what makes the track-cut systematics
// (nHitsFit 15/17/20, chi2 3/5, DCA 1/2/3, pT 0.15/0.2) reproducible from the tree.
Bool_t PassDeuteronNominal(const femto_phi_tree::TrackRowV2& t) {
  const UInt_t need = femto_phi_tree::kSelTrackQualityNom | femto_phi_tree::kSelNominalPid |
                      femto_phi_tree::kSelNominalFemto;
  return (t.selFlags & need) == need;
}

Bool_t PassDeuteronVariation(const femto_phi_tree::TrackRowV2& t, Double_t nSigmaMax,
                             Double_t dcaMax, Int_t minNHitsDedx) {
  if (t.speciesCode != femto_phi_tree::kSpeciesDeuteron) return kFALSE;
  if (!PassNuclearDedxHits(t, minNHitsDedx)) return kFALSE;
  if (!(t.selFlags & femto_phi_tree::kSelNominalPid)) return kFALSE;
  if (nSigmaMax < 0 && dcaMax < 0) return PassDeuteronNominal(t);
  const Bool_t tightenOnly = (nSigmaMax < 0 || nSigmaMax <= 2.0) && (dcaMax < 0 || dcaMax <= 1.0);
  if (tightenOnly) {
    if (!PassDeuteronNominal(t)) return kFALSE;
    if (nSigmaMax > 0 && TMath::Abs(t.NSigmaDeuteron()) > nSigmaMax) return kFALSE;
    if (dcaMax > 0 && t.Dca() >= dcaMax) return kFALSE;
    return kTRUE;
  }
  if (nSigmaMax > 0 && TMath::Abs(t.NSigmaDeuteron()) > nSigmaMax) return kFALSE;
  if (dcaMax > 0 && t.Dca() >= dcaMax) return kFALSE;
  return (t.selFlags & femto_phi_tree::kSelTrackQualityNom) ? kTRUE : kFALSE;
}

Bool_t SharedTrack(const Cand& phi, const Cand& d) {
  return d.trackIndex >= 0 && (d.trackIndex == phi.dau1 || d.trackIndex == phi.dau2);
}
}  // namespace

void anaFemtoPhiTreeDownstreamV3(const Char_t* treeFile, const Char_t* outFile,
                                 const Char_t* mainconf, Int_t bufferSizeOverride = -1,
                                 const Char_t* mixingModeOverride = "bufferAll",
                                 Double_t nSigmaDeuteronMax = -1.0, Double_t dcaDeuteronMax = -1.0,
                                 Double_t signalMin = -1.0, Double_t signalMax = -1.0,
                                 Int_t minNHitsDedxNuclear = -1,
                                 Bool_t recomputeDaughterPid = kFALSE,
                                 Bool_t recomputeMixBin = kFALSE) {
  using namespace femto_phi_tree;
  if (!ConfigManager::GetInstance().LoadConfig(mainconf)) {
    std::cerr << "ERROR: failed to load " << mainconf << std::endl;
    return;
  }
  PhiCutConfig& phiCfg = ConfigManager::GetInstance().GetPhiCuts();
  if (!phiCfg.FinalizeRapidityFrame(ConfigManager::GetInstance().GetCentralityCuts())) {
    std::cerr << "ERROR: FinalizeRapidityFrame failed" << std::endl;
    return;
  }
  // Schema 3 does not carry the track helix origin, so the KK decay DCA cannot be formed.
  if (phiCfg.maxDCAKK < kDcaKKDisabledAbove) {
    std::cerr << "ERROR: maxDCAKK = " << phiCfg.maxDCAKK << " is an active cut, but schema 3 "
              << "does not store originX/Y/Z or the track B field, so the KK DCA cannot be "
              << "computed. Either leave the cut disabled or produce the tree with the origin "
              << "companion fields (implementation-plan-20260913.md sec 10.5)." << std::endl;
    return;
  }

  MixingConfig& mix = ConfigManager::GetInstance().GetMixingConfig();
  const FemtoConfig& femtoCfg = ConfigManager::GetInstance().GetFemtoConfig();
  const Int_t bufferSize = (bufferSizeOverride > 0) ? bufferSizeOverride : mix.bufferSize;
  TString mode = mixingModeOverride ? mixingModeOverride : mix.mixingMode.c_str();
  Double_t sigMin = 1.012, sigMax = 1.026;
  const FemtoConfig::ChannelDef* ch = femtoCfg.FindChannel("phi_deuteron_signal");
  if (ch) { sigMin = ch->signalMin; sigMax = ch->signalMax; }
  if (signalMin > 0) sigMin = signalMin;
  if (signalMax > 0) sigMax = signalMax;

  TFile* fin = TFile::Open(treeFile, "READ");
  if (!fin || fin->IsZombie()) { std::cerr << "ERROR: cannot open " << treeFile << std::endl; return; }
  TTree* ev = (TTree*)fin->Get("FemtoEventTree");
  TTree* tr = (TTree*)fin->Get("FemtoTrackTree");
  if (!ev || !tr) { std::cerr << "ERROR: missing trees" << std::endl; return; }

  EventRowV3 e;
  TrackRowV2 t;
  ev->SetBranchAddress("schemaVersion", &e.schemaVersion);
  ev->SetBranchAddress("eventUID", &e.eventUID);
  ev->SetBranchAddress("cent9", &e.cent9);
  ev->SetBranchAddress("vz", &e.vz);
  ev->SetBranchAddress("psi2", &e.psi2);
  ev->SetBranchAddress("mixBin", &e.mixBin);

  tr->SetBranchAddress("eventUID", &t.eventUID);
  tr->SetBranchAddress("trackIndex", &t.trackIndex);
  tr->SetBranchAddress("speciesCode", &t.speciesCode);
  tr->SetBranchAddress("pT", &t.pT);
  tr->SetBranchAddress("eta", &t.eta);
  tr->SetBranchAddress("phi", &t.phi);
  tr->SetBranchAddress("nSigmaDeuteron", &t.nSigmaDeuteron);
  tr->SetBranchAddress("tofBeta", &t.tofBeta);
  tr->SetBranchAddress("dca", &t.dca);
  tr->SetBranchAddress("nHitsFit", &t.nHitsFit);
  tr->SetBranchAddress("nHitsMax", &t.nHitsMax);
  tr->SetBranchAddress("nHitsDedx", &t.nHitsDedx);
  tr->SetBranchAddress("selFlags", &t.selFlags);

  ev->GetEntry(0);
  if (e.schemaVersion != kSchemaVersionV3) {
    std::cerr << "ERROR: tree is schema " << e.schemaVersion << ", this reader needs "
              << kSchemaVersionV3 << ". Use anaFemtoPhiTreeDownstream.C for schema 1."
              << std::endl;
    return;
  }

  std::map<ULong64_t, std::vector<Long64_t> > tracksByEvent;
  const Long64_t nTr = tr->GetEntries();
  for (Long64_t i = 0; i < nTr; ++i) {
    tr->GetEntry(i);
    tracksByEvent[t.eventUID].push_back(i);
  }

  TFile* fout = new TFile(outFile, "RECREATE");
  TH1D* hMkk = new TH1D("hMkk", "M_{KK};M_{KK} (GeV/c^{2});counts", 80, 0.98, 1.06);
  TH1D* hKstarSE = new TH1D("hKstarSE_phi_deuteron_signal", "SE k*;k* (GeV/c);counts", 50, 0.0, 1.0);
  TH1D* hKstarME = new TH1D("hKstarME_phi_deuteron_signal", "ME k*;k* (GeV/c);counts", 50, 0.0, 1.0);
  TH1D* hNPhi = new TH1D("hNPhi", "nominal #phi / event", 21, -0.5, 20.5);
  TH1D* hND = new TH1D("hND", "nominal d / event", 21, -0.5, 20.5);
  TH1D* hRejectShared = new TH1D("hRejectShared", "shared-track rejects", 2, -0.5, 1.5);
  TH1D* hSameEventME = new TH1D("hSameEventME", "ME pairs from same eventUID (must be 0)", 2, -0.5, 1.5);

  std::map<Int_t, std::deque<MixEvent> > pool;
  Long64_t nSE = 0, nME = 0, nShared = 0, nPhiTot = 0, nDTot = 0, nSameEventME = 0;
  const Long64_t nEv = ev->GetEntries();
  const Int_t maxMixed = mix.maxMixedPairsPerEvent;
  TRandom3 rng(1);

  for (Long64_t ie = 0; ie < nEv; ++ie) {
    ev->GetEntry(ie);
    std::vector<Cand> kp, km, deuterons, phis;
    std::map<Int_t, TrackRowV2> kmap;
    const std::vector<Long64_t>& idxs = tracksByEvent[e.eventUID];
    for (size_t k = 0; k < idxs.size(); ++k) {
      tr->GetEntry(idxs[k]);
      if (t.speciesCode == kSpeciesKp || t.speciesCode == kSpeciesKm) {
        if (!PassKaonForPhi(t)) continue;
        kmap[t.trackIndex] = t;
        Cand c;
        c.trackIndex = t.trackIndex;
        c.mKK = 0;
        c.dau1 = c.dau2 = -1;
        TVector3 p = Momentum(t);
        c.p4 = TLorentzVector(p, TMath::Sqrt(kKaonMass * kKaonMass + p.Mag2()));
        if (t.speciesCode == kSpeciesKp) kp.push_back(c); else km.push_back(c);
      } else if (t.speciesCode == kSpeciesDeuteron) {
        if (!PassDeuteronVariation(t, nSigmaDeuteronMax, dcaDeuteronMax, minNHitsDedxNuclear)) continue;
        Cand c;
        c.trackIndex = t.trackIndex;
        c.mKK = 0;
        c.dau1 = c.dau2 = -1;
        TVector3 p = Momentum(t);
        c.p4 = TLorentzVector(p, TMath::Sqrt(kDeuteronMass * kDeuteronMass + p.Mag2()));
        deuterons.push_back(c);
      }
    }

    for (size_t i = 0; i < kp.size(); ++i) {
      for (size_t j = 0; j < km.size(); ++j) {
        const TrackRowV2& a = kmap[kp[i].trackIndex];
        const TrackRowV2& b = kmap[km[j].trackIndex];
        if (!PassDaughterPid(a, recomputeDaughterPid) || !PassDaughterPid(b, recomputeDaughterPid))
          continue;
        TVector3 pa = Momentum(a), pb = Momentum(b);
        Double_t ea = TMath::Sqrt(kKaonMass * kKaonMass + pa.Mag2());
        Double_t eb = TMath::Sqrt(kKaonMass * kKaonMass + pb.Mag2());
        TVector3 phiMom = pa + pb;
        Double_t etot = ea + eb;
        Double_t m2 = etot * etot - phiMom.Mag2();
        if (m2 <= 0) continue;
        Double_t invMass = TMath::Sqrt(m2);
        Double_t cosT = pa.Dot(pb) / (pa.Mag() * pb.Mag());
        if (cosT > 1) cosT = 1;
        if (cosT < -1) cosT = -1;
        Double_t opening = TMath::ACos(cosT);
        Double_t pzPhi = phiMom.Z();
        Double_t yLab = 0.5 * TMath::Log((etot + pzPhi) / (etot - pzPhi));
        Double_t yPair = phiCfg.ApplyAnalysisRapidity(yLab);
        if (opening < phiCfg.minOpeningAngle || opening > phiCfg.maxOpeningAngle) continue;
        if (yPair < phiCfg.minPairRapidity || yPair > phiCfg.maxPairRapidity) continue;
        hMkk->Fill(invMass);
        Cand phi;
        phi.trackIndex = -1;
        phi.dau1 = a.trackIndex;
        phi.dau2 = b.trackIndex;
        phi.mKK = (Float_t)invMass;
        phi.p4 = TLorentzVector(phiMom, etot);
        phis.push_back(phi);
        nPhiTot++;
      }
    }
    nDTot += deuterons.size();
    hNPhi->Fill(phis.size());
    hND->Fill(deuterons.size());

    for (size_t i = 0; i < phis.size(); ++i) {
      if (phis[i].mKK < sigMin || phis[i].mKK > sigMax) continue;
      for (size_t j = 0; j < deuterons.size(); ++j) {
        if (SharedTrack(phis[i], deuterons[j])) { nShared++; hRejectShared->Fill(1); continue; }
        hKstarSE->Fill(KStar(phis[i].p4, deuterons[j].p4));
        nSE++;
      }
    }

    const Int_t mixBin =
        recomputeMixBin ? MixBinOf(e.Vz(), (Int_t)e.cent9, e.Psi2()) : (Int_t)e.mixBin;
    std::deque<MixEvent>& binPool = pool[mixBin];
    std::vector<MixRef> refs;
    for (size_t ib = 0; ib < binPool.size(); ++ib) {
      const MixEvent& buf = binPool[ib];
      for (size_t i = 0; i < phis.size(); ++i) {
        if (phis[i].mKK < sigMin || phis[i].mKK > sigMax) continue;
        for (size_t j = 0; j < buf.deuterons.size(); ++j) {
          MixRef r; r.ib = ib; r.reverse = 0; r.i = i; r.j = j; refs.push_back(r);
        }
      }
      if (mix.mixBothDirections) {
        for (size_t i = 0; i < buf.phis.size(); ++i) {
          if (buf.phis[i].mKK < sigMin || buf.phis[i].mKK > sigMax) continue;
          for (size_t j = 0; j < deuterons.size(); ++j) {
            MixRef r; r.ib = ib; r.reverse = 1; r.i = i; r.j = j; refs.push_back(r);
          }
        }
      }
    }

    std::vector<size_t> pick;
    const Bool_t isRandom = (mode == "randomSample");
    if (!isRandom || maxMixed <= 0 || (Int_t)refs.size() <= maxMixed) {
      for (size_t ir = 0; ir < refs.size(); ++ir) pick.push_back(ir);
    } else {
      std::vector<size_t> all;
      for (size_t ir = 0; ir < refs.size(); ++ir) all.push_back(ir);
      for (Int_t n = 0; n < maxMixed && !all.empty(); ++n) {
        size_t k = (size_t)(rng.Rndm() * all.size());
        if (k >= all.size()) k = all.size() - 1;
        pick.push_back(all[k]);
        all.erase(all.begin() + k);
      }
    }
    for (size_t ip = 0; ip < pick.size(); ++ip) {
      const MixRef& r = refs[pick[ip]];
      const MixEvent& buf = binPool[r.ib];
      if (buf.eventUID == e.eventUID) { nSameEventME++; hSameEventME->Fill(1); continue; }
      const Cand& A = r.reverse ? buf.phis[r.i] : phis[r.i];
      const Cand& B = r.reverse ? deuterons[r.j] : buf.deuterons[r.j];
      hKstarME->Fill(KStar(A.p4, B.p4));
      nME++;
    }

    if (!phis.empty() || !deuterons.empty()) {
      MixEvent me;
      me.eventUID = e.eventUID;
      me.phis = phis;
      me.deuterons = deuterons;
      binPool.push_back(me);
      if ((Int_t)binPool.size() > bufferSize) binPool.pop_front();
    }
  }

  std::cout << "[downstreamV3] events=" << nEv << " phi=" << nPhiTot << " d=" << nDTot
            << " SE=" << nSE << " ME=" << nME << " shared=" << nShared
            << " sameEventME=" << nSameEventME << std::endl;

  fout->cd();
  TNamed("treeFile", treeFile).Write();
  TNamed("mixingMode", mode.Data()).Write();
  TNamed("schemaVersion", TString::Format("%u", e.schemaVersion).Data()).Write();
  TNamed("nSE", TString::Format("%lld", nSE).Data()).Write();
  TNamed("nME", TString::Format("%lld", nME).Data()).Write();
  TNamed("nPhi", TString::Format("%lld", nPhiTot).Data()).Write();
  TNamed("nD", TString::Format("%lld", nDTot).Data()).Write();
  TNamed("daughterPid", recomputeDaughterPid ? "recomputed" : "kSelNominalPid").Write();
  TNamed("mixBinSource", recomputeMixBin ? "recomputed" : "stored").Write();
  TNamed("minNHitsDedxNuclear",
         TString::Format("%d", (minNHitsDedxNuclear > 0)
                                   ? minNHitsDedxNuclear
                                   : (Int_t)ConfigManager::GetInstance()
                                         .GetNuclearIdCuts()
                                         .minNHitsDedxNuclear)
             .Data())
      .Write();
  TNamed("nSameEventME", TString::Format("%lld", nSameEventME).Data()).Write();
  hMkk->Write(); hKstarSE->Write(); hKstarME->Write(); hNPhi->Write(); hND->Write();
  hRejectShared->Write(); hSameEventME->Write();
  fout->Close();
  fin->Close();
}
