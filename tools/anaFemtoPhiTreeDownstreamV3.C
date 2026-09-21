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
//   originX/Y/Z        -> carried for KAONS ONLY, in the FemtoKaonOriginTree companion, relative
//                        to the primary vertex. That is all the KK decay DCA needs beyond the
//                        track row, so maxDCAKK stays available. A tree written without the
//                        companion still makes this reader refuse an active maxDCAKK rather than
//                        feed the cut a garbage number.
//
// C++98 / ROOT 5 ACLiC.
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TH3F.h"
#include "TProfile.h"
#include "TParameter.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TMath.h"
#include "TString.h"
#include "TNamed.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TRandom3.h"
#include "ConfigManager.h"
#include "cuts/PhiCutConfig.h"
#include "cuts/MixingConfig.h"
#include "cuts/FemtoConfig.h"
#include "cuts/EventCutConfig.h"
#include "cuts/PIDCutConfig.h"
#include "cuts/NuclearIdCutConfig.h"
#include "cuts/TrackCutConfig.h"
#include "FemtoPhiTreeSchema.h"
#include "FemtoFlagConfigSnapshot.h"
#include "FemtoMixingSampler.h"
#include "FemtoPhiMixSampler.h"
#include "Data006Config.h"
#include "StPhiKKReconstruction.h"
#include "StNuclearIdHelper.h"

#include <deque>
#include <iostream>
#include <limits>
#include <map>
#include <vector>

namespace {
const Double_t kKaonMass = 0.493677;
const Double_t kProtonMass = 0.938272;

// Bachelor species paired against the phi. Index order is fixed and used for the histogram
// and counter arrays below.
enum { kPartDeuteron = 0, kPartProton = 1, kPartTriton = 2, kPartHe3 = 3, kPartHe4 = 4,
       kNPart = 5 };
const Char_t* const kPartName[kNPart] = {"deuteron", "proton", "triton", "he3", "he4"};

// The nuclear-ID species behind each bachelor index, or -1 for the proton, which does not go
// through StNuclearIdHelper at all. The four-vector and the rapidity are built with the maker's
// own helper rather than a local mass table: NuclearP4 also applies the Z = 2 rigidity scaling
// for 3He and 4He, which a mass-only treatment would silently drop.
Int_t NucSpeciesOf(Int_t part) {
  switch (part) {
    case kPartDeuteron: return (Int_t)kNucDeuteron;
    case kPartTriton:   return (Int_t)kNucTriton;
    case kPartHe3:      return (Int_t)kNucHe3;
    case kPartHe4:      return (Int_t)kNucHe4;
    default:            return -1;
  }
}

Int_t PartOfSpeciesCode(UChar_t code) {
  switch (code) {
    case femto_phi_tree::kSpeciesDeuteron: return kPartDeuteron;
    case femto_phi_tree::kSpeciesProton:   return kPartProton;
    case femto_phi_tree::kSpeciesTriton:   return kPartTriton;
    case femto_phi_tree::kSpeciesHe3:      return kPartHe3;
    case femto_phi_tree::kSpeciesHe4:      return kPartHe4;
    default:                               return -1;
  }
}
// Species that can sit on the A side of a channel. The legacy ROT/MIX stores stay present while
// the validation-only charge legs retain which kaon was rotated/current.
enum { kAPhi = 0, kAPhiRot = 1, kAPhiMix = 2, kAPhiRotKp = 3, kAPhiRotKm = 4,
       kAPhiMixCurKp = 5, kAPhiMixCurKm = 6, kAKPlus = 7, kAKMinus = 8, kNA = 9 };
const Char_t* const kAName[kNA] = {"phi", "phi_rot", "phi_mix", "phi_rot_kp", "phi_rot_km",
                                   "phi_mix_curkp", "phi_mix_curkm", "phikaon_plus",
                                   "phikaon_minus"};

Int_t AOfKey(const std::string& key) {
  for (Int_t i = 0; i < kNA; ++i)
    if (key == kAName[i]) return i;
  return -1;
}

Int_t BOfKey(const std::string& key) {
  for (Int_t i = 0; i < kNPart; ++i)
    if (key == kPartName[i]) return i;
  return -1;
}

Int_t PartOfKey(const std::string& key) { return BOfKey(key); }

// A physical KK decay DCA is a fraction of a cm; anything at or above this means the cut was
// deliberately parked in the "off" position.
const Double_t kDcaKKDisabledAbove = 50.0;

struct Cand {
  Int_t trackIndex;
  TLorentzVector p4;
  Float_t mKK;
  // Constituent tracks, as (event, index) pairs: one for a single track, two for a KK pair. The
  // event has to travel with the index because a fully-mixed phi takes its daughters from two
  // different events, and there an index alone is not an identity (FemtoCandidatesShareTrack).
  Int_t nCon;
  ULong64_t conEv[2];
  Int_t conIdx[2];
  // daughter kinematics, kept so the close-pair cut can work on the tracks that actually
  // interfere in the TPC rather than on the reconstructed phi
  Double_t dEta[2], dPhi[2], dPt[2];
  Short_t dCharge[2];
  // single-track kinematics, for a bachelor
  Double_t tEta, tPhi, tPt;
  Short_t tCharge;

  void Reset() {
    trackIndex = -1;
    mKK = 0;
    nCon = 0;
    for (Int_t i = 0; i < 2; ++i) {
      conEv[i] = 0; conIdx[i] = -1;
      dEta[i] = 0; dPhi[i] = 0; dPt[i] = -1; dCharge[i] = 0;
    }
    tEta = 0; tPhi = 0; tPt = -1; tCharge = 0;
  }
  void AddCon(ULong64_t ev, Int_t idx) {
    if (nCon < 2) { conEv[nCon] = ev; conIdx[nCon] = idx; ++nCon; }
  }
};

// One channel of the maker's channel list, with the histograms it owns. Driving the downstream
// from femtoCfg.channels rather than from a second hard-coded list is what keeps the histogram
// names, the mass windows and the A x B combinations identical to the maker's by construction.
struct ChanHists {
  TString name;
  Int_t ia, ib, aPart;
  Bool_t resonance;
  Bool_t trackTrack;
  Bool_t identical;
  Bool_t doMixing;
  // Extended pair QA: the phi-p / phi-d family carries the same per-pair information as the
  // track-track DATA-006 channels, so the two analyses can be compared observable by observable.
  Bool_t extendedQa;
  // Background templates (phi_rot, phi_mix) run a wider channel mass window than the signal
  // channel. Their QA is restricted to the signal window so a template subtracts a like for like.
  Bool_t qaSignalGate;
  Double_t mLo, mHi;
  Double_t normQMin, normQMax;
  Double_t closePairDEta, closePairDPhiStar;
  TH1D* se;
  TH1D* me;
  TH2D* seCent;
  TH2D* meCent;
  TH3F* mkkSE;
  TH3F* mkkME;
  TH3F* mkkSEWide;
  TH3F* mkkMEWide;
  TH3F* pairMtKstarSECent;
  TH3F* pairMtKstarMECent;
  TH2D* pairYSECent;
  TH2D* pairYMECent;
  TH2D* pairYNormSECent;
  TH2D* pairYNormMECent;
  TH3F* deltaEtaPhiStarSECent;
  TH3F* deltaEtaPhiStarMECent;
  TH2D* constituentSECent[8];
  TH2D* constituentMECent[8];
  TH1D* cutFlow;
  TH1D* closeSE;
  TH1D* closeME;
  Long64_t nSE, nME, nCloseSE, nCloseME, nShared;
};

// The mass window the extended QA is filled in. A signal, sideband or wide channel has already
// cut on its own window by the time this is called. The rot / mix background templates run the
// wider nominal window, so their QA is restricted to the signal window: a template is only
// subtractable from the signal distribution if both were taken over the same M(KK) region.
Bool_t QaMassPass(const ChanHists& c, const Cand& a, const FemtoConfig::ChannelDef* sig) {
  if (!c.qaSignalGate) return kTRUE;
  if (!sig) return kFALSE;
  return a.mKK >= sig->signalMin && a.mKK <= sig->signalMax;
}

struct ChannelInput {
  std::string name;
  std::string partA;
  std::string partB;
  Bool_t enabled;
  Bool_t doMixing;
  Double_t signalMin;
  Double_t signalMax;
  Double_t normQMin;
  Double_t normQMax;
  Double_t closePairDEta;
  Double_t closePairDPhiStar;
};

struct MixEvent {
  ULong64_t eventUID;
  Double_t bField;
  Int_t cent9;
  std::vector<Cand> a[kNA];
  std::vector<Cand> part[kNPart];
};

struct MixRef { size_t ib; Int_t reverse; size_t i; size_t j; };

// One contiguous block in the flat d(E0) x K+(E1) x K-(E2) population. E0 is the
// current event; E1 and E2 are distinct events already stored in the same mixing bin.
// Keeping only O(bufferSize^2) blocks lets the diagnostic sample the full triplet
// population without materialising every combination.
struct ThreeEventBlock {
  size_t kpEvent;
  size_t kmEvent;
  femto_phi_mix::PairCount begin;
  femto_phi_mix::PairCount end;
  femto_phi_mix::PairCount nD;
  femto_phi_mix::PairCount nKp;
  femto_phi_mix::PairCount nKm;
};

Bool_t CheckedTripletCount(femto_phi_mix::PairCount a, femto_phi_mix::PairCount b,
                           femto_phi_mix::PairCount c,
                           femto_phi_mix::PairCount& result) {
  const femto_phi_mix::PairCount maxv =
      std::numeric_limits<femto_phi_mix::PairCount>::max();
  if ((a != 0 && b > maxv / a) || (a * b != 0 && c > maxv / (a * b))) return kFALSE;
  result = a * b * c;
  return kTRUE;
}

Bool_t ResolveThreeEventIndex(const std::vector<ThreeEventBlock>& blocks,
                              femto_phi_mix::PairCount flat, size_t& blockIndex,
                              size_t& iD, size_t& iKp, size_t& iKm) {
  size_t lo = 0, hi = blocks.size();
  while (lo < hi) {
    const size_t mid = lo + (hi - lo) / 2;
    if (flat < blocks[mid].end)
      hi = mid;
    else
      lo = mid + 1;
  }
  if (lo >= blocks.size() || flat < blocks[lo].begin) return kFALSE;
  const ThreeEventBlock& b = blocks[lo];
  femto_phi_mix::PairCount local = flat - b.begin;
  iD = (size_t)(local % b.nD);
  local /= b.nD;
  iKm = (size_t)(local % b.nKm);
  local /= b.nKm;
  iKp = (size_t)local;
  if (iKp >= b.nKp) return kFALSE;
  blockIndex = lo;
  return kTRUE;
}

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

// ---------------------------------------------------------------------------
// Systematic variations (Step 5).
//
// Spec string, "key=value,key=value"; empty means nominal everywhere. The rule throughout is
// that a nominal decision is READ from the flag the maker stored, and only a varied one is
// recomputed from the packed row -- a variation moves a cut by far more than the quantization,
// while the nominal must reproduce the maker exactly.
//
//   nHitsFit  trkDca  trkPt   track cuts; chi2 is held at nominal through kSelTrackChi2Nom
//   chi2                      track chi2; only the storage-envelope value is representable,
//                             because chi2 itself is not stored (see the warning below)
//   kNSigma                   phi-daughter kaon |nSigma_K|
//   kM2Lo  kM2Hi              phi-daughter kaon m2 window
//   kmRequireTof              0/1, phiDaughterKaonMinusRequireTof
//   dNSigma  dDca  dNHitsDedx deuteron
//   pNSigma  pDca             proton
//   dTofPThr  pTofPThr        bachelor TofMomentumThreshold [GeV/c]: TOF is required at and above
//                             this momentum, so RAISING it loosens. 99 reproduces a dE/dx-only
//                             selection, which is what maker_auau3p85fxt_anaFemtoPhi_dedxOnly.yaml
//                             uses (deuteron 99.0, proton 2.0) against production's 0.0.
//                             Proton is the harder of the two: the maker folded the TOF rule into
//                             kSelNominalPid itself (PassTofProtonPid && |nSigma| <= cut), so a
//                             TOF-unmatched proton carries no nominal PID bit and this variation
//                             has to rebuild the PID as well as the femto cut. Both inputs are
//                             stored, so that is exact rather than approximate. The deuteron
//                             needs no such care: its kSelNominalPid is StNuclearIdHelper's
//                             dE/dx identification, which never looked at the TOF.
// ---------------------------------------------------------------------------
struct Var {
  Double_t nHitsFit, trkDca, trkPt, chi2;
  Double_t kNSigma, kM2Lo, kM2Hi;
  Int_t kmRequireTof;
  Double_t dNSigma, dDca;
  Int_t dNHitsDedx;
  Double_t pNSigma, pDca;
  Double_t dTofPThr, pTofPThr;
  Bool_t anyTrack, anyKaon, anyPid;
  TString spec;
  Var()
      : nHitsFit(-1), trkDca(-1), trkPt(-1), chi2(-1), kNSigma(-1), kM2Lo(-1), kM2Hi(-1),
        kmRequireTof(-1), dNSigma(-1), dDca(-1), dNHitsDedx(-1), pNSigma(-1), pDca(-1),
        dTofPThr(-1), pTofPThr(-1),
        anyTrack(kFALSE), anyKaon(kFALSE), anyPid(kFALSE), spec("nominal") {}
};

Bool_t ParseVariation(const Char_t* text, Var& v) {
  TString in(text ? text : "");
  in.ReplaceAll(" ", "");
  if (in.Length() == 0) return kTRUE;
  v.spec = in;
  TObjArray* items = in.Tokenize(",");
  Bool_t ok = kTRUE;
  for (Int_t i = 0; i < items->GetEntries(); ++i) {
    TString kv = ((TObjString*)items->At(i))->GetString();
    Int_t eq = kv.Index("=");
    if (eq <= 0) { std::cerr << "ERROR: bad variation item '" << kv << "'" << std::endl; ok = kFALSE; continue; }
    TString k = kv(0, eq);
    Double_t val = TString(kv(eq + 1, kv.Length())).Atof();
    if (k == "nHitsFit") { v.nHitsFit = val; v.anyTrack = kTRUE; }
    else if (k == "trkDca") { v.trkDca = val; v.anyTrack = kTRUE; }
    else if (k == "trkPt") { v.trkPt = val; v.anyTrack = kTRUE; }
    else if (k == "chi2") { v.chi2 = val; v.anyTrack = kTRUE; }
    else if (k == "kNSigma") { v.kNSigma = val; v.anyKaon = kTRUE; }
    else if (k == "kM2Lo") { v.kM2Lo = val; v.anyPid = kTRUE; }
    else if (k == "kM2Hi") { v.kM2Hi = val; v.anyPid = kTRUE; }
    else if (k == "kmRequireTof") { v.kmRequireTof = (Int_t)val; v.anyPid = kTRUE; }
    else if (k == "dNSigma") v.dNSigma = val;
    else if (k == "dDca") v.dDca = val;
    else if (k == "dNHitsDedx") v.dNHitsDedx = (Int_t)val;
    else if (k == "pNSigma") v.pNSigma = val;
    else if (k == "pDca") v.pDca = val;
    else if (k == "dTofPThr") v.dTofPThr = val;
    else if (k == "pTofPThr") v.pTofPThr = val;
    else { std::cerr << "ERROR: unknown variation key '" << k << "'" << std::endl; ok = kFALSE; }
  }
  delete items;
  return ok;
}

// The nominal track-quality decision is the kSelTrackQualityNom bit. Any track-cut variation has
// to be rebuilt from stored fields instead -- and that is only possible because kSelTrackChi2Nom
// carries the one input that is not stored. Without that bit a row failing kSelTrackQualityNom
// could have failed on nHitsFit, on chi2, or on both, and the two are indistinguishable.
Bool_t PassTrackQuality(const femto_phi_tree::TrackRowV2& t, const Var& v) {
  if (!v.anyTrack) return (t.selFlags & femto_phi_tree::kSelTrackQualityNom) != 0;
  const TrackCutConfig& tc = ConfigManager::GetInstance().GetTrackCuts();
  if (v.chi2 < 0 || v.chi2 <= tc.maxChi2) {
    if (!(t.selFlags & femto_phi_tree::kSelTrackChi2Nom)) return kFALSE;
  }
  const Int_t nhf = (v.nHitsFit > 0) ? (Int_t)v.nHitsFit : (Int_t)tc.minNHitsFit;
  const Double_t dcaMax = (v.trkDca > 0) ? v.trkDca : tc.maxDCA;
  const Double_t ptMin = (v.trkPt > 0) ? v.trkPt : tc.minPt;
  if (t.NHitsFit() < nhf) return kFALSE;
  if (t.nHitsMax <= 0) return kFALSE;
  if ((Double_t)t.NHitsFit() / (Double_t)t.nHitsMax < tc.minNHitsRatio) return kFALSE;
  if ((Int_t)t.nHitsDedx < tc.minNHitsDedx) return kFALSE;
  if (t.Pt() < ptMin || t.Pt() > tc.maxPt) return kFALSE;
  if (t.Eta() < tc.minEta || t.Eta() > tc.maxEta) return kFALSE;
  if (t.Dca() > dcaMax) return kFALSE;
  return kTRUE;
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
Bool_t RecomputeDaughterPid(const femto_phi_tree::TrackRowV2& t, const Var& v) {
  const PIDCutConfig& pid = ConfigManager::GetInstance().GetPIDCuts();
  const Double_t p = t.P();
  const Bool_t tof = t.HasTof();
  const Short_t q = (Short_t)t.Charge();
  const Bool_t kmTof =
      (v.kmRequireTof >= 0) ? (v.kmRequireTof != 0) : pid.phiDaughterKaonMinusRequireTof;
  const Double_t m2Lo = (v.kM2Lo > 0) ? v.kM2Lo : pid.minMass2Kaon;
  const Double_t m2Hi = (v.kM2Hi > 0) ? v.kM2Hi : pid.maxMass2Kaon;
  if (!tof) {
    if (q < 0 && kmTof) return kFALSE;
    return p <= pid.pMomKaonPID;
  }
  Bool_t pass = kTRUE;
  if (pid.tofUseMass2Cut) {
    const Double_t m2 = t.Mass2();
    pass = pass && (m2 >= m2Lo && m2 <= m2Hi);
  }
  if (pid.tofUseDeltaInvBetaCut) {
    const Double_t d = t.DeltaOneOverBeta(kKaonMass);
    pass = pass && (TMath::Abs(d) <= pid.maxAbsDeltaOneOverBetaKaon);
  }
  return pass;
}

Bool_t PassDaughterPid(const femto_phi_tree::TrackRowV2& t, const Var& v, Bool_t recompute) {
  if (!recompute && !v.anyPid) return (t.selFlags & femto_phi_tree::kSelNominalPid) != 0;
  return RecomputeDaughterPid(t, v);
}

// StPhiKKReconstruction::PassTofKaonPid, the loose-collection filter the maker stored as
// kSelLoosePid. It reads the SAME pid.minMass2Kaon / maxMass2Kaon the production daughter PID
// uses, so an m2 systematic moves both and the flag cannot be reused: the loose filter has to be
// rebuilt here too, or the widened window admits nothing. How far it can be widened is bounded by
// envKaonMass2Lo / envKaonMass2Hi, the storage gate.
Bool_t RecomputeLoosePid(const femto_phi_tree::TrackRowV2& t, const Var& v) {
  const PIDCutConfig& pid = ConfigManager::GetInstance().GetPIDCuts();
  if (!pid.requireTOF) return kTRUE;
  TString fb(pid.tofFallbackMode.c_str());
  fb.ToLower();
  if (fb.IsNull()) fb = "acceptlowpt";
  if (fb == "acceptlowpt" && t.Pt() <= pid.pTofFallbackMax) return kTRUE;
  if (!t.HasTof()) return (fb == "tpconly");
  Bool_t pass = kTRUE;
  if (pid.tofUseMass2Cut) {
    const Double_t m2Lo = (v.kM2Lo > 0) ? v.kM2Lo : pid.minMass2Kaon;
    const Double_t m2Hi = (v.kM2Hi > 0) ? v.kM2Hi : pid.maxMass2Kaon;
    pass = pass && (t.Mass2() >= m2Lo && t.Mass2() <= m2Hi);
  }
  if (pid.tofUseDeltaInvBetaCut) {
    pass = pass && (TMath::Abs(t.DeltaOneOverBeta(kKaonMass)) <= pid.maxAbsDeltaOneOverBetaKaon);
  }
  return pass;
}

// kSelKaonCutsNom is PassNominalKaonCuts = nominal track cuts + phi.maxDCAKaon + phi.nSigmaKaon.
// A track-cut or kaon-nSigma variation has to rebuild it; everything it needs is stored.
Bool_t PassKaonForPhi(const femto_phi_tree::TrackRowV2& t, const Var& v) {
  if (v.anyPid) {
    if (!RecomputeLoosePid(t, v)) return kFALSE;
  } else if (!(t.selFlags & femto_phi_tree::kSelLoosePid)) {
    return kFALSE;
  }
  if (!v.anyTrack && !v.anyKaon) return (t.selFlags & femto_phi_tree::kSelKaonCutsNom) != 0;
  const PhiCutConfig& phi = ConfigManager::GetInstance().GetPhiCuts();
  if (!PassTrackQuality(t, v)) return kFALSE;
  if (t.Dca() > phi.maxDCAKaon) return kFALSE;
  const Double_t ns = (v.kNSigma > 0) ? v.kNSigma : phi.nSigmaKaon;
  if (TMath::Abs(t.NSigmaKaon()) > ns) return kFALSE;
  return kTRUE;
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
Bool_t PassDeuteronNominal(const femto_phi_tree::TrackRowV2& t, const Var& v) {
  const UInt_t need = femto_phi_tree::kSelNominalPid | femto_phi_tree::kSelNominalFemto;
  if ((t.selFlags & need) != need) return kFALSE;
  return PassTrackQuality(t, v);
}

// Rebuild the FULL species femto cut from stored fields, with only the varied parameters moved.
//
// A loosening variation cannot ride on kSelNominalFemto -- the flag has already applied the
// nominal value -- so the whole predicate has to be reconstructed. An earlier version rebuilt only
// nSigma and DCA and returned true, silently dropping the TOF / mass2 rule, the pT windows, eta,
// the hit counts and the rapidity window: dNSigma=3.0 then returned 685,895 deuterons against a
// nominal 240,995, a factor 2.8 that looked like a working variation. Every field this needs is
// stored, so there is no reason to approximate it.
Bool_t PassFemtoSpeciesRebuild(const femto_phi_tree::TrackRowV2& t, Bool_t isProton,
                               Double_t nsMax, Double_t dcaMax, Double_t tofPThrVar = -1.0) {
  const FemtoConfig& fc = ConfigManager::GetInstance().GetFemtoConfig();
  PhiCutConfig& phiCfg = ConfigManager::GetInstance().GetPhiCuts();
  // The maker takes the deuteron mass from StNuclearIdHelper (1.87561), not from a local
  // constant, and the rapidity window has to be evaluated against the same number.
  const Double_t mass = isProton ? kProtonMass : StNuclearIdHelper::SpeciesMass(kNucDeuteron);
  const Double_t maxDca = (dcaMax > 0) ? dcaMax
                                       : (isProton ? fc.protonMaxDca : fc.deuteronMaxDca);
  const Double_t maxNSigma =
      (nsMax > 0) ? nsMax : (isProton ? fc.protonMaxAbsNSigma : fc.deuteronMaxAbsNSigma);
  const Double_t minPtPre = isProton ? fc.protonMinPtPre : fc.deuteronMinPtPre;
  const Double_t maxPtPre = isProton ? 1e9 : fc.deuteronMaxPtPre;
  const Double_t maxAbsEta = isProton ? fc.protonMaxAbsEta : fc.deuteronMaxAbsEta;
  const Short_t minNHitsFit = isProton ? fc.protonMinNHitsFit : fc.deuteronMinNHitsFit;
  const Double_t minNHitsRatio = isProton ? fc.protonMinNHitsRatio : fc.deuteronMinNHitsRatio;
  // A threshold of 0.0 is meaningful (TOF required everywhere), so the "not varied" sentinel is
  // negative rather than zero.
  const Double_t tofPThr =
      (tofPThrVar >= 0.0) ? tofPThrVar
                          : (isProton ? fc.protonTofMomentumThreshold
                                      : fc.deuteronTofMomentumThreshold);
  const Double_t minMass2 = isProton ? fc.protonMinMass2 : fc.deuteronMinMass2;
  const Double_t maxMass2 = isProton ? fc.protonMaxMass2 : fc.deuteronMaxMass2;
  const Double_t minPtPair = isProton ? fc.protonMinPtPair : fc.deuteronMinPtPair;
  const Double_t maxPtPair = isProton ? fc.protonMaxPtPair : fc.deuteronMaxPtPair;
  const Double_t minYCm = isProton ? fc.protonMinRapidityCm : fc.deuteronMinRapidityCm;
  const Double_t maxYCm = isProton ? fc.protonMaxRapidityCm : fc.deuteronMaxRapidityCm;

  if (isProton) {
    if (fc.protonChargeMode == "positive" && t.Charge() <= 0) return kFALSE;
    if (fc.protonChargeMode == "negative" && t.Charge() >= 0) return kFALSE;
  } else {
    if (t.Charge() <= 0) return kFALSE;
  }
  if (t.Dca() >= maxDca) return kFALSE;
  const Double_t pmom = t.P();
  if (!isProton && (pmom < fc.deuteronMinPMom || pmom > fc.deuteronMaxPMom)) return kFALSE;
  if (t.Pt() < minPtPre || t.Pt() > maxPtPre) return kFALSE;
  if (TMath::Abs(t.Eta()) >= maxAbsEta) return kFALSE;
  const Double_t ns = isProton ? t.NSigmaProton() : t.NSigmaNuclear();
  if (TMath::Abs(ns) >= maxNSigma) return kFALSE;
  if (t.NHitsFit() < minNHitsFit) return kFALSE;
  if (t.nHitsMax <= 0) return kFALSE;
  if ((Double_t)t.NHitsFit() / (Double_t)t.nHitsMax < minNHitsRatio) return kFALSE;
  const Bool_t passTofRule =
      (pmom < tofPThr) ||
      (pmom > tofPThr && t.HasTof() && t.Mass2() >= minMass2 && t.Mass2() <= maxMass2);
  if (!passTofRule) return kFALSE;
  if (t.Pt() < minPtPair || t.Pt() > maxPtPair) return kFALSE;
  TVector3 p = Momentum(t);
  TLorentzVector lv(p, TMath::Sqrt(mass * mass + p.Mag2()));
  const Double_t yCm = phiCfg.ApplyAnalysisRapidity(lv.Rapidity());
  if (yCm < minYCm || yCm > maxYCm) return kFALSE;
  return kTRUE;
}

Bool_t PassDeuteronVariation(const femto_phi_tree::TrackRowV2& t, const Var& v) {
  if (t.speciesCode != femto_phi_tree::kSpeciesDeuteron) return kFALSE;
  if (!PassNuclearDedxHits(t, v.dNHitsDedx)) return kFALSE;
  if (!(t.selFlags & femto_phi_tree::kSelNominalPid)) return kFALSE;
  if (!PassTrackQuality(t, v)) return kFALSE;
  if (v.dNSigma < 0 && v.dDca < 0 && v.dTofPThr < 0) {
    return (t.selFlags & femto_phi_tree::kSelNominalFemto) != 0;
  }
  const FemtoConfig& fc = ConfigManager::GetInstance().GetFemtoConfig();
  // Tightening can ride on top of kSelNominalFemto; loosening cannot, because the flag has
  // already applied the nominal value, so the whole femto cut has to be rebuilt. A TOF threshold
  // ABOVE the nominal is a loosening -- it stops requiring the match over a range of momenta --
  // which is why it is compared the other way round from the nSigma and DCA windows.
  const Bool_t tightenOnly = (v.dNSigma < 0 || v.dNSigma <= fc.deuteronMaxAbsNSigma) &&
                             (v.dDca < 0 || v.dDca <= fc.deuteronMaxDca) &&
                             (v.dTofPThr < 0 || v.dTofPThr <= fc.deuteronTofMomentumThreshold);
  if (tightenOnly && !(t.selFlags & femto_phi_tree::kSelNominalFemto)) return kFALSE;
  return PassFemtoSpeciesRebuild(t, kFALSE, v.dNSigma, v.dDca, v.dTofPThr);
}

// The maker's proton PID, rebuilt from stored fields with the TOF threshold moved.
//
// StFemtoPhiTreeMaker sets kSelNominalPid for a proton as
//   PassTofProtonPid(trk) && |nSigmaProton| <= pid.nSigmaProton,
// and PassTofProtonPid is "below the threshold, accept; at or above it, require a TOF match and
// mass2 inside [protonMinMass2, protonMaxMass2]". With production's threshold of 0.0 that means
// every accepted proton has a TOF match, so the stored bit cannot express a selection that does
// not require one: on one merged block, 107,828,359 proton rows carry the PID bit with a TOF
// match and only 3,864 without. Rebuilding is therefore not an optimisation here but the only
// way to reach the 215,931,886 TOF-unmatched proton rows the envelope kept.
Bool_t RecomputeProtonPid(const femto_phi_tree::TrackRowV2& t, Double_t tofPThrVar) {
  const FemtoConfig& fc = ConfigManager::GetInstance().GetFemtoConfig();
  const PIDCutConfig& pid = ConfigManager::GetInstance().GetPIDCuts();
  const Double_t thr = (tofPThrVar >= 0.0) ? tofPThrVar : fc.protonTofMomentumThreshold;
  const Double_t pmom = t.P();
  Bool_t tofOk;
  if (pmom < thr) tofOk = kTRUE;
  else if (!t.HasTof()) tofOk = kFALSE;
  else tofOk = (t.Mass2() >= fc.protonMinMass2 && t.Mass2() <= fc.protonMaxMass2);
  if (!tofOk) return kFALSE;
  return TMath::Abs(t.NSigmaProton()) <= pid.nSigmaProton;
}

Bool_t PassProtonVariation(const femto_phi_tree::TrackRowV2& t, const Var& v) {
  if (t.speciesCode != femto_phi_tree::kSpeciesProton) return kFALSE;
  // Only a TOF variation may bypass the stored PID bit, and then it rebuilds the same predicate.
  if (v.pTofPThr < 0) {
    if (!(t.selFlags & femto_phi_tree::kSelNominalPid)) return kFALSE;
  } else {
    if (!RecomputeProtonPid(t, v.pTofPThr)) return kFALSE;
  }
  if (!PassTrackQuality(t, v)) return kFALSE;
  if (v.pNSigma < 0 && v.pDca < 0 && v.pTofPThr < 0) {
    return (t.selFlags & femto_phi_tree::kSelNominalFemto) != 0;
  }
  const FemtoConfig& fc = ConfigManager::GetInstance().GetFemtoConfig();
  const Bool_t tightenOnly = (v.pNSigma < 0 || v.pNSigma <= fc.protonMaxAbsNSigma) &&
                             (v.pDca < 0 || v.pDca <= fc.protonMaxDca) &&
                             (v.pTofPThr < 0 || v.pTofPThr <= fc.protonTofMomentumThreshold);
  if (tightenOnly && !(t.selFlags & femto_phi_tree::kSelNominalFemto)) return kFALSE;
  return PassFemtoSpeciesRebuild(t, kTRUE, v.pNSigma, v.pDca, v.pTofPThr);
}

// StFemtoMaker reaches the proton collection through PassProtonCuts, which begins with the same
// PassTrackCuts gate the deuteron branch inherits, then |nSigmaProton| <= pid.nSigmaProton and
// IsProton (= PassTofProtonPid); the charge mode and the pT / rapidity windows live inside
// PassFemtoProtonCuts. The tree maker records exactly that split across the three flags below, so
// the nominal proton is the same triple as the nominal deuteron. There is no nuclear-ID dE/dx hit
// requirement on this path.

// triton / 3He / 4He: nominal only. No planned systematic varies them and `Var` carries no knobs
// for them, so the selection is exactly the flag triple the maker stored plus the two cuts the
// tree deliberately leaves to the reader -- the analysis-level nuclear dE/dx hit count (the
// storage envelope is looser, so that the cut stays variable) and the nominal track quality.
// StFemtoPhiTreeMaker::StoreNuclearSpecies sets kSelNominalPid from IsTriton / IsHe3 / IsHe4 and
// kSelNominalFemto from PassFemtoNuclearCuts, which is field-for-field the maker's
// PassFemtoTritonCuts / PassFemtoHe3Cuts / PassFemtoHe4Cuts including the pre-cuts the maker
// applies in its own track loop.
Bool_t PassNuclearNominal(const femto_phi_tree::TrackRowV2& t, const Var& v, Int_t part) {
  if (PartOfSpeciesCode(t.speciesCode) != part) return kFALSE;
  if (!PassNuclearDedxHits(t, v.dNHitsDedx)) return kFALSE;
  const UInt_t need = femto_phi_tree::kSelNominalPid | femto_phi_tree::kSelNominalFemto;
  if ((t.selFlags & need) != need) return kFALSE;
  return PassTrackQuality(t, v);
}

Bool_t PassBachelor(const femto_phi_tree::TrackRowV2& t, const Var& v, Int_t part) {
  if (part == kPartDeuteron) return PassDeuteronVariation(t, v);
  if (part == kPartProton) return PassProtonVariation(t, v);
  return PassNuclearNominal(t, v, part);
}

// The nuclear four-vectors come from the maker's own StNuclearIdHelper::NuclearP4, which applies
// the species mass AND the Z = 2 rigidity scaling for 3He / 4He. The stored pT, eta, phi are the
// PicoDst pMom, i.e. rigidity: for a Z = 2 nucleus the momentum is twice that, and a mass-only
// treatment would put 3He and 4He at half their true momentum.
TLorentzVector BachelorP4(const femto_phi_tree::TrackRowV2& t, Int_t part) {
  TVector3 p = Momentum(t);
  const Int_t nuc = NucSpeciesOf(part);
  if (nuc < 0) return TLorentzVector(p, TMath::Sqrt(kProtonMass * kProtonMass + p.Mag2()));
  return StNuclearIdHelper::NuclearP4(p, (NuclearSpecies)nuc);
}

// ---------------------------------------------------------------------------
// TPC two-track (close-pair) cut, Step 4b.
//
// phi* is the azimuth of a helix extrapolated to transverse radius r; two tracks that stay close
// in (eta, phi*) over the TPC volume are candidates for merging or splitting. The cut acts on
// (daughter kaon, bachelor) track pairs -- the phi itself is not a track and cannot merge with
// anything. It must be applied identically in same-event and mixed-event, or it sculpts the
// correlation function.
// ---------------------------------------------------------------------------
Double_t PhiStar(Double_t phi, Double_t pt, Short_t q, Double_t bT, Double_t r) {
  Double_t arg = -0.3 * q * bT * r / (2.0 * pt);
  if (arg > 1.0) arg = 1.0;
  if (arg < -1.0) arg = -1.0;
  return phi + TMath::ASin(arg);
}

Double_t DPhiStarMin(Double_t phiA, Double_t ptA, Short_t qA, Double_t phiB, Double_t ptB,
                     Short_t qB, Double_t bT, const FemtoConfig& fc) {
  Double_t best = 1e9;
  for (Double_t r = fc.closePairRadiusMin; r <= fc.closePairRadiusMax + 1e-9;
       r += fc.closePairRadiusStep) {
    Double_t d = PhiStar(phiA, ptA, qA, bT, r) - PhiStar(phiB, ptB, qB, bT, r);
    while (d > TMath::Pi()) d -= TMath::TwoPi();
    while (d < -TMath::Pi()) d += TMath::TwoPi();
    if (TMath::Abs(d) < TMath::Abs(best)) best = d;
  }
  return best;
}

// kTRUE means the pair is rejected
Bool_t ClosePairReject(const Cand& phi, const Cand& bach, Double_t bachEta, Double_t bachPhi,
                       Double_t bachPt, Short_t bachCharge, Double_t bT, Int_t species,
                       const FemtoConfig& fc) {
  if (!fc.closePairEnabled) return kFALSE;
  // Measured windows exist for the deuteron and the proton only. The heavier nuclei are closer
  // to the deuteron in curvature than to the proton, so they take the deuteron window until they
  // have statistics of their own (closepair-step4b-20260914.md).
  const Bool_t isProtonBach = (species == kPartProton);
  const Double_t wEta =
      isProtonBach ? fc.closePairDEtaProton : fc.closePairDEtaDeuteron;
  const Double_t wPhi =
      isProtonBach ? fc.closePairDPhiStarProton : fc.closePairDPhiStarDeuteron;
  if (wEta <= 0 && wPhi <= 0) return kFALSE;
  const Bool_t ellipse = (fc.closePairShape == "ellipse");
  for (Int_t i = 0; i < 2; ++i) {
    if (phi.dPt[i] <= 0) continue;
    const Double_t dEta = phi.dEta[i] - bachEta;
    const Double_t dPhiS = DPhiStarMin(phi.dPhi[i], phi.dPt[i], phi.dCharge[i], bachPhi, bachPt,
                                       bachCharge, bT, fc);
    if (ellipse) {
      Double_t a = (wEta > 0) ? dEta / wEta : 0.0;
      Double_t b = (wPhi > 0) ? dPhiS / wPhi : 0.0;
      if (a * a + b * b < 1.0) return kTRUE;
    } else {
      const Bool_t inEta = (wEta <= 0) || (TMath::Abs(dEta) < wEta);
      const Bool_t inPhi = (wPhi <= 0) || (TMath::Abs(dPhiS) < wPhi);
      if (inEta && inPhi) return kTRUE;
    }
  }
  (void)bach;
  return kFALSE;
}

// Uncut two-track coordinates of a phi-bachelor pair, one entry per phi daughter. This is the
// resonance-side counterpart of the DATA-006 track-track QA below: it is filled whether or not
// the close-pair veto is enabled, so the veto can be studied after the fact.
const Double_t kCloseQaDEtaMax = 0.16;   // the histogram's own |delta eta| range

void FillPhiBachelorCloseQA(TH3F* h, const Cand& phi, Double_t bachEta, Double_t bachPhi,
                            Double_t bachPt, Short_t bachCharge, Double_t bT,
                            const FemtoConfig& fc, Double_t centX) {
  if (!h) return;
  for (Int_t i = 0; i < 2; ++i) {
    if (phi.dPt[i] <= 0) continue;
    const Double_t dEta = phi.dEta[i] - bachEta;
    // The radius scan is the expensive part of the whole reader, and a daughter this far away
    // in eta can neither be vetoed (the widest window is 0.04) nor land anywhere but the
    // overflow bin. The QA therefore covers |delta eta| < kCloseQaDEtaMax by construction.
    if (TMath::Abs(dEta) >= kCloseQaDEtaMax) continue;
    const Double_t dPhiS = DPhiStarMin(phi.dPhi[i], phi.dPt[i], phi.dCharge[i], bachPhi, bachPt,
                                       bachCharge, bT, fc);
    h->Fill(dEta, dPhiS, centX);
  }
}

// Track-track two-track QA/cut for DATA-006. Unlike the phi-bachelor cut above, both candidates
// are single tracks. The uncut coordinates are always histogrammed, including in the off baseline.
Bool_t Data006ClosePair(const Cand& a, const Cand& b, Double_t bT,
                        const Data006Config& cfg, const ChanHists& ch,
                        Double_t& dEta, Double_t& dPhiStar) {
  dEta = a.tEta - b.tEta;
  Double_t best = 1e9;
  for (Double_t r = cfg.closePairRadiusMin; r <= cfg.closePairRadiusMax + 1e-9;
       r += cfg.closePairRadiusStep) {
    Double_t d = PhiStar(a.tPhi, a.tPt, a.tCharge, bT, r) -
                 PhiStar(b.tPhi, b.tPt, b.tCharge, bT, r);
    while (d > TMath::Pi()) d -= TMath::TwoPi();
    while (d < -TMath::Pi()) d += TMath::TwoPi();
    if (TMath::Abs(d) < TMath::Abs(best)) best = d;
  }
  dPhiStar = best;
  if (!cfg.closePairEnabled) return kFALSE;

  const Double_t wEta = ch.closePairDEta;
  const Double_t wPhi = ch.closePairDPhiStar;
  if (wEta <= 0.0 && wPhi <= 0.0) return kFALSE;
  if (cfg.closePairShape == "ellipse") {
    const Double_t x = (wEta > 0.0) ? dEta / wEta : 0.0;
    const Double_t y = (wPhi > 0.0) ? dPhiStar / wPhi : 0.0;
    return x * x + y * y < 1.0;
  }
  return ((wEta <= 0.0) || TMath::Abs(dEta) < wEta) &&
         ((wPhi <= 0.0) || TMath::Abs(dPhiStar) < wPhi);
}

// Build the state StPhiKKReconstruction needs from a packed row plus the companion origin. The
// origin is stored relative to the primary vertex and is used as-is: translating both helices by
// the same vector leaves their distance of closest approach unchanged, so the vertex is never
// added back and its 0.01 cm quantisation cannot reach the KK DCA.
PhiKkTrackState MakeKkState(const femto_phi_tree::TrackRowV2& t,
                                                   const femto_phi_tree::KaonOriginRow& o,
                                                   Float_t bField) {
  PhiKkTrackState s;
  s.pT = (Float_t)t.Pt();
  s.eta = (Float_t)t.Eta();
  s.phi = (Float_t)t.Phi();
  s.charge = (Short_t)t.Charge();
  s.originX = (Float_t)o.Dx();
  s.originY = (Float_t)o.Dy();
  s.originZ = (Float_t)o.Dz();
  s.momentumX = (Float_t)t.Px();
  s.momentumY = (Float_t)t.Py();
  s.momentumZ = (Float_t)t.Pz();
  s.BField = bField;
  s.tofMatch = t.HasTof();
  s.mass2 = (Float_t)t.Mass2();
  s.deltaOneOverBeta = (Float_t)t.DeltaOneOverBeta(kKaonMass);
  return s;
}

// Compare this job's configuration against the one the tree was produced with.
//
// Every selFlags bit is a decision the maker already made, so it carries the maker's config with
// it. A downstream configured differently reads the flag and gets the producer's answer silently
// -- measured: a tree built with deuteronTofMomentumThreshold 99.0, read with the production
// value 0.0, returned 554,931 deuterons instead of ~240,912. Refusing is the only safe default.
//
// Returns the number of mismatches; 0 means the tree and this job agree.
Int_t CheckFlagConfig(TFile* fin) {
  TTree* t = (TTree*)fin->Get("FemtoFlagConfig");
  if (!t) {
    std::cerr << "ERROR: this tree carries no FemtoFlagConfig snapshot, so the configuration it "
              << "was produced with cannot be checked against this job's. Every selFlags bit "
              << "depends on that configuration. Re-produce with a maker that writes the "
              << "snapshot." << std::endl;
    return -1;
  }
  TString* key = 0;
  TString* val = 0;
  t->SetBranchAddress("key", &key);
  t->SetBranchAddress("value", &val);

  // After hadd the rows repeat once per subjob. Collapse them, and treat a key that carries two
  // different values as its own failure: that means trees from two configurations were merged.
  std::map<TString, TString> stored;
  std::vector<TString> inconsistent;
  for (Long64_t i = 0; i < t->GetEntries(); ++i) {
    t->GetEntry(i);
    std::map<TString, TString>::iterator it = stored.find(*key);
    if (it == stored.end()) {
      stored[*key] = *val;
    } else if (it->second != *val) {
      inconsistent.push_back(TString::Format("%s: merged trees disagree (%s vs %s)", key->Data(),
                                             it->second.Data(), val->Data()));
    }
  }

  std::vector<femto_flag_config::Entry> mine = femto_flag_config::Collect();
  std::vector<TString> diff;
  for (size_t i = 0; i < mine.size(); ++i) {
    std::map<TString, TString>::const_iterator it = stored.find(mine[i].first);
    if (it == stored.end()) {
      diff.push_back(TString::Format("%-40s tree: (absent)      this job: %s",
                                     mine[i].first.Data(), mine[i].second.Data()));
    } else if (it->second != mine[i].second) {
      diff.push_back(TString::Format("%-40s tree: %-14s this job: %s", mine[i].first.Data(),
                                     it->second.Data(), mine[i].second.Data()));
    }
  }

  if (inconsistent.empty() && diff.empty()) {
    std::cout << "[downstreamV3] flag-config snapshot: " << stored.size() << " keys, all match"
              << std::endl;
    return 0;
  }
  std::cerr << "\nERROR: this job's configuration does not match the one the tree was produced "
            << "with.\n       The stored selFlags encode the producer's cuts, so running anyway "
            << "would\n       silently give you the producer's selection, not yours.\n"
            << std::endl;
  for (size_t i = 0; i < inconsistent.size(); ++i)
    std::cerr << "  MERGE  " << inconsistent[i] << std::endl;
  for (size_t i = 0; i < diff.size(); ++i) std::cerr << "  DIFF   " << diff[i] << std::endl;
  std::cerr << "\n       Use the config the tree was produced with, or re-produce the tree."
            << std::endl;
  return (Int_t)(inconsistent.size() + diff.size());
}

// The maker's FemtoCandidatesShareTrack: two candidates overlap when they have a constituent
// track in common, compared as (event, index).
Bool_t SharedTrack(const Cand& a, const Cand& b) {
  for (Int_t i = 0; i < a.nCon; ++i)
    for (Int_t j = 0; j < b.nCon; ++j)
      if (a.conIdx[i] >= 0 && a.conIdx[i] == b.conIdx[j] && a.conEv[i] == b.conEv[j]) return kTRUE;
  return kFALSE;
}

// phi_rot / phi_mix / phikaon_* all reach the close-pair cut through the same daughter slots, so
// the window choice keys on the bachelor, not on the A species.
const std::string WideMkkSuffix(const std::string& channelName) {
  const std::string sig = "_signal";
  if (channelName.size() > sig.size() &&
      channelName.compare(channelName.size() - sig.size(), sig.size(), sig) == 0)
    return channelName.substr(0, channelName.size() - sig.size()) + "_wide";
  if (channelName.compare(0, 8, "phi_rot_") == 0 || channelName.compare(0, 8, "phi_mix_") == 0)
    return channelName + "_wide";
  return "";
}

// Construct one rotated charge leg. In postRotation mode the pair topology is evaluated after
// the rotation; preRotationLegacy deliberately keeps the historical pre-rotation decision.
Bool_t BuildRotatedCandidate(const femto_phi_tree::TrackRowV2& kp,
                             const femto_phi_tree::TrackRowV2& km, ULong64_t eventUID,
                             Double_t deltaPhi, Bool_t rotatePlus, Bool_t applyPostCut,
                             const PhiCutConfig& phiCfg, Cand& out) {
  TVector3 pPlus = Momentum(kp);
  TVector3 pMinus = Momentum(km);
  Double_t phiRot = (rotatePlus ? kp.Phi() : km.Phi()) + deltaPhi;
  while (phiRot > TMath::Pi()) phiRot -= TMath::TwoPi();
  while (phiRot < -TMath::Pi()) phiRot += TMath::TwoPi();
  if (rotatePlus)
    pPlus.SetPtEtaPhi(kp.Pt(), kp.Eta(), phiRot);
  else
    pMinus.SetPtEtaPhi(km.Pt(), km.Eta(), phiRot);

  const Double_t ePlus = TMath::Sqrt(kKaonMass * kKaonMass + pPlus.Mag2());
  const Double_t eMinus = TMath::Sqrt(kKaonMass * kKaonMass + pMinus.Mag2());
  const TVector3 phiMom = pPlus + pMinus;
  const Double_t etot = ePlus + eMinus;
  const Double_t m2 = etot * etot - phiMom.Mag2();
  if (m2 <= 0.0) return kFALSE;
  const Double_t invMass = TMath::Sqrt(m2);
  const Double_t opening = pPlus.Angle(pMinus);
  const Double_t yLab = 0.5 * TMath::Log((etot + phiMom.Z()) / (etot - phiMom.Z()));
  const Double_t yPair = phiCfg.ApplyAnalysisRapidity(yLab);
  if (applyPostCut &&
      (opening < phiCfg.minOpeningAngle || opening > phiCfg.maxOpeningAngle ||
       yPair < phiCfg.minPairRapidity || yPair > phiCfg.maxPairRapidity))
    return kFALSE;

  out.Reset();
  out.AddCon(eventUID, kp.trackIndex);
  out.AddCon(eventUID, km.trackIndex);
  out.mKK = (Float_t)invMass;
  out.p4 = TLorentzVector(phiMom, etot);
  out.dEta[0] = kp.Eta(); out.dPhi[0] = pPlus.Phi(); out.dPt[0] = kp.Pt();
  out.dCharge[0] = (Short_t)kp.Charge();
  out.dEta[1] = km.Eta(); out.dPhi[1] = pMinus.Phi(); out.dPt[1] = km.Pt();
  out.dCharge[1] = (Short_t)km.Charge();
  return kTRUE;
}

// Evaluator for the fully-mixed KK background, handed to the maker's own
// femto_phi_mix::SampleEligiblePairs so that the cap, the without-replacement draw and the flat
// index numbering are the maker's and not a second implementation. A functor rather than a
// lambda, to keep the file compilable as C++98.
struct MixKKEval {
  const std::vector<Cand>* curKp;
  const std::vector<Cand>* curKm;
  const std::vector<Cand>* bufKp;   // indexed by pool event
  const std::vector<Cand>* bufKm;
  const ULong64_t* bufUID;
  size_t nPool;
  ULong64_t curUID;
  PhiCutConfig* phiCfg;
  std::vector<Cand>* out;
  std::vector<Bool_t>* reverseTags;

  femto_phi_mix::EvalStatus operator()(femto_phi_mix::PairCount,
                                       const femto_mixing::PairReference& ref) const {
    if (ref.poolEventIndex >= nPool) return femto_phi_mix::kEvalResolveFail;
    const Cand* kpC = 0;
    const Cand* kmC = 0;
    ULong64_t evP = 0, evM = 0;
    if (!ref.reverse) {
      if (ref.firstIndex >= curKp->size() || ref.secondIndex >= bufKm[ref.poolEventIndex].size())
        return femto_phi_mix::kEvalResolveFail;
      kpC = &(*curKp)[ref.firstIndex];
      kmC = &bufKm[ref.poolEventIndex][ref.secondIndex];
      evP = curUID;
      evM = bufUID[ref.poolEventIndex];
    } else {
      if (ref.firstIndex >= bufKp[ref.poolEventIndex].size() || ref.secondIndex >= curKm->size())
        return femto_phi_mix::kEvalResolveFail;
      kpC = &bufKp[ref.poolEventIndex][ref.firstIndex];
      kmC = &(*curKm)[ref.secondIndex];
      evP = bufUID[ref.poolEventIndex];
      evM = curUID;
    }
    const TLorentzVector pKK = kpC->p4 + kmC->p4;
    const Double_t invMass = pKK.M();
    if (invMass <= 0.0) return femto_phi_mix::kEvalReject;
    const TVector3 pPlus = kpC->p4.Vect();
    const TVector3 pMinus = kmC->p4.Vect();
    const Double_t opening = pPlus.Angle(pMinus);
    const TVector3 phiMom = pKK.Vect();
    const Double_t etot = kpC->p4.E() + kmC->p4.E();
    const Double_t pz = phiMom.Z();
    const Double_t yLab = 0.5 * TMath::Log((etot + pz) / (etot - pz));
    const Double_t yPair = phiCfg->ApplyAnalysisRapidity(yLab);
    if (opening < phiCfg->minOpeningAngle || opening > phiCfg->maxOpeningAngle)
      return femto_phi_mix::kEvalReject;
    if (yPair < phiCfg->minPairRapidity || yPair > phiCfg->maxPairRapidity)
      return femto_phi_mix::kEvalReject;
    Cand c;
    c.Reset();
    c.AddCon(evP, kpC->trackIndex);
    c.AddCon(evM, kmC->trackIndex);
    c.mKK = (Float_t)invMass;
    c.p4 = TLorentzVector(phiMom, TMath::Sqrt(invMass * invMass + phiMom.Mag2()));
    c.dEta[0] = kpC->tEta; c.dPhi[0] = kpC->tPhi; c.dPt[0] = kpC->tPt; c.dCharge[0] = kpC->tCharge;
    c.dEta[1] = kmC->tEta; c.dPhi[1] = kmC->tPhi; c.dPt[1] = kmC->tPt; c.dCharge[1] = kmC->tCharge;
    out->push_back(c);
    reverseTags->push_back(ref.reverse ? kTRUE : kFALSE);
    return femto_phi_mix::kEvalAccept;
  }
};

}  // namespace

void anaFemtoPhiTreeDownstreamV3(const Char_t* treeFile, const Char_t* outFile,
                                 const Char_t* mainconf, Int_t bufferSizeOverride = -1,
                                 const Char_t* mixingModeOverride = "bufferAll",
                                 const Char_t* variation = "", Double_t signalMin = -1.0,
                                 Double_t signalMax = -1.0,
                                 Bool_t recomputeDaughterPid = kFALSE,
                                 Bool_t recomputeMixBin = kFALSE,
                                 Long64_t maxEvents = -1, Int_t pairMtBins = 240,
                                 Double_t pairMtMin = 0.8, Double_t pairMtMax = 3.2,
                                 Bool_t data005DiagnosticsOnly = kFALSE,
                                 const Char_t* data006ConfigPath = "",
                                 Long64_t firstEvent = 0,
                                 Bool_t threeEventEnabled = kFALSE,
                                 Long64_t threeEventMaxRawSamplesPerEvent = 0,
                                 Long64_t threeEventSamplingSeed = 0) {
  if (threeEventEnabled &&
      (threeEventMaxRawSamplesPerEvent <= 0 || threeEventSamplingSeed <= 0)) {
    std::cerr << "ERROR: three-event mixing requires a positive raw-sample cap and seed"
              << std::endl;
    return;
  }
  if (threeEventEnabled && (data005DiagnosticsOnly || (data006ConfigPath && data006ConfigPath[0]))) {
    std::cerr << "ERROR: three-event phi-background mixing is only defined for the standard "
              << "phi analysis" << std::endl;
    return;
  }
  Data006Config data006Cfg;
  const Bool_t data006Enabled = data006ConfigPath && data006ConfigPath[0];
  if (data006Enabled) {
    if (!data006Cfg.Load(data006ConfigPath)) {
      std::cerr << "ERROR: failed to load DATA-006 config " << data006ConfigPath << std::endl;
      return;
    }
    pairMtBins = data006Cfg.pairMtBins;
    pairMtMin = data006Cfg.pairMtMin;
    pairMtMax = data006Cfg.pairMtMax;
  }
  using namespace femto_phi_tree;
  if (pairMtBins <= 0 || !(pairMtMax > pairMtMin)) {
    std::cerr << "ERROR: invalid pair-mT axis: bins=" << pairMtBins << " range=["
              << pairMtMin << "," << pairMtMax << "]" << std::endl;
    return;
  }
  Var var;
  if (!ParseVariation(variation, var)) {
    std::cerr << "ERROR: could not parse variation '" << variation << "'" << std::endl;
    return;
  }
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
  const FemtoConfig::ChannelDef* pairMtSignal = femtoCfg.FindChannel("phi_proton_signal");
  if (!data006Enabled && !pairMtSignal) {
    std::cerr << "ERROR: missing phi_proton_signal channel required for pair-mT diagnostics"
              << std::endl;
    return;
  }
  const Int_t bufferSize = (bufferSizeOverride > 0) ? bufferSizeOverride : mix.bufferSize;
  TString mode = mixingModeOverride ? mixingModeOverride : mix.mixingMode.c_str();
  // The mass window now comes from each channel's own definition; signalMin / signalMax on the
  // command line override it for the phi channels only (see the channel table below).

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
  ev->SetBranchAddress("rawMult", &e.rawMult);
  ev->SetBranchAddress("vz", &e.vz);
  ev->SetBranchAddress("psi2", &e.psi2);
  ev->SetBranchAddress("mixBin", &e.mixBin);
  ev->SetBranchAddress("bField", &e.bField);

  tr->SetBranchAddress("eventUID", &t.eventUID);
  tr->SetBranchAddress("trackIndex", &t.trackIndex);
  tr->SetBranchAddress("speciesCode", &t.speciesCode);
  tr->SetBranchAddress("pT", &t.pT);
  tr->SetBranchAddress("eta", &t.eta);
  tr->SetBranchAddress("phi", &t.phi);
  tr->SetBranchAddress("nSigmaKaon", &t.nSigmaKaon);
  tr->SetBranchAddress("nSigmaProton", &t.nSigmaProton);
  tr->SetBranchAddress("nSigmaDeuteron", &t.nSigmaDeuteron);
  tr->SetBranchAddress("tofBeta", &t.tofBeta);
  tr->SetBranchAddress("dca", &t.dca);
  tr->SetBranchAddress("nHitsFit", &t.nHitsFit);
  tr->SetBranchAddress("nHitsMax", &t.nHitsMax);
  tr->SetBranchAddress("nHitsDedx", &t.nHitsDedx);
  tr->SetBranchAddress("selFlags", &t.selFlags);

  if (CheckFlagConfig(fin) != 0) return;

  // Companion tree with the phi-daughter kaon helix origin. Absent in trees written before
  // 2026-09-14; without it the KK DCA cannot be formed and an active maxDCAKK must be refused.
  TTree* ko = (TTree*)fin->Get("FemtoKaonOriginTree");
  const Bool_t dcaKKActive = (phiCfg.maxDCAKK < kDcaKKDisabledAbove);
  std::map<ULong64_t, KaonOriginRow> originOf;
  if (ko && dcaKKActive) {
    KaonOriginRow o;
    ko->SetBranchAddress("eventUID", &o.eventUID);
    ko->SetBranchAddress("trackIndex", &o.trackIndex);
    ko->SetBranchAddress("originDx", &o.originDx);
    ko->SetBranchAddress("originDy", &o.originDy);
    ko->SetBranchAddress("originDz", &o.originDz);
    for (Long64_t i = 0; i < ko->GetEntries(); ++i) {
      ko->GetEntry(i);
      originOf[o.eventUID * 100000ull + (ULong64_t)o.trackIndex] = o;
    }
  }
  if (dcaKKActive && !ko) {
    std::cerr << "ERROR: maxDCAKK = " << phiCfg.maxDCAKK << " is an active cut, but this tree has "
              << "no FemtoKaonOriginTree, so the KK DCA cannot be computed. Re-produce with "
              << "treeWriteKaonOrigin: true, or leave the cut disabled." << std::endl;
    return;
  }
  std::cout << "[downstreamV3] kaon origin rows=" << (ko ? ko->GetEntries() : 0)
            << "  maxDCAKK=" << phiCfg.maxDCAKK << (dcaKKActive ? " (ACTIVE)" : " (off)")
            << std::endl;

  ev->GetEntry(0);
  if (e.schemaVersion != kSchemaVersionV3) {
    std::cerr << "ERROR: tree is schema " << e.schemaVersion << ", this reader needs "
              << kSchemaVersionV3 << ". Use anaFemtoPhiTreeDownstream.C for schema 1."
              << std::endl;
    return;
  }

  // The track rows are walked with a cursor rather than indexed into a map. The maker fills all of
  // an event's tracks before moving to the next event, and hadd concatenates both trees in the
  // same input-file order, so the rows of one event are contiguous and arrive in the same order as
  // the event rows. The map this replaces cost a measured 1,003 B per event -- 5.0 GB for an
  // average run, 11.4 GB for the largest -- which would have forced the whole production to be
  // chunked around the reader's memory instead of around anything physical.
  const Long64_t nTr = tr->GetEntries();
  Long64_t trackCursor = 0;
  Long64_t nTracksSeen = 0;

  TFile* fout = new TFile(outFile, "RECREATE");
  TH1D* hMkk = new TH1D("hMkk", "M_{KK};M_{KK} (GeV/c^{2});counts", 80, 0.98, 1.06);
  // Same name and binning as the maker's own rotated-background QA, so the two M(KK) shapes can
  // be compared directly. The rotated background cannot be reproduced pair for pair -- the maker
  // draws its angles from gRandom and a reader cannot replay that stream -- so this distribution
  // and the candidate count are what the closure has to rest on.
  TH1D* hPhiRotMkk = new TH1D("hPhiRot_MKK", "rotated #phi M_{KK};M_{KK} (GeV/c^{2});counts",
                              200, 0.98, 1.18);
  TH1D* hNPart[kNPart];
  for (Int_t sp = 0; sp < kNPart; ++sp) {
    // 0-100, not 0-20: the proton multiplicity runs past 20 per event and a saturated axis
    // makes the histogram unusable as a candidate-count cross-check.
    hNPart[sp] = new TH1D(TString::Format("hN%s", kPartName[sp]).Data(),
                          TString::Format("nominal %s / event", kPartName[sp]).Data(), 101, -0.5,
                          100.5);
  }
  TH1D* hNA[kNA];
  for (Int_t ia = 0; ia < kNA; ++ia) {
    hNA[ia] = new TH1D(TString::Format("hNA_%s", kAName[ia]).Data(),
                       TString::Format("%s candidates / event", kAName[ia]).Data(), 201, -0.5,
                       200.5);
  }
  TH1D* hNPhi = new TH1D("hNPhi", "nominal #phi / event", 21, -0.5, 20.5);
  TH1D* hDcaKK = new TH1D("hDcaKK", "KK decay DCA of accepted daughter pairs;DCA_{KK} (cm);pairs",
                          200, 0.0, 20.0);
  TH1D* hRejectShared = new TH1D("hRejectShared", "shared-track rejects", 2, -0.5, 1.5);
  TH1D* hSameEventME = new TH1D("hSameEventME", "ME pairs from same eventUID (must be 0)", 2, -0.5, 1.5);
  TH1D* hAcceptedEventsVsCent9 = new TH1D(
      "hAcceptedEventsVsCent9", "accepted reduced-tree events;cent9;events", 9, -0.5, 8.5);
  TProfile* pFxtMultVsCent9 = new TProfile(
      "pFxtMultVsCent9", "raw fxtMult for accepted events;cent9;raw fxtMult", 9, -0.5, 8.5);
  TH1D* hFxtMultSumVsCent9 = new TH1D(
      "hFxtMultSumVsCent9", "sum raw fxtMult;cent9;sum raw fxtMult", 9, -0.5, 8.5);
  TH1D* hFxtMultSum2VsCent9 = new TH1D(
      "hFxtMultSum2VsCent9", "sum raw fxtMult squared;cent9;sum raw fxtMult^{2}", 9, -0.5, 8.5);

  // Channel list: normal operation follows FemtoConfig; DATA-006 supplies a downstream-only
  // overlay while the producer mainconf remains loaded for event/PID cuts and snapshot checking.
  std::vector<ChannelInput> channelInputs;
  if (data006Enabled) {
    for (size_t i = 0; i < data006Cfg.channels.size(); ++i) {
      const Data006ChannelConfig& in = data006Cfg.channels[i];
      ChannelInput c;
      c.name = in.name; c.partA = in.partA; c.partB = in.partB;
      c.enabled = in.enabled; c.doMixing = in.doMixing;
      c.signalMin = c.signalMax = 0.0;
      c.normQMin = in.normQMin; c.normQMax = in.normQMax;
      c.closePairDEta = in.closePairDEta;
      c.closePairDPhiStar = in.closePairDPhiStar;
      channelInputs.push_back(c);
    }
  } else {
    for (size_t i = 0; i < femtoCfg.channels.size(); ++i) {
      const FemtoConfig::ChannelDef& in = femtoCfg.channels[i];
      ChannelInput c;
      c.name = in.name; c.partA = in.partA; c.partB = in.partB;
      c.enabled = in.enabled; c.doMixing = in.doMixing;
      c.signalMin = in.signalMin; c.signalMax = in.signalMax;
      c.normQMin = in.normQMin; c.normQMax = in.normQMax;
      c.closePairDEta = c.closePairDPhiStar = 0.0;
      channelInputs.push_back(c);
    }
  }

  const Int_t kNKstar = data006Enabled ? data006Cfg.kstarBins : 300;
  const Double_t kKstarLo = data006Enabled ? data006Cfg.kstarMin : 0.0;
  const Double_t kKstarHi = data006Enabled ? data006Cfg.kstarMax : 3.0;
  const Int_t kNMkk = 200;
  const Double_t kMkkLo = 0.98, kMkkHi = 1.18;
  const Int_t kNCent = 9;
  const Double_t kCentLo = -0.5, kCentHi = 8.5;
  TH3D* hThreeEvent = 0;
  if (threeEventEnabled) {
    hThreeEvent = new TH3D(
        "hPhiMKK_vs_Kstar3E_deuteron_wide",
        "three-event d(E0), K+(E1), K-(E2);M_{KK} (GeV/c^{2});k* (GeV/c);cent9",
        kNMkk, kMkkLo, kMkkHi, kNKstar, kKstarLo, kKstarHi,
        kNCent, kCentLo, kCentHi);
    hThreeEvent->Sumw2();
  }
  // The pair-rapidity axis is deliberately the DATA-006 one in both runs: the phi-p / phi-d and
  // the p-p / p-d / d-d acceptance checks are only comparable if the bin edges are the same.
  const Int_t kNPairY = data006Enabled ? data006Cfg.pairRapidityBins : 160;
  const Double_t kPairYLo = data006Enabled ? data006Cfg.pairRapidityMin : -1.5;
  const Double_t kPairYHi = data006Enabled ? data006Cfg.pairRapidityMax : 0.5;

  std::vector<ChanHists> chans;
  for (size_t ic = 0; ic < channelInputs.size(); ++ic) {
    const ChannelInput& cd = channelInputs[ic];
    if (!cd.enabled) continue;
    if (data005DiagnosticsOnly && cd.name != "phi_proton_signal" &&
        cd.name != "phi_rot_proton") continue;
    ChanHists c;
    c.ia = AOfKey(cd.partA);
    c.aPart = PartOfKey(cd.partA);
    c.ib = BOfKey(cd.partB);
    c.trackTrack = data006Enabled && c.aPart >= 0 && c.ib >= 0;
    if ((!c.trackTrack && c.ia < 0) || c.ib < 0) {
      std::cout << "[downstreamV3] channel " << cd.name << " (" << cd.partA << " x " << cd.partB
                << ") skipped: species not built by this reader" << std::endl;
      continue;
    }
    c.name = cd.name.c_str();
    c.resonance = !c.trackTrack && (c.ia >= kAPhi && c.ia <= kAPhiMixCurKm);
    c.identical = c.trackTrack && c.aPart == c.ib;
    // Only the proton and deuteron bachelors get the extended QA. The triton / He3 / He4
    // channels would double the histogram memory for channels no correlation analysis reads yet.
    c.extendedQa = c.trackTrack ||
                   (c.resonance && (c.ib == kPartProton || c.ib == kPartDeuteron));
    c.qaSignalGate = c.resonance && c.ia != kAPhi;
    c.doMixing = cd.doMixing;
    c.mLo = cd.signalMin;
    c.mHi = cd.signalMax;
    c.normQMin = cd.normQMin;
    c.normQMax = cd.normQMax;
    c.closePairDEta = cd.closePairDEta;
    c.closePairDPhiStar = cd.closePairDPhiStar;
    if (c.ia == kAPhi) {
      if (signalMin > 0) c.mLo = signalMin;
      if (signalMax > 0) c.mHi = signalMax;
    }
    c.se = new TH1D(TString::Format("hKstarSE_%s", c.name.Data()).Data(),
                    TString::Format("SE k* (%s);k* (GeV/c);pairs", c.name.Data()).Data(),
                    kNKstar, kKstarLo, kKstarHi);
    c.me = new TH1D(TString::Format("hKstarME_%s", c.name.Data()).Data(),
                    TString::Format("ME k* (%s);k* (GeV/c);pairs", c.name.Data()).Data(),
                    kNKstar, kKstarLo, kKstarHi);
    c.seCent = new TH2D(TString::Format("hKstarSEVsCent_%s", c.name.Data()).Data(),
                        TString::Format("SE k* vs cent9 (%s);k* (GeV/c);cent9", c.name.Data()).Data(),
                        kNKstar, kKstarLo, kKstarHi, kNCent, kCentLo, kCentHi);
    c.meCent = new TH2D(TString::Format("hKstarMEVsCent_%s", c.name.Data()).Data(),
                        TString::Format("ME k* vs cent9 (%s);k* (GeV/c);cent9", c.name.Data()).Data(),
                        kNKstar, kKstarLo, kKstarHi, kNCent, kCentLo, kCentHi);
    c.mkkSE = c.mkkME = c.mkkSEWide = c.mkkMEWide = 0;
    c.pairMtKstarSECent = c.pairMtKstarMECent = 0;
    c.pairYSECent = c.pairYMECent = 0;
    c.pairYNormSECent = c.pairYNormMECent = 0;
    c.deltaEtaPhiStarSECent = c.deltaEtaPhiStarMECent = 0;
    c.cutFlow = 0;
    for (Int_t ik = 0; ik < 8; ++ik) {
      c.constituentSECent[ik] = 0;
      c.constituentMECent[ik] = 0;
    }
    if (c.resonance) {
      c.mkkSE = new TH3F(TString::Format("hPhiMKK_vs_KstarSE_%s", c.name.Data()).Data(),
                         TString::Format("M_{KK} vs k* SE vs cent9 (%s);M_{KK} (GeV/c^{2});k* (GeV/c);cent9",
                                         c.name.Data()).Data(),
                         kNMkk, kMkkLo, kMkkHi, kNKstar, kKstarLo, kKstarHi, kNCent, kCentLo, kCentHi);
      c.mkkME = new TH3F(TString::Format("hPhiMKK_vs_KstarME_%s", c.name.Data()).Data(),
                         TString::Format("M_{KK} vs k* ME vs cent9 (%s);M_{KK} (GeV/c^{2});k* (GeV/c);cent9",
                                         c.name.Data()).Data(),
                         kNMkk, kMkkLo, kMkkHi, kNKstar, kKstarLo, kKstarHi, kNCent, kCentLo, kCentHi);
      const std::string wide = WideMkkSuffix(cd.name);
      if (!wide.empty()) {
        c.mkkSEWide = new TH3F(TString::Format("hPhiMKK_vs_KstarSE_%s", wide.c_str()).Data(),
                               TString::Format("M_{KK} vs k* SE vs cent9 (%s);M_{KK} (GeV/c^{2});k* (GeV/c);cent9",
                                               wide.c_str()).Data(),
                               kNMkk, kMkkLo, kMkkHi, kNKstar, kKstarLo, kKstarHi, kNCent, kCentLo, kCentHi);
        c.mkkMEWide = new TH3F(TString::Format("hPhiMKK_vs_KstarME_%s", wide.c_str()).Data(),
                               TString::Format("M_{KK} vs k* ME vs cent9 (%s);M_{KK} (GeV/c^{2});k* (GeV/c);cent9",
                                               wide.c_str()).Data(),
                               kNMkk, kMkkLo, kMkkHi, kNKstar, kKstarLo, kKstarHi, kNCent, kCentLo, kCentHi);
      }
    }
    if (c.extendedQa) {
      c.pairMtKstarSECent = new TH3F(
          TString::Format("hPairMtHalf_vs_KstarSEVsCent_%s", c.name.Data()).Data(),
          TString::Format("SE pair m_{T}/2 vs k* vs cent9 (%s);m_{T,pair}/2 (GeV/c^{2});k* (GeV/c);cent9",
                          c.name.Data()).Data(),
          pairMtBins, pairMtMin, pairMtMax, kNKstar, kKstarLo, kKstarHi,
          kNCent, kCentLo, kCentHi);
    }
    if (c.extendedQa) {
      c.pairMtKstarMECent = new TH3F(
          TString::Format("hPairMtHalf_vs_KstarMEVsCent_%s", c.name.Data()).Data(),
          TString::Format("ME pair m_{T}/2 vs k* vs cent9 (%s);m_{T,pair}/2 (GeV/c^{2});k* (GeV/c);cent9",
                          c.name.Data()).Data(),
          pairMtBins, pairMtMin, pairMtMax, kNKstar, kKstarLo, kKstarHi,
          kNCent, kCentLo, kCentHi);
      c.pairYSECent = new TH2D(TString::Format("hPairRapiditySEVsCent_%s", c.name.Data()).Data(),
                               "SE pair y_{CM};y_{CM};cent9",
                               kNPairY, kPairYLo, kPairYHi, kNCent, kCentLo, kCentHi);
      c.pairYNormSECent = new TH2D(
          TString::Format("hPairRapidityNormSEVsCent_%s", c.name.Data()).Data(),
          "SE pair y_{CM} in normalization k*;y_{CM};cent9",
          kNPairY, kPairYLo, kPairYHi, kNCent, kCentLo, kCentHi);
      c.pairYNormMECent = new TH2D(
          TString::Format("hPairRapidityNormMEVsCent_%s", c.name.Data()).Data(),
          "ME pair y_{CM} in normalization k*;y_{CM};cent9",
          kNPairY, kPairYLo, kPairYHi, kNCent, kCentLo, kCentHi);
      c.pairYMECent = new TH2D(TString::Format("hPairRapidityMEVsCent_%s", c.name.Data()).Data(),
                               "ME pair y_{CM};y_{CM};cent9",
                               kNPairY, kPairYLo, kPairYHi, kNCent, kCentLo, kCentHi);
      c.deltaEtaPhiStarSECent = new TH3F(
          TString::Format("hDeltaEtaDeltaPhiStarSEVsCent_%s", c.name.Data()).Data(),
          "SE two-track QA;#Delta#eta;min #Delta#phi* (rad);cent9",
          160, -0.16, 0.16, 160, -0.16, 0.16, kNCent, kCentLo, kCentHi);
      c.deltaEtaPhiStarMECent = new TH3F(
          TString::Format("hDeltaEtaDeltaPhiStarMEVsCent_%s", c.name.Data()).Data(),
          "ME two-track QA;#Delta#eta;min #Delta#phi* (rad);cent9",
          160, -0.16, 0.16, 160, -0.16, 0.16, kNCent, kCentLo, kCentHi);
      // Ten stages, not eight: the resonance channels also drop pairs on the channel mass
      // window, and a cut flow that hides that stage cannot be compared with the track-track
      // one. The mass bins stay empty for p-p / p-d / d-d.
      c.cutFlow = new TH1D(TString::Format("hPairCutFlow_%s", c.name.Data()).Data(),
                           "pair cut flow;stage;pairs", 10, 0.5, 10.5);
      const Char_t* labels[10] = {"SE eligible","SE shared","SE mass","SE close","SE accepted",
                                  "ME eligible","ME shared","ME mass","ME close","ME accepted"};
      for (Int_t ibin = 1; ibin <= 10; ++ibin) c.cutFlow->GetXaxis()->SetBinLabel(ibin, labels[ibin-1]);
      const Char_t* kinName[8] = {"PartAPt","PartAP","PartARapidity","PartAEta",
                                  "PartBPt","PartBP","PartBRapidity","PartBEta"};
      const Char_t* kinTitle[8] = {"part A p_{T};p_{T} (GeV/c);cent9",
                                   "part A p;p (GeV/c);cent9",
                                   "part A y_{CM};y_{CM};cent9",
                                   "part A #eta;#eta;cent9",
                                   "part B p_{T};p_{T} (GeV/c);cent9",
                                   "part B p;p (GeV/c);cent9",
                                   "part B y_{CM};y_{CM};cent9",
                                   "part B #eta;#eta;cent9"};
      for (Int_t ik = 0; ik < 8; ++ik) {
        const Bool_t momentum = (ik == 0 || ik == 1 || ik == 4 || ik == 5);
        const Bool_t rapidity = (ik == 2 || ik == 6);
        const Int_t nb = momentum ? 250 : (rapidity ? kNPairY : 200);
        const Double_t lo = momentum ? 0.0 : (rapidity ? kPairYLo : -2.5);
        const Double_t hi = momentum ? 5.0 : (rapidity ? kPairYHi : 2.5);
        c.constituentSECent[ik] = new TH2D(
            TString::Format("h%sSEVsCent_%s", kinName[ik], c.name.Data()).Data(),
            kinTitle[ik], nb, lo, hi, kNCent, kCentLo, kCentHi);
        c.constituentMECent[ik] = new TH2D(
            TString::Format("h%sMEVsCent_%s", kinName[ik], c.name.Data()).Data(),
            kinTitle[ik], nb, lo, hi, kNCent, kCentLo, kCentHi);
      }
    }
    c.closeSE = new TH1D(TString::Format("hClosePairRejectSE_%s", c.name.Data()).Data(),
                         "SE pairs removed by the close-pair cut;k* (GeV/c);pairs",
                         kNKstar, kKstarLo, kKstarHi);
    c.closeME = new TH1D(TString::Format("hClosePairRejectME_%s", c.name.Data()).Data(),
                         "ME pairs removed by the close-pair cut;k* (GeV/c);pairs",
                         kNKstar, kKstarLo, kKstarHi);
    c.nSE = c.nME = c.nCloseSE = c.nCloseME = c.nShared = 0;
    chans.push_back(c);
  }
  std::cout << "[downstreamV3] channels: " << chans.size() << " of " << channelInputs.size()
            << " enabled in the selected channel config" << std::endl;

  std::map<Int_t, std::deque<MixEvent> > pool;
  Long64_t nPartTot[kNPart] = {0};
  Long64_t nATot[kNA] = {0};
  Long64_t nShared = 0, nPhiTot = 0, nSameEventME = 0, nDcaKKReject = 0, nNoOrigin = 0;
  ULong64_t nMixAttemptedForward = 0, nMixAttemptedReverse = 0;
  ULong64_t nMixStoredForward = 0, nMixStoredReverse = 0, nMixCapHitEvents = 0;
  ULong64_t nThreeEventPopulation = 0, nThreeEventRawSampled = 0;
  ULong64_t nThreeEventAccepted = 0, nThreeEventPairCutRejected = 0;
  ULong64_t nThreeEventCapHitEvents = 0, nThreeEventPopulationEvents = 0;
  ULong64_t nThreeEventUidRejectedBlocks = 0, nThreeEventResolveErrors = 0;
  ULong64_t nThreeEventOverflowEvents = 0;
  const Long64_t nEvTotal = ev->GetEntries();
  // One job may take a slice [evBegin, evEnd) of a block instead of the whole thing. A p-p pass
  // over a full 11 M-event block is a day of wall clock, which no batch slot will hold, so the
  // block is split and the pieces are added afterwards. The only difference this makes to the
  // physics is that event mixing does not cross a slice boundary: for a slice of millions of
  // events that is a handful of events at one seam, and the mixing buffer is a few events deep.
  const Long64_t evBegin = (firstEvent > 0)
                               ? ((firstEvent < nEvTotal) ? firstEvent : nEvTotal)
                               : 0;
  const Long64_t evEnd = (maxEvents > 0 && evBegin + maxEvents < nEvTotal) ? evBegin + maxEvents
                                                                           : nEvTotal;
  const Long64_t nEv = evEnd - evBegin;
  if (evBegin > 0) {
    // Walk the track cursor to the first row of the first event of this slice. The rows are
    // grouped by event in the event tree's order, which the completeness check below verifies
    // whenever a job reads a whole block.
    ev->GetEntry(evBegin);
    const ULong64_t uidFirst = e.eventUID;
    while (trackCursor < nTr) {
      tr->GetEntry(trackCursor);
      if (t.eventUID == uidFirst) break;
      ++trackCursor;
    }
    std::cout << "[downstreamV3] slice [" << evBegin << "," << evEnd << ") of " << nEvTotal
              << " events; track cursor starts at row " << trackCursor << std::endl;
  }
  const Int_t maxMixed = mix.maxMixedPairsPerEvent;
  TRandom3 rng(1);
  // The rotated background is a random construct: the maker draws its angles from gRandom and a
  // reader cannot replay that stream, so this is an independent draw. The candidate COUNT is
  // deterministic (accepted pairs x rotationN) and must match exactly; the M(KK) shape can only
  // agree statistically.
  TRandom3 rotRng(femtoCfg.rotationSeed != 0 ? (UInt_t)femtoCfg.rotationSeed : 20260915u);
  // The fully-mixed sampler's RNG is seeded once and consumed across the whole run, exactly as
  // StFemtoMaker seeds m_phiMixRng in Init(). Reproducing phi_mix pair for pair therefore needs
  // the events visited in the same order and the sampler called on every event, including the
  // ones that yield nothing.
  femto_phi_mix::SplitMix64 mixKKRng(
      femtoCfg.fullyMixedSamplingSeed > 0 ? (UInt_t)femtoCfg.fullyMixedSamplingSeed : 314159u);

  for (Long64_t ie = evBegin; ie < evEnd; ++ie) {
    ev->GetEntry(ie);
    std::vector<Cand> kp, km;
    std::vector<Cand> A[kNA];
    std::vector<Cand> B[kNPart];
    std::map<Int_t, TrackRowV2> kmap;
    while (trackCursor < nTr) {
      tr->GetEntry(trackCursor);
      if (t.eventUID != e.eventUID) break;   // next event's rows start here
      ++trackCursor;
      ++nTracksSeen;
      if (t.speciesCode == kSpeciesKp || t.speciesCode == kSpeciesKm) {
        if (data006Enabled) continue;
        if (!PassKaonForPhi(t, var)) continue;
        kmap[t.trackIndex] = t;
        Cand c;
        c.Reset();
        c.trackIndex = t.trackIndex;
        c.AddCon(e.eventUID, t.trackIndex);
        TVector3 p = Momentum(t);
        c.p4 = TLorentzVector(p, TMath::Sqrt(kKaonMass * kKaonMass + p.Mag2()));
        c.tEta = t.Eta(); c.tPhi = t.Phi(); c.tPt = t.Pt(); c.tCharge = (Short_t)t.Charge();
        c.dEta[0] = t.Eta(); c.dPhi[0] = t.Phi(); c.dPt[0] = t.Pt();
        c.dCharge[0] = (Short_t)t.Charge();
        if (t.speciesCode == kSpeciesKp) kp.push_back(c); else km.push_back(c);
      } else {
        // One branch for every bachelor species. A row that fails its own selection must not
        // stop the others from being read, which is exactly the coupling that had to be removed
        // from StFemtoMaker (closure-proton-20260914.md); here the species are separate rows, so
        // the only thing to get right is that `continue` skips a row, not a species.
        const Int_t sp = PartOfSpeciesCode(t.speciesCode);
        if (sp < 0) continue;
        if (!PassBachelor(t, var, sp)) continue;
        Cand c;
        c.Reset();
        c.trackIndex = t.trackIndex;
        c.AddCon(e.eventUID, t.trackIndex);
        c.p4 = BachelorP4(t, sp);
        c.tEta = t.Eta(); c.tPhi = t.Phi(); c.tPt = t.Pt(); c.tCharge = (Short_t)t.Charge();
        B[sp].push_back(c);
      }
    }

    // h-K species: the phi-daughter kaons under the production PID, which is what
    // StFemtoMaker::MakePhiDaughterKaonCandidate is given. Not the loose collection.
    for (size_t i = 0; i < kp.size(); ++i)
      if (PassDaughterPid(kmap[kp[i].trackIndex], var, recomputeDaughterPid))
        A[kAKPlus].push_back(kp[i]);
    for (size_t i = 0; i < km.size(); ++i)
      if (PassDaughterPid(kmap[km[i].trackIndex], var, recomputeDaughterPid))
        A[kAKMinus].push_back(km[i]);

    for (size_t i = 0; i < kp.size(); ++i) {
      for (size_t j = 0; j < km.size(); ++j) {
        const TrackRowV2& a = kmap[kp[i].trackIndex];
        const TrackRowV2& b = kmap[km[j].trackIndex];
        if (!PassDaughterPid(a, var, recomputeDaughterPid) ||
            !PassDaughterPid(b, var, recomputeDaughterPid))
          continue;
        if (dcaKKActive) {
          std::map<ULong64_t, KaonOriginRow>::const_iterator ia =
              originOf.find(e.eventUID * 100000ull + (ULong64_t)a.trackIndex);
          std::map<ULong64_t, KaonOriginRow>::const_iterator ib =
              originOf.find(e.eventUID * 100000ull + (ULong64_t)b.trackIndex);
          if (ia == originOf.end() || ib == originOf.end()) { ++nNoOrigin; continue; }
          TVector3 dcaPosA, dcaPosB;
          const Double_t dcaKK = StPhiKKReconstruction::CalculateDCA(
              MakeKkState(a, ia->second, e.bField), MakeKkState(b, ib->second, e.bField), dcaPosA,
              dcaPosB);
          hDcaKK->Fill(dcaKK);
          if (dcaKK > phiCfg.maxDCAKK) { ++nDcaKKReject; continue; }
        }
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
        const Bool_t passUnrotatedPair =
            (opening >= phiCfg.minOpeningAngle && opening <= phiCfg.maxOpeningAngle &&
             yPair >= phiCfg.minPairRapidity && yPair <= phiCfg.maxPairRapidity);
        if (passUnrotatedPair) {
          hMkk->Fill(invMass);
          Cand phi;
          phi.Reset();
          phi.AddCon(e.eventUID, a.trackIndex);
          phi.AddCon(e.eventUID, b.trackIndex);
          phi.mKK = (Float_t)invMass;
          phi.p4 = TLorentzVector(phiMom, etot);
          phi.dEta[0] = a.Eta(); phi.dPhi[0] = a.Phi(); phi.dPt[0] = a.Pt();
          phi.dCharge[0] = (Short_t)a.Charge();
          phi.dEta[1] = b.Eta(); phi.dPhi[1] = b.Phi(); phi.dPt[1] = b.Pt();
          phi.dCharge[1] = (Short_t)b.Charge();
          A[kAPhi].push_back(phi);
          nPhiTot++;
        }

        // preRotationLegacy reproduces the historical K+ rotation exactly. postRotation starts
        // from every PID/DCA-qualified pair and applies opening-angle/rapidity cuts to each leg's
        // rotated kinematics. Both charge legs share the same deltaPhi draw.
        const Bool_t postRotation = (femtoCfg.rotationPairCutMode == "postRotation");
        if (femtoCfg.rotationEnabled && (postRotation || passUnrotatedPair)) {
          for (Int_t irot = 0; irot < femtoCfg.rotationN; ++irot) {
            const Double_t dPhiRot =
                rotRng.Uniform(femtoCfg.rotationMinAngle, femtoCfg.rotationMaxAngle);
            if (femtoCfg.rotationChargeMode != "legacyKm") {
              Cand rotKp;
              if (BuildRotatedCandidate(a, b, e.eventUID, dPhiRot, kTRUE, postRotation,
                                        phiCfg, rotKp)) {
                // The legacy store remains the K+ leg in bothSeparated mode, making regression
                // checks exact rather than merely statistical.
                A[kAPhiRot].push_back(rotKp);
                hPhiRotMkk->Fill(rotKp.mKK);
                if (femtoCfg.rotationChargeMode == "bothSeparated")
                  A[kAPhiRotKp].push_back(rotKp);
              }
            }
            if (femtoCfg.rotationChargeMode != "legacyKp") {
              Cand rotKm;
              if (BuildRotatedCandidate(a, b, e.eventUID, dPhiRot, kFALSE, postRotation,
                                        phiCfg, rotKm)) {
                if (femtoCfg.rotationChargeMode == "legacyKm") {
                  A[kAPhiRot].push_back(rotKm);
                  hPhiRotMkk->Fill(rotKm.mKK);
                } else {
                  A[kAPhiRotKm].push_back(rotKm);
                }
              }
            }
          }
        }
      }
    }

    // ---- fully-mixed KK background ----
    // StFemtoMaker::BuildFullyMixedPhiCandidates reads the pool before StoreEventForMixing, so
    // the current event is not in it, and draws with the shared sampler under the
    // fullyMixedMaxCandidates cap. The population is the production-PID kaons on both sides,
    // which is exactly what A[kAKPlus] / A[kAKMinus] hold.
    const Int_t mixBinNow =
        recomputeMixBin ? MixBinOf(e.Vz(), (Int_t)e.cent9, e.Psi2()) : (Int_t)e.mixBin;
    if (!data005DiagnosticsOnly && femtoCfg.fullyMixedEnabled &&
        !(A[kAKPlus].empty() && A[kAKMinus].empty())) {
      std::deque<MixEvent>& poolNow = pool[mixBinNow];
      if (!poolNow.empty()) {
        std::vector<std::vector<Cand> > bufKp(poolNow.size()), bufKm(poolNow.size());
        std::vector<ULong64_t> bufUID(poolNow.size());
        std::vector<femto_mixing::EventCandidateCounts> counts;
        counts.reserve(poolNow.size());
        for (size_t ib = 0; ib < poolNow.size(); ++ib) {
          bufKp[ib] = poolNow[ib].a[kAKPlus];
          bufKm[ib] = poolNow[ib].a[kAKMinus];
          bufUID[ib] = poolNow[ib].eventUID;
          counts.push_back(femto_mixing::EventCandidateCounts(bufKp[ib].size(), bufKm[ib].size()));
        }
        const femto_mixing::SamplingPlan plan = femto_mixing::BuildSamplingPlan(
            A[kAKPlus].size(), A[kAKMinus].size(), counts, true);
        const femto_phi_mix::PairCount maxCand =
            (femtoCfg.fullyMixedMaxCandidates > 0)
                ? (femto_phi_mix::PairCount)femtoCfg.fullyMixedMaxCandidates
                : 0;
        MixKKEval evalKK;
        evalKK.curKp = &A[kAKPlus];
        evalKK.curKm = &A[kAKMinus];
        evalKK.bufKp = bufKp.empty() ? 0 : &bufKp[0];
        evalKK.bufKm = bufKm.empty() ? 0 : &bufKm[0];
        evalKK.bufUID = bufUID.empty() ? 0 : &bufUID[0];
        evalKK.nPool = poolNow.size();
        evalKK.curUID = e.eventUID;
        evalKK.phiCfg = &phiCfg;
        evalKK.out = &A[kAPhiMix];
        std::vector<Bool_t> reverseTags;
        evalKK.reverseTags = &reverseTags;
        std::vector<femto_phi_mix::PairCount> stored;
        femto_phi_mix::CapSampleStats stats;
        femto_phi_mix::SampleEligiblePairs(plan, maxCand, mixKKRng, evalKK, stored, stats);
        // The sampler may evaluate the cap+1 eligible pair to learn whether the cap was hit;
        // that pair is not stored, so drop the overflow exactly as the maker does.
        if (A[kAPhiMix].size() > stored.size()) A[kAPhiMix].resize(stored.size());
        if (reverseTags.size() > stored.size()) reverseTags.resize(stored.size());
        nMixAttemptedForward += (ULong64_t)stats.attemptedForward;
        nMixAttemptedReverse += (ULong64_t)stats.attemptedReverse;
        nMixStoredForward += (ULong64_t)stats.storedForward;
        nMixStoredReverse += (ULong64_t)stats.storedReverse;
        if (stats.capHit) ++nMixCapHitEvents;
        if (femtoCfg.mixChargeMode == "bothSeparated") {
          for (size_t im = 0; im < A[kAPhiMix].size(); ++im) {
            if (reverseTags[im])
              A[kAPhiMixCurKm].push_back(A[kAPhiMix][im]);
            else
              A[kAPhiMixCurKp].push_back(A[kAPhiMix][im]);
          }
        }
      }
    }

    hNPhi->Fill(A[kAPhi].size());
    for (Int_t ia = 0; ia < kNA; ++ia) {
      nATot[ia] += A[ia].size();
      hNA[ia]->Fill(A[ia].size());
    }
    for (Int_t sp = 0; sp < kNPart; ++sp) {
      nPartTot[sp] += B[sp].size();
      hNPart[sp]->Fill(B[sp].size());
    }

    const Double_t centX = (e.cent9 >= 0) ? (Double_t)e.cent9 : -0.5;
    hAcceptedEventsVsCent9->Fill(centX);
    pFxtMultVsCent9->Fill(centX, (Double_t)e.rawMult);
    hFxtMultSumVsCent9->Fill(centX, (Double_t)e.rawMult);
    hFxtMultSum2VsCent9->Fill(centX, (Double_t)e.rawMult * (Double_t)e.rawMult);

    // ---- fully factorised three-event background ----
    // Target distribution: every d(E0) x K+(E1) x K-(E2) triplet with all three
    // eventUIDs distinct and all events in the same mixing bin. A fixed number of raw
    // triplets is sampled without replacement. Each accepted triplet carries the inverse
    // sampling fraction, so a cap controls CPU without changing the expected distribution.
    if (threeEventEnabled && hThreeEvent && !B[kPartDeuteron].empty()) {
      std::deque<MixEvent>& pool3 = pool[mixBinNow];
      std::vector<ThreeEventBlock> blocks;
      femto_phi_mix::PairCount population = 0;
      Bool_t overflow = kFALSE;
      for (size_t ip = 0; ip < pool3.size() && !overflow; ++ip) {
        if (pool3[ip].eventUID == e.eventUID) {
          ++nThreeEventUidRejectedBlocks;
          continue;
        }
        const femto_phi_mix::PairCount nKp = pool3[ip].a[kAKPlus].size();
        if (nKp == 0) continue;
        for (size_t im = 0; im < pool3.size(); ++im) {
          if (im == ip) continue;
          if (pool3[im].eventUID == e.eventUID ||
              pool3[im].eventUID == pool3[ip].eventUID) {
            ++nThreeEventUidRejectedBlocks;
            continue;
          }
          const femto_phi_mix::PairCount nKm = pool3[im].a[kAKMinus].size();
          if (nKm == 0) continue;
          femto_phi_mix::PairCount blockSize = 0;
          if (!CheckedTripletCount(B[kPartDeuteron].size(), nKp, nKm, blockSize) ||
              blockSize > std::numeric_limits<femto_phi_mix::PairCount>::max() - population) {
            overflow = kTRUE;
            break;
          }
          ThreeEventBlock block;
          block.kpEvent = ip;
          block.kmEvent = im;
          block.begin = population;
          block.end = population + blockSize;
          block.nD = B[kPartDeuteron].size();
          block.nKp = nKp;
          block.nKm = nKm;
          blocks.push_back(block);
          population = block.end;
        }
      }
      if (overflow) {
        ++nThreeEventOverflowEvents;
      } else if (population > 0) {
        ++nThreeEventPopulationEvents;
        nThreeEventPopulation += population;
        const femto_phi_mix::PairCount requested =
            (femto_phi_mix::PairCount)threeEventMaxRawSamplesPerEvent;
        const femto_phi_mix::PairCount nDraw = (requested < population) ? requested : population;
        if (nDraw < population) ++nThreeEventCapHitEvents;
        const Double_t samplingWeight = (Double_t)population / (Double_t)nDraw;
        // Per-event seeding makes a cap scan nested (cap N is a prefix of cap 2N) and
        // keeps a block reproducible even when it is processed in slices.
        femto_phi_mix::SplitMix64 threeEventRng(
            (femto_phi_mix::PairCount)threeEventSamplingSeed ^
            ((femto_phi_mix::PairCount)e.eventUID + 0x9E3779B97F4A7C15ULL));
        femto_phi_mix::LazyPermutationSampler sampler(population, threeEventRng);
        for (femto_phi_mix::PairCount idraw = 0; idraw < nDraw; ++idraw) {
          femto_phi_mix::PairCount flat = 0;
          if (!sampler.Next(flat)) {
            ++nThreeEventResolveErrors;
            break;
          }
          size_t iblock = 0, iD = 0, iKp = 0, iKm = 0;
          if (!ResolveThreeEventIndex(blocks, flat, iblock, iD, iKp, iKm)) {
            ++nThreeEventResolveErrors;
            continue;
          }
          ++nThreeEventRawSampled;
          const ThreeEventBlock& block = blocks[iblock];
          const Cand& d = B[kPartDeuteron][iD];
          const Cand& kp3 = pool3[block.kpEvent].a[kAKPlus][iKp];
          const Cand& km3 = pool3[block.kmEvent].a[kAKMinus][iKm];
          const TLorentzVector pKK = kp3.p4 + km3.p4;
          const Double_t invMass = pKK.M();
          if (!(invMass > 0.0)) {
            ++nThreeEventPairCutRejected;
            continue;
          }
          const Double_t opening = kp3.p4.Vect().Angle(km3.p4.Vect());
          const Double_t yPair = phiCfg.ApplyAnalysisRapidity(pKK.Rapidity());
          if (opening < phiCfg.minOpeningAngle || opening > phiCfg.maxOpeningAngle ||
              yPair < phiCfg.minPairRapidity || yPair > phiCfg.maxPairRapidity) {
            ++nThreeEventPairCutRejected;
            continue;
          }
          hThreeEvent->Fill(invMass, KStar(pKK, d.p4), centX, samplingWeight);
          ++nThreeEventAccepted;
        }
      }
    }

    // ---- same event ----
    for (size_t icha = 0; icha < chans.size(); ++icha) {
      ChanHists& c = chans[icha];
      const std::vector<Cand>& av = c.trackTrack ? B[c.aPart] : A[c.ia];
      const std::vector<Cand>& bv = B[c.ib];
      for (size_t i = 0; i < av.size(); ++i) {
        const size_t jBegin = c.identical ? i + 1 : 0;
        for (size_t j = jBegin; j < bv.size(); ++j) {
          if (c.cutFlow) c.cutFlow->Fill(1);
          if (femtoCfg.closePairVetoSameTrack && SharedTrack(av[i], bv[j])) {
            c.nShared++;
            nShared++;
            hRejectShared->Fill(1);
            if (c.cutFlow) c.cutFlow->Fill(2);
            continue;
          }
          const Double_t ks = KStar(av[i].p4, bv[j].p4);
          if (c.resonance && c.mkkSEWide) c.mkkSEWide->Fill(av[i].mKK, ks, centX);
          if (c.resonance && (av[i].mKK < c.mLo || av[i].mKK > c.mHi)) {
            if (c.cutFlow) c.cutFlow->Fill(3);
            continue;
          }

          Bool_t rejectClose = kFALSE;
          if (c.trackTrack) {
            Double_t dEta = 0.0, dPhiStar = 0.0;
            rejectClose = Data006ClosePair(av[i], bv[j], e.bField / 10.0, data006Cfg, c,
                                           dEta, dPhiStar);
            c.deltaEtaPhiStarSECent->Fill(dEta, dPhiStar, centX);
          } else {
            FillPhiBachelorCloseQA(c.deltaEtaPhiStarSECent, av[i], bv[j].tEta, bv[j].tPhi,
                                   bv[j].tPt, bv[j].tCharge, e.bField / 10.0, femtoCfg, centX);
            rejectClose = ClosePairReject(av[i], bv[j], bv[j].tEta, bv[j].tPhi, bv[j].tPt,
                                          bv[j].tCharge, e.bField / 10.0, c.ib, femtoCfg);
          }
          if (rejectClose) {
            c.nCloseSE++;
            c.closeSE->Fill(ks);
            if (c.cutFlow) c.cutFlow->Fill(4);
            continue;
          }

          c.se->Fill(ks);
          c.seCent->Fill(ks, centX);
          if (c.cutFlow) c.cutFlow->Fill(5);
          if (c.extendedQa && QaMassPass(c, av[i], pairMtSignal)) {
            const TLorentzVector pair = av[i].p4 + bv[j].p4;
            c.pairMtKstarSECent->Fill(0.5 * pair.Mt(), ks, centX);
            const Double_t pairY = phiCfg.ApplyAnalysisRapidity(pair.Rapidity());
            c.pairYSECent->Fill(pairY, centX);
            if (ks >= c.normQMin && ks < c.normQMax)
              c.pairYNormSECent->Fill(pairY, centX);
            const Double_t kin[8] = {
                av[i].p4.Pt(), av[i].p4.P(), phiCfg.ApplyAnalysisRapidity(av[i].p4.Rapidity()),
                c.trackTrack ? av[i].tEta : av[i].p4.Eta(), bv[j].p4.Pt(), bv[j].p4.P(),
                phiCfg.ApplyAnalysisRapidity(bv[j].p4.Rapidity()), bv[j].tEta};
            for (Int_t ik = 0; ik < 8; ++ik) c.constituentSECent[ik]->Fill(kin[ik], centX);
          }
          if (c.mkkSE) c.mkkSE->Fill(av[i].mKK, ks, centX);
          c.nSE++;
        }
      }
    }

    // ---- mixed event ----
    std::deque<MixEvent>& binPool = pool[mixBinNow];
    // Channels share the event pool but never share a pair, so the sampling is done per channel:
    // with maxMixedPairsPerEvent set, each channel gets its own cap, exactly as the maker draws
    // separately per channel.
    for (size_t icha = 0; icha < chans.size(); ++icha) {
      ChanHists& c = chans[icha];
      if (!c.doMixing) continue;
      std::vector<MixRef> refs;
      const std::vector<Cand>& curA = c.trackTrack ? B[c.aPart] : A[c.ia];
      for (size_t ib = 0; ib < binPool.size(); ++ib) {
        const MixEvent& buf = binPool[ib];
        const std::vector<Cand>& bufA = c.trackTrack ? buf.part[c.aPart] : buf.a[c.ia];
        for (size_t i = 0; i < curA.size(); ++i)
          for (size_t j = 0; j < buf.part[c.ib].size(); ++j) {
            MixRef r; r.ib = ib; r.reverse = 0; r.i = i; r.j = j; refs.push_back(r);
          }
        // For identical species, current x buffer already contains every unordered cross-event
        // pair once. Adding the reverse direction would duplicate the exact same physical pairs.
        if (mix.mixBothDirections && !c.identical) {
          for (size_t i = 0; i < bufA.size(); ++i)
            for (size_t j = 0; j < B[c.ib].size(); ++j) {
              MixRef r; r.ib = ib; r.reverse = 1; r.i = i; r.j = j; refs.push_back(r);
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
        const std::vector<Cand>& bufA = c.trackTrack ? buf.part[c.aPart] : buf.a[c.ia];
        const Cand& AA = r.reverse ? bufA[r.i] : curA[r.i];
        const Cand& BB = r.reverse ? B[c.ib][r.j] : buf.part[c.ib][r.j];
        if (c.cutFlow) c.cutFlow->Fill(6);
        if (femtoCfg.closePairVetoSameTrack && SharedTrack(AA, BB)) {
          c.nShared++;
          nShared++;
          hRejectShared->Fill(1);
          if (c.cutFlow) c.cutFlow->Fill(7);
          continue;
        }
        const Double_t ks = KStar(AA.p4, BB.p4);
        if (c.resonance && c.mkkMEWide) c.mkkMEWide->Fill(AA.mKK, ks, centX);
        if (c.resonance && (AA.mKK < c.mLo || AA.mKK > c.mHi)) {
          if (c.cutFlow) c.cutFlow->Fill(8);
          continue;
        }
        const Double_t bT = (r.reverse ? buf.bField : e.bField) / 10.0;
        Bool_t rejectClose = kFALSE;
        if (c.trackTrack) {
          Double_t dEta = 0.0, dPhiStar = 0.0;
          rejectClose = Data006ClosePair(AA, BB, bT, data006Cfg, c, dEta, dPhiStar);
          c.deltaEtaPhiStarMECent->Fill(dEta, dPhiStar, centX);
        } else {
          FillPhiBachelorCloseQA(c.deltaEtaPhiStarMECent, AA, BB.tEta, BB.tPhi, BB.tPt,
                                 BB.tCharge, bT, femtoCfg, centX);
          rejectClose = ClosePairReject(AA, BB, BB.tEta, BB.tPhi, BB.tPt, BB.tCharge,
                                        bT, c.ib, femtoCfg);
        }
        if (rejectClose) {
          c.nCloseME++;
          c.closeME->Fill(ks);
          if (c.cutFlow) c.cutFlow->Fill(9);
          continue;
        }
        c.me->Fill(ks);
        c.meCent->Fill(ks, centX);
        if (c.cutFlow) c.cutFlow->Fill(10);
        if (c.extendedQa && QaMassPass(c, AA, pairMtSignal)) {
          const TLorentzVector pair = AA.p4 + BB.p4;
          c.pairMtKstarMECent->Fill(0.5 * pair.Mt(), ks, centX);
          const Double_t pairY = phiCfg.ApplyAnalysisRapidity(pair.Rapidity());
          c.pairYMECent->Fill(pairY, centX);
          if (ks >= c.normQMin && ks < c.normQMax)
            c.pairYNormMECent->Fill(pairY, centX);
          const Double_t kin[8] = {
              AA.p4.Pt(), AA.p4.P(), phiCfg.ApplyAnalysisRapidity(AA.p4.Rapidity()),
              c.trackTrack ? AA.tEta : AA.p4.Eta(), BB.p4.Pt(), BB.p4.P(),
              phiCfg.ApplyAnalysisRapidity(BB.p4.Rapidity()), BB.tEta};
          for (Int_t ik = 0; ik < 8; ++ik) c.constituentMECent[ik]->Fill(kin[ik], centX);
        }
        if (c.mkkME) c.mkkME->Fill(AA.mKK, ks, centX);
        c.nME++;
      }
    }

    // Every accepted event enters the pool, including ones with no phi and no bachelor.
    // StFemtoMaker::StoreEventForMixing guards only on m_eventCandidates.empty(), and
    // m_eventCandidates[speciesKey] inserts a key for every enabled species on every event, so
    // that guard never fires. An empty event contributes no pairs but still occupies a buffer
    // slot and evicts an older one, and 18.7% of events here have no nominal deuteron, so
    // skipping them would give the downstream a different buffer than the maker.
    if (!data005DiagnosticsOnly) {
      MixEvent me;
      me.eventUID = e.eventUID;
      me.bField = e.bField;
      me.cent9 = (Int_t)e.cent9;
      for (Int_t ia = 0; ia < kNA; ++ia) me.a[ia] = A[ia];
      for (Int_t sp = 0; sp < kNPart; ++sp) me.part[sp] = B[sp];
      binPool.push_back(me);
      if ((Int_t)binPool.size() > bufferSize) binPool.pop_front();
    }
  }

  // The cursor walk assumes the track rows are grouped by event and in the same order as the event
  // rows. If that ever stops being true -- a differently ordered producer, a merge that interleaves
  // -- rows would be silently skipped, so say so rather than quietly analysing a subset.
  if (evBegin == 0 && evEnd == nEvTotal && nTracksSeen != nTr) {
    std::cerr << "ERROR: read " << nTracksSeen << " of " << nTr << " track rows. The track tree is "
              << "not grouped by event in the event tree's order, which this reader requires.\n"
              << "       Refusing to report results from a subset of the tree." << std::endl;
    fout->Close();
    fin->Close();
    return;
  }
  std::cout << "[downstreamV3] events=" << nEv << "/" << nEvTotal << " [" << evBegin << ","
            << evEnd << ") phi=" << nPhiTot
            << std::endl;
  for (Int_t ia = 0; ia < kNA; ++ia)
    std::cout << "[downstreamV3]   A " << kAName[ia] << ": n=" << nATot[ia] << std::endl;
  for (Int_t sp = 0; sp < kNPart; ++sp)
    std::cout << "[downstreamV3]   B " << kPartName[sp] << ": n=" << nPartTot[sp] << std::endl;
  for (size_t icha = 0; icha < chans.size(); ++icha) {
    const ChanHists& c = chans[icha];
    std::cout << "[downstreamV3]   " << c.name << ": SE=" << c.nSE << " ME=" << c.nME;
    if ((data006Enabled && data006Cfg.closePairEnabled) ||
        (!data006Enabled && femtoCfg.closePairEnabled))
      std::cout << " closeSE=" << c.nCloseSE << " closeME=" << c.nCloseME;
    std::cout << std::endl;
  }
  std::cout << "[downstreamV3] shared=" << nShared << " sameEventME=" << nSameEventME
            << " dcaKKReject=" << nDcaKKReject << " noOrigin=" << nNoOrigin << std::endl;
  if (threeEventEnabled) {
    std::cout << "[downstreamV3] threeEvent population=" << nThreeEventPopulation
              << " sampled=" << nThreeEventRawSampled << " accepted=" << nThreeEventAccepted
              << " pairCutRejected=" << nThreeEventPairCutRejected
              << " capHitEvents=" << nThreeEventCapHitEvents
              << " uidRejectedBlocks=" << nThreeEventUidRejectedBlocks
              << " resolveErrors=" << nThreeEventResolveErrors
              << " overflowEvents=" << nThreeEventOverflowEvents << std::endl;
  }

  fout->cd();
  TNamed("treeFile", treeFile).Write();
  TNamed("nEventsProcessed", TString::Format("%lld", nEv).Data()).Write();
  TNamed("eventSliceBegin", TString::Format("%lld", evBegin).Data()).Write();
  TNamed("eventSliceEnd", TString::Format("%lld", evEnd).Data()).Write();
  TNamed("nEventsTotal", TString::Format("%lld", nEvTotal).Data()).Write();
  TNamed("maxEvents", TString::Format("%lld", maxEvents).Data()).Write();
  TNamed("data006Enabled", data006Enabled ? "true" : "false").Write();
  TNamed("data006Config", data006Enabled ? data006ConfigPath : "NOT_AVAILABLE").Write();
  TNamed("mixingMode", mode.Data()).Write();
  TNamed("rotationChargeMode", femtoCfg.rotationChargeMode.c_str()).Write();
  TNamed("rotationPairCutMode", femtoCfg.rotationPairCutMode.c_str()).Write();
  TNamed("mixChargeMode", femtoCfg.mixChargeMode.c_str()).Write();
  TNamed("schemaVersion", TString::Format("%u", e.schemaVersion).Data()).Write();
  TNamed("nP", TString::Format("%lld", nPartTot[kPartProton]).Data()).Write();
  for (size_t icha = 0; icha < chans.size(); ++icha) {
    const ChanHists& c = chans[icha];
    TNamed(TString::Format("nSE_%s", c.name.Data()).Data(),
           TString::Format("%lld", c.nSE).Data()).Write();
    TNamed(TString::Format("nME_%s", c.name.Data()).Data(),
           TString::Format("%lld", c.nME).Data()).Write();
  }
  for (Int_t ia = 0; ia < kNA; ++ia)
    TNamed(TString::Format("nCandA_%s", kAName[ia]).Data(),
           TString::Format("%lld", nATot[ia]).Data()).Write();
  // Exact per-species candidate totals. hN<species> saturates -- one bin per candidate -- so a
  // reader must take the totals from here, not from the histogram. The SE/ME totals are written
  // per CHANNEL above, since a species can feed several channels.
  for (Int_t sp = 0; sp < kNPart; ++sp)
    TNamed(TString::Format("nCand_%s", kPartName[sp]).Data(),
           TString::Format("%lld", nPartTot[sp]).Data()).Write();
  TNamed("nPhi", TString::Format("%lld", nPhiTot).Data()).Write();
  TNamed("nD", TString::Format("%lld", nPartTot[kPartDeuteron]).Data()).Write();
  TNamed("variation", var.spec.Data()).Write();
  TNamed("daughterPid",
         (recomputeDaughterPid || var.anyPid) ? "recomputed" : "kSelNominalPid").Write();
  TNamed("mixBinSource", recomputeMixBin ? "recomputed" : "stored").Write();
  TNamed("maxDCAKK", TString::Format("%.3f%s", phiCfg.maxDCAKK,
                                     dcaKKActive ? " (active)" : " (off)").Data()).Write();

  TString closePairSummary("disabled");
  if (data006Enabled && data006Cfg.closePairEnabled) {
    closePairSummary = TString::Format(
        "%s r=%.2f-%.2f/%.2f; per-channel windows in DATA-006 config",
        data006Cfg.closePairShape.c_str(), data006Cfg.closePairRadiusMin,
        data006Cfg.closePairRadiusMax, data006Cfg.closePairRadiusStep);
  } else if (!data006Enabled && femtoCfg.closePairEnabled) {
    closePairSummary = TString::Format(
        "%s dEtaD=%.3f dPhiD=%.3f dEtaP=%.3f dPhiP=%.3f r=%.2f-%.2f/%.2f",
        femtoCfg.closePairShape.c_str(), femtoCfg.closePairDEtaDeuteron,
        femtoCfg.closePairDPhiStarDeuteron, femtoCfg.closePairDEtaProton,
        femtoCfg.closePairDPhiStarProton, femtoCfg.closePairRadiusMin,
        femtoCfg.closePairRadiusMax, femtoCfg.closePairRadiusStep);
  }
  TNamed("closePair", closePairSummary.Data()).Write();
  TNamed("minNHitsDedxNuclear",
         TString::Format("%d", (var.dNHitsDedx > 0)
                                   ? var.dNHitsDedx
                                   : (Int_t)ConfigManager::GetInstance()
                                         .GetNuclearIdCuts()
                                         .minNHitsDedxNuclear).Data()).Write();
  TNamed("nSameEventME", TString::Format("%lld", nSameEventME).Data()).Write();
  TNamed("nMixAttemptedForward", TString::Format("%llu", nMixAttemptedForward).Data()).Write();
  TNamed("nMixAttemptedReverse", TString::Format("%llu", nMixAttemptedReverse).Data()).Write();
  TNamed("nMixStoredForward", TString::Format("%llu", nMixStoredForward).Data()).Write();
  TNamed("nMixStoredReverse", TString::Format("%llu", nMixStoredReverse).Data()).Write();
  TNamed("nMixCapHitEvents", TString::Format("%llu", nMixCapHitEvents).Data()).Write();
  TNamed("threeEventEnabled", threeEventEnabled ? "true" : "false").Write();
  TNamed("threeEventTopology",
         "d=current E0; K+=buffer E1; K-=buffer E2; E0,E1,E2 eventUID-distinct; same mix bin").Write();
  TNamed("threeEventSampling",
         "uniform without replacement over raw triplets per E0; accepted fills weighted by population/nDraw").Write();
  TNamed("threeEventMaxRawSamplesPerEvent",
         TString::Format("%lld", threeEventMaxRawSamplesPerEvent).Data()).Write();
  TNamed("threeEventSamplingSeed",
         TString::Format("%lld", threeEventSamplingSeed).Data()).Write();
  TNamed("nThreeEventPopulation", TString::Format("%llu", nThreeEventPopulation).Data()).Write();
  TNamed("nThreeEventPopulationEvents",
         TString::Format("%llu", nThreeEventPopulationEvents).Data()).Write();
  TNamed("nThreeEventRawSampled", TString::Format("%llu", nThreeEventRawSampled).Data()).Write();
  TNamed("nThreeEventAccepted", TString::Format("%llu", nThreeEventAccepted).Data()).Write();
  TNamed("nThreeEventPairCutRejected",
         TString::Format("%llu", nThreeEventPairCutRejected).Data()).Write();
  TNamed("nThreeEventCapHitEvents",
         TString::Format("%llu", nThreeEventCapHitEvents).Data()).Write();
  TNamed("nThreeEventUidRejectedBlocks",
         TString::Format("%llu", nThreeEventUidRejectedBlocks).Data()).Write();
  TNamed("nThreeEventResolveErrors",
         TString::Format("%llu", nThreeEventResolveErrors).Data()).Write();
  TNamed("nThreeEventOverflowEvents",
         TString::Format("%llu", nThreeEventOverflowEvents).Data()).Write();
  TNamed("multiplicityEstimator",
         "raw StPicoEvent::fxtMult() = numberOfPrimaryTracks(); no efficiency correction").Write();
  TNamed("pairMtDefinition",
         data006Enabled
             ? "0.5*sqrt((E_A+E_B)^2-(pz_A+pz_B)^2), accepted SE/ME track-track pairs"
             : "0.5*sqrt((E_phi+E_p)^2-(pz_phi+pz_p)^2), raw SE pairs in phi_proton_signal mass window after pair cuts").Write();
  TNamed("pairMtUnit", "GeV/c^2").Write();
  TNamed("data005DiagnosticsOnly", data005DiagnosticsOnly ? "true" : "false").Write();
  TNamed("identicalPairConvention",
         data006Enabled
             ? "SE unordered i<j; ME current x buffer once; no duplicated reverse direction"
             : "NOT_APPLICABLE").Write();
  TNamed("sharedHitInformation",
         data006Enabled ? "NOT_AVAILABLE_IN_REDUCED_TREE" : "NOT_APPLICABLE").Write();
  TNamed("duplicateTrackPolicy",
         data006Enabled
             ? "one TrackRow per species/index; same (eventUID,trackIndex) rejected; identical SE uses i<j"
             : "FemtoCandidatesShareTrack").Write();
  TParameter<Int_t>("pairMtBins", pairMtBins).Write();
  TParameter<Double_t>("pairMtMin", pairMtMin).Write();
  TParameter<Double_t>("pairMtMax", pairMtMax).Write();
  if (data006Enabled) {
    TParameter<Int_t>("data006KstarBins", data006Cfg.kstarBins).Write();
    TParameter<Double_t>("data006KstarMin", data006Cfg.kstarMin).Write();
    TParameter<Double_t>("data006KstarMax", data006Cfg.kstarMax).Write();
    TParameter<Double_t>("data006DeliveryKstarMax", data006Cfg.deliveryKstarMax).Write();
  }

  hMkk->Write();
  hPhiRotMkk->Write();
  hDcaKK->Write();
  hNPhi->Write();
  for (Int_t sp = 0; sp < kNPart; ++sp) hNPart[sp]->Write();
  for (Int_t ia = 0; ia < kNA; ++ia) hNA[ia]->Write();
  for (size_t icha = 0; icha < chans.size(); ++icha) {
    ChanHists& c = chans[icha];
    c.se->Write();
    c.me->Write();
    c.seCent->Write();
    c.meCent->Write();
    if (c.mkkSE) c.mkkSE->Write();
    if (c.mkkME) c.mkkME->Write();
    if (c.mkkSEWide) c.mkkSEWide->Write();
    if (c.mkkMEWide) c.mkkMEWide->Write();
    if (c.pairMtKstarSECent) c.pairMtKstarSECent->Write();
    if (c.pairMtKstarMECent) c.pairMtKstarMECent->Write();
    if (c.pairYSECent) c.pairYSECent->Write();
    if (c.pairYMECent) c.pairYMECent->Write();
    if (c.pairYNormSECent) c.pairYNormSECent->Write();
    if (c.pairYNormMECent) c.pairYNormMECent->Write();
    if (c.deltaEtaPhiStarSECent) c.deltaEtaPhiStarSECent->Write();
    if (c.deltaEtaPhiStarMECent) c.deltaEtaPhiStarMECent->Write();
    if (c.cutFlow) c.cutFlow->Write();
    for (Int_t ik = 0; ik < 8; ++ik)
      if (c.constituentSECent[ik]) c.constituentSECent[ik]->Write();
    for (Int_t ik = 0; ik < 8; ++ik)
      if (c.constituentMECent[ik]) c.constituentMECent[ik]->Write();
    c.closeSE->Write();
    c.closeME->Write();
    TParameter<Double_t>(TString::Format("normQMin_%s", c.name.Data()).Data(),
                         c.normQMin).Write();
    TParameter<Double_t>(TString::Format("normQMax_%s", c.name.Data()).Data(),
                         c.normQMax).Write();
  }
  hAcceptedEventsVsCent9->Write();
  pFxtMultVsCent9->Write();
  hFxtMultSumVsCent9->Write();
  hFxtMultSum2VsCent9->Write();
  hRejectShared->Write();
  hSameEventME->Write();
  if (hThreeEvent) hThreeEvent->Write();
  fout->Close();
  fin->Close();
}
