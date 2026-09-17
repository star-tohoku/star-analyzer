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
#include "TH3F.h"
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
#include "StPhiKKReconstruction.h"
#include "StNuclearIdHelper.h"

#include <deque>
#include <iostream>
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
// Species that can sit on the A side of a channel: the phi and its two backgrounds, and the two
// phi-daughter kaon species that carry the h-K correlations. The maker builds each of these into
// m_eventCandidates under exactly these keys, and a channel names one of them as partA.
enum { kAPhi = 0, kAPhiRot = 1, kAPhiMix = 2, kAKPlus = 3, kAKMinus = 4, kNA = 5 };
const Char_t* const kAName[kNA] = {"phi", "phi_rot", "phi_mix", "phikaon_plus", "phikaon_minus"};

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
  Int_t ia, ib;
  Bool_t resonance;
  Double_t mLo, mHi;
  TH1D* se;
  TH1D* me;
  TH2D* seCent;
  TH2D* meCent;
  TH3F* mkkSE;
  TH3F* mkkME;
  TH3F* mkkSEWide;
  TH3F* mkkMEWide;
  TH1D* closeSE;
  TH1D* closeME;
  Long64_t nSE, nME, nCloseSE, nCloseME, nShared;
};

struct MixEvent {
  ULong64_t eventUID;
  Double_t bField;
  Int_t cent9;
  std::vector<Cand> a[kNA];
  std::vector<Cand> part[kNPart];
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
                                 Bool_t recomputeMixBin = kFALSE) {
  using namespace femto_phi_tree;
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
  std::map<ULong64_t, KaonOriginRow> originOf;
  if (ko) {
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
  const Bool_t dcaKKActive = (phiCfg.maxDCAKK < kDcaKKDisabledAbove);
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

  // ------------------------------------------------------------------------------------------
  // Channel table, built from the same YAML the maker reads. Names, mass windows and the A x B
  // combinations therefore cannot drift from the maker's, and a histogram of a given name can be
  // compared with the maker's bin for bin. The binning is the one in
  // config/hist/hist_anaFemtoPhi.yaml, so no rebinning is needed on either side.
  // ------------------------------------------------------------------------------------------
  const Int_t kNKstar = 300;
  const Double_t kKstarHi = 3.0;
  const Int_t kNMkk = 200;
  const Double_t kMkkLo = 0.98, kMkkHi = 1.18;
  const Int_t kNCent = 9;
  const Double_t kCentLo = -0.5, kCentHi = 8.5;

  std::vector<ChanHists> chans;
  for (size_t ic = 0; ic < femtoCfg.channels.size(); ++ic) {
    const FemtoConfig::ChannelDef& cd = femtoCfg.channels[ic];
    if (!cd.enabled) continue;
    ChanHists c;
    c.ia = AOfKey(cd.partA);
    c.ib = BOfKey(cd.partB);
    if (c.ia < 0 || c.ib < 0) {
      std::cout << "[downstreamV3] channel " << cd.name << " (" << cd.partA << " x " << cd.partB
                << ") skipped: species not built by this reader" << std::endl;
      continue;
    }
    c.name = cd.name.c_str();
    c.resonance = (c.ia == kAPhi || c.ia == kAPhiRot || c.ia == kAPhiMix);
    c.mLo = cd.signalMin;
    c.mHi = cd.signalMax;
    // The command-line window overrides the phi channels only; the rotated and fully-mixed
    // backgrounds keep the window their channel declares, exactly as the maker applies it.
    if (c.ia == kAPhi) {
      if (signalMin > 0) c.mLo = signalMin;
      if (signalMax > 0) c.mHi = signalMax;
    }
    c.se = new TH1D(TString::Format("hKstarSE_%s", c.name.Data()).Data(),
                    TString::Format("SE k* (%s);k* (GeV/c);pairs", c.name.Data()).Data(),
                    kNKstar, 0.0, kKstarHi);
    c.me = new TH1D(TString::Format("hKstarME_%s", c.name.Data()).Data(),
                    TString::Format("ME k* (%s);k* (GeV/c);pairs", c.name.Data()).Data(),
                    kNKstar, 0.0, kKstarHi);
    c.seCent = new TH2D(TString::Format("hKstarSEVsCent_%s", c.name.Data()).Data(),
                        TString::Format("SE k* vs cent9 (%s);k* (GeV/c);cent9", c.name.Data()).Data(),
                        kNKstar, 0.0, kKstarHi, kNCent, kCentLo, kCentHi);
    c.meCent = new TH2D(TString::Format("hKstarMEVsCent_%s", c.name.Data()).Data(),
                        TString::Format("ME k* vs cent9 (%s);k* (GeV/c);cent9", c.name.Data()).Data(),
                        kNKstar, 0.0, kKstarHi, kNCent, kCentLo, kCentHi);
    c.mkkSE = c.mkkME = c.mkkSEWide = c.mkkMEWide = 0;
    if (c.resonance) {
      c.mkkSE = new TH3F(TString::Format("hPhiMKK_vs_KstarSE_%s", c.name.Data()).Data(),
                         TString::Format("M_{KK} vs k* SE vs cent9 (%s);M_{KK} (GeV/c^{2});k* (GeV/c);cent9",
                                         c.name.Data()).Data(),
                         kNMkk, kMkkLo, kMkkHi, kNKstar, 0.0, kKstarHi, kNCent, kCentLo, kCentHi);
      c.mkkME = new TH3F(TString::Format("hPhiMKK_vs_KstarME_%s", c.name.Data()).Data(),
                         TString::Format("M_{KK} vs k* ME vs cent9 (%s);M_{KK} (GeV/c^{2});k* (GeV/c);cent9",
                                         c.name.Data()).Data(),
                         kNMkk, kMkkLo, kMkkHi, kNKstar, 0.0, kKstarHi, kNCent, kCentLo, kCentHi);
      // The "wide" TH3 is filled before the mass window, which is what the kstarMassFitCF
      // background template needs: the full M(KK) shape at each k*.
      const std::string wide = WideMkkSuffix(cd.name);
      if (!wide.empty()) {
        c.mkkSEWide = new TH3F(TString::Format("hPhiMKK_vs_KstarSE_%s", wide.c_str()).Data(),
                               TString::Format("M_{KK} vs k* SE vs cent9 (%s);M_{KK} (GeV/c^{2});k* (GeV/c);cent9",
                                               wide.c_str()).Data(),
                               kNMkk, kMkkLo, kMkkHi, kNKstar, 0.0, kKstarHi, kNCent, kCentLo, kCentHi);
        c.mkkMEWide = new TH3F(TString::Format("hPhiMKK_vs_KstarME_%s", wide.c_str()).Data(),
                               TString::Format("M_{KK} vs k* ME vs cent9 (%s);M_{KK} (GeV/c^{2});k* (GeV/c);cent9",
                                               wide.c_str()).Data(),
                               kNMkk, kMkkLo, kMkkHi, kNKstar, 0.0, kKstarHi, kNCent, kCentLo, kCentHi);
      }
    }
    c.closeSE = new TH1D(TString::Format("hClosePairRejectSE_%s", c.name.Data()).Data(),
                         "SE pairs removed by the close-pair cut;k* (GeV/c);pairs", kNKstar, 0.0,
                         kKstarHi);
    c.closeME = new TH1D(TString::Format("hClosePairRejectME_%s", c.name.Data()).Data(),
                         "ME pairs removed by the close-pair cut;k* (GeV/c);pairs", kNKstar, 0.0,
                         kKstarHi);
    c.nSE = c.nME = c.nCloseSE = c.nCloseME = c.nShared = 0;
    chans.push_back(c);
  }
  std::cout << "[downstreamV3] channels: " << chans.size() << " of " << femtoCfg.channels.size()
            << " enabled in the config" << std::endl;

  std::map<Int_t, std::deque<MixEvent> > pool;
  Long64_t nPartTot[kNPart] = {0};
  Long64_t nATot[kNA] = {0};
  Long64_t nShared = 0, nPhiTot = 0, nSameEventME = 0, nDcaKKReject = 0, nNoOrigin = 0;
  const Long64_t nEv = ev->GetEntries();
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

  for (Long64_t ie = 0; ie < nEv; ++ie) {
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
        if (opening < phiCfg.minOpeningAngle || opening > phiCfg.maxOpeningAngle) continue;
        if (yPair < phiCfg.minPairRapidity || yPair > phiCfg.maxPairRapidity) continue;
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

        // Rotated background. StFemtoMaker::BuildRotatedPhiCandidates accepts the pair on the
        // UNROTATED kinematics -- opening angle and pair rapidity are evaluated before the
        // rotation and not re-applied after it -- then turns the K+ azimuth by a uniform angle
        // and recomputes only the mass and the momentum. Reproduced here in that order.
        if (femtoCfg.rotationEnabled) {
          for (Int_t irot = 0; irot < femtoCfg.rotationN; ++irot) {
            const Double_t dPhiRot =
                rotRng.Uniform(femtoCfg.rotationMinAngle, femtoCfg.rotationMaxAngle);
            Double_t phiRot = a.Phi() + dPhiRot;
            while (phiRot > TMath::Pi()) phiRot -= TMath::TwoPi();
            while (phiRot < -TMath::Pi()) phiRot += TMath::TwoPi();
            TVector3 paRot;
            paRot.SetPtEtaPhi(a.Pt(), a.Eta(), phiRot);
            const Double_t eaRot = TMath::Sqrt(kKaonMass * kKaonMass + paRot.Mag2());
            TVector3 phiMomRot = paRot + pb;
            const Double_t m2Rot = (eaRot + eb) * (eaRot + eb) - phiMomRot.Mag2();
            if (m2Rot <= 0) continue;
            Cand rot;
            rot.Reset();
            rot.AddCon(e.eventUID, a.trackIndex);
            rot.AddCon(e.eventUID, b.trackIndex);
            rot.mKK = (Float_t)TMath::Sqrt(m2Rot);
            rot.p4 = TLorentzVector(phiMomRot, eaRot + eb);
            rot.dEta[0] = a.Eta(); rot.dPhi[0] = phiRot; rot.dPt[0] = a.Pt();
            rot.dCharge[0] = (Short_t)a.Charge();
            rot.dEta[1] = b.Eta(); rot.dPhi[1] = b.Phi(); rot.dPt[1] = b.Pt();
            rot.dCharge[1] = (Short_t)b.Charge();
            hPhiRotMkk->Fill(rot.mKK);
            A[kAPhiRot].push_back(rot);
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
    if (femtoCfg.fullyMixedEnabled && !(A[kAKPlus].empty() && A[kAKMinus].empty())) {
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
        std::vector<femto_phi_mix::PairCount> stored;
        femto_phi_mix::CapSampleStats stats;
        femto_phi_mix::SampleEligiblePairs(plan, maxCand, mixKKRng, evalKK, stored, stats);
        // The sampler may evaluate the cap+1 eligible pair to learn whether the cap was hit;
        // that pair is not stored, so drop the overflow exactly as the maker does.
        if (A[kAPhiMix].size() > stored.size()) A[kAPhiMix].resize(stored.size());
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

    // ---- same event ----
    for (size_t icha = 0; icha < chans.size(); ++icha) {
      ChanHists& c = chans[icha];
      const std::vector<Cand>& av = A[c.ia];
      const std::vector<Cand>& bv = B[c.ib];
      for (size_t i = 0; i < av.size(); ++i) {
        for (size_t j = 0; j < bv.size(); ++j) {
          if (femtoCfg.closePairVetoSameTrack && SharedTrack(av[i], bv[j])) {
            c.nShared++;
            nShared++;
            hRejectShared->Fill(1);
            continue;
          }
          const Double_t ks = KStar(av[i].p4, bv[j].p4);
          if (c.resonance && c.mkkSEWide) c.mkkSEWide->Fill(av[i].mKK, ks, centX);
          if (c.resonance && (av[i].mKK < c.mLo || av[i].mKK > c.mHi)) continue;
          if (ClosePairReject(av[i], bv[j], bv[j].tEta, bv[j].tPhi, bv[j].tPt, bv[j].tCharge,
                              e.bField / 10.0, c.ib, femtoCfg)) {
            c.nCloseSE++;
            c.closeSE->Fill(ks);
            continue;
          }
          c.se->Fill(ks);
          c.seCent->Fill(ks, centX);
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
      std::vector<MixRef> refs;
      for (size_t ib = 0; ib < binPool.size(); ++ib) {
        const MixEvent& buf = binPool[ib];
        for (size_t i = 0; i < A[c.ia].size(); ++i)
          for (size_t j = 0; j < buf.part[c.ib].size(); ++j) {
            MixRef r; r.ib = ib; r.reverse = 0; r.i = i; r.j = j; refs.push_back(r);
          }
        if (mix.mixBothDirections) {
          for (size_t i = 0; i < buf.a[c.ia].size(); ++i)
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
        const Cand& AA = r.reverse ? buf.a[c.ia][r.i] : A[c.ia][r.i];
        const Cand& BB = r.reverse ? B[c.ib][r.j] : buf.part[c.ib][r.j];
        if (femtoCfg.closePairVetoSameTrack && SharedTrack(AA, BB)) {
          c.nShared++;
          nShared++;
          hRejectShared->Fill(1);
          continue;
        }
        const Double_t ks = KStar(AA.p4, BB.p4);
        // The centrality of a mixed pair is the current event's, as in the maker: cent is a
        // property of the fill, and the mixing bin already ties the two events to the same class.
        if (c.resonance && c.mkkMEWide) c.mkkMEWide->Fill(AA.mKK, ks, centX);
        if (c.resonance && (AA.mKK < c.mLo || AA.mKK > c.mHi)) continue;
        // Mixed pairs get the identical cut. The B field is taken from the event the A candidate
        // came from; within a run it is constant, and the mixing bins do not span runs in practice.
        const Double_t bT = (r.reverse ? buf.bField : e.bField) / 10.0;
        if (ClosePairReject(AA, BB, BB.tEta, BB.tPhi, BB.tPt, BB.tCharge, bT, c.ib, femtoCfg)) {
          c.nCloseME++;
          c.closeME->Fill(ks);
          continue;
        }
        c.me->Fill(ks);
        c.meCent->Fill(ks, centX);
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
    {
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
  if (nTracksSeen != nTr) {
    std::cerr << "ERROR: read " << nTracksSeen << " of " << nTr << " track rows. The track tree is "
              << "not grouped by event in the event tree's order, which this reader requires.\n"
              << "       Refusing to report results from a subset of the tree." << std::endl;
    fout->Close();
    fin->Close();
    return;
  }
  std::cout << "[downstreamV3] events=" << nEv << " phi=" << nPhiTot << std::endl;
  for (Int_t ia = 0; ia < kNA; ++ia)
    std::cout << "[downstreamV3]   A " << kAName[ia] << ": n=" << nATot[ia] << std::endl;
  for (Int_t sp = 0; sp < kNPart; ++sp)
    std::cout << "[downstreamV3]   B " << kPartName[sp] << ": n=" << nPartTot[sp] << std::endl;
  for (size_t icha = 0; icha < chans.size(); ++icha) {
    const ChanHists& c = chans[icha];
    std::cout << "[downstreamV3]   " << c.name << ": SE=" << c.nSE << " ME=" << c.nME;
    if (femtoCfg.closePairEnabled)
      std::cout << " closeSE=" << c.nCloseSE << " closeME=" << c.nCloseME;
    std::cout << std::endl;
  }
  std::cout << "[downstreamV3] shared=" << nShared << " sameEventME=" << nSameEventME
            << " dcaKKReject=" << nDcaKKReject << " noOrigin=" << nNoOrigin << std::endl;

  fout->cd();
  TNamed("treeFile", treeFile).Write();
  TNamed("mixingMode", mode.Data()).Write();
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
         (recomputeDaughterPid || var.anyPid) ? "recomputed" : "kSelNominalPid")
      .Write();
  TNamed("mixBinSource", recomputeMixBin ? "recomputed" : "stored").Write();
  TNamed("maxDCAKK", TString::Format("%.3f%s", phiCfg.maxDCAKK,
                                     dcaKKActive ? " (active)" : " (off)")
                         .Data())
      .Write();
  TNamed("closePair",
         femtoCfg.closePairEnabled
             ? TString::Format("%s dEtaD=%.3f dPhiD=%.3f dEtaP=%.3f dPhiP=%.3f r=%.2f-%.2f/%.2f",
                               femtoCfg.closePairShape.c_str(), femtoCfg.closePairDEtaDeuteron,
                               femtoCfg.closePairDPhiStarDeuteron, femtoCfg.closePairDEtaProton,
                               femtoCfg.closePairDPhiStarProton, femtoCfg.closePairRadiusMin,
                               femtoCfg.closePairRadiusMax, femtoCfg.closePairRadiusStep)
                   .Data()
             : "disabled")
      .Write();
  TNamed("minNHitsDedxNuclear",
         TString::Format("%d", (var.dNHitsDedx > 0)
                                   ? var.dNHitsDedx
                                   : (Int_t)ConfigManager::GetInstance()
                                         .GetNuclearIdCuts()
                                         .minNHitsDedxNuclear)
             .Data())
      .Write();
  TNamed("nSameEventME", TString::Format("%lld", nSameEventME).Data()).Write();
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
    c.closeSE->Write();
    c.closeME->Write();
  }
  hRejectShared->Write(); hSameEventME->Write();
  fout->Close();
  fin->Close();
}
