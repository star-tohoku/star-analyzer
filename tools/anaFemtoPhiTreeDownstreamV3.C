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

#include <deque>
#include <iostream>
#include <map>
#include <vector>

namespace {
const Double_t kKaonMass = 0.493677;
const Double_t kDeuteronMass = 1.875612;
const Double_t kProtonMass = 0.938272;

// Bachelor species paired against the phi. Index order is fixed and used for the histogram
// and counter arrays below.
enum { kPartDeuteron = 0, kPartProton = 1, kNPart = 2 };
const Char_t* const kPartName[kNPart] = {"deuteron", "proton"};
// A physical KK decay DCA is a fraction of a cm; anything at or above this means the cut was
// deliberately parked in the "off" position.
const Double_t kDcaKKDisabledAbove = 50.0;

struct Cand {
  Int_t trackIndex;
  TLorentzVector p4;
  Float_t mKK;
  Int_t dau1;
  Int_t dau2;
  // daughter kinematics, kept so the close-pair cut can work on the tracks that actually
  // interfere in the TPC rather than on the reconstructed phi
  Double_t dEta[2], dPhi[2], dPt[2];
  Short_t dCharge[2];
  // single-track kinematics, for a bachelor
  Double_t tEta, tPhi, tPt;
  Short_t tCharge;
};

struct MixEvent {
  ULong64_t eventUID;
  Double_t bField;
  std::vector<Cand> phis;
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
// ---------------------------------------------------------------------------
struct Var {
  Double_t nHitsFit, trkDca, trkPt, chi2;
  Double_t kNSigma, kM2Lo, kM2Hi;
  Int_t kmRequireTof;
  Double_t dNSigma, dDca;
  Int_t dNHitsDedx;
  Double_t pNSigma, pDca;
  Bool_t anyTrack, anyKaon, anyPid;
  TString spec;
  Var()
      : nHitsFit(-1), trkDca(-1), trkPt(-1), chi2(-1), kNSigma(-1), kM2Lo(-1), kM2Hi(-1),
        kmRequireTof(-1), dNSigma(-1), dDca(-1), dNHitsDedx(-1), pNSigma(-1), pDca(-1),
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

Bool_t PassDeuteronVariation(const femto_phi_tree::TrackRowV2& t, const Var& v) {
  if (t.speciesCode != femto_phi_tree::kSpeciesDeuteron) return kFALSE;
  if (!PassNuclearDedxHits(t, v.dNHitsDedx)) return kFALSE;
  if (!(t.selFlags & femto_phi_tree::kSelNominalPid)) return kFALSE;
  if (!PassTrackQuality(t, v)) return kFALSE;
  if (v.dNSigma < 0 && v.dDca < 0) {
    return (t.selFlags & femto_phi_tree::kSelNominalFemto) != 0;
  }
  const FemtoConfig& fc = ConfigManager::GetInstance().GetFemtoConfig();
  // Tightening can ride on top of kSelNominalFemto; loosening cannot, because the flag has
  // already applied the nominal value, so the whole femto cut has to be rebuilt.
  const Bool_t tightenOnly = (v.dNSigma < 0 || v.dNSigma <= fc.deuteronMaxAbsNSigma) &&
                             (v.dDca < 0 || v.dDca <= fc.deuteronMaxDca);
  if (tightenOnly && !(t.selFlags & femto_phi_tree::kSelNominalFemto)) return kFALSE;
  const Double_t ns = (v.dNSigma > 0) ? v.dNSigma : fc.deuteronMaxAbsNSigma;
  const Double_t dc = (v.dDca > 0) ? v.dDca : fc.deuteronMaxDca;
  if (TMath::Abs(t.NSigmaDeuteron()) >= ns) return kFALSE;
  if (t.Dca() >= dc) return kFALSE;
  return kTRUE;
}

Bool_t PassProtonVariation(const femto_phi_tree::TrackRowV2& t, const Var& v) {
  if (t.speciesCode != femto_phi_tree::kSpeciesProton) return kFALSE;
  if (!(t.selFlags & femto_phi_tree::kSelNominalPid)) return kFALSE;
  if (!PassTrackQuality(t, v)) return kFALSE;
  if (v.pNSigma < 0 && v.pDca < 0) {
    return (t.selFlags & femto_phi_tree::kSelNominalFemto) != 0;
  }
  const FemtoConfig& fc = ConfigManager::GetInstance().GetFemtoConfig();
  const Bool_t tightenOnly = (v.pNSigma < 0 || v.pNSigma <= fc.protonMaxAbsNSigma) &&
                             (v.pDca < 0 || v.pDca <= fc.protonMaxDca);
  if (tightenOnly && !(t.selFlags & femto_phi_tree::kSelNominalFemto)) return kFALSE;
  const Double_t ns = (v.pNSigma > 0) ? v.pNSigma : fc.protonMaxAbsNSigma;
  const Double_t dc = (v.pDca > 0) ? v.pDca : fc.protonMaxDca;
  if (TMath::Abs(t.NSigmaProton()) >= ns) return kFALSE;
  if (t.Dca() >= dc) return kFALSE;
  return kTRUE;
}

// StFemtoMaker reaches the proton collection through PassProtonCuts, which begins with the same
// PassTrackCuts gate the deuteron branch inherits, then |nSigmaProton| <= pid.nSigmaProton and
// IsProton (= PassTofProtonPid); the charge mode and the pT / rapidity windows live inside
// PassFemtoProtonCuts. The tree maker records exactly that split across the three flags below, so
// the nominal proton is the same triple as the nominal deuteron. There is no nuclear-ID dE/dx hit
// requirement on this path.
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
  const Double_t wEta =
      (species == 0) ? fc.closePairDEtaDeuteron : fc.closePairDEtaProton;
  const Double_t wPhi =
      (species == 0) ? fc.closePairDPhiStarDeuteron : fc.closePairDPhiStarProton;
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

Bool_t SharedTrack(const Cand& phi, const Cand& d) {
  return d.trackIndex >= 0 && (d.trackIndex == phi.dau1 || d.trackIndex == phi.dau2);
}
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
  TH1D* hKstarSE[kNPart];
  TH1D* hKstarME[kNPart];
  TH1D* hNPart[kNPart];
  for (Int_t sp = 0; sp < kNPart; ++sp) {
    hKstarSE[sp] = new TH1D(TString::Format("hKstarSE_phi_%s_signal", kPartName[sp]),
                            "SE k*;k* (GeV/c);counts", 50, 0.0, 1.0);
    hKstarME[sp] = new TH1D(TString::Format("hKstarME_phi_%s_signal", kPartName[sp]),
                            "ME k*;k* (GeV/c);counts", 50, 0.0, 1.0);
    hNPart[sp] = new TH1D(TString::Format("hN%s", kPartName[sp]),
                          TString::Format("nominal %s / event", kPartName[sp]), 21, -0.5, 20.5);
  }
  TH1D* hNPhi = new TH1D("hNPhi", "nominal #phi / event", 21, -0.5, 20.5);
  // Removal rate of the close-pair cut, binned in k*: the quantity that says whether the cut
  // sculpts the correlation function or merely trims it.
  TH1D* hCloseSE[kNPart];
  TH1D* hCloseME[kNPart];
  for (Int_t sp = 0; sp < kNPart; ++sp) {
    hCloseSE[sp] = new TH1D(TString::Format("hClosePairRejectSE_phi_%s", kPartName[sp]),
                            "SE pairs removed by the close-pair cut;k* (GeV/c);pairs", 50, 0.0,
                            1.0);
    hCloseME[sp] = new TH1D(TString::Format("hClosePairRejectME_phi_%s", kPartName[sp]),
                            "ME pairs removed by the close-pair cut;k* (GeV/c);pairs", 50, 0.0,
                            1.0);
  }
  Long64_t nCloseSE[kNPart] = {0, 0}, nCloseME[kNPart] = {0, 0};
  TH1D* hRejectShared = new TH1D("hRejectShared", "shared-track rejects", 2, -0.5, 1.5);
  TH1D* hSameEventME = new TH1D("hSameEventME", "ME pairs from same eventUID (must be 0)", 2, -0.5, 1.5);

  std::map<Int_t, std::deque<MixEvent> > pool;
  Long64_t nSE[kNPart] = {0, 0}, nME[kNPart] = {0, 0}, nPartTot[kNPart] = {0, 0};
  Long64_t nShared = 0, nPhiTot = 0, nSameEventME = 0;
  const Long64_t nEv = ev->GetEntries();
  const Int_t maxMixed = mix.maxMixedPairsPerEvent;
  TRandom3 rng(1);

  for (Long64_t ie = 0; ie < nEv; ++ie) {
    ev->GetEntry(ie);
    std::vector<Cand> kp, km, phis;
    std::vector<Cand> part[kNPart];
    std::map<Int_t, TrackRowV2> kmap;
    const std::vector<Long64_t>& idxs = tracksByEvent[e.eventUID];
    for (size_t k = 0; k < idxs.size(); ++k) {
      tr->GetEntry(idxs[k]);
      if (t.speciesCode == kSpeciesKp || t.speciesCode == kSpeciesKm) {
        if (!PassKaonForPhi(t, var)) continue;
        kmap[t.trackIndex] = t;
        Cand c;
        c.trackIndex = t.trackIndex;
        c.mKK = 0;
        c.dau1 = c.dau2 = -1;
        TVector3 p = Momentum(t);
        c.p4 = TLorentzVector(p, TMath::Sqrt(kKaonMass * kKaonMass + p.Mag2()));
        if (t.speciesCode == kSpeciesKp) kp.push_back(c); else km.push_back(c);
      } else if (t.speciesCode == kSpeciesDeuteron) {
        if (!PassDeuteronVariation(t, var)) continue;
        Cand c;
        c.trackIndex = t.trackIndex;
        c.mKK = 0;
        c.dau1 = c.dau2 = -1;
        TVector3 p = Momentum(t);
        c.p4 = TLorentzVector(p, TMath::Sqrt(kDeuteronMass * kDeuteronMass + p.Mag2()));
        c.tEta = t.Eta(); c.tPhi = t.Phi(); c.tPt = t.Pt(); c.tCharge = (Short_t)t.Charge();
        part[kPartDeuteron].push_back(c);
      } else if (t.speciesCode == kSpeciesProton) {
        if (!PassProtonVariation(t, var)) continue;
        Cand c;
        c.trackIndex = t.trackIndex;
        c.mKK = 0;
        c.dau1 = c.dau2 = -1;
        TVector3 p = Momentum(t);
        c.p4 = TLorentzVector(p, TMath::Sqrt(kProtonMass * kProtonMass + p.Mag2()));
        c.tEta = t.Eta(); c.tPhi = t.Phi(); c.tPt = t.Pt(); c.tCharge = (Short_t)t.Charge();
        part[kPartProton].push_back(c);
      }
    }

    for (size_t i = 0; i < kp.size(); ++i) {
      for (size_t j = 0; j < km.size(); ++j) {
        const TrackRowV2& a = kmap[kp[i].trackIndex];
        const TrackRowV2& b = kmap[km[j].trackIndex];
        if (!PassDaughterPid(a, var, recomputeDaughterPid) ||
            !PassDaughterPid(b, var, recomputeDaughterPid))
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
        phi.dEta[0] = a.Eta(); phi.dPhi[0] = a.Phi(); phi.dPt[0] = a.Pt();
        phi.dCharge[0] = (Short_t)a.Charge();
        phi.dEta[1] = b.Eta(); phi.dPhi[1] = b.Phi(); phi.dPt[1] = b.Pt();
        phi.dCharge[1] = (Short_t)b.Charge();
        phis.push_back(phi);
        nPhiTot++;
      }
    }
    hNPhi->Fill(phis.size());
    for (Int_t sp = 0; sp < kNPart; ++sp) {
      nPartTot[sp] += part[sp].size();
      hNPart[sp]->Fill(part[sp].size());
    }

    for (Int_t sp = 0; sp < kNPart; ++sp) {
      for (size_t i = 0; i < phis.size(); ++i) {
        if (phis[i].mKK < sigMin || phis[i].mKK > sigMax) continue;
        for (size_t j = 0; j < part[sp].size(); ++j) {
          if (femtoCfg.closePairVetoSameTrack && SharedTrack(phis[i], part[sp][j])) {
            nShared++;
            hRejectShared->Fill(1);
            continue;
          }
          const Double_t ks = KStar(phis[i].p4, part[sp][j].p4);
          if (ClosePairReject(phis[i], part[sp][j], part[sp][j].tEta, part[sp][j].tPhi,
                              part[sp][j].tPt, part[sp][j].tCharge, e.bField / 10.0, sp,
                              femtoCfg)) {
            nCloseSE[sp]++;
            hCloseSE[sp]->Fill(ks);
            continue;
          }
          hKstarSE[sp]->Fill(ks);
          nSE[sp]++;
        }
      }
    }

    const Int_t mixBin =
        recomputeMixBin ? MixBinOf(e.Vz(), (Int_t)e.cent9, e.Psi2()) : (Int_t)e.mixBin;
    std::deque<MixEvent>& binPool = pool[mixBin];
    // The bachelor species share the event pool but never share a pair, so the sampling is done
    // per species: with maxMixedPairsPerEvent set, each species gets its own cap, exactly as the
    // maker draws separately per channel.
    for (Int_t sp = 0; sp < kNPart; ++sp) {
      std::vector<MixRef> refs;
      for (size_t ib = 0; ib < binPool.size(); ++ib) {
        const MixEvent& buf = binPool[ib];
        for (size_t i = 0; i < phis.size(); ++i) {
          if (phis[i].mKK < sigMin || phis[i].mKK > sigMax) continue;
          for (size_t j = 0; j < buf.part[sp].size(); ++j) {
            MixRef r; r.ib = ib; r.reverse = 0; r.i = i; r.j = j; refs.push_back(r);
          }
        }
        if (mix.mixBothDirections) {
          for (size_t i = 0; i < buf.phis.size(); ++i) {
            if (buf.phis[i].mKK < sigMin || buf.phis[i].mKK > sigMax) continue;
            for (size_t j = 0; j < part[sp].size(); ++j) {
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
        const Cand& B = r.reverse ? part[sp][r.j] : buf.part[sp][r.j];
        const Double_t ks = KStar(A.p4, B.p4);
        // Mixed pairs get the identical cut. The B field is taken from the event the phi came
        // from; within a run it is constant, and the mixing bins do not span runs in practice.
        const Double_t bT = (r.reverse ? buf.bField : e.bField) / 10.0;
        if (ClosePairReject(A, B, B.tEta, B.tPhi, B.tPt, B.tCharge, bT, sp, femtoCfg)) {
          nCloseME[sp]++;
          hCloseME[sp]->Fill(ks);
          continue;
        }
        hKstarME[sp]->Fill(ks);
        nME[sp]++;
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
      me.phis = phis;
      for (Int_t sp = 0; sp < kNPart; ++sp) me.part[sp] = part[sp];
      binPool.push_back(me);
      if ((Int_t)binPool.size() > bufferSize) binPool.pop_front();
    }
  }

  std::cout << "[downstreamV3] events=" << nEv << " phi=" << nPhiTot << std::endl;
  for (Int_t sp = 0; sp < kNPart; ++sp) {
    std::cout << "[downstreamV3]   phi-" << kPartName[sp] << ": n=" << nPartTot[sp]
              << " SE=" << nSE[sp] << " ME=" << nME[sp] << std::endl;
  }
  std::cout << "[downstreamV3] shared=" << nShared << " sameEventME=" << nSameEventME << std::endl;
  if (femtoCfg.closePairEnabled) {
    for (Int_t sp = 0; sp < kNPart; ++sp) {
      const Double_t fse = (nSE[sp] + nCloseSE[sp]) > 0
                               ? 100.0 * nCloseSE[sp] / (nSE[sp] + nCloseSE[sp]) : 0.0;
      const Double_t fme = (nME[sp] + nCloseME[sp]) > 0
                               ? 100.0 * nCloseME[sp] / (nME[sp] + nCloseME[sp]) : 0.0;
      std::cout << "[downstreamV3]   close-pair phi-" << kPartName[sp] << ": SE removed "
                << nCloseSE[sp] << " (" << fse << "%)  ME removed " << nCloseME[sp] << " ("
                << fme << "%)" << std::endl;
    }
  }

  fout->cd();
  TNamed("treeFile", treeFile).Write();
  TNamed("mixingMode", mode.Data()).Write();
  TNamed("schemaVersion", TString::Format("%u", e.schemaVersion).Data()).Write();
  TNamed("nSE", TString::Format("%lld", nSE[kPartDeuteron]).Data()).Write();
  TNamed("nME", TString::Format("%lld", nME[kPartDeuteron]).Data()).Write();
  TNamed("nSEProton", TString::Format("%lld", nSE[kPartProton]).Data()).Write();
  TNamed("nMEProton", TString::Format("%lld", nME[kPartProton]).Data()).Write();
  TNamed("nP", TString::Format("%lld", nPartTot[kPartProton]).Data()).Write();
  TNamed("nPhi", TString::Format("%lld", nPhiTot).Data()).Write();
  TNamed("nD", TString::Format("%lld", nPartTot[kPartDeuteron]).Data()).Write();
  TNamed("variation", var.spec.Data()).Write();
  TNamed("daughterPid",
         (recomputeDaughterPid || var.anyPid) ? "recomputed" : "kSelNominalPid")
      .Write();
  TNamed("mixBinSource", recomputeMixBin ? "recomputed" : "stored").Write();
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
  hNPhi->Write();
  for (Int_t sp = 0; sp < kNPart; ++sp) {
    hKstarSE[sp]->Write();
    hKstarME[sp]->Write();
    hNPart[sp]->Write();
    hCloseSE[sp]->Write();
    hCloseME[sp]->Write();
  }
  hRejectShared->Write(); hSameEventME->Write();
  fout->Close();
  fin->Close();
}
