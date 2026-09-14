#ifndef FEMTO_PHI_TREE_SCHEMA_H
#define FEMTO_PHI_TREE_SCHEMA_H

#include "Rtypes.h"
#include <cmath>

namespace femto_phi_tree {

// Schema 1 : Float_t track row, written up to 2026-09-13.
// Schema 2 : scaled-integer ("packed") track row. ROOT 5.34 ignores the leaflist
//            Float16/Double32 range syntax, so narrowing must be done with the C++
//            type itself. See implementation-plan-20260913.md §3.2 and §10.
const UInt_t kSchemaVersionV1 = 1;
const UInt_t kSchemaVersionV2 = 2;  // packed track row, Float_t event row (pilot only)
const UInt_t kSchemaVersionV3 = 3;  // packed event row as well
const UInt_t kSchemaVersion = kSchemaVersionV3;

enum SpeciesCode {
  kSpeciesUnknown = 0,
  kSpeciesKp = 1,
  kSpeciesKm = 2,
  kSpeciesDeuteron = 3,
  kSpeciesProton = 4,
  // Light nuclei beyond the deuteron. Their nuclear nSigma rides in the nSigmaDeuteron slot --
  // see NSigmaNuclear() -- because a row carries exactly one nuclear hypothesis and three more
  // Short fields would cost ~6 B on EVERY row, roughly +20% of the tree, to serve 3% of it.
  // For He3 and He4 the stored pT/eta/phi are the TRACK (rigidity) values; the physical momentum
  // is ChargeNumber x p, which the reader applies from speciesCode alone.
  kSpeciesTriton = 5,
  kSpeciesHe3 = 6,
  kSpeciesHe4 = 7
};

enum SelFlag {
  kSelTrackQualityNom = 1u << 0,
  kSelLoosePid = 1u << 1,
  kSelNominalPid = 1u << 2,
  kSelNominalFemto = 1u << 3,
  kSelTofMatch = 1u << 4,
  kSelKaonCutsNom = 1u << 5,
  // chi2 is the only nominal track cut whose input is not stored in the row, so the AND-ed
  // kSelTrackQualityNom bit cannot be decomposed: a row that fails it might have failed on
  // nHitsFit, on chi2, or on both, and no track-cut variation can be reproduced without knowing
  // which. This bit records the chi2 decision on its own. Every other nominal track cut
  // (nHitsFit, hit ratio, nHitsDedx, pT, eta, global DCA) is re-evaluable from stored fields,
  // and requirePrimaryTrack is already guaranteed by the storage envelope, so one bit is enough
  // to make the whole track-cut family reproducible -- at no cost, since selFlags is a UShort
  // with ten bits still free.
  kSelTrackChi2Nom = 1u << 6
};

enum EventFlag {
  kEvtPassEventCuts = 1u << 0,
  kEvtPassPileup = 1u << 1,
  kEvtPassCent = 1u << 2,
  kEvtPassMaxNTr = 1u << 3
};

enum PairFlag {
  kPairPassTof = 1u << 0,
  kPairPassKinematics = 1u << 1
};

// Schema 1 rows. Kept so that trees written before 2026-09-13 remain readable.
struct EventRow {
  UInt_t schemaVersion;
  ULong64_t eventUID;
  Int_t runId;
  Int_t eventId;
  UInt_t sourceFileHash;
  Long64_t sourceEntry;
  UInt_t subjobId;
  Float_t vx;
  Float_t vy;
  Float_t vz;
  Float_t vr;
  Float_t vzVpd;
  Float_t bField;
  Int_t refMult;
  Int_t rawMult;
  Int_t nBTOFMatch;
  Int_t nTracks;
  Float_t refMultCorr;
  Float_t centWeight;
  Float_t centralityPercent;
  Int_t cent9;
  Int_t cent16;
  Float_t qx;
  Float_t qy;
  Float_t psi2;
  Int_t mixVzBin;
  Int_t mixCentBin;
  Int_t mixEpBin;
  Int_t mixBin;
  UInt_t eventFlags;
  Int_t nKp;
  Int_t nKm;
  Int_t nDeuteron;
  Int_t nProton;

  EventRow() { Reset(); }

  void Reset() {
    schemaVersion = kSchemaVersionV1;
    eventUID = 0;
    runId = 0;
    eventId = 0;
    sourceFileHash = 0;
    sourceEntry = -1;
    subjobId = 0;
    vx = vy = vz = vr = vzVpd = bField = 0.0f;
    refMult = rawMult = nBTOFMatch = nTracks = 0;
    refMultCorr = centWeight = centralityPercent = 0.0f;
    cent9 = cent16 = -1;
    qx = qy = 0.0f;
    psi2 = -1.0f;
    mixVzBin = mixCentBin = mixEpBin = mixBin = 0;
    eventFlags = 0;
    nKp = nKm = nDeuteron = nProton = 0;
  }
};

struct TrackRow {
  ULong64_t eventUID;
  Int_t trackIndex;
  UChar_t speciesCode;
  Short_t charge;
  Float_t px;
  Float_t py;
  Float_t pz;
  Float_t dEdx;
  Float_t nSigmaKaon;
  Float_t nSigmaDeuteron;
  Float_t nSigmaProton;
  Char_t tofMatch;
  Float_t tofBeta;
  Float_t mass2;
  Float_t deltaOneOverBeta;
  Float_t dca;
  Short_t nHitsFit;
  Short_t nHitsMax;
  Short_t nHitsDedx;
  Float_t chi2;
  Float_t originX;
  Float_t originY;
  Float_t originZ;
  Float_t bField;
  UInt_t selFlags;

  TrackRow() { Reset(); }

  void Reset() {
    eventUID = 0;
    trackIndex = -1;
    speciesCode = kSpeciesUnknown;
    charge = 0;
    px = py = pz = 0.0f;
    dEdx = 0.0f;
    nSigmaKaon = nSigmaDeuteron = nSigmaProton = 0.0f;
    tofMatch = 0;
    tofBeta = -1.0f;
    mass2 = -999.0f;
    deltaOneOverBeta = 999.0f;
    dca = 0.0f;
    nHitsFit = nHitsMax = nHitsDedx = 0;
    chi2 = 0.0f;
    originX = originY = originZ = 0.0f;
    bField = 0.0f;
    selFlags = 0;
  }
};

struct PhiPairRow {
  ULong64_t eventUID;
  Int_t dauPlusIndex;
  Int_t dauMinusIndex;
  Float_t px;
  Float_t py;
  Float_t pz;
  Float_t mKK;
  Float_t dcaDaughters;
  Float_t openingAngle;
  Float_t rapidity;
  UInt_t pairFlags;

  PhiPairRow() { Reset(); }

  void Reset() {
    eventUID = 0;
    dauPlusIndex = dauMinusIndex = -1;
    px = py = pz = 0.0f;
    mKK = 0.0f;
    dcaDaughters = 0.0f;
    openingAngle = 0.0f;
    rapidity = 0.0f;
    pairFlags = 0;
  }
};


// ---------------------------------------------------------------------------
// Schema 2 packed track row
//
// Step 0 specification table. Every stored quantity, its C++ type, the integer
// scale factor, the representable range and the resulting resolution.
//
//   branch        type       scale   range                    resolution
//   eventUID      ULong64_t    -     (runId<<32)|eventId       exact
//   trackIndex    UShort_t     -     0 .. 65535                exact
//   speciesCode   UChar_t      -     0 .. 255                  exact
//   pT            UShort_t   6e3     0 .. 10.9225 GeV/c        1.67e-4 GeV/c
//   eta           Short_t    1e4     -3.2768 .. 3.2767         1e-4
//   phi           Short_t    1e4     -3.2768 .. 3.2767 rad     1e-4 rad
//   dEdx          UShort_t   600     0 .. 109.225 keV/cm       1.67e-3 keV/cm
//   nSigmaKaon    Short_t    100     -327.68 .. 327.66         0.01
//                                    32767 = "not computable" sentinel
//   nSigmaPion    Short_t    100     same                      0.01
//   nSigmaProton  Short_t    100     same                      0.01
//   nSigmaDeutron Short_t    100     same                      0.01
//   tofBeta       Short_t    8e3     -4.096 .. 4.0959          1.25e-4
//   dca           UShort_t   1e4     0 .. 6.5535 cm            1e-4 cm
//   nHitsFit      Char_t       -     -127 .. 127, sign = charge  exact
//   nHitsMax      UChar_t      -     0 .. 255                  exact
//   nHitsDedx     UChar_t      -     0 .. 255                  exact
//   selFlags      UShort_t     -     16 bits                   exact
//
// Dropped relative to schema 1, with the reason:
//   px, py, pz            -> replaced by pT/eta/phi (same information, half the bytes)
//   charge                -> sign of nHitsFit
//   tofMatch              -> kSelTofMatch bit of selFlags
//   mass2                 -> recomputable: p^2 (1/beta^2 - 1)
//   deltaOneOverBeta      -> recomputable: 1/beta - sqrt(m^2 + p^2)/p
//   chi2                  -> not used by any planned cut
//   originX/Y/Z, bField   -> only needed for the maxDCAKK cut, which is disabled
//                            (maxDCAKK: 200). See plan §10.5: if it is ever enabled
//                            these go into a companion tree, not the primary tree.
//   TPC hit topology      -> measured to carry no two-track information (plan §10).
// ---------------------------------------------------------------------------

namespace scale {
// Every range below was set from the measured extremes of the loose storage envelope,
// not from a guess. Step 2 asserts that no value is ever clipped.
const Double_t kPt = 6000.0;       // 0 .. 10.9225 GeV/c   (envMaxPt = 10.0)
const Double_t kEta = 10000.0;     // +-3.2767             (envMinEta = -2.0)
const Double_t kPhi = 10000.0;     // +-3.2767 rad
const Double_t kDedx = 250.0;      // 0 .. 262.1 keV/cm. Singly-charged tracks reach 92, but
                                  // Z=2 nuclei ionise ~4x and clipped at the old 600 (max 109.2).      // 0 .. 109.2 keV/cm    (observed max 91.98)
const Double_t kNSigma = 100.0;    // +-327.67             (PicoDst caps K/pi/p at +-32.767)
const Double_t kTofBeta = 8000.0;  // +-4.096              (PicoDst caps beta at 3.2768)
const Double_t kDca = 10000.0;     // 0 .. 6.5535 cm       (envMaxDca = 3.0)
}  // namespace scale

// NuclearIdDeDxVsMom::GetNSigma returns 999.0 to mean "not computable here" (momentum
// outside the parameterisation, or sigma <= 0). That is a sentinel, not a large nSigma,
// so it is preserved as its own packed code instead of being clipped to +327.67.
const Double_t kNSigmaSentinel = 999.0;
// A track that reached the BTOF cannot physically have beta this small at the envelope
// minimum momentum of 0.15 GeV/c; such values are timing failures. Recomputing m^2 from
// them produces numbers of order 1e8, so they are reported as "not measured" instead.
const Double_t kMinUsableBeta = 0.1;
const Double_t kNSigmaSentinelThreshold = 900.0;
const Short_t kNSigmaPackedInvalid = 32767;

// Saturation bookkeeping. Every packing call that clips increments a counter, so a
// production job can assert that no physics value was silently truncated.
struct PackStats {
  Long64_t pT, eta, phi, dEdx, tofBeta, dca, trackIndex, nHits;
  Long64_t nSigmaKaon, nSigmaPion, nSigmaProton, nSigmaDeuteron;
  Long64_t nSigmaSentinel;  // deliberate sentinel encodings, not a failure
  // The event row gets its own counters. Sharing a track counter hides which field clipped.
  Long64_t evVertex, evQ, evMult, evCent;
  Long64_t kaonOrigin;
  PackStats() { Reset(); }
  void Reset() {
    pT = eta = phi = dEdx = tofBeta = dca = trackIndex = nHits = 0;
    nSigmaKaon = nSigmaPion = nSigmaProton = nSigmaDeuteron = 0;
    nSigmaSentinel = 0;
    evVertex = evQ = evMult = evCent = 0;
    kaonOrigin = 0;
  }
  // Sentinel encodings are excluded: they are intended, not truncation.
  Long64_t Total() const {
    return pT + eta + phi + dEdx + tofBeta + dca + trackIndex + nHits + nSigmaKaon +
           nSigmaPion + nSigmaProton + nSigmaDeuteron + evVertex + evQ + evMult + evCent +
           kaonOrigin;
  }
};

inline UShort_t PackU16(Double_t v, Double_t s, Long64_t& sat) {
  Double_t r = (v >= 0) ? std::floor(v * s + 0.5) : std::ceil(v * s - 0.5);
  if (r > 65535.0) { r = 65535.0; ++sat; }
  if (r < 0.0) { r = 0.0; ++sat; }
  return (UShort_t)r;
}
inline Short_t PackI16(Double_t v, Double_t s, Long64_t& sat) {
  Double_t r = (v >= 0) ? std::floor(v * s + 0.5) : std::ceil(v * s - 0.5);
  if (r > 32767.0) { r = 32767.0; ++sat; }
  if (r < -32768.0) { r = -32768.0; ++sat; }
  return (Short_t)r;
}
inline Double_t UnpackU16(UShort_t v, Double_t s) { return (Double_t)v / s; }
inline Double_t UnpackI16(Short_t v, Double_t s) { return (Double_t)v / s; }

inline Short_t PackNSigma(Double_t v, Long64_t& sat, Long64_t& sentinel) {
  if (v >= kNSigmaSentinelThreshold || v <= -kNSigmaSentinelThreshold) {
    ++sentinel;
    return kNSigmaPackedInvalid;
  }
  return PackI16(v, scale::kNSigma, sat);
}
inline Double_t UnpackNSigma(Short_t v) {
  if (v == kNSigmaPackedInvalid) return kNSigmaSentinel;
  return UnpackI16(v, scale::kNSigma);
}
inline Bool_t NSigmaValid(Short_t v) { return v != kNSigmaPackedInvalid; }

struct TrackRowV2 {
  ULong64_t eventUID;
  UShort_t trackIndex;
  UChar_t speciesCode;
  UShort_t pT;
  Short_t eta;
  Short_t phi;
  UShort_t dEdx;
  Short_t nSigmaKaon;
  Short_t nSigmaPion;
  Short_t nSigmaProton;
  Short_t nSigmaDeuteron;
  Short_t tofBeta;
  UShort_t dca;
  Char_t nHitsFit;  // sign carries the charge
  UChar_t nHitsMax;
  UChar_t nHitsDedx;
  UShort_t selFlags;

  TrackRowV2() { Reset(); }

  void Reset() {
    eventUID = 0;
    trackIndex = 0;
    speciesCode = kSpeciesUnknown;
    pT = 0;
    eta = phi = 0;
    dEdx = 0;
    nSigmaKaon = nSigmaPion = nSigmaProton = nSigmaDeuteron = 0;
    tofBeta = 0;
    dca = 0;
    nHitsFit = 0;
    nHitsMax = nHitsDedx = 0;
    selFlags = 0;
  }

  // ---- decoded accessors, the only supported way to read a packed row ----
  Double_t Pt() const { return UnpackU16(pT, scale::kPt); }
  Double_t Eta() const { return UnpackI16(eta, scale::kEta); }
  Double_t Phi() const { return UnpackI16(phi, scale::kPhi); }
  Double_t Px() const { return Pt() * std::cos(Phi()); }
  Double_t Py() const { return Pt() * std::sin(Phi()); }
  Double_t Pz() const { return Pt() * std::sinh(Eta()); }
  Double_t P() const { return Pt() * std::cosh(Eta()); }
  Double_t Dedx() const { return UnpackU16(dEdx, scale::kDedx); }
  Double_t NSigmaKaon() const { return UnpackNSigma(nSigmaKaon); }
  Bool_t NSigmaKaonValid() const { return NSigmaValid(nSigmaKaon); }
  Double_t NSigmaPion() const { return UnpackNSigma(nSigmaPion); }
  Bool_t NSigmaPionValid() const { return NSigmaValid(nSigmaPion); }
  Double_t NSigmaProton() const { return UnpackNSigma(nSigmaProton); }
  Bool_t NSigmaProtonValid() const { return NSigmaValid(nSigmaProton); }
  Double_t NSigmaDeuteron() const { return UnpackNSigma(nSigmaDeuteron); }
  // The nuclear hypothesis this row was identified under: deuteron for speciesCode 3, triton for
  // 5, He3 for 6, He4 for 7. Prefer this to NSigmaDeuteron() on any nuclear row.
  Double_t NSigmaNuclear() const { return UnpackNSigma(nSigmaDeuteron); }
  // Charge number of the stored species. For He3/He4 the stored pT/eta/phi are the TRACK
  // (rigidity) values and the physical momentum is this factor times them.
  Int_t NuclearChargeNumber() const {
    return (speciesCode == kSpeciesHe3 || speciesCode == kSpeciesHe4) ? 2 : 1;
  }
  Bool_t NSigmaDeuteronValid() const { return NSigmaValid(nSigmaDeuteron); }
  Double_t TofBeta() const { return UnpackI16(tofBeta, scale::kTofBeta); }
  Double_t Dca() const { return UnpackU16(dca, scale::kDca); }
  Int_t NHitsFit() const { return (nHitsFit >= 0) ? nHitsFit : -nHitsFit; }
  Int_t Charge() const { return (nHitsFit >= 0) ? 1 : -1; }
  Bool_t HasTof() const { return (selFlags & kSelTofMatch) != 0; }

  // Derived quantities that schema 1 stored explicitly.
  //
  // Precision note. m^2 = p^2 (1/beta^2 - 1), so d(m^2)/d(beta) = -2 p^2 / beta^3: the
  // 1.25e-4 quantization step of tofBeta is amplified as 1/beta^3. Measured against the
  // schema 1 stored value on 2.97e6 TOF tracks:
  //   whole sample                       RMS 6.1e-4 GeV^2/c^4
  //   -0.2 < m^2 < 6 (99.3% of tracks,   RMS 4.5e-4, i.e. 1/442 of the narrowest planned
  //   where every planned PID window lives)      PID window, kaon m^2 0.16-0.36
  // Outside that region the amplification grows; those are slow or mistimed tracks that no
  // planned selection uses. Below kMinUsableBeta the result is not a measurement at all, so
  // the sentinel is returned rather than a number that looks like one.
  Double_t Mass2() const {
    if (!HasTof()) return -999.0;
    Double_t b = TofBeta();
    if (b <= kMinUsableBeta) return -999.0;
    Double_t p = P();
    return p * p * (1.0 / (b * b) - 1.0);
  }
  Double_t DeltaOneOverBeta(Double_t mass) const {
    if (!HasTof()) return 999.0;
    Double_t b = TofBeta();
    if (b <= kMinUsableBeta) return 999.0;
    Double_t p = P();
    if (p <= 0) return 999.0;
    return 1.0 / b - std::sqrt(mass * mass + p * p) / p;
  }
};

struct EventRowV2 {
  UInt_t schemaVersion;
  ULong64_t eventUID;
  Int_t runId;
  Int_t eventId;
  UInt_t sourceFileHash;    // 32-bit FNV-1a, kept for backward comparison only
  UShort_t sourceFileIndex; // index into the per-file path table; the invertible key
  Long64_t sourceEntry;
  UInt_t subjobId;
  UInt_t triggerId;    // first configured trigger this event fired, 0 if none
  UShort_t triggerBits;  // bitmask over the configured trigger list
  Float_t vx, vy, vz, vr, vzVpd;
  Float_t bField;
  Int_t refMult, rawMult, nBTOFMatch, nTracks;
  Float_t refMultCorr, centWeight, centralityPercent;
  Int_t cent9, cent16;
  Float_t qx, qy, psi2;
  Int_t mixVzBin, mixCentBin, mixEpBin, mixBin;
  UInt_t eventFlags;
  Int_t nKp, nKm, nDeuteron, nProton;

  EventRowV2() { Reset(); }

  void Reset() {
    schemaVersion = kSchemaVersionV2;
    eventUID = 0;
    runId = eventId = 0;
    sourceFileHash = 0;
    sourceFileIndex = 0;
    sourceEntry = -1;
    subjobId = 0;
    triggerId = 0;
    triggerBits = 0;
    vx = vy = vz = vr = vzVpd = bField = 0.0f;
    refMult = rawMult = nBTOFMatch = nTracks = 0;
    refMultCorr = centWeight = centralityPercent = 0.0f;
    cent9 = cent16 = -1;
    qx = qy = 0.0f;
    psi2 = -1.0f;
    mixVzBin = mixCentBin = mixEpBin = mixBin = 0;
    eventFlags = 0;
    nKp = nKm = nDeuteron = nProton = 0;
  }
};


// ---------------------------------------------------------------------------
// Schema 3 packed event row
//
// Schema 2 packed the track row but left the event row as Float_t; measured at
// 49.33 B/event against 20.25 B/event packed, i.e. 66 GB over the full dataset.
//
//   branch             type       scale  range                  resolution
//   eventUID           ULong64_t    -    (runId<<32)|eventId     exact
//   sourceFileHash     UInt_t       -    FNV-1a of the path      exact
//   sourceFileIndex    UShort_t     -    row of SourceFileTable  exact
//   sourceEntry        Int_t        -    entry in that file      exact
//   subjobId           UInt_t       -    FNV-1a of the jobid     exact
//   triggerId          UInt_t       -    first configured match  exact
//   triggerBits        UShort_t     -    mask over the list      exact
//   vx, vy             Short_t    1e3    +-32.767 cm             10 um
//   vz                 Short_t    1e2    +-327.67 cm             100 um
//   vzVpd              Short_t    1e1    +-3276.7 cm             1 mm
//                                  the FXT VPD is not usable: measured range -1725..+1289 cm,
//                                  and the event cut ignores it unless |vzVpd| < maxAbsVzVpd
//   vr                 UShort_t   1e3    0 .. 65.535 cm          10 um
//   bField             Float_t      -    kilogauss               exact
//   refMult ..nTracks  UShort_t     -    0 .. 65535              exact
//   refMultCorr        UShort_t   1e2    0 .. 655.35             0.01
//   centWeight         UShort_t   1e4    0 .. 6.5535             1e-4
//   centralityPercent  UShort_t   1e2    0 .. 655.35 %           0.01 %
//   cent9, cent16      Char_t       -    -1 .. 15                exact
//   qx, qy             Short_t    1e2    +-327.67                0.01
//   psi2               Short_t    1e4    +-3.2767 rad            1e-4
//   mixBin             UShort_t     -    0 .. 65535              exact
//   eventFlags         UChar_t      -    8 bits                  exact
//   nKp..nProton       UShort_t     -    0 .. 65535              exact
//
// Dropped relative to schema 2, with the reason:
//   runId, eventId   -> both recoverable from eventUID
//   mixVzBin, mixCentBin, mixEpBin
//                    -> derived from vz / cent9 / psi2 and the mixing YAML, and each is a
//                       factor of the composite mixBin that is stored.
// mixBin itself is kept (2 B/event, 0.24% of the row) after the Step 4 closure showed that
// recomputing it from the packed vz reassigns 1.4% of events to a neighbouring bin: the vz
// packing step is 0.01 cm and the bins are 0.4 cm wide. Keeping it does not freeze the
// binning -- the downstream still recomputes on request for a mixing-binning systematic --
// it only makes the nominal reproduce the maker exactly.
// ---------------------------------------------------------------------------

namespace escale {
const Double_t kVxy = 1000.0;
const Double_t kVz = 100.0;
const Double_t kVzVpd = 10.0;
const Double_t kVr = 1000.0;
const Double_t kRefMultCorr = 100.0;
const Double_t kCentWeight = 10000.0;
const Double_t kCentPercent = 100.0;
const Double_t kQ = 100.0;
const Double_t kPsi2 = 10000.0;
// Helix origin of a phi-daughter kaon, stored as (origin - primary vertex). Measured worst
// component on real data is 1.996 cm, so Short x 1e4 (+-3.2767 cm, 1e-4 cm steps) is ample.
const Double_t kOrigin = 10000.0;
}  // namespace escale

// ---------------------------------------------------------------------------
// Companion row: the helix origin of a phi-daughter kaon.
//
// The KK decay DCA needs each daughter's helix -- origin, momentum, charge, B field. Momentum and
// charge come from the packed track row and the B field from the event row, so the origin is the
// only missing piece, and it is only ever needed for kaons. Carrying it on every row would cost
// +32% of the tree; carrying it for kaons alone costs +0.2%, because kaons are 0.66% of the
// stored rows once the daughter-PID reach gate is applied.
//
// The origin is stored RELATIVE TO THE PRIMARY VERTEX. That is not only a packing trick: two
// helices translated by the same vector have the same distance of closest approach, so a reader
// can build both helices in the vertex frame and never add the vertex back. The 0.01 cm
// quantisation of the stored vz therefore cannot reach the KK DCA at all.
// ---------------------------------------------------------------------------
struct KaonOriginRow {
  ULong64_t eventUID;
  UShort_t trackIndex;
  Short_t originDx, originDy, originDz;

  KaonOriginRow() { Reset(); }
  void Reset() {
    eventUID = 0;
    trackIndex = 0;
    originDx = originDy = originDz = 0;
  }
  Double_t Dx() const { return originDx / escale::kOrigin; }
  Double_t Dy() const { return originDy / escale::kOrigin; }
  Double_t Dz() const { return originDz / escale::kOrigin; }
};

struct EventRowV3 {
  UInt_t schemaVersion;
  ULong64_t eventUID;
  UInt_t sourceFileHash;
  UShort_t sourceFileIndex;
  Int_t sourceEntry;
  UInt_t subjobId;
  UInt_t triggerId;
  UShort_t triggerBits;
  Short_t vx, vy, vz, vzVpd;
  UShort_t vr;
  Float_t bField;
  UShort_t refMult, rawMult, nBTOFMatch, nTracks;
  UShort_t refMultCorr, centWeight, centralityPercent;
  Char_t cent9, cent16;
  Short_t qx, qy, psi2;
  // The mixing-bin index the maker assigned, stored rather than recomputed. vz is packed at
  // 0.01 cm while the mixing bins are 0.4 cm wide, so recomputing moves 1.4% of events to a
  // neighbouring bin (measured on the 5-file pilot) and the mixed-event pair set stops matching
  // the maker-direct chain exactly. Recomputation from vz/cent9/psi2 stays available and is the
  // right thing to do for a mixing-binning systematic, where a 0.01 cm boundary reshuffle is
  // far below the size of the variation itself.
  UShort_t mixBin;
  UChar_t eventFlags;
  UShort_t nKp, nKm, nDeuteron, nProton;

  EventRowV3() { Reset(); }

  void Reset() {
    schemaVersion = kSchemaVersionV3;
    eventUID = 0;
    sourceFileHash = 0;
    sourceFileIndex = 0;
    sourceEntry = -1;
    subjobId = 0;
    triggerId = 0;
    triggerBits = 0;
    vx = vy = vz = vzVpd = 0;
    vr = 0;
    bField = 0.0f;
    refMult = rawMult = nBTOFMatch = nTracks = 0;
    refMultCorr = centWeight = centralityPercent = 0;
    cent9 = cent16 = -1;
    qx = qy = 0;
    psi2 = 0;
    mixBin = 0;
    eventFlags = 0;
    nKp = nKm = nDeuteron = nProton = 0;
  }

  Int_t RunId() const { return (Int_t)(UInt_t)(eventUID >> 32); }
  Int_t EventId() const { return (Int_t)(UInt_t)(eventUID & 0xFFFFFFFFull); }
  Double_t Vx() const { return UnpackI16(vx, escale::kVxy); }
  Double_t Vy() const { return UnpackI16(vy, escale::kVxy); }
  Double_t Vz() const { return UnpackI16(vz, escale::kVz); }
  Double_t VzVpd() const { return UnpackI16(vzVpd, escale::kVzVpd); }
  Double_t Vr() const { return UnpackU16(vr, escale::kVr); }
  Double_t RefMultCorr() const { return UnpackU16(refMultCorr, escale::kRefMultCorr); }
  Double_t CentWeight() const { return UnpackU16(centWeight, escale::kCentWeight); }
  Double_t CentralityPercent() const { return UnpackU16(centralityPercent, escale::kCentPercent); }
  Double_t Qx() const { return UnpackI16(qx, escale::kQ); }
  Double_t Qy() const { return UnpackI16(qy, escale::kQ); }
  Double_t Psi2() const { return UnpackI16(psi2, escale::kPsi2); }
};

inline ULong64_t MakeEventUID(Int_t runId, Int_t eventId) {
  return (static_cast<ULong64_t>(static_cast<UInt_t>(runId)) << 32) |
         static_cast<ULong64_t>(static_cast<UInt_t>(eventId));
}

inline UInt_t HashStringFnv(const char* s) {
  UInt_t h = 2166136261u;
  if (!s) return 0;
  while (*s) {
    h ^= static_cast<UInt_t>(static_cast<unsigned char>(*s++));
    h *= 16777619u;
  }
  return h;
}

}  // namespace femto_phi_tree

#endif
