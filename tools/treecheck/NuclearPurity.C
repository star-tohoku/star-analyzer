// How pure is a TPC-only selection for each bachelor species?
//
// Generalises ProtonPurity.C to deuteron, triton, 3He and 4He. The question is the same: the
// production analysis requires a TOF match on every bachelor (all the TofMomentumThreshold keys
// are 0.0), and dropping that requirement roughly doubles the sample. Whether that is a good
// trade depends on how pure the dE/dx selection is on its own, which the TOF-matched subsample
// can measure: take every track the analysis would accept WITHOUT the TOF rule, keep those that
// happen to carry a TOF match, and see where their m2 lands.
//
// Two numbers are reported per momentum slice, as for the proton:
//   purity(win)  fraction inside the species m2 window. Simple, but it calls the peak's own tails
//                contamination, so it is a lower bound.
//   purity(sb)   every track minus what lands outside the window, minus a flat mismatch background
//                taken from a sideband where no species peak sits. No peak shape is assumed; a fit
//                was tried for the proton and described the data far worse (chi2/ndf 30-116).
//
// Rigidity. The stored pT/eta/phi are track quantities, so for 3He and 4He (Z = 2) they are the
// rigidity and the physical momentum is twice them. TrackRowV2::Mass2() uses the stored momentum,
// so what it returns for a Z = 2 nucleus is m2/Z^2 -- and that is exactly the convention the YAML
// windows are written in (he3MinMass2 0.8 / MaxMass2 2.8 brackets 7.887/4 = 1.97). So Mass2() is
// compared against the configured window directly, as the maker does, and the momentum axis is
// rigidity for every species.
//
// Cut values are those of maker_auau3p85fxt_anaFemtoPhiTree_prod.yaml, printed at run time so a
// drift from the YAML is visible rather than silent.
//
// Usage: NuclearPurity.C("<merged block>.root", "deuteron", stride, "<out.root>")
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TString.h"
#include <cstdio>
#include "FemtoPhiTreeSchema.h"

namespace {
struct SpeciesCut {
  const char* name;
  UChar_t code;
  Double_t m2Peak;      // rigidity-based, i.e. what Mass2() returns for a real one
  Double_t m2Lo, m2Hi;  // the YAML window
  Double_t sbLo, sbHi;  // flat-background sideband: no species peak may sit here
  Double_t m2HistLo, m2HistHi;
  Double_t nSigmaMax, dcaMax, ptMin, ptMax, etaMax, pMin, pMax;
  Int_t nHitsFitMin, nHitsDedxMin;
  Double_t ratioMin;
};

// The sidebands are chosen to avoid every other species' peak in rigidity m2:
// proton 0.88, He3 1.97, deuteron 3.52, He4 3.47, triton 7.89.
const SpeciesCut kSpecies[5] = {
  {"proton",   femto_phi_tree::kSpeciesProton,   0.880, 0.6, 1.2,  1.4,  2.2, -0.5,  2.5,
   2.0, 1.0, 0.4, 2.0, 2.0, 0.0, 1e9, 15,  0, 0.52},
  {"deuteron", femto_phi_tree::kSpeciesDeuteron, 3.518, 2.5, 4.5,  5.0,  6.5, -0.5,  9.0,
   2.0, 1.0, 0.4, 3.0, 2.0, 0.2, 3.0, 15, 10, 0.52},
  {"triton",   femto_phi_tree::kSpeciesTriton,   7.890, 5.8, 9.8, 11.0, 14.0, -0.5, 16.0,
   2.0, 1.0, 0.4, 3.0, 2.0, 0.2, 3.0, 15, 10, 0.52},
  {"he3",      femto_phi_tree::kSpeciesHe3,      1.972, 0.8, 2.8,  4.5,  6.0, -0.5,  9.0,
   2.0, 1.0, 0.4, 3.0, 2.0, 0.2, 3.0, 15, 10, 0.52},
  {"he4",      femto_phi_tree::kSpeciesHe4,      3.473, 2.5, 4.5,  5.5,  7.0, -0.5,  9.0,
   2.0, 1.0, 0.4, 3.0, 2.0, 0.2, 3.0, 15, 10, 0.52}
};
}  // namespace

void NuclearPurity(const char* treeFile, const char* species = "deuteron", Long64_t stride = 20,
                   const char* outFile = "") {
  const SpeciesCut* sc = 0;
  for (Int_t i = 0; i < 5; ++i) if (TString(species) == kSpecies[i].name) sc = &kSpecies[i];
  if (!sc) { printf("unknown species '%s'\n", species); return; }

  TFile* f = TFile::Open(treeFile);
  if (!f || f->IsZombie()) { printf("cannot open %s\n", treeFile); return; }
  TTree* t = (TTree*)f->Get("FemtoTrackTree");
  if (!t) { printf("no FemtoTrackTree\n"); return; }

  femto_phi_tree::TrackRowV2 r;
  const char* need[] = {"speciesCode", "pT", "eta", "phi", "nSigmaProton", "nSigmaDeuteron",
                        "tofBeta", "dca", "nHitsFit", "nHitsMax", "nHitsDedx", "selFlags"};
  t->SetBranchStatus("*", 0);
  for (Int_t i = 0; i < 12; ++i) t->SetBranchStatus(need[i], 1);
  t->SetBranchAddress("speciesCode", &r.speciesCode);
  t->SetBranchAddress("pT", &r.pT);
  t->SetBranchAddress("eta", &r.eta);
  t->SetBranchAddress("phi", &r.phi);
  t->SetBranchAddress("nSigmaProton", &r.nSigmaProton);
  t->SetBranchAddress("nSigmaDeuteron", &r.nSigmaDeuteron);
  t->SetBranchAddress("tofBeta", &r.tofBeta);
  t->SetBranchAddress("dca", &r.dca);
  t->SetBranchAddress("nHitsFit", &r.nHitsFit);
  t->SetBranchAddress("nHitsMax", &r.nHitsMax);
  t->SetBranchAddress("nHitsDedx", &r.nHitsDedx);
  t->SetBranchAddress("selFlags", &r.selFlags);

  printf("\n=== %s ===\n", sc->name);
  printf("cuts (from maker_auau3p85fxt_anaFemtoPhiTree_prod.yaml, TOF rule removed):\n"
         "  |nSigma| < %.1f, DCA < %.1f, pT %.1f-%.1f, |eta| < %.1f, nHitsFit >= %d, ratio >= %.2f",
         sc->nSigmaMax, sc->dcaMax, sc->ptMin, sc->ptMax, sc->etaMax, sc->nHitsFitMin, sc->ratioMin);
  if (sc->nHitsDedxMin > 0) printf(", nHitsDedx >= %d", sc->nHitsDedxMin);
  if (sc->pMax < 1e8) printf(", p %.1f-%.1f", sc->pMin, sc->pMax);
  printf("\n  m2 window [%.1f, %.1f], peak at %.3f, flat-background sideband [%.1f, %.1f]\n",
         sc->m2Lo, sc->m2Hi, sc->m2Peak, sc->sbLo, sc->sbHi);

  const Int_t kNP = 7;
  const Double_t pEdge[kNP + 1] = {0.4, 0.7, 1.0, 1.3, 1.6, 2.0, 2.4, 3.0};
  Long64_t nAll[kNP], nTof[kNP], nIn[kNP];
  Double_t nSb[kNP];
  for (Int_t i = 0; i < kNP; ++i) { nAll[i] = nTof[i] = nIn[i] = 0; nSb[i] = 0; }

  TH2D* hM2vsP = new TH2D(Form("hM2vsP_%s", sc->name),
                          Form("%s: m^{2} vs rigidity;p/Z (GeV/c);m^{2}/Z^{2} (GeV^{2})", sc->name),
                          52, 0.0, 3.2, 300, sc->m2HistLo, sc->m2HistHi);

  const Long64_t n = t->GetEntries();
  if (stride < 1) stride = 1;
  Long64_t seen = 0;
  for (Long64_t i = 0; i < n; i += stride) {
    t->GetEntry(i);
    if (r.speciesCode != sc->code) continue;
    if (sc->code == femto_phi_tree::kSpeciesProton) { if (r.Charge() <= 0) continue; }
    else if (r.Charge() <= 0) continue;
    const Double_t ns = (sc->code == femto_phi_tree::kSpeciesProton) ? r.NSigmaProton()
                                                                    : r.NSigmaNuclear();
    if (TMath::Abs(ns) >= sc->nSigmaMax) continue;
    if (r.Dca() >= sc->dcaMax) continue;
    if (r.Pt() < sc->ptMin || r.Pt() > sc->ptMax) continue;
    if (TMath::Abs(r.Eta()) >= sc->etaMax) continue;
    if (r.NHitsFit() < sc->nHitsFitMin) continue;
    if (sc->nHitsDedxMin > 0 && (Int_t)r.nHitsDedx < sc->nHitsDedxMin) continue;
    if (r.nHitsMax <= 0 || (Double_t)r.NHitsFit() / (Double_t)r.nHitsMax < sc->ratioMin) continue;
    const Double_t p = r.P();
    if (sc->pMax < 1e8 && (p < sc->pMin || p > sc->pMax)) continue;
    ++seen;

    Int_t slice = -1;
    for (Int_t s = 0; s < kNP; ++s) if (p >= pEdge[s] && p < pEdge[s + 1]) { slice = s; break; }
    if (slice < 0) continue;
    ++nAll[slice];
    if (!r.HasTof()) continue;
    const Double_t m2 = r.Mass2();
    if (m2 < -900.0) continue;
    ++nTof[slice];
    hM2vsP->Fill(p, m2);
    if (m2 >= sc->m2Lo && m2 <= sc->m2Hi) ++nIn[slice];
    if (m2 >= sc->sbLo && m2 <= sc->sbHi) nSb[slice] += 1.0;
  }

  // The sideband gives a flat level per unit m2; scale it over the whole histogram range to get
  // the mismatch background that hides under every part of the distribution, the window included.
  const Double_t sbWidth = sc->sbHi - sc->sbLo;
  const Double_t fullWidth = sc->m2HistHi - sc->m2HistLo;

  printf("\n  %-12s %10s %10s %10s %10s %11s %11s\n", "p/Z (GeV/c)", "selected", "TOF", "in window",
         "flat bkg", "purity(win)", "purity(sb)");
  Long64_t sAll = 0, sTof = 0, sIn = 0;
  Double_t sFlat = 0;
  for (Int_t s = 0; s < kNP; ++s) {
    if (nTof[s] < 20) continue;
    const Double_t flat = (sbWidth > 0) ? nSb[s] / sbWidth * fullWidth : 0.0;
    const Double_t pw = (Double_t)nIn[s] / (Double_t)nTof[s];
    const Double_t ps = ((Double_t)nTof[s] - flat > 0)
                            ? ((Double_t)nTof[s] - ((Double_t)nTof[s] - (Double_t)nIn[s]) - 0.0)
                                  / (Double_t)nTof[s]
                            : 0.0;
    // purity(sb): everything in the window, minus the flat background that also sits inside it.
    const Double_t flatInWin = (sbWidth > 0) ? nSb[s] / sbWidth * (sc->m2Hi - sc->m2Lo) : 0.0;
    const Double_t purSb = ((Double_t)nIn[s] - flatInWin) / (Double_t)nTof[s];
    (void)ps;
    printf("  %4.1f - %4.1f  %10lld %10lld %10lld %10.0f %11.4f %11.4f\n", pEdge[s], pEdge[s + 1],
           nAll[s], nTof[s], nIn[s], flatInWin, pw, purSb);
    sAll += nAll[s]; sTof += nTof[s]; sIn += nIn[s]; sFlat += flatInWin;
  }
  if (sTof) {
    printf("  %-12s %10lld %10lld %10lld %10.0f %11.4f %11.4f\n", "all", sAll, sTof, sIn, sFlat,
           (Double_t)sIn / (Double_t)sTof, ((Double_t)sIn - sFlat) / (Double_t)sTof);
    printf("\n  sampled every %lld of %lld rows; %lld tracks pass without the TOF rule\n",
           stride, n, seen);
    printf("  TOF matching efficiency over the selected sample: %.3f\n",
           sAll ? (Double_t)sTof / (Double_t)sAll : 0.0);
  } else {
    printf("  no slice had enough TOF-matched tracks to quote a purity\n");
  }

  if (TString(outFile).Length()) {
    TFile* o = TFile::Open(outFile, "RECREATE");
    hM2vsP->Write();
    o->Close();
    printf("  m2 vs p written to %s\n", outFile);
  }
  f->Close();
}
