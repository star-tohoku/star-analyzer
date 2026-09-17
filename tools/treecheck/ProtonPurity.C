// How pure is a TPC-only proton selection, and where does it stop being pure?
//
// The production analysis requires a TOF match on every proton (protonTofMomentumThreshold 0.0).
// Relaxing that to the dE/dx-only study's 2.0 GeV/c raises the phi-p correlation significance from
// 1.50 to 3.59 on the full data set, and the correlation amplitude does not visibly fall -- but
// that is an indirect argument. This measures the thing directly.
//
// The method is the standard one: the TOF-matched subsample tells you what the TPC-only selection
// would have contained. Take every track the analysis would accept as a proton WITHOUT the TOF
// rule, keep those that happen to carry a TOF match, and ask what fraction of them lands inside
// the proton m2 window. That fraction is the purity of the TPC-only selection, measured in slices
// of momentum so the answer says not just "how pure" but "up to what momentum".
//
// The one assumption is that TOF matching is unbiased with respect to species. It is not exactly
// true -- the match probability depends on momentum and on where a track lands on the tray -- so
// the number is an estimate, not a calibration. What it is good for is finding the momentum at
// which contamination stops being negligible, which is what the cut has to be set against.
//
// Usage: ProtonPurity.C("<merged block>.root", stride)
//   stride samples every Nth row; the full block holds 4.3e8 rows and the answer needs ~1e6.
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TString.h"
#include <cstdio>
#include "FemtoPhiTreeSchema.h"

void ProtonPurity(const char* treeFile, Long64_t stride = 100, const char* outFile = "") {
  TFile* f = TFile::Open(treeFile);
  if (!f || f->IsZombie()) { printf("cannot open %s\n", treeFile); return; }
  TTree* t = (TTree*)f->Get("FemtoTrackTree");
  if (!t) { printf("no FemtoTrackTree\n"); return; }

  femto_phi_tree::TrackRowV2 r;
  const char* need[] = {"speciesCode", "pT", "eta", "phi", "nSigmaProton", "nSigmaPion",
                        "nSigmaKaon", "tofBeta", "dca", "nHitsFit", "nHitsMax", "selFlags"};
  t->SetBranchStatus("*", 0);
  for (Int_t i = 0; i < 12; ++i) t->SetBranchStatus(need[i], 1);
  t->SetBranchAddress("speciesCode", &r.speciesCode);
  t->SetBranchAddress("pT", &r.pT);
  t->SetBranchAddress("eta", &r.eta);
  t->SetBranchAddress("phi", &r.phi);
  t->SetBranchAddress("nSigmaProton", &r.nSigmaProton);
  t->SetBranchAddress("nSigmaPion", &r.nSigmaPion);
  t->SetBranchAddress("nSigmaKaon", &r.nSigmaKaon);
  t->SetBranchAddress("tofBeta", &r.tofBeta);
  t->SetBranchAddress("dca", &r.dca);
  t->SetBranchAddress("nHitsFit", &r.nHitsFit);
  t->SetBranchAddress("nHitsMax", &r.nHitsMax);
  t->SetBranchAddress("selFlags", &r.selFlags);

  // The production proton cuts, minus the TOF rule. Values from
  // maker_auau3p85fxt_anaFemtoPhiTree_prod.yaml and pid_auau3p85fxt_anaPhi.yaml.
  const Double_t kNSigmaMax = 2.0, kDcaMax = 1.0, kPtMin = 0.4, kPtMax = 2.0, kEtaMax = 2.0;
  const Short_t kNHitsFitMin = 15;
  const Double_t kRatioMin = 0.52;
  const Double_t kM2Lo = 0.6, kM2Hi = 1.2;

  const Int_t kNP = 10;
  const Double_t pEdge[kNP + 1] = {0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.4, 3.0};
  Long64_t nTof[kNP], nIn[kNP], nAll[kNP];
  for (Int_t i = 0; i < kNP; ++i) { nTof[i] = nIn[i] = nAll[i] = 0; }

  TH2D* hM2vsP = new TH2D("hM2vsP", "m^{2} vs p for TPC-selected protons;p (GeV/c);m^{2} (GeV^{2})",
                          60, 0.0, 3.0, 200, -0.5, 2.5);

  const Long64_t n = t->GetEntries();
  if (stride < 1) stride = 1;
  Long64_t seen = 0;
  for (Long64_t i = 0; i < n; i += stride) {
    t->GetEntry(i);
    if (r.speciesCode != femto_phi_tree::kSpeciesProton) continue;
    // Everything the analysis asks of a proton except the TOF rule.
    if (r.Charge() <= 0) continue;
    if (TMath::Abs(r.NSigmaProton()) >= kNSigmaMax) continue;
    if (r.Dca() >= kDcaMax) continue;
    if (r.Pt() < kPtMin || r.Pt() > kPtMax) continue;
    if (TMath::Abs(r.Eta()) >= kEtaMax) continue;
    if (r.NHitsFit() < kNHitsFitMin) continue;
    if (r.nHitsMax <= 0 || (Double_t)r.NHitsFit() / (Double_t)r.nHitsMax < kRatioMin) continue;
    ++seen;

    const Double_t p = r.P();
    Int_t slice = -1;
    for (Int_t s = 0; s < kNP; ++s) if (p >= pEdge[s] && p < pEdge[s + 1]) { slice = s; break; }
    if (slice < 0) continue;
    ++nAll[slice];
    if (!r.HasTof()) continue;
    const Double_t m2 = r.Mass2();
    if (m2 < -900.0) continue;          // TOF matched but the timing is unusable
    ++nTof[slice];
    hM2vsP->Fill(p, m2);
    if (m2 >= kM2Lo && m2 <= kM2Hi) ++nIn[slice];
  }

  printf("\nsampled every %lld rows of %lld; %lld tracks pass the proton cuts without the TOF rule\n",
         stride, n, seen);
  printf("\n  %-14s %12s %12s %12s %10s %10s\n", "p (GeV/c)", "selected", "TOF matched",
         "in m2 window", "purity", "stat err");
  Long64_t sAll = 0, sTof = 0, sIn = 0;
  for (Int_t s = 0; s < kNP; ++s) {
    if (nTof[s] == 0) continue;
    const Double_t pur = (Double_t)nIn[s] / (Double_t)nTof[s];
    const Double_t err = TMath::Sqrt(pur * (1.0 - pur) / (Double_t)nTof[s]);
    printf("  %5.1f - %5.1f  %12lld %12lld %12lld %10.4f %10.4f\n", pEdge[s], pEdge[s + 1],
           nAll[s], nTof[s], nIn[s], pur, err);
    sAll += nAll[s]; sTof += nTof[s]; sIn += nIn[s];
  }
  if (sTof) {
    const Double_t pur = (Double_t)sIn / (Double_t)sTof;
    printf("  %-14s %12lld %12lld %12lld %10.4f %10.4f\n", "all", sAll, sTof, sIn, pur,
           TMath::Sqrt(pur * (1.0 - pur) / (Double_t)sTof));
    printf("\nTOF matching efficiency over the selected sample: %.3f\n",
           sAll ? (Double_t)sTof / (Double_t)sAll : 0.0);
  }
  printf("\nPurity here is the fraction of TOF-matched tracks inside m2 [%.1f, %.1f]; it estimates\n"
         "what a TPC-only selection would contain, assuming TOF matching does not prefer a species.\n",
         kM2Lo, kM2Hi);

  if (TString(outFile).Length()) {
    TFile* o = TFile::Open(outFile, "RECREATE");
    hM2vsP->Write();
    o->Close();
    printf("m2 vs p written to %s\n", outFile);
  }
  f->Close();
}
