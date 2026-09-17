// What the tree can and cannot say about phi-daughter kaon purity.
//
// The baryon measurement in NuclearPurity.C works because the storage envelope for those species
// puts no condition on m2: the whole distribution survives, so the contamination outside the
// species window can be counted. The kaon envelope is different. With
// envKaonRequireDaughterPidReach: true a TOF-matched kaon row is stored only if its m2 is inside
// [envKaonMass2Lo, envKaonMass2Hi] = [0.05, 0.50], and the pion (m2 = 0.019) and proton (0.88)
// peaks fall outside that. The contamination was discarded at production time.
//
// So an ABSOLUTE kaon purity cannot be recovered from the tree; it needs the PicoDst. What can be
// measured is the fraction of stored candidates that land inside the analysis window
// [minMass2Kaon, maxMass2Kaon] = [0.16, 0.36], separately for the two charges. That is a lower
// bound on the contamination, and comparing the charges needs no absolute scale at all.
//
// The macro prints the observed m2 range as evidence of the truncation, so the limitation is
// visible in the output rather than having to be taken on trust.
//
// Usage: KaonPurity.C("<merged block>.root", stride, "<out.root>")
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TString.h"
#include <cstdio>
#include "FemtoPhiTreeSchema.h"

void KaonPurity(const char* treeFile, Long64_t stride = 20, const char* outFile = "") {
  TFile* f = TFile::Open(treeFile);
  if (!f || f->IsZombie()) { printf("cannot open %s\n", treeFile); return; }
  TTree* t = (TTree*)f->Get("FemtoTrackTree");
  if (!t) { printf("no FemtoTrackTree\n"); return; }

  femto_phi_tree::TrackRowV2 r;
  const char* need[] = {"speciesCode", "pT", "eta", "phi", "nSigmaKaon", "tofBeta", "dca",
                        "nHitsFit", "nHitsMax", "selFlags"};
  t->SetBranchStatus("*", 0);
  for (Int_t i = 0; i < 10; ++i) t->SetBranchStatus(need[i], 1);
  t->SetBranchAddress("speciesCode", &r.speciesCode);
  t->SetBranchAddress("pT", &r.pT);
  t->SetBranchAddress("eta", &r.eta);
  t->SetBranchAddress("phi", &r.phi);
  t->SetBranchAddress("nSigmaKaon", &r.nSigmaKaon);
  t->SetBranchAddress("tofBeta", &r.tofBeta);
  t->SetBranchAddress("dca", &r.dca);
  t->SetBranchAddress("nHitsFit", &r.nHitsFit);
  t->SetBranchAddress("nHitsMax", &r.nHitsMax);
  t->SetBranchAddress("selFlags", &r.selFlags);

  const Double_t kNSigma = 3.0;          // pid.nSigmaKaon
  const Double_t kWinLo = 0.16, kWinHi = 0.36;   // pid.minMass2Kaon / maxMass2Kaon
  const Double_t kEnvLo = 0.05, kEnvHi = 0.50;   // envKaonMass2Lo / Hi
  const Int_t kNP = 6;
  const Double_t pEdge[kNP + 1] = {0.0, 0.3, 0.5, 0.8, 1.1, 1.5, 2.5};

  Long64_t nTof[2][kNP], nIn[2][kNP], nNoTof[2][kNP];
  Double_t mn[2] = {1e9, 1e9}, mx[2] = {-1e9, -1e9};
  for (Int_t c = 0; c < 2; ++c) for (Int_t s = 0; s < kNP; ++s)
    { nTof[c][s] = nIn[c][s] = nNoTof[c][s] = 0; }

  TH2D* h2[2];
  h2[0] = new TH2D("hKpM2vsP", "K^{+} candidates: m^{2} vs p;p (GeV/c);m^{2} (GeV^{2})",
                   50, 0.0, 2.5, 120, 0.0, 0.55);
  h2[1] = new TH2D("hKmM2vsP", "K^{-} candidates: m^{2} vs p;p (GeV/c);m^{2} (GeV^{2})",
                   50, 0.0, 2.5, 120, 0.0, 0.55);

  const Long64_t n = t->GetEntries();
  if (stride < 1) stride = 1;
  for (Long64_t i = 0; i < n; i += stride) {
    t->GetEntry(i);
    Int_t c = -1;
    if (r.speciesCode == femto_phi_tree::kSpeciesKp) c = 0;
    else if (r.speciesCode == femto_phi_tree::kSpeciesKm) c = 1;
    if (c < 0) continue;
    if (TMath::Abs(r.NSigmaKaon()) >= kNSigma) continue;
    const Double_t p = r.P();
    Int_t s = -1;
    for (Int_t k = 0; k < kNP; ++k) if (p >= pEdge[k] && p < pEdge[k + 1]) { s = k; break; }
    if (s < 0) continue;
    if (!r.HasTof()) { ++nNoTof[c][s]; continue; }
    const Double_t m2 = r.Mass2();
    if (m2 < -900.0) continue;
    ++nTof[c][s];
    h2[c]->Fill(p, m2);
    if (m2 < mn[c]) mn[c] = m2;
    if (m2 > mx[c]) mx[c] = m2;
    if (m2 >= kWinLo && m2 <= kWinHi) ++nIn[c][s];
  }

  printf("\nkaon m2 envelope [%.2f, %.2f]; observed range K+ [%.4f, %.4f], K- [%.4f, %.4f]\n",
         kEnvLo, kEnvHi, mn[0], mx[0], mn[1], mx[1]);
  printf("The observed range stopping at the envelope is the evidence that everything outside it\n"
         "was discarded at production time: no ABSOLUTE purity can be quoted from this tree.\n");
  const char* lab[2] = {"K+", "K-"};
  for (Int_t c = 0; c < 2; ++c) {
    printf("\n  === %s (|nSigma_K| < %.1f) ===\n", lab[c], kNSigma);
    printf("  %-12s %12s %12s %12s %12s\n", "p (GeV/c)", "TOF matched", "no TOF",
           "in [0.16,0.36]", "fraction");
    Long64_t sT = 0, sI = 0, sN = 0;
    for (Int_t s = 0; s < kNP; ++s) {
      if (nTof[c][s] + nNoTof[c][s] == 0) continue;
      printf("  %4.1f - %4.1f  %12lld %12lld %12lld %12.4f\n", pEdge[s], pEdge[s + 1],
             nTof[c][s], nNoTof[c][s], nIn[c][s],
             nTof[c][s] ? (Double_t)nIn[c][s] / (Double_t)nTof[c][s] : 0.0);
      sT += nTof[c][s]; sI += nIn[c][s]; sN += nNoTof[c][s];
    }
    printf("  %-12s %12lld %12lld %12lld %12.4f\n", "all", sT, sN, sI,
           sT ? (Double_t)sI / (Double_t)sT : 0.0);
  }

  if (TString(outFile).Length()) {
    TFile* o = TFile::Open(outFile, "RECREATE");
    h2[0]->Write(); h2[1]->Write();
    o->Close();
    printf("\n  m2 vs p written to %s\n", outFile);
  }
  f->Close();
}
