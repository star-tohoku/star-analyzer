// Scan the phi-daughter kaon PID and judge it on the phi signal, not on the correlation function.
//
// Why not the correlation function. The k* bins where the correlation lives hold 10-15 same-event
// pairs in the whole 2.0e9-event data set, and section 16 of the 2026-09-16 record shows that any
// significance built on them is unstable -- a bin went from 9 to 13 pairs and its error got worse.
// The phi signal itself is a different matter: tens of thousands of candidates per block, so S, B
// and S/sqrt(S+B) are stable to well under a percent and can rank cut choices.
//
// Why one pass. Every variation is a different kaon PID applied to the same tracks, so all of them
// can be filled in a single sweep of the tree. That replaces one farm pass per point (141 jobs)
// with one local run of a few minutes.
//
// What is reproduced exactly, and what is not. The kaon track selection is read from the flag the
// maker stored (kSelKaonCutsNom), so it is the maker's own. The daughter PID is reimplemented here
// because varying it is the whole point, and it follows StPhiKKReconstruction::PassPhiDaughterTofPid:
//   K+ : no TOF match -> accept if p <= pMomKaonPID; TOF match -> m2 inside the window
//   K- : TOF match always required, then m2 inside the window
// The pair-level selection applied is the lab-frame rapidity window [minPairRapidity,
// maxPairRapidity] = [-1.5, 0.0]. No mixed-event or rotated background is built: the comparison is
// between variations of the same quantity, and a same-event fit separates signal from combinatorial
// background well enough to rank them.
//
// Usage: PhiYieldScan.C("<merged block>.root", nEventsMax, "<out.root>")
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TF1.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TString.h"
#include <cstdio>
#include <vector>
#include "FemtoPhiTreeSchema.h"

namespace {
const Double_t kKaonMass = 0.493677;
const Double_t kPMomKaonPID = 0.5;       // pid.pMomKaonPID
const Double_t kMinPairY = -1.5, kMaxPairY = 0.0;   // minPairRapidity / maxPairRapidity, lab frame

struct Variation {
  const char* name;
  Double_t nSigma;   // |nSigma_K| limit
  Double_t m2Lo, m2Hi;
  Bool_t kmRequireTof;
};

// Nominal is nSigma 3.0, window [0.16, 0.36], K- requires TOF. The window scan goes both ways:
// section 14e found that only 17% of stored K- candidates land inside the window, which suggests
// narrowing may drop more background than signal.
const Variation kVar[] = {
  {"nominal",        3.0, 0.16, 0.36, kTRUE},
  {"nSigma2.0",      2.0, 0.16, 0.36, kTRUE},
  {"nSigma2.5",      2.5, 0.16, 0.36, kTRUE},
  {"nSigma4.0",      4.0, 0.16, 0.36, kTRUE},
  {"nSigma5.0",      5.0, 0.16, 0.36, kTRUE},
  {"m2_0.20-0.32",   3.0, 0.20, 0.32, kTRUE},
  {"m2_0.18-0.34",   3.0, 0.18, 0.34, kTRUE},
  {"m2_0.14-0.40",   3.0, 0.14, 0.40, kTRUE},
  {"m2_0.10-0.45",   3.0, 0.10, 0.45, kTRUE},
  {"kmNoTof",        3.0, 0.16, 0.36, kFALSE}
};
const Int_t kNVar = 10;

Bool_t PassDaughterPid(const femto_phi_tree::TrackRowV2& r, const Variation& v) {
  if (TMath::Abs(r.NSigmaKaon()) >= v.nSigma) return kFALSE;
  const Bool_t neg = (r.Charge() < 0);
  if (!r.HasTof()) {
    if (neg && v.kmRequireTof) return kFALSE;
    return r.P() <= kPMomKaonPID;
  }
  const Double_t m2 = r.Mass2();
  return (m2 >= v.m2Lo && m2 <= v.m2Hi);
}
}  // namespace

void FillPairs(TH1D** h, const std::vector<femto_phi_tree::TrackRowV2>& kp,
               const std::vector<femto_phi_tree::TrackRowV2>& km) {
  for (Int_t v = 0; v < kNVar; ++v) {
    for (size_t a = 0; a < kp.size(); ++a) {
      if (!PassDaughterPid(kp[a], kVar[v])) continue;
      TLorentzVector p1;
      p1.SetPtEtaPhiM(kp[a].Pt(), kp[a].Eta(), kp[a].Phi(), kKaonMass);
      for (size_t b = 0; b < km.size(); ++b) {
        if (!PassDaughterPid(km[b], kVar[v])) continue;
        TLorentzVector p2;
        p2.SetPtEtaPhiM(km[b].Pt(), km[b].Eta(), km[b].Phi(), kKaonMass);
        const TLorentzVector pair = p1 + p2;
        const Double_t y = pair.Rapidity();
        if (y < kMinPairY || y > kMaxPairY) continue;
        h[v]->Fill(pair.M());
      }
    }
  }
}

void PhiYieldScan(const char* treeFile, Long64_t nEventsMax = -1, const char* outFile = "") {
  TFile* f = TFile::Open(treeFile);
  if (!f || f->IsZombie()) { printf("cannot open %s\n", treeFile); return; }
  TTree* t = (TTree*)f->Get("FemtoTrackTree");
  if (!t) { printf("no FemtoTrackTree\n"); return; }

  femto_phi_tree::TrackRowV2 r;
  ULong64_t uid = 0;
  const char* need[] = {"eventUID", "speciesCode", "pT", "eta", "phi", "nSigmaKaon", "tofBeta",
                        "nHitsFit", "selFlags"};
  t->SetBranchStatus("*", 0);
  for (Int_t i = 0; i < 9; ++i) t->SetBranchStatus(need[i], 1);
  t->SetBranchAddress("eventUID", &uid);
  t->SetBranchAddress("speciesCode", &r.speciesCode);
  t->SetBranchAddress("pT", &r.pT);
  t->SetBranchAddress("eta", &r.eta);
  t->SetBranchAddress("phi", &r.phi);
  t->SetBranchAddress("nSigmaKaon", &r.nSigmaKaon);
  t->SetBranchAddress("tofBeta", &r.tofBeta);
  t->SetBranchAddress("nHitsFit", &r.nHitsFit);
  t->SetBranchAddress("selFlags", &r.selFlags);

  TH1D* h[kNVar];
  for (Int_t v = 0; v < kNVar; ++v)
    h[v] = new TH1D(Form("hMkk_%d", v), Form("M(KK) %s;M_{KK} (GeV/c^{2});pairs", kVar[v].name),
                    160, 0.98, 1.10);

  std::vector<femto_phi_tree::TrackRowV2> kp, km;
  ULong64_t prevUid = 0;
  Bool_t first = kTRUE;
  Long64_t nEvt = 0, nPairEvt = 0;
  const Long64_t n = t->GetEntries();

  // The rows are written event by event, so a change of eventUID closes an event. An event with
  // no kaon at all still has to advance prevUid: an earlier version only advanced it when the
  // buffers were non-empty, which left a stale id behind and made every later row look like a new
  // event, so no two kaons ever shared one. It showed up as a phi signal of one pair in 200,000
  // events.
  for (Long64_t i = 0; i < n; ++i) {
    t->GetEntry(i);
    if (first) { prevUid = uid; first = kFALSE; }
    else if (uid != prevUid) {
      if (kp.size() && km.size()) { ++nPairEvt; FillPairs(h, kp, km); }
      kp.clear(); km.clear();
      ++nEvt;
      prevUid = uid;
      if (nEventsMax > 0 && nEvt >= nEventsMax) break;
    }
    if (!(r.selFlags & femto_phi_tree::kSelKaonCutsNom)) continue;
    if (r.speciesCode == femto_phi_tree::kSpeciesKp) kp.push_back(r);
    else if (r.speciesCode == femto_phi_tree::kSpeciesKm) km.push_back(r);
  }
  if (kp.size() && km.size()) { ++nPairEvt; FillPairs(h, kp, km); }
  if (!first) ++nEvt;

  printf("\n%lld events read, %lld with both a K+ and a K-, from %s\n", nEvt, nPairEvt, treeFile);
  printf("\n  %-16s %10s %10s %8s %10s %10s\n", "variation", "signal", "background", "S/B",
         "S/sqrt(S+B)", "rel. to nom");
  Double_t nomFom = 0;
  for (Int_t v = 0; v < kNVar; ++v) {
    // A Gaussian for the phi on a quadratic combinatorial background, over the usual fit range.
    TF1* fit = new TF1(Form("f%d", v), "gaus(0)+pol2(3)", 0.99, 1.06);
    fit->SetParameters(h[v]->GetMaximum() * 0.3, 1.0195, 0.005, h[v]->GetMaximum(), 0.0, 0.0);
    fit->SetParLimits(1, 1.014, 1.025);
    fit->SetParLimits(2, 0.0015, 0.015);
    h[v]->Fit(fit, "QNR");
    const Double_t bw = h[v]->GetBinWidth(1);
    const Double_t sig = fit->GetParameter(0) * fit->GetParameter(2) * TMath::Sqrt(2 * TMath::Pi()) / bw;
    TF1* bkg = new TF1(Form("b%d", v), "pol2", 0.99, 1.06);
    bkg->SetParameters(fit->GetParameter(3), fit->GetParameter(4), fit->GetParameter(5));
    const Double_t lo = fit->GetParameter(1) - 2 * fit->GetParameter(2);
    const Double_t hi = fit->GetParameter(1) + 2 * fit->GetParameter(2);
    const Double_t bkgIn = bkg->Integral(lo, hi) / bw;
    const Double_t sigIn = sig * 0.9545;          // a 2-sigma window holds 95.45% of a Gaussian
    const Double_t fom = (sigIn + bkgIn > 0) ? sigIn / TMath::Sqrt(sigIn + bkgIn) : 0.0;
    if (v == 0) nomFom = fom;
    printf("  %-16s %10.0f %10.0f %8.3f %10.1f %10.3f\n", kVar[v].name, sigIn, bkgIn,
           bkgIn > 0 ? sigIn / bkgIn : 0.0, fom, nomFom > 0 ? fom / nomFom : 0.0);
  }
  printf("\nS and B are counted inside +-2 sigma of the fitted phi peak. S/sqrt(S+B) is the figure\n"
         "of merit: it is what the statistical error on an extracted phi yield goes as.\n");

  if (TString(outFile).Length()) {
    TFile* o = TFile::Open(outFile, "RECREATE");
    for (Int_t v = 0; v < kNVar; ++v) { h[v]->SetName(Form("hMkk_%s", kVar[v].name)); h[v]->Write(); }
    o->Close();
    printf("M(KK) histograms written to %s\n", outFile);
  }
  f->Close();
}
