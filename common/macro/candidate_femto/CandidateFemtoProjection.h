#ifndef CANDIDATE_FEMTO_PROJECTION_H
#define CANDIDATE_FEMTO_PROJECTION_H

// Projections of the StCandidateFemtoMaker pair histograms (THnSparseF hPairSE_<channel> /
// hPairME_<channel>, axes: 0 k*, 1 cent9, 2 mass class A, 3 mass class B, 4 mass A, 5 mass B,
// 6 BDT score A, 7 BDT score B) for the analysis macros in common/macro/candidate_femto/.
// Selections: a cent9 list (empty = all) and BDT score thresholds per species (a score bin is
// kept when its lower edge is >= the threshold; the threshold must be a score-axis bin edge).
// Mass classes: 0 none, 1 signal, 2 left sideband, 3 right sideband (CandidateFemtoConfig.h).

#include "THnSparse.h"
#include "TAxis.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TString.h"

#include <iostream>
#include <vector>

struct CandidatePairSelection {
  CandidatePairSelection() : scoreMinA(-1e9), scoreMinB(-1e9) {}
  std::vector<Int_t> cent9;  // empty = all centralities
  Double_t scoreMinA, scoreMinB;
  bool AcceptCent(Int_t c) const {
    if (cent9.empty()) return true;
    for (size_t i = 0; i < cent9.size(); ++i) if (cent9[i] == c) return true;
    return false;
  }
  // "4,5,6,7,8" -> {4, 5, 6, 7, 8}; "" -> all
  void SetCentList(const char* list) {
    cent9.clear();
    TObjArray* tok = TString(list ? list : "").Tokenize(",");
    for (Int_t i = 0; i < tok->GetEntries(); ++i) cent9.push_back(((TObjString*)tok->At(i))->GetString().Atoi());
    delete tok;
  }
  TString Describe() const {
    TString s = "cent9 ";
    if (cent9.empty()) s += "all";
    for (size_t i = 0; i < cent9.size(); ++i) s += TString::Format("%s%d", i ? "," : "", cent9[i]);
    if (scoreMinA > -1e8) s += TString::Format(", score A >= %g", scoreMinA);
    if (scoreMinB > -1e8) s += TString::Format(", score B >= %g", scoreMinB);
    return s;
  }
};

// kFALSE (and a message) if a score threshold is not a bin edge of the score axis.
inline Bool_t CandidateScoreThresholdIsEdge(THnSparse* h, Int_t axis, Double_t threshold) {
  if (threshold < -1e8) return kTRUE;
  const TAxis* a = h->GetAxis(axis);
  for (Int_t b = 1; b <= a->GetNbins() + 1; ++b)
    if (TMath::Abs(a->GetBinLowEdge(b) - threshold) < 1e-9) return kTRUE;
  std::cerr << "ERROR: score threshold " << threshold << " is not a bin edge of axis " << a->GetName() << " of "
            << h->GetName() << std::endl;
  return kFALSE;
}

// Loop over the filled bins that pass the selection; f(coords, content, error2).
template <class F>
void CandidateForEachSelectedBin(THnSparse* h, const CandidatePairSelection& sel, F& f) {
  std::vector<Int_t> idx(h->GetNdimensions());
  const TAxis* centAxis = h->GetAxis(1);
  const TAxis* scoreA = h->GetAxis(6);
  const TAxis* scoreB = h->GetAxis(7);
  for (Long64_t i = 0; i < h->GetNbins(); ++i) {
    const Double_t w = h->GetBinContent(i, &idx[0]);
    if (w == 0) continue;
    if (!sel.AcceptCent(TMath::Nint(centAxis->GetBinCenter(idx[1])))) continue;
    if (scoreA->GetBinLowEdge(idx[6]) < sel.scoreMinA - 1e-9) continue;
    if (scoreB->GetBinLowEdge(idx[7]) < sel.scoreMinB - 1e-9) continue;
    f(&idx[0], w, h->GetBinError2(i));
  }
}

// (k*, mass A, mass B) for every mass class, restricted to the mass ranges [aLo, aHi) x [bLo, bHi)
// (whole bins: the ranges should be bin edges of the mass axes).
struct CandidateFillKstarMassMass {
  TH3D* out;
  THnSparse* in;
  void operator()(const Int_t* idx, Double_t w, Double_t e2) {
    const Int_t b = out->FindBin(in->GetAxis(0)->GetBinCenter(idx[0]), in->GetAxis(4)->GetBinCenter(idx[4]),
                                 in->GetAxis(5)->GetBinCenter(idx[5]));
    if (out->IsBinUnderflow(b) || out->IsBinOverflow(b)) return;
    out->AddBinContent(b, w);
    out->GetSumw2()->fArray[b] += e2;
  }
};

inline TH3D* CandidateProjectKstarMassMass(THnSparse* h, const CandidatePairSelection& sel, Double_t aLo, Double_t aHi,
                                           Double_t bLo, Double_t bHi, const char* name) {
  TAxis* k = h->GetAxis(0);
  TAxis* a = h->GetAxis(4);
  TAxis* b = h->GetAxis(5);
  const Int_t a1 = a->FindBin(aLo + 1e-9), a2 = a->FindBin(aHi - 1e-9);
  const Int_t b1 = b->FindBin(bLo + 1e-9), b2 = b->FindBin(bHi - 1e-9);
  TH3D* out = new TH3D(name, TString::Format("%s;k* [GeV/c];mass A;mass B", h->GetTitle()), k->GetNbins(), k->GetXmin(),
                       k->GetXmax(), a2 - a1 + 1, a->GetBinLowEdge(a1), a->GetBinUpEdge(a2), b2 - b1 + 1,
                       b->GetBinLowEdge(b1), b->GetBinUpEdge(b2));
  out->SetDirectory(0);
  out->Sumw2();
  CandidateFillKstarMassMass f;
  f.out = out;
  f.in = h;
  CandidateForEachSelectedBin(h, sel, f);
  out->SetEntries(out->GetSumOfWeights());
  return out;
}

// k* (and k* vs cent9) for one mass-class combination (classA, classB); -1 = any class.
struct CandidateFillKstar {
  TH1D* out;
  TH2D* outCent;
  THnSparse* in;
  Int_t clsA, clsB;
  void operator()(const Int_t* idx, Double_t w, Double_t e2) {
    const Int_t ca = TMath::Nint(in->GetAxis(2)->GetBinCenter(idx[2]));
    const Int_t cb = TMath::Nint(in->GetAxis(3)->GetBinCenter(idx[3]));
    if ((clsA >= 0 && ca != clsA) || (clsB >= 0 && cb != clsB)) return;
    const Double_t ks = in->GetAxis(0)->GetBinCenter(idx[0]);
    const Int_t b = out->FindBin(ks);
    out->AddBinContent(b, w);
    out->GetSumw2()->fArray[b] += e2;
    if (outCent) {
      const Int_t bc = outCent->FindBin(ks, in->GetAxis(1)->GetBinCenter(idx[1]));
      outCent->AddBinContent(bc, w);
      outCent->GetSumw2()->fArray[bc] += e2;
    }
  }
};

inline TH1D* CandidateProjectKstar(THnSparse* h, const CandidatePairSelection& sel, Int_t clsA, Int_t clsB,
                                   const char* name, TH2D** kstarVsCent = 0) {
  const TAxis* k = h->GetAxis(0);
  TH1D* out = new TH1D(name, TString::Format("%s;k* [GeV/c];pairs", h->GetTitle()), k->GetNbins(), k->GetXmin(), k->GetXmax());
  out->SetDirectory(0);
  out->Sumw2();
  CandidateFillKstar f;
  f.out = out;
  f.outCent = 0;
  if (kstarVsCent) {
    f.outCent = new TH2D(TString(name) + "_vsCent", ";k* [GeV/c];cent9", k->GetNbins(), k->GetXmin(), k->GetXmax(), 9, -0.5, 8.5);
    f.outCent->SetDirectory(0);
    f.outCent->Sumw2();
    *kstarVsCent = f.outCent;
  }
  f.in = h;
  f.clsA = clsA;
  f.clsB = clsB;
  CandidateForEachSelectedBin(h, sel, f);
  out->SetEntries(out->GetSumOfWeights());
  return out;
}

#endif
