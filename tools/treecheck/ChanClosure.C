// Full channel-set closure: StFemtoMaker direct vs schema 3 tree + anaFemtoPhiTreeDownstreamV3.
// Both sides now write hKstarSE_<channel> / hKstarME_<channel> with the same binning, so the
// comparison is bin for bin with no rebinning.
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TKey.h"
#include "TString.h"
#include "TMath.h"
#include <cstdio>
#include <set>
#include <string>
#include <vector>

namespace {
void Cmp(const char* label, TH1* mk, TH1* tr, Int_t& nExact, Int_t& nTotals, Int_t& nBad) {
  if (!mk && !tr) return;
  if (!mk || !tr) {
    printf("  %-34s %s\n", label, !mk ? "MAKER MISSING" : "TREE MISSING");
    ++nBad;
    return;
  }
  const Double_t im = mk->Integral(0, -1), it = tr->Integral(0, -1);
  Int_t nDiff = 0;
  Double_t sumAbs = 0, maxAbs = 0;
  // config/hist/hist_anaFemtoPhi.yaml gives the h-K channels with a nuclear bachelor a coarse
  // 150-bin k* axis while everything else gets 300 bins over the same range. The downstream
  // books the fine axis for every channel, so the finer side is rebinned here before comparing.
  TH1* trUse = tr;
  if (tr->GetNbinsX() > mk->GetNbinsX() && tr->GetNbinsX() % mk->GetNbinsX() == 0) {
    trUse = (TH1*)tr->Clone(Form("%s_reb", tr->GetName()));
    trUse->Rebin(tr->GetNbinsX() / mk->GetNbinsX());
  }
  const Int_t nb = TMath::Min(mk->GetNbinsX(), trUse->GetNbinsX());
  for (Int_t i = 1; i <= nb; ++i) {
    const Double_t d = trUse->GetBinContent(i) - mk->GetBinContent(i);
    if (d != 0) { ++nDiff; sumAbs += TMath::Abs(d); if (TMath::Abs(d) > maxAbs) maxAbs = TMath::Abs(d); }
  }
  const char* verdict = (im == it && nDiff == 0) ? "EXACT"
                        : (im == it)            ? "totals match"
                                                : "DIFFER";
  if (im == it && nDiff == 0) ++nExact; else if (im == it) ++nTotals; else ++nBad;
  printf("  %-34s maker=%10.0f tree=%10.0f  diff=%+8.0f  bins=%4d sum|d|=%8.0f  %s\n", label, im,
         it, it - im, nDiff, sumAbs, verdict);
}
}  // namespace

void ChanClosure(const char* makerFile, const char* downFile) {
  TFile* fm = TFile::Open(makerFile);
  TFile* fd = TFile::Open(downFile);
  if (!fm || !fd) { printf("cannot open inputs\n"); return; }
  printf("\n=== full channel-set closure ===\nmaker: %s\ntree : %s\n\n", makerFile, downFile);

  // channel list = every hKstarSE_* the downstream wrote
  std::vector<std::string> chans;
  std::set<std::string> seen;
  TIter next(fd->GetListOfKeys());
  TKey* k;
  while ((k = (TKey*)next())) {
    TString n(k->GetName());
    if (!n.BeginsWith("hKstarSE_")) continue;
    std::string ch(n.Data() + 9);
    if (seen.insert(ch).second) chans.push_back(ch);
  }

  Int_t nExact = 0, nTotals = 0, nBad = 0;
  for (size_t i = 0; i < chans.size(); ++i) {
    const char* c = chans[i].c_str();
    printf("--- %s ---\n", c);
    Cmp(Form("SE k*"), (TH1*)fm->Get(Form("hKstarSE_%s", c)), (TH1*)fd->Get(Form("hKstarSE_%s", c)),
        nExact, nTotals, nBad);
    Cmp(Form("ME k*"), (TH1*)fm->Get(Form("hKstarME_%s", c)), (TH1*)fd->Get(Form("hKstarME_%s", c)),
        nExact, nTotals, nBad);
    TH2* m2 = (TH2*)fm->Get(Form("hKstarSEVsCent_%s", c));
    TH2* t2 = (TH2*)fd->Get(Form("hKstarSEVsCent_%s", c));
    if (m2 && t2)
      printf("  %-34s maker=%10.0f tree=%10.0f  %s\n", "SE k* vs cent (TH2)", m2->Integral(0, -1, 0, -1),
             t2->Integral(0, -1, 0, -1),
             m2->Integral(0, -1, 0, -1) == t2->Integral(0, -1, 0, -1) ? "totals match" : "DIFFER");
    TH3* m3 = (TH3*)fm->Get(Form("hPhiMKK_vs_KstarSE_%s", c));
    TH3* t3 = (TH3*)fd->Get(Form("hPhiMKK_vs_KstarSE_%s", c));
    if (m3 && t3)
      printf("  %-34s maker=%10.0f tree=%10.0f  %s\n", "M(KK) vs k* vs cent (TH3)",
             m3->Integral(0, -1, 0, -1, 0, -1), t3->Integral(0, -1, 0, -1, 0, -1),
             m3->Integral(0, -1, 0, -1, 0, -1) == t3->Integral(0, -1, 0, -1, 0, -1) ? "totals match"
                                                                                   : "DIFFER");
  }
  printf("\nsummary: %d exact, %d totals-only, %d differing (over %d channels)\n", nExact, nTotals,
         nBad, (Int_t)chans.size());
  fm->Close(); fd->Close();
}
