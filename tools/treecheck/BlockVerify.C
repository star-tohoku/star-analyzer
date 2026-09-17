// Post-merge gate: does a merged block hold exactly what its inputs held?
//
// This is the check that has to pass before the subjob outputs that went into a block can be
// removed. It is deliberately stronger than "the file opens": for every tree in the block it
// compares the merged entry count against the sum over the inputs, and it does the same for every
// bin of hCounters. Entry counts come from the tree headers, so the cost is one file open per
// input and no event data is read.
//
// Two ROOT details the comparison has to allow for. The tree maker calls AutoSave, so each tree
// is present in an input file under more than one cycle; hadd keeps only the highest cycle of a
// given name, and TFile::Get returns that same highest cycle, so both sides count the same thing.
// FemtoFlagConfig is itself a TTree of configuration rows rather than an event tree, so hadd
// concatenates it: its merged length is the per-file length times the number of inputs, and it is
// reported separately rather than treated as a mismatch.
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TString.h"
#include "TSystem.h"
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

namespace {
const char* kTrees[] = {"FemtoEventTree", "FemtoTrackTree", "FemtoKaonOriginTree",
                        "SourceFileTable"};
const Int_t kNTrees = 4;
}  // namespace

// Returns 0 when the block matches its inputs, non-zero otherwise, so that a caller can gate a
// deletion on the return value rather than on parsing the printout.
Int_t BlockVerify(const char* mergedFile, const char* inputList) {
  std::ifstream in(inputList);
  if (!in.is_open()) { printf("FAIL cannot open input list %s\n", inputList); return 2; }

  TFile* m = TFile::Open(mergedFile);
  if (!m || m->IsZombie()) { printf("FAIL cannot open merged file %s\n", mergedFile); return 2; }

  Long64_t mergedN[kNTrees];
  for (Int_t t = 0; t < kNTrees; ++t) {
    TTree* tr = (TTree*)m->Get(kTrees[t]);
    mergedN[t] = tr ? tr->GetEntries() : -1;
  }
  TH1* mh = (TH1*)m->Get("hCounters");
  TTree* mflag = (TTree*)m->Get("FemtoFlagConfig");
  const Long64_t mergedFlag = mflag ? mflag->GetEntries() : -1;

  Long64_t sumN[kNTrees];
  for (Int_t t = 0; t < kNTrees; ++t) sumN[t] = 0;
  std::vector<Double_t> sumCounters;
  Long64_t sumFlag = 0;
  Int_t nInputs = 0, nUnreadable = 0;
  std::string path;
  while (std::getline(in, path)) {
    if (path.empty()) continue;
    ++nInputs;
    TFile* f = TFile::Open(path.c_str());
    if (!f || f->IsZombie()) {
      printf("  cannot open input %s\n", path.c_str());
      ++nUnreadable;
      if (f) f->Close();
      continue;
    }
    for (Int_t t = 0; t < kNTrees; ++t) {
      TTree* tr = (TTree*)f->Get(kTrees[t]);
      if (tr) sumN[t] += tr->GetEntries();
    }
    TTree* fl = (TTree*)f->Get("FemtoFlagConfig");
    if (fl) sumFlag += fl->GetEntries();
    TH1* h = (TH1*)f->Get("hCounters");
    if (h) {
      if (sumCounters.empty()) sumCounters.assign(h->GetNbinsX() + 1, 0.0);
      for (Int_t b = 1; b <= h->GetNbinsX() && b < (Int_t)sumCounters.size(); ++b)
        sumCounters[b] += h->GetBinContent(b);
    }
    f->Close();
  }

  Int_t bad = nUnreadable;
  printf("%s\n  %d inputs", gSystem->BaseName(mergedFile), nInputs);
  if (nUnreadable) printf("  (%d UNREADABLE)", nUnreadable);
  printf("\n");
  for (Int_t t = 0; t < kNTrees; ++t) {
    const Bool_t okTree = (mergedN[t] == sumN[t]);
    if (!okTree) ++bad;
    printf("  %-20s merged %12lld  inputs %12lld  %s\n", kTrees[t], mergedN[t], sumN[t],
           okTree ? "ok" : "MISMATCH");
  }
  if (mh) {
    Int_t badBins = 0;
    for (Int_t b = 1; b <= mh->GetNbinsX() && b < (Int_t)sumCounters.size(); ++b) {
      if (mh->GetBinContent(b) != sumCounters[b]) {
        printf("  hCounters bin %d (%s): merged %.0f, inputs %.0f  MISMATCH\n", b,
               mh->GetXaxis()->GetBinLabel(b), mh->GetBinContent(b), sumCounters[b]);
        ++badBins;
      }
    }
    bad += badBins;
    printf("  %-20s %d bins, %d mismatched\n", "hCounters", mh->GetNbinsX(), badBins);
  } else {
    printf("  hCounters            MISSING\n");
    ++bad;
  }
  // Informational: the configuration snapshot is concatenated, not summed over events.
  printf("  %-20s merged %12lld  inputs %12lld  (concatenated by hadd)\n", "FemtoFlagConfig",
         mergedFlag, sumFlag);
  m->Close();

  printf("%s %s\n", bad ? "FAIL" : "PASS", mergedFile);
  return bad;
}
