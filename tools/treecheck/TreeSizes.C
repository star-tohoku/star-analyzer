// Per-tree compressed size of a reduced tree file, to settle what a quoted B/event includes.
#include "TFile.h"
#include "TTree.h"
#include "TKey.h"
#include "TList.h"
#include <cstdio>
#include <set>
#include <string>

void TreeSizes(const char* file) {
  TFile* f = TFile::Open(file);
  if (!f) { printf("cannot open %s\n", file); return; }
  printf("\n%s   file size %lld B\n", file, f->GetSize());
  Long64_t tot = 0;
  std::set<std::string> seen;
  TIter next(f->GetListOfKeys());
  TKey* k;
  while ((k = (TKey*)next())) {
    if (TString(k->GetClassName()) != "TTree") continue;
    if (!seen.insert(k->GetName()).second) continue;   // keys are cycled; count once
    TTree* t = (TTree*)f->Get(k->GetName());
    if (!t) continue;
    printf("  %-24s %12lld entries  %12lld B zip\n", t->GetName(), t->GetEntries(),
           (Long64_t)t->GetZipBytes());
    tot += (Long64_t)t->GetZipBytes();
  }
  printf("  %-24s %12s  %12lld B zip\n", "TOTAL trees", "", tot);
  f->Close();
}
