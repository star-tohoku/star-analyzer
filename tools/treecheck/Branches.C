#include "TFile.h"
#include "TTree.h"
#include "TObjArray.h"
#include "TBranch.h"
#include <cstdio>
void Branches(const char* f, const char* tree) {
  TFile* fin = TFile::Open(f);
  TTree* t = (TTree*)fin->Get(tree);
  if (!t) { printf("no tree %s\n", tree); return; }
  TObjArray* br = t->GetListOfBranches();
  printf("%s: %d branches, %lld entries\n", tree, br->GetEntries(), t->GetEntries());
  for (Int_t i = 0; i < br->GetEntries(); ++i)
    printf("  %s\n", ((TBranch*)br->At(i))->GetName());
  fin->Close();
}
