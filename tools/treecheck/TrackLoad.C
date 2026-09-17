#include "TFile.h"
#include "TTree.h"
#include <cstdio>
void TrackLoad(const char* f) {
  TFile* fin = TFile::Open(f);
  TTree* ev = (TTree*)fin->Get("FemtoEventTree");
  TTree* tr = (TTree*)fin->Get("FemtoTrackTree");
  UShort_t nTracks = 0; ev->SetBranchAddress("nTracks", &nTracks);  // schema 3 packs it as UShort_t
  Double_t sum = 0;
  const Long64_t n = ev->GetEntries();
  for (Long64_t i = 0; i < n; ++i) { ev->GetEntry(i); sum += nTracks; }
  printf("events %lld\n", n);
  printf("PicoDst tracks per event : %.1f\n", sum / n);
  printf("stored rows per event    : %.1f\n", (Double_t)tr->GetEntries() / n);
  printf("reduction                : %.1f x\n", sum / tr->GetEntries());
  fin->Close();
}
