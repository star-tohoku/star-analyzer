#include "TFile.h"
#include "TH1.h"
#include <cstdio>
void CounterCheck(const char* f) {
  TFile* fin = TFile::Open(f);
  TH1* h = (TH1*)fin->Get("hCounters");
  if (!h) { printf("no hCounters\n"); return; }
  printf("hCounters: %d bins\n", h->GetNbinsX());
  for (Int_t i = 1; i <= h->GetNbinsX(); ++i)
    printf("  %-18s %12.0f\n", h->GetXaxis()->GetBinLabel(i), h->GetBinContent(i));
  fin->Close();
}
