// Bytes per input event for one produced tree, together with the input file it came from, so the
// cost can be split by input stream (st_physics vs st_physics_adc).
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TString.h"
#include "TSystem.h"
#include <cstdio>
void PerFileCost(const char* f) {
  TFile* fin = TFile::Open(f);
  if (!fin || fin->IsZombie()) { printf("BAD %s\n", f); return; }
  TH1* h = (TH1*)fin->Get("hCounters");
  TTree* tbl = (TTree*)fin->Get("SourceFileTable");
  if (!h || !tbl) { printf("INCOMPLETE %s\n", f); return; }
  TString* path = 0;
  tbl->SetBranchAddress("path", &path);
  tbl->GetEntry(0);
  TString base = path ? gSystem->BaseName(path->Data()) : "?";
  const Double_t nIn = h->GetBinContent(1), nAcc = h->GetBinContent(2);
  const Double_t size = (Double_t)fin->GetSize();
  printf("%-42s nIn=%8.0f  nAcc=%8.0f  %9.0f B  %7.1f B/ev  %s\n", base.Data(), nIn, nAcc, size,
         nIn > 0 ? size / nIn : 0.0, base.Contains("_adc_") ? "ADC" : "plain");
  fin->Close();
}
