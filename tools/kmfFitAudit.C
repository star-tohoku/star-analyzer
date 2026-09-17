// T1 of the Step 10 plan: audit every per-k* mass fit behind a kstarMassFitCF sidecar.
//
// The correlation function has one independent M(KK) fit behind every point. checkHist records a
// status for each (0 OK, 1 fit failed, 2 negative yield, 3 low statistics, 4 normalisation failed,
// 5 non-finite) and the extracted SE/ME yields, but nothing reads them back: a CF point built on a
// failed fit looks exactly like one built on a good fit. This macro reads them back.
//
// Usage: kmfFitAudit.C("<...CFkmf.root>", "out.tsv")
//   Prints every non-OK point, and a per-graph summary with the relative yield error, so that a
//   point whose error does not shrink with statistics (Step 8 found one at k* = 0.125) is visible.
#include "TFile.h"
#include "TGraphErrors.h"
#include "TKey.h"
#include "TList.h"
#include "TString.h"
#include <cstdio>
#include <map>
#include <string>
#include <vector>

namespace {
const char* kStatusName[6] = {"OK", "FIT_FAIL", "NEG_YIELD", "LOW_STAT", "NORM_FAIL", "NON_FINITE"};

TString yieldNameFor(const TString& statusName, const char* which) {
  // slice_<id>_kmf_fitstatus_<rest>  ->  slice_<id>_kmf_Y_<which>_<rest>
  TString out = statusName;
  out.ReplaceAll("kmf_fitstatus_", Form("kmf_Y_%s_", which));
  return out;
}
}  // namespace

void kmfFitAudit(const char* sidecar, const char* outTsv) {
  TFile* f = TFile::Open(sidecar);
  if (!f || f->IsZombie()) { printf("cannot open %s\n", sidecar); return; }
  FILE* out = fopen(outTsv, "w");
  if (!out) { printf("cannot write %s\n", outTsv); return; }
  fprintf(out, "graph\tkstar\tstatus\tstatusName\tY_SE\teY_SE\trelErr_SE\n");

  std::vector<TString> statusGraphs;
  TIter next(f->GetListOfKeys());
  TKey* k;
  while ((k = (TKey*)next())) {
    TString n = k->GetName();
    if (n.Contains("kmf_fitstatus_")) statusGraphs.push_back(n);
  }
  printf("fit-status graphs: %d\n", (Int_t)statusGraphs.size());

  std::map<std::string, Int_t> tally;
  Int_t nPointsTotal = 0, nBad = 0;
  for (size_t i = 0; i < statusGraphs.size(); ++i) {
    TGraphErrors* gs = (TGraphErrors*)f->Get(statusGraphs[i]);
    if (!gs) continue;
    TGraphErrors* gy = (TGraphErrors*)f->Get(yieldNameFor(statusGraphs[i], "SE"));
    for (Int_t p = 0; p < gs->GetN(); ++p) {
      Double_t x, st;
      gs->GetPoint(p, x, st);
      const Int_t s = (Int_t)(st + 0.5);
      ++nPointsTotal;
      Double_t y = 0, ey = 0;
      if (gy && p < gy->GetN()) { Double_t xx; gy->GetPoint(p, xx, y); ey = gy->GetErrorY(p); }
      const Double_t rel = (y > 0) ? ey / y : -1.0;
      if (s != 0) {
        ++nBad;
        tally[(s >= 0 && s < 6) ? kStatusName[s] : "UNKNOWN"]++;
        printf("  %-62s k*=%.3f  %s  Y=%.1f +- %.1f\n", statusGraphs[i].Data(), x,
               (s >= 0 && s < 6) ? kStatusName[s] : "UNKNOWN", y, ey);
      }
      fprintf(out, "%s\t%.4f\t%d\t%s\t%.3f\t%.3f\t%.5f\n", statusGraphs[i].Data(), x, s,
              (s >= 0 && s < 6) ? kStatusName[s] : "UNKNOWN", y, ey, rel);
    }
  }
  fclose(out);
  printf("\npoints %d   non-OK %d\n", nPointsTotal, nBad);
  for (std::map<std::string, Int_t>::const_iterator it = tally.begin(); it != tally.end(); ++it)
    printf("  %-12s %d\n", it->first.c_str(), it->second);
  printf("written: %s\n", outTsv);
}
