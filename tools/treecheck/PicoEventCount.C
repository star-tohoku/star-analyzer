// Count the events actually present in a PicoDst, for the Step 7 ledger.
//
// The catalog's "events" column is what SUMS wrote into every per-process .list, and it is the
// only expected-event number available offline. Where the production's own counters disagree with
// it, the ledger has to know whether the file was short or the job was. This macro answers the
// first half: open each PicoDst and report its entry count next to the catalog number.
//
// Needs xrootd, so it runs under root4star (script/singularity_treecheck.sh), not the container's
// bare ROOT.
//
// Usage: PicoEventCount.C("urls.txt", "out.tsv")
//   urls.txt: one "<url> <catalogEvents>" per line, exactly as the SUMS .list files are written.
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>

void PicoEventCount(const char* urlList, const char* outFile) {
  std::ifstream in(urlList);
  if (!in) { printf("cannot open %s\n", urlList); return; }
  FILE* out = fopen(outFile, "w");
  if (!out) { printf("cannot write %s\n", outFile); return; }
  fprintf(out, "fileName\tcatalogEvents\tentries\tstatus\n");

  std::string line;
  Long64_t nOk = 0, nFail = 0, sumCat = 0, sumEnt = 0;
  while (std::getline(in, line)) {
    if (line.empty()) continue;
    std::istringstream is(line);
    std::string url; Long64_t cat = -1;
    is >> url >> cat;
    const TString base = gSystem->BaseName(url.c_str());
    TFile* f = TFile::Open(url.c_str());
    if (!f || f->IsZombie()) {
      fprintf(out, "%s\t%lld\t-1\topen_failed\n", base.Data(), cat);
      ++nFail;
      if (f) delete f;
      continue;
    }
    TTree* t = (TTree*)f->Get("PicoDst");
    if (!t) {
      fprintf(out, "%s\t%lld\t-1\tno_PicoDst\n", base.Data(), cat);
      ++nFail;
    } else {
      const Long64_t n = t->GetEntries();
      fprintf(out, "%s\t%lld\t%lld\t%s\n", base.Data(), cat, n, (cat == n) ? "equal" : "differs");
      ++nOk; sumCat += cat; sumEnt += n;
    }
    f->Close();
    delete f;
  }
  fclose(out);
  printf("read %lld  failed %lld  catalog %lld  entries %lld  diff %lld\n",
         nOk, nFail, sumCat, sumEnt, sumEnt - sumCat);
}
