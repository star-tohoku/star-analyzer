// Step 7 completeness ledger: dump, for every merged block, the SourceFileTable rows and the
// summed hCounters. Nothing here reads event data -- the tables are small and the counters are
// histogram bins -- so a 11 GB block costs a header read, not a scan.
//
// Usage: SourceLedger.C("blocklist.txt", "out_prefix")
//   writes <out_prefix>.src.tsv  : block, subjobId, sourceFileIndex, fileName
//          <out_prefix>.cnt.tsv  : block, nInputEvents ... nJobs (17 labelled counters)
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TString.h"
#include "TSystem.h"
#include <cstdio>
#include <fstream>
#include <string>

void SourceLedger(const char* blockList, const char* outPrefix) {
  std::ifstream in(blockList);
  if (!in) { printf("cannot open %s\n", blockList); return; }

  FILE* fsrc = fopen(Form("%s.src.tsv", outPrefix), "w");
  FILE* fcnt = fopen(Form("%s.cnt.tsv", outPrefix), "w");
  if (!fsrc || !fcnt) { printf("cannot write %s.*\n", outPrefix); return; }
  fprintf(fsrc, "block\tsubjobId\tsourceFileIndex\tfileName\n");

  std::string line;
  Int_t nBlock = 0, nBad = 0;
  Long64_t nRow = 0;
  Bool_t header = kFALSE;
  while (std::getline(in, line)) {
    if (line.empty()) continue;
    const TString block = gSystem->BaseName(line.c_str());
    TFile* f = TFile::Open(line.c_str());
    if (!f || f->IsZombie()) { printf("BAD %s\n", block.Data()); ++nBad; continue; }

    TTree* tbl = (TTree*)f->Get("SourceFileTable");
    if (tbl) {
      UInt_t sub = 0; UShort_t idx = 0; TString* name = 0;
      tbl->SetBranchAddress("subjobId", &sub);
      tbl->SetBranchAddress("sourceFileIndex", &idx);
      tbl->SetBranchAddress("fileName", &name);
      for (Long64_t i = 0; i < tbl->GetEntries(); ++i) {
        tbl->GetEntry(i);
        fprintf(fsrc, "%s\t%u\t%u\t%s\n", block.Data(), sub, (UInt_t)idx,
                name ? name->Data() : "");
        ++nRow;
      }
    } else {
      printf("NO SourceFileTable %s\n", block.Data());
      ++nBad;
    }

    TH1* h = (TH1*)f->Get("hCounters");
    if (h) {
      if (!header) {
        fprintf(fcnt, "block");
        for (Int_t b = 1; b <= h->GetNbinsX(); ++b)
          fprintf(fcnt, "\t%s", h->GetXaxis()->GetBinLabel(b));
        fprintf(fcnt, "\n");
        header = kTRUE;
      }
      fprintf(fcnt, "%s", block.Data());
      for (Int_t b = 1; b <= h->GetNbinsX(); ++b)
        fprintf(fcnt, "\t%lld", (Long64_t)(h->GetBinContent(b) + 0.5));
      fprintf(fcnt, "\n");
    } else {
      printf("NO hCounters %s\n", block.Data());
      ++nBad;
    }

    f->Close();
    delete f;
    ++nBlock;
    if (nBlock % 10 == 0) printf("... %d blocks, %lld source rows\n", nBlock, nRow);
  }
  fclose(fsrc);
  fclose(fcnt);
  printf("blocks %d  source rows %lld  problems %d\n", nBlock, nRow, nBad);
}
