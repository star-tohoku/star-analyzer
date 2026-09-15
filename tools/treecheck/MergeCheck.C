// What only a merged, multi-subjob tree can be asked:
//   - is subjobId actually unique per subjob (it is a 32-bit hash of the SUMS $JOBID string)
//   - does (subjobId, sourceFileIndex) still name one PicoDst after hadd
//   - are there duplicate events
//   - did hCounters sum
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TString.h"
#include <cstdio>
#include <map>
#include <set>
#include <string>

void MergeCheck(const char* mergedFile) {
  TFile* f = TFile::Open(mergedFile);
  if (!f) { printf("cannot open %s\n", mergedFile); return; }

  TH1* h = (TH1*)f->Get("hCounters");
  if (h) {
    printf("=== hCounters (summed) ===\n");
    for (Int_t i = 1; i <= h->GetNbinsX(); ++i)
      printf("  %-18s %14.0f\n", h->GetXaxis()->GetBinLabel(i), h->GetBinContent(i));
  } else printf("hCounters MISSING\n");

  TTree* tbl = (TTree*)f->Get("SourceFileTable");
  std::set<UInt_t> subjobs;
  std::map<std::string, std::string> pathOf;
  Int_t ambiguous = 0;
  if (tbl) {
    UInt_t sub = 0; UShort_t idx = 0; TString* path = 0;
    tbl->SetBranchAddress("subjobId", &sub);
    tbl->SetBranchAddress("sourceFileIndex", &idx);
    tbl->SetBranchAddress("path", &path);
    for (Long64_t i = 0; i < tbl->GetEntries(); ++i) {
      tbl->GetEntry(i);
      subjobs.insert(sub);
      const std::string key = Form("%u:%u", sub, (UInt_t)idx);
      const std::string val = path ? path->Data() : "";
      std::map<std::string, std::string>::iterator it = pathOf.find(key);
      if (it == pathOf.end()) pathOf[key] = val;
      else if (it->second != val) {
        if (ambiguous < 3)
          printf("  AMBIGUOUS key %s -> '%s' and '%s'\n", key.c_str(), it->second.c_str(),
                 val.c_str());
        ++ambiguous;
      }
    }
    printf("\n=== SourceFileTable ===\n");
    printf("  rows                  %lld\n", tbl->GetEntries());
    printf("  distinct subjobId     %d\n", (Int_t)subjobs.size());
    printf("  distinct (subjob,idx) %d\n", (Int_t)pathOf.size());
    printf("  ambiguous keys        %d   %s\n", ambiguous,
           ambiguous == 0 ? "(one PicoDst per key)" : "*** COLLISION ***");
  } else printf("SourceFileTable MISSING\n");

  TTree* ev = (TTree*)f->Get("FemtoEventTree");
  if (ev) {
    ULong64_t uid = 0; UInt_t sub = 0;
    ev->SetBranchAddress("eventUID", &uid);
    ev->SetBranchAddress("subjobId", &sub);
    std::set<ULong64_t> seen;
    std::set<UInt_t> evSubjobs;
    Long64_t dup = 0;
    const Long64_t n = ev->GetEntries();
    for (Long64_t i = 0; i < n; ++i) {
      ev->GetEntry(i);
      if (!seen.insert(uid).second) ++dup;
      evSubjobs.insert(sub);
    }
    printf("\n=== FemtoEventTree ===\n");
    printf("  events                %lld\n", n);
    printf("  distinct eventUID     %d\n", (Int_t)seen.size());
    printf("  duplicate events      %lld   %s\n", dup, dup == 0 ? "(none)" : "*** DUPLICATES ***");
    printf("  distinct subjobId     %d\n", (Int_t)evSubjobs.size());
  } else printf("FemtoEventTree MISSING\n");
  f->Close();
}
