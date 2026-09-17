// Pre-merge gate: which subjob outputs are complete enough to hadd.
//
// A job whose output copy was interrupted still leaves a file, and that file can still carry an
// hCounters histogram while its trees are short -- the 2026-09-15 pilot left 19 such stubs among
// 218 outputs, and merging them would have inflated the counters against the tree contents.
// Every file must satisfy: openable, has the trees, event rows == nAcceptedEvents, no pack
// saturation.
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TNamed.h"
#include "TString.h"
#include "TSystem.h"
#include <cstdio>
#include <fstream>
#include <string>

// The good-list is written as "<runNumber> <path>" so that `sort | awk` groups it for the
// per-run merge without opening every file a second time.
void FileIntegrity(const char* listFile, const char* goodListOut = "") {
  std::ifstream in(listFile);
  if (!in.is_open()) { printf("cannot open list %s\n", listFile); return; }
  FILE* good = (TString(goodListOut).Length() > 0) ? fopen(goodListOut, "w") : 0;
  std::string path;
  Int_t nOk = 0, nBad = 0, nMultiRun = 0;
  Long64_t sumAccepted = 0, sumInput = 0, sumBytes = 0;
  while (std::getline(in, path)) {
    if (path.empty()) continue;
    TFile* f = TFile::Open(path.c_str());
    TString why;
    if (!f || f->IsZombie()) why = "cannot open";
    else {
      TTree* ev = (TTree*)f->Get("FemtoEventTree");
      TTree* tbl = (TTree*)f->Get("SourceFileTable");
      TH1* h = (TH1*)f->Get("hCounters");
      TNamed* sat = (TNamed*)f->Get("packSaturationTotal");
      if (!ev) why = "no FemtoEventTree";
      else if (!tbl) why = "no SourceFileTable";
      else if (!h) why = "no hCounters";
      else if (ev->GetEntries() != (Long64_t)h->GetBinContent(2))
        why = TString::Format("event rows %lld != nAcceptedEvents %.0f", ev->GetEntries(),
                              h->GetBinContent(2));
      else if (sat && TString(sat->GetTitle()) != "0")
        why = TString::Format("packSaturationTotal %s", sat->GetTitle());
      if (why.Length() == 0) {
        ++nOk;
        sumAccepted += ev->GetEntries();
        sumInput += (Long64_t)h->GetBinContent(1);
        sumBytes += f->GetSize();
        // Emit the run number with the path, so that one pass over the outputs serves both the
        // integrity gate and the grouping the per-run merge needs. The run is taken from the
        // PicoDst basename the job recorded, e.g. st_physics_22166015_raw_4000002.picoDst.root;
        // a job that read files from more than one run is reported rather than guessed at.
        TString* srcName = 0;
        tbl->SetBranchAddress("fileName", &srcName);
        std::string run, mixed;
        for (Long64_t r = 0; r < tbl->GetEntries(); ++r) {
          tbl->GetEntry(r);
          TString n = srcName ? *srcName : TString("");
          Ssiz_t p1 = n.Index("_2");                       // run numbers here start with 2
          TString cand = (p1 > 0) ? TString(n(p1 + 1, 8)) : TString("");
          if (!cand.IsDigit()) continue;
          if (run.empty()) run = cand.Data();
          else if (run != cand.Data()) mixed = cand.Data();
        }
        if (!mixed.empty()) ++nMultiRun;
        if (good) fprintf(good, "%s %s\n", run.empty() ? "unknown" : run.c_str(), path.c_str());
      }
    }
    if (why.Length() > 0) {
      ++nBad;
      printf("  BAD  %-90s %s\n", gSystem->BaseName(path.c_str()), why.Data());
    }
    if (f) f->Close();
  }
  if (good) fclose(good);
  printf("\n%d complete, %d incomplete\n", nOk, nBad);
  // SUMS packs whatever files it picked into a subjob, without regard to run boundaries, so a
  // large fraction of the outputs hold events from two runs. That is why the merge is by block of
  // subjobs and not by run: an output cannot be assigned to one run.
  printf("%d of the complete files span more than one run (%.0f%%)\n", nMultiRun,
         nOk ? 100.0 * nMultiRun / nOk : 0.0);
  printf("complete files: %lld input events, %lld accepted, %.3f GB, %.1f B/input event\n",
         sumInput, sumAccepted, sumBytes / 1e9, sumInput ? (Double_t)sumBytes / sumInput : 0.0);
  if (nBad) printf("\nRe-run the %d incomplete subjobs before merging; do not hadd them.\n", nBad);
}
