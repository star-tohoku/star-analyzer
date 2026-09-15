// The effective provenance test the plan's Step 2 was missing: not "can the table be read", but
// "does (subjobId, sourceFileIndex, sourceEntry) reopen the same event in the PicoDst".
#include "TFile.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TString.h"
#include "TRandom3.h"
#include <cstdio>
#include <map>
#include <string>

void ProvCheck(const char* treeFile, Int_t nSample = 100) {
  TFile* f = TFile::Open(treeFile);
  if (!f) { printf("cannot open %s\n", treeFile); return; }
  TTree* tbl = (TTree*)f->Get("SourceFileTable");
  TTree* ev = (TTree*)f->Get("FemtoEventTree");
  if (!tbl || !ev) { printf("missing SourceFileTable or FemtoEventTree\n"); return; }

  UInt_t tSub = 0; UShort_t tIdx = 0; TString* tPath = 0;
  const Bool_t hasSubjob = (tbl->GetBranch("subjobId") != 0);
  if (hasSubjob) tbl->SetBranchAddress("subjobId", &tSub);
  tbl->SetBranchAddress("sourceFileIndex", &tIdx);
  tbl->SetBranchAddress("path", &tPath);
  std::map<std::string, std::string> pathOf;   // "subjob:index" -> path
  for (Long64_t i = 0; i < tbl->GetEntries(); ++i) {
    tbl->GetEntry(i);
    pathOf[Form("%u:%u", hasSubjob ? tSub : 0u, (UInt_t)tIdx)] = tPath ? tPath->Data() : "";
  }
  printf("SourceFileTable: %lld rows, subjobId branch %s\n", tbl->GetEntries(),
         hasSubjob ? "PRESENT" : "ABSENT");

  // schema 3 packs the identity: there is no runId / eventId branch, only
  // eventUID = (runId << 32) | eventId.
  // schema 3 writes sourceEntry as Int_t ("sourceEntry/I"); asking ROOT for a Long64_t
  // makes SetBranchAddress refuse the address and the variable stays 0.
  ULong64_t uid = 0; UInt_t subjob = 0; UShort_t fidx = 0; Int_t sentry = 0;
  ev->SetBranchAddress("eventUID", &uid);
  ev->SetBranchAddress("subjobId", &subjob);
  ev->SetBranchAddress("sourceFileIndex", &fidx);
  ev->SetBranchAddress("sourceEntry", &sentry);

  TRandom3 rng(20260915);
  const Long64_t n = ev->GetEntries();
  Int_t ok = 0, bad = 0, noFile = 0;
  Long64_t maxEntry = 0;
  std::map<std::string, TFile*> open;
  for (Int_t s = 0; s < nSample; ++s) {
    ev->GetEntry((Long64_t)(rng.Rndm() * n));
    const Int_t runId = (Int_t)(uid >> 32);
    const Int_t eventId = (Int_t)(uid & 0xFFFFFFFFull);
    if ((Long64_t)sentry > maxEntry) maxEntry = sentry;
    std::map<std::string, std::string>::const_iterator it =
        pathOf.find(Form("%u:%u", subjob, (UInt_t)fidx));
    if (it == pathOf.end()) { ++noFile; continue; }
    TFile*& pf = open[it->second];
    if (!pf) pf = TFile::Open(it->second.c_str());
    if (!pf || pf->IsZombie()) { ++noFile; continue; }
    TTree* pico = (TTree*)pf->Get("PicoDst");
    if (!pico) { ++noFile; continue; }
    if (sentry < 0 || sentry >= pico->GetEntries()) { ++bad; continue; }
    pico->GetEntry(sentry);
    TLeaf* lr = pico->GetLeaf("Event.mRunId");
    TLeaf* le = pico->GetLeaf("Event.mEventId");
    if (!lr || !le) { printf("  PicoDst leaves not found\n"); return; }
    const Int_t pRun = (Int_t)lr->GetValue(0), pEvt = (Int_t)le->GetValue(0);
    if (pRun == runId && pEvt == eventId) ++ok;
    else {
      if (bad < 3)
        printf("  MISMATCH tree(run %d, evt %d) vs pico entry %d (run %d, evt %d) in %s\n", runId,
               eventId, sentry, pRun, pEvt, it->second.c_str());
      ++bad;
    }
  }
  printf("sampled %d events: %d reopened correctly, %d mismatched, %d unresolvable\n", nSample, ok,
         bad, noFile);
  printf("max sourceEntry seen in the sample = %lld\n", (Long64_t)maxEntry);
  f->Close();
}
