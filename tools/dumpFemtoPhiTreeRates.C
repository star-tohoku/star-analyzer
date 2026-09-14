// Per-run / centrality size and candidate rates from a reduced tree. CINT-safe (ROOT 5).
void dumpFemtoPhiTreeRates(const Char_t* rootPath, const Char_t* outPath) {
  TFile* f = TFile::Open(rootPath, "READ");
  if (!f || f->IsZombie()) {
    printf("ERROR: cannot open %s\n", rootPath);
    return;
  }
  TTree* ev = (TTree*)f->Get("FemtoEventTree");
  TTree* tr = (TTree*)f->Get("FemtoTrackTree");
  if (!ev || !tr) {
    printf("ERROR: missing trees\n");
    return;
  }
  Int_t runId = 0, cent9 = 0;
  ev->SetBranchAddress("runId", &runId);
  ev->SetBranchAddress("cent9", &cent9);

  UChar_t speciesCode = 0;
  UInt_t selFlags = 0;
  tr->SetBranchAddress("speciesCode", &speciesCode);
  tr->SetBranchAddress("selFlags", &selFlags);

  const Int_t kMaxRun = 64;
  Int_t runIds[64];
  Long64_t runEv[64];
  Int_t nRun = 0;
  Long64_t nCent[16];
  Int_t ic;
  for (ic = 0; ic < 16; ++ic) nCent[ic] = 0;
  for (ic = 0; ic < kMaxRun; ++ic) {
    runIds[ic] = 0;
    runEv[ic] = 0;
  }

  Long64_t nEv = ev->GetEntries();
  Long64_t nTr = tr->GetEntries();
  Long64_t i;
  for (i = 0; i < nEv; ++i) {
    ev->GetEntry(i);
    if (cent9 >= 0 && cent9 < 16) nCent[cent9]++;
    Int_t found = -1;
    Int_t ir;
    for (ir = 0; ir < nRun; ++ir) {
      if (runIds[ir] == runId) {
        found = ir;
        break;
      }
    }
    if (found < 0) {
      if (nRun >= kMaxRun) continue;
      found = nRun;
      runIds[nRun] = runId;
      nRun++;
    }
    runEv[found]++;
  }

  Long64_t nSp[5] = {0, 0, 0, 0, 0};
  Long64_t nNomPid[5] = {0, 0, 0, 0, 0};
  Long64_t nNomFemto[5] = {0, 0, 0, 0, 0};
  Long64_t nNomBoth[5] = {0, 0, 0, 0, 0};
  for (i = 0; i < nTr; ++i) {
    tr->GetEntry(i);
    Int_t s = (Int_t)speciesCode;
    if (s < 0 || s > 4) s = 0;
    nSp[s]++;
    if (selFlags & (1u << 2)) nNomPid[s]++;
    if (selFlags & (1u << 3)) nNomFemto[s]++;
    if ((selFlags & 0xcu) == 0xcu) nNomBoth[s]++;
  }

  FILE* out = fopen(outPath, "w");
  if (!out) {
    printf("ERROR: cannot write %s\n", outPath);
    return;
  }
  Long_t id = 0, flags = 0, modtime = 0;
  Long64_t sz = 0;
  gSystem->GetPathInfo(rootPath, &id, &sz, &flags, &modtime);
  fprintf(out, "path\t%s\n", rootPath);
  fprintf(out, "size_bytes\t%lld\n", (long long)sz);
  fprintf(out, "nEvent\t%lld\n", (long long)nEv);
  fprintf(out, "nTrack\t%lld\n", (long long)nTr);
  fprintf(out, "bytes_per_accepted\t%.4f\n", nEv > 0 ? (Double_t)sz / (Double_t)nEv : 0.0);
  const char* names[5] = {"unknown", "Kp", "Km", "d", "p"};
  fprintf(out, "species\tnStore\tnomPid\tnomFemto\tnomBoth\n");
  for (ic = 0; ic < 5; ++ic) {
    fprintf(out, "%s\t%lld\t%lld\t%lld\t%lld\n", names[ic], (long long)nSp[ic], (long long)nNomPid[ic],
            (long long)nNomFemto[ic], (long long)nNomBoth[ic]);
  }
  fprintf(out, "cent9\tnEvent\n");
  for (ic = 0; ic < 16; ++ic) {
    if (nCent[ic] > 0) fprintf(out, "%d\t%lld\n", ic, (long long)nCent[ic]);
  }
  fprintf(out, "runId\tnEvent\n");
  for (ic = 0; ic < nRun; ++ic) {
    fprintf(out, "%d\t%lld\n", runIds[ic], (long long)runEv[ic]);
  }
  fclose(out);
  printf("wrote %s nEv=%lld size=%lld nP=%lld nRun=%d\n", outPath, (long long)nEv, (long long)sz,
         (long long)nSp[4], nRun);
  f->Close();
}
