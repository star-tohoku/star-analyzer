void dumpFemtoPhiTreeSelFlags(const Char_t* rootPath, const Char_t* outPath) {
  TFile* f = TFile::Open(rootPath, "READ");
  if (!f || f->IsZombie()) {
    printf("ERROR: cannot open %s\n", rootPath);
    return;
  }
  TTree* tr = (TTree*)f->Get("FemtoTrackTree");
  if (!tr) {
    printf("ERROR: no FemtoTrackTree\n");
    return;
  }
  UChar_t speciesCode = 0;
  UInt_t selFlags = 0;
  tr->SetBranchAddress("speciesCode", &speciesCode);
  tr->SetBranchAddress("selFlags", &selFlags);
  Long64_t n[5] = {0, 0, 0, 0, 0};
  Long64_t nomPid[5] = {0, 0, 0, 0, 0};
  Long64_t nomFemto[5] = {0, 0, 0, 0, 0};
  Long64_t kaonCuts[5] = {0, 0, 0, 0, 0};
  Long64_t loosePid[5] = {0, 0, 0, 0, 0};
  Long64_t trackQ[5] = {0, 0, 0, 0, 0};
  const Long64_t nent = tr->GetEntries();
  for (Long64_t i = 0; i < nent; ++i) {
    tr->GetEntry(i);
    Int_t s = (Int_t)speciesCode;
    if (s < 0 || s > 4) s = 0;
    n[s]++;
    if (selFlags & (1u << 0)) trackQ[s]++;
    if (selFlags & (1u << 1)) loosePid[s]++;
    if (selFlags & (1u << 2)) nomPid[s]++;
    if (selFlags & (1u << 3)) nomFemto[s]++;
    if (selFlags & (1u << 5)) kaonCuts[s]++;
  }
  FILE* out = fopen(outPath, "w");
  if (!out) out = stdout;
  fprintf(out, "species\tn\ttrackQ\tloosePid\tnomPid\tnomFemto\tkaonCutsNom\n");
  const char* names[5] = {"unknown", "Kp", "Km", "d", "p"};
  for (int s = 0; s < 5; ++s) {
    fprintf(out, "%s\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\n", names[s], (long long)n[s], (long long)trackQ[s],
            (long long)loosePid[s], (long long)nomPid[s], (long long)nomFemto[s], (long long)kaonCuts[s]);
  }
  if (out != stdout) fclose(out);
  printf("wrote %s nTrack=%lld\n", outPath, (long long)nent);
  f->Close();
}
