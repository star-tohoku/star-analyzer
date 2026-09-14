void dumpFemtoPhiTreeInfo(const Char_t* rootPath, const Char_t* outDir) {
  TString dir(outDir);
  gSystem->mkdir(dir, kTRUE);
  TFile* f = TFile::Open(rootPath, "READ");
  if (!f || f->IsZombie()) {
    printf("ERROR: cannot open %s\n", rootPath);
    return;
  }

  FILE* schema = fopen(Form("%s/tree-schema.txt", outDir), "w");
  FILE* sizeTsv = fopen(Form("%s/branch-size.tsv", outDir), "w");
  FILE* files = fopen(Form("%s/root-files.tsv", outDir), "w");
  fprintf(sizeTsv, "file\ttree\tbranch\ttype\ttotBytes\tzipBytes\tentries\n");
  fprintf(files, "path\tsize_bytes\tkeys\tevent_entries\ttrack_entries\tpair_entries\n");
  fprintf(schema, "file %s\n", rootPath);

  const char* tnames[3] = {"FemtoEventTree", "FemtoTrackTree", "FemtoPhiPairTree"};
  Long64_t nEnt[3] = {0, 0, 0};
  for (int it = 0; it < 3; ++it) {
    TTree* t = (TTree*)f->Get(tnames[it]);
    if (!t) continue;
    nEnt[it] = t->GetEntries();
    fprintf(schema, "tree %s title=\"%s\" entries=%lld totBytes=%lld zipBytes=%lld\n",
            t->GetName(), t->GetTitle(), (long long)t->GetEntries(),
            (long long)t->GetTotBytes(), (long long)t->GetZipBytes());
    if (it == 0) {
      fprintf(schema, "sourceEntry_min=%.0f sourceEntry_max=%.0f eventId_min=%.0f eventId_max=%.0f\n",
              t->GetMinimum("sourceEntry"), t->GetMaximum("sourceEntry"),
              t->GetMinimum("eventId"), t->GetMaximum("eventId"));
    }
    TObjArray* branches = t->GetListOfBranches();
    if (!branches) continue;
    for (int i = 0; i < branches->GetEntries(); ++i) {
      TBranch* br = (TBranch*)branches->At(i);
      if (!br) continue;
      const char* typ = "?";
      TObjArray* leaves = br->GetListOfLeaves();
      if (leaves && leaves->GetEntries() > 0) {
        TLeaf* leaf = (TLeaf*)leaves->At(0);
        if (leaf) typ = leaf->GetTypeName();
      }
      fprintf(schema, "  branch %s type=%s title=\"%s\"\n", br->GetName(), typ, br->GetTitle());
      fprintf(sizeTsv, "%s\t%s\t%s\t%s\t%lld\t%lld\t%lld\n", rootPath, t->GetName(), br->GetName(), typ,
              (long long)br->GetTotBytes(), (long long)br->GetZipBytes(), (long long)t->GetEntries());
    }
  }

  TList* keys = f->GetListOfKeys();
  if (keys) {
    for (int i = 0; i < keys->GetEntries(); ++i) {
      TKey* key = (TKey*)keys->At(i);
      if (!key) continue;
      fprintf(schema, "key %s class=%s title=\"%s\"\n", key->GetName(), key->GetClassName(), key->GetTitle());
    }
  }

  Long_t id = 0, flags = 0, modtime = 0;
  Long64_t sz = 0;
  gSystem->GetPathInfo(rootPath, &id, &sz, &flags, &modtime);
  fprintf(files, "%s\t%lld\t%d\t%lld\t%lld\t%lld\n", rootPath, (long long)sz,
          keys ? keys->GetEntries() : 0, (long long)nEnt[0], (long long)nEnt[1], (long long)nEnt[2]);
  fclose(schema);
  fclose(sizeTsv);
  fclose(files);
  printf("size_bytes=%lld nEvent=%lld nTrack=%lld nPair=%lld\n", (long long)sz, (long long)nEnt[0],
         (long long)nEnt[1], (long long)nEnt[2]);
  f->Close();
}
