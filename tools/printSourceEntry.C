void printSourceEntry(const char* path) {
  TFile* f = TFile::Open(path, "READ");
  if (!f || f->IsZombie()) { printf("fail %s\n", path); return; }
  TTree* t = (TTree*)f->Get("FemtoEventTree");
  printf("%s entries=%lld sourceEntry min=%.0f max=%.0f eventId min=%.0f max=%.0f\n",
         path, t->GetEntries(), t->GetMinimum("sourceEntry"), t->GetMaximum("sourceEntry"),
         t->GetMinimum("eventId"), t->GetMaximum("eventId"));
  f->Close();
}
