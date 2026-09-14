void printFemtoPhiTreeClosureStats(const Char_t* path, const Char_t* tag) {
  TFile* f = TFile::Open(path, "READ");
  if (!f || f->IsZombie()) {
    printf("fail %s\n", path);
    return;
  }
  printf("=== %s %s ===\n", tag, path);
  const char* names[] = {
    "hN", "hNKaonPlus", "hNKaonMinus", "hDeuteron_NCand", "hPhiPair_NPairs_tofStrict",
    "hPhiPair_Mass_tofStrict", "hPhiPair_Mass_stage0", "hMkk",
    "hKstarSE_phi_deuteron_signal", "hKstarME_phi_deuteron_signal",
    "hKstarSE_phi_deuteron_leftSB", "hKstarME_phi_deuteron_leftSB",
    "hNPhi", "hND", "hSameEventME", "hRejectShared", "hFlow", 0
  };
  for (int i = 0; names[i]; ++i) {
    TH1* h = (TH1*)f->Get(names[i]);
    if (!h) {
      printf("MISSING %s\n", names[i]);
      continue;
    }
    printf("%s entries=%.0f integral=%.4f mean=%.6f\n", names[i], h->GetEntries(), h->Integral(), h->GetMean());
  }
  TH1* hM = (TH1*)f->Get("hMkk");
  if (!hM) hM = (TH1*)f->Get("hPhiPair_Mass_tofStrict");
  if (hM) {
    const Int_t b1 = hM->FindBin(1.012001);
    const Int_t b2 = hM->FindBin(1.025999);
    printf("MKK_signal_1.012_1.026 integral=%.4f\n", hM->Integral(b1, b2));
  }
  TNamed* ns = (TNamed*)f->Get("nSE");
  TNamed* nm = (TNamed*)f->Get("nME");
  TNamed* np = (TNamed*)f->Get("nPhi");
  TNamed* nd = (TNamed*)f->Get("nD");
  TNamed* nSame = (TNamed*)f->Get("nSameEventME");
  if (ns) printf("meta nSE=%s nME=%s nPhi=%s nD=%s nSameEventME=%s\n", ns->GetTitle(),
                 nm ? nm->GetTitle() : "-", np ? np->GetTitle() : "-", nd ? nd->GetTitle() : "-",
                 nSame ? nSame->GetTitle() : "-");
  f->Close();
}
