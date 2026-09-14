// Classify tree nom-d tracks that StFemtoMaker would not store.
// Does not modify StFemtoMaker. CINT-safe (ROOT 5).
// Classes:
//   ndedx_lt_nuclear  — maker requires nHitsDedx >= nuclearid minNHitsDedxNuclear before IsDeuteron
//   proton_continue   — PassProtonCuts && !IsProton then continue (p>=TOF threshold failing m2)
//   both              — both of the above
// Cut numbers are written in the TSV header and must match YAML.
void diagFemtoPhiTreeDeuteronDiff(const Char_t* rootPath, const Char_t* outPath,
                                  Double_t nSigmaProtonMax = 2.0,
                                  Double_t protonTofP = 2.0,
                                  Double_t protonMinMass2 = 0.6,
                                  Double_t protonMaxMass2 = 1.2,
                                  Int_t minNHitsDedxNuclear = 15) {
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
  ULong64_t eventUID = 0;
  Int_t trackIndex = -1;
  UChar_t speciesCode = 0;
  UInt_t selFlags = 0;
  Char_t tofMatch = 0;
  Short_t nHitsDedx = 0;
  Float_t px = 0, py = 0, pz = 0, mass2 = 0, nSigmaProton = 0, nSigmaDeuteron = 0, dca = 0;
  tr->SetBranchAddress("eventUID", &eventUID);
  tr->SetBranchAddress("trackIndex", &trackIndex);
  tr->SetBranchAddress("speciesCode", &speciesCode);
  tr->SetBranchAddress("selFlags", &selFlags);
  tr->SetBranchAddress("tofMatch", &tofMatch);
  tr->SetBranchAddress("nHitsDedx", &nHitsDedx);
  tr->SetBranchAddress("px", &px);
  tr->SetBranchAddress("py", &py);
  tr->SetBranchAddress("pz", &pz);
  tr->SetBranchAddress("mass2", &mass2);
  tr->SetBranchAddress("nSigmaProton", &nSigmaProton);
  tr->SetBranchAddress("nSigmaDeuteron", &nSigmaDeuteron);
  tr->SetBranchAddress("dca", &dca);

  const UInt_t kTrackQ = 1u << 0;
  const UInt_t kNeedD = (1u << 2) | (1u << 3);

  FILE* out = fopen(outPath, "w");
  if (!out) {
    printf("ERROR: cannot write %s\n", outPath);
    return;
  }
  fprintf(out, "# source=%s\n", rootPath);
  fprintf(out, "# nSigmaProtonMax=%.3f protonTofP=%.3f mass2=[%.3f,%.3f] minNHitsDedxNuclear=%d\n",
          nSigmaProtonMax, protonTofP, protonMinMass2, protonMaxMass2, minNHitsDedxNuclear);
  fprintf(out, "eventUID\ttrackIndex\tnHitsDedx\tpmom\tnSigmaProton\tnSigmaDeuteron\tdca\ttofMatch\tmass2\tclass\n");

  Long64_t nNomD = 0, nNdedx = 0, nPcont = 0, nBoth = 0, nOther = 0, nNoQ = 0;
  const Long64_t nent = tr->GetEntries();
  Long64_t i;
  for (i = 0; i < nent; ++i) {
    tr->GetEntry(i);
    if (speciesCode != 3) continue;
    if ((selFlags & kNeedD) != kNeedD) continue;
    nNomD++;
    const Double_t pmom = TMath::Sqrt((Double_t)px * px + (Double_t)py * py + (Double_t)pz * pz);
    const Bool_t passTrackQ = (selFlags & kTrackQ) != 0;
    const Bool_t passProtonCuts = passTrackQ && (TMath::Abs(nSigmaProton) <= nSigmaProtonMax);
    Bool_t isProton = kTRUE;
    if (pmom >= protonTofP) {
      isProton = tofMatch && mass2 >= protonMinMass2 && mass2 <= protonMaxMass2;
    }
    const Bool_t noTrackQ = (selFlags & kTrackQ) == 0;
    const Bool_t ndedx = (nHitsDedx < minNHitsDedxNuclear);
    const Bool_t pcont = passProtonCuts && !isProton;
    const char* cls = "unclassified";
    if (noTrackQ) {
      cls = "no_trackQ";
      nNoQ++;
    } else if (ndedx && pcont) {
      cls = "both";
      nBoth++;
    } else if (ndedx) {
      cls = "ndedx_lt_nuclear";
      nNdedx++;
    } else if (pcont) {
      cls = "proton_continue";
      nPcont++;
    } else {
      cls = "unclassified";
      nOther++;
    }
    if (strcmp(cls, "unclassified") != 0) {
      fprintf(out, "%llu\t%d\t%d\t%.6f\t%.4f\t%.4f\t%.4f\t%d\t%.4f\t%s\n",
              (unsigned long long)eventUID, trackIndex, (int)nHitsDedx, pmom, nSigmaProton, nSigmaDeuteron,
              dca, (int)tofMatch, mass2, cls);
    }
  }
  const Long64_t nExpl = nNoQ + nNdedx + nPcont + nBoth;
  fprintf(out, "# nNomD=%lld nNoTrackQ=%lld nNdedx=%lld nPcont=%lld nBoth=%lld nUnclassified=%lld nExplained=%lld\n",
          (long long)nNomD, (long long)nNoQ, (long long)nNdedx, (long long)nPcont, (long long)nBoth,
          (long long)nOther, (long long)nExpl);
  fclose(out);
  printf("nNomD=%lld noQ=%lld ndedx=%lld pcont=%lld both=%lld other=%lld explained=%lld wrote %s\n",
         (long long)nNomD, (long long)nNoQ, (long long)nNdedx, (long long)nPcont, (long long)nBoth,
         (long long)nOther, (long long)nExpl, outPath);
  f->Close();
}
