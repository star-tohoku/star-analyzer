void run_anaFemtoPhiTreeDownstreamV3(const Char_t* treeFile, const Char_t* outFile, const Char_t* configPath,
                                   Int_t bufferSize = -1, const Char_t* mixingMode = "bufferAll",
                                   Double_t nSigmaDeuteronMax = -1.0, Double_t dcaDeuteronMax = -1.0,
                                   Double_t signalMin = -1.0, Double_t signalMax = -1.0,
                                   Int_t minNHitsDedxNuclear = -1,
                                   Bool_t recomputeDaughterPid = kFALSE,
                                   Bool_t recomputeMixBin = kFALSE) {
  const char* pwd = gSystem->Getenv("PWD");
  if (!pwd) pwd = ".";
  gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
  loadSharedLibraries();
  if (gSystem->Load(TString(pwd) + "/lib/libStarAnaConfig.so") < 0) return;
  if (gSystem->Load(TString(pwd) + "/lib/libStRefMultCorr.so") < 0) return;
  if (gSystem->Load(TString(pwd) + "/lib/libStCommon.so") < 0) return;
  gInterpreter->AddIncludePath(pwd);
  gInterpreter->AddIncludePath(TString::Format("%s/include", pwd));
  gInterpreter->AddIncludePath(TString::Format("%s/StMaker/common", pwd));
  gSystem->AddLinkedLibs(TString::Format("-L%s/lib -lStarAnaConfig -lStRefMultCorr -lStCommon -Wl,-rpath,%s/lib", pwd, pwd));
  gROOT->ProcessLine(TString::Format(".L %s/tools/anaFemtoPhiTreeDownstreamV3.C+", pwd));
  anaFemtoPhiTreeDownstreamV3(treeFile, outFile, configPath, bufferSize, mixingMode, nSigmaDeuteronMax,
                            dcaDeuteronMax, signalMin, signalMax, minNHitsDedxNuclear,
                            recomputeDaughterPid, recomputeMixBin);
}
