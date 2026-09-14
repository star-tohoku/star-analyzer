void run_anaFemtoPhiTreeDownstream(const Char_t* treeFile, const Char_t* outFile, const Char_t* configPath,
                                   Int_t bufferSize = -1, const Char_t* mixingMode = "bufferAll",
                                   Double_t nSigmaDeuteronMax = -1.0, Double_t dcaDeuteronMax = -1.0,
                                   Double_t signalMin = -1.0, Double_t signalMax = -1.0) {
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
  gROOT->ProcessLine(TString::Format(".L %s/tools/anaFemtoPhiTreeDownstream.C+", pwd));
  anaFemtoPhiTreeDownstream(treeFile, outFile, configPath, bufferSize, mixingMode, nSigmaDeuteronMax,
                            dcaDeuteronMax, signalMin, signalMax);
}
