void run_anaFemtoPhiTreeDownstream(const Char_t* treeFile, const Char_t* outFile, const Char_t* configPath,
                                   Int_t bufferSize = -1, const Char_t* mixingMode = "bufferAll",
                                   Double_t nSigmaDeuteronMax = -1.0, Double_t dcaDeuteronMax = -1.0,
                                   Double_t signalMin = -1.0, Double_t signalMax = -1.0) {
  // The process's real working directory, not $PWD: csh (every SUMS job) does not update
  // the PWD environment variable on cd, so on the farm it named the wrong directory.
  // See results/pilot-farm-20260915.md.
  TString cwd = gSystem->WorkingDirectory();
  const char* pwd = cwd.Data();
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
