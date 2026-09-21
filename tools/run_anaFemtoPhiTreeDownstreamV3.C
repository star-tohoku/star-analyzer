void run_anaFemtoPhiTreeDownstreamV3(const Char_t* treeFile, const Char_t* outFile, const Char_t* configPath,
                                   Int_t bufferSize = -1, const Char_t* mixingMode = "bufferAll",
                                   const Char_t* variation = "", Double_t signalMin = -1.0,
                                   Double_t signalMax = -1.0,
                                   Bool_t recomputeDaughterPid = kFALSE,
                                   Bool_t recomputeMixBin = kFALSE,
                                   Long64_t maxEvents = -1, Int_t pairMtBins = 240,
                                   Double_t pairMtMin = 0.8, Double_t pairMtMax = 3.2,
                                   Bool_t data005DiagnosticsOnly = kFALSE,
                                   const Char_t* data006ConfigPath = "",
                                   Long64_t firstEvent = 0,
                                   Bool_t threeEventEnabled = kFALSE,
                                   Long64_t threeEventMaxRawSamplesPerEvent = 0,
                                   Long64_t threeEventSamplingSeed = 0) {
  // The process's real working directory, not $PWD: a SUMS job is a csh script that cd's
  // into its runtime bundle, and csh does not update the PWD environment variable, so the
  // libraries and macros would be looked for in the wrong place on the farm.
  // See results/pilot-farm-20260915.md.
  TString cwdStr = gSystem->WorkingDirectory();
  const char* pwd = cwdStr.Data();
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
  anaFemtoPhiTreeDownstreamV3(treeFile, outFile, configPath, bufferSize, mixingMode, variation,
                             signalMin, signalMax, recomputeDaughterPid, recomputeMixBin, maxEvents,
                             pairMtBins, pairMtMin, pairMtMax, data005DiagnosticsOnly,
                             data006ConfigPath, firstEvent, threeEventEnabled,
                             threeEventMaxRawSamplesPerEvent, threeEventSamplingSeed);
}
