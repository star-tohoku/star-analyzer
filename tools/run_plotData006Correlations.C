void run_plotData006Correlations(const Char_t* packageDirectory, const Char_t* data006Config,
                                 const Char_t* outputDirectory, Int_t displayRebin = 4,
                                 Double_t cfYMin = 0.0, Double_t cfYMax = 2.5,
                                 const Char_t* runLabel = "DATA-006 TEST",
                                 const Char_t* filePrefix = "data006_test") {
  TString cwdStr = gSystem->WorkingDirectory();
  const char* pwd = cwdStr.Data();
  gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
  loadSharedLibraries();
  if (gSystem->Load(TString(pwd) + "/lib/libStarAnaConfig.so") < 0) return;
  gInterpreter->AddIncludePath(pwd);
  gInterpreter->AddIncludePath(TString::Format("%s/include", pwd));
  gSystem->AddLinkedLibs(TString::Format("-L%s/lib -lStarAnaConfig -Wl,-rpath,%s/lib", pwd, pwd));
  gROOT->ProcessLine(TString::Format(".L %s/tools/plotData006Correlations.C+", pwd));
  plotData006Correlations(packageDirectory, data006Config, outputDirectory,
                          displayRebin, cfYMin, cfYMax, runLabel, filePrefix);
}
