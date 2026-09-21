void run_finalizeData006(const Char_t* mergedInput, const Char_t* data006Config,
                         const Char_t* outputDirectory, const Char_t* productionTag,
                         const Char_t* gitCommit) {
  TString cwdStr = gSystem->WorkingDirectory();
  const char* pwd = cwdStr.Data();
  gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
  loadSharedLibraries();
  if (gSystem->Load(TString(pwd) + "/lib/libStarAnaConfig.so") < 0) return;
  gInterpreter->AddIncludePath(pwd);
  gInterpreter->AddIncludePath(TString::Format("%s/include", pwd));
  gSystem->AddLinkedLibs(TString::Format("-L%s/lib -lStarAnaConfig -Wl,-rpath,%s/lib", pwd, pwd));
  gROOT->ProcessLine(TString::Format(".L %s/tools/finalizeData006.C+", pwd));
  finalizeData006(mergedInput, data006Config, outputDirectory, productionTag, gitCommit);
  gROOT->ProcessLine(TString::Format(".L %s/tools/validateData006Package.C+", pwd));
  if (!validateData006Package(outputDirectory)) gSystem->Exit(1);
}
