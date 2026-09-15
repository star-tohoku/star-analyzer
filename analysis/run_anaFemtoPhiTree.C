// run_anaFemtoPhiTree.C - load libs and call anaFemtoPhiTree
// Usage: root4star -b -q 'run_anaFemtoPhiTree.C("input.list","output.root","0",100)'

void run_anaFemtoPhiTree(const Char_t* inputFile,
                         const Char_t* outputFile,
                         const Char_t* jobid = "0",
                         Long64_t nEventsMax = -1,
                         const Char_t* configPath = 0,
                         Long64_t nSkip = 0)
{
  // The process's real working directory, not $PWD: a SUMS job is a csh script that cd's
  // into its runtime bundle, and csh does not update the PWD environment variable, so the
  // libraries and macros would be looked for in the wrong place on the farm.
  // See results/pilot-farm-20260915.md.
  TString cwdStr = gSystem->WorkingDirectory();
  const char* pwd = cwdStr.Data();

  gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
  loadSharedLibraries();
  gSystem->Load("StPicoEvent");
  gSystem->Load("StPicoDstMaker");

  if (gSystem->Load(TString(pwd) + "/lib/libStarAnaConfig.so") < 0) {
    std::cerr << "ERROR: failed to load libStarAnaConfig.so" << std::endl;
    return;
  }
  if (gSystem->Load(TString(pwd) + "/lib/libStRefMultCorr.so") < 0) {
    std::cerr << "ERROR: failed to load libStRefMultCorr.so" << std::endl;
    return;
  }
  if (gSystem->Load(TString(pwd) + "/lib/libStCommon.so") < 0) {
    std::cerr << "ERROR: failed to load libStCommon.so" << std::endl;
    return;
  }
  if (gSystem->Load(TString(pwd) + "/lib/libStFemtoPhiTreeMaker.so") < 0) {
    std::cerr << "ERROR: failed to load libStFemtoPhiTreeMaker.so" << std::endl;
    return;
  }

  gInterpreter->AddIncludePath(pwd);
  gInterpreter->AddIncludePath(TString::Format("%s/include", pwd));
  gInterpreter->AddIncludePath(TString::Format("%s/StMaker/common", pwd));
  gInterpreter->AddIncludePath(TString::Format("%s/StMaker/StFemtoPhiTreeMaker", pwd));
  gInterpreter->AddIncludePath("$STAR/StRoot");
  gSystem->AddLinkedLibs(TString::Format(
      "-L%s/lib -lStarAnaConfig -lStRefMultCorr -lStCommon -lStFemtoPhiTreeMaker -Wl,-rpath,%s/lib", pwd, pwd));

  gROOT->ProcessLine(TString::Format(".L %s/analysis/anaFemtoPhiTree.C+", pwd));
  anaFemtoPhiTree(inputFile, outputFile, jobid, nEventsMax, configPath, nSkip);
}
