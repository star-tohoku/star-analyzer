void run_anaPhiMesicNucleus(const Char_t* inputFile,
                            const Char_t* outputFile,
                            const Char_t* jobid = "0",
                            Long64_t nEventsMax = -1,
                            const Char_t* configPath = 0) {
  // The process's real working directory, not $PWD: csh (every SUMS job) does not update
  // the PWD environment variable on cd, so on the farm it named the wrong directory.
  // See results/pilot-farm-20260915.md.
  TString cwd = gSystem->WorkingDirectory();
  const char* pwd = cwd.Data();

  gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
  loadSharedLibraries();
  gSystem->Load("StPicoEvent");
  gSystem->Load("StPicoDstMaker");

  if (gSystem->Load(TString(pwd) + "/lib/libStarAnaConfig.so") < 0) return;
  if (gSystem->Load(TString(pwd) + "/lib/libStRefMultCorr.so") < 0) return;
  if (gSystem->Load(TString(pwd) + "/lib/libStPhiMesicNucleusMaker.so") < 0) return;

  gInterpreter->AddIncludePath(pwd);
  gInterpreter->AddIncludePath(TString::Format("%s/include", pwd));
  gInterpreter->AddIncludePath(TString::Format("%s/StMaker/common", pwd));
  gInterpreter->AddIncludePath("$STAR/StRoot");
  gSystem->AddLinkedLibs(
      TString::Format("-L%s/lib -lStarAnaConfig -lStRefMultCorr -lStPhiMesicNucleusMaker -Wl,-rpath,%s/lib", pwd, pwd));

  gROOT->ProcessLine(TString::Format(".L %s/analysis/anaPhiMesicNucleus.C+", pwd));
  anaPhiMesicNucleus(inputFile, outputFile, jobid, nEventsMax, configPath);
}
