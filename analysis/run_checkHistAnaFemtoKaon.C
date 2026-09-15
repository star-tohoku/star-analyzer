// run_checkHistAnaFemtoKaon.C - Load lib and compile checkHist macro
// Usage: root4star -b -q 'run_checkHistAnaFemtoKaon.C("rootfile/...","anaName","config/mainconf/...")'

void run_checkHistAnaFemtoKaon(const Char_t* rootFile,
                               const Char_t* anaName,
                               const Char_t* mainconfPath = 0)
{
  // The process's real working directory, not $PWD (see results/pilot-farm-20260915.md).
  TString cwd = gSystem->WorkingDirectory();
  const char* pwd = cwd.Data();

  if (gSystem->Load(TString(pwd) + "/lib/libStarAnaConfig.so") < 0) {
    std::cerr << "ERROR: failed to load libStarAnaConfig.so" << std::endl;
    return;
  }

  gInterpreter->AddIncludePath(pwd);
  gInterpreter->AddIncludePath(TString::Format("%s/include", pwd));
  gSystem->AddLinkedLibs(TString::Format("-L%s/lib -lStarAnaConfig -Wl,-rpath,%s/lib", pwd, pwd));

  gROOT->ProcessLine(TString::Format(".L %s/common/macro/checkHistAnaFemtoKaon.C+", pwd));
  checkHistAnaFemtoKaon(rootFile, anaName, mainconfPath);
}
