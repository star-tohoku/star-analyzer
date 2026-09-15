// run_checkHistAnaPhi.C - Load libStarAnaConfig, compile checkHistAnaPhi.C+, and call it.
// ROOT's CINT does not resolve symbols from dynamically loaded libs; we must compile with ACLiC and link.
// Usage: root4star -b -q 'run_checkHistAnaPhi.C("rootfile/...","anaName","config/mainconf/main_anaName.yaml")'

void run_checkHistAnaPhi(const Char_t* rootFile,
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

  gROOT->ProcessLine(TString::Format(".L %s/common/macro/checkHistAnaPhi.C+", pwd));
  checkHistAnaPhi(rootFile, anaName, mainconfPath);
}
