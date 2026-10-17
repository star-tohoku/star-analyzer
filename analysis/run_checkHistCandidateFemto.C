// run_checkHistCandidateFemto.C - Load libStarAnaConfig (YamlParser), compile
// common/macro/candidate_femto/checkHistCandidateFemto.C with ACLiC and call it.
// Usage (from the project root, inside the STAR container):
//   root4star -b -q 'analysis/run_checkHistCandidateFemto.C("out.root","config/mainconf/main_....yaml","qa.pdf","8;6,7;4,5")'
// The ACLiC build goes to tmp/kf-aclic/<host>-<pid> (per process, nothing next to the sources).

void run_checkHistCandidateFemto(const Char_t* rootFile, const Char_t* mainconfPath, const Char_t* outPdf,
                                 const Char_t* centGroups)
{
  const TString pwd = gSystem->WorkingDirectory();
  if (!rootFile || !*rootFile || !mainconfPath || !*mainconfPath || !outPdf || !*outPdf || !centGroups || !*centGroups) {
    std::cerr << "ERROR: usage run_checkHistCandidateFemto(rootFile, mainconf, outPdf, centGroups)" << std::endl;
    gSystem->Exit(1);
    return;
  }
  if (gSystem->Load(pwd + "/lib/libStarAnaConfig.so") < 0) {
    std::cerr << "ERROR: failed to load libStarAnaConfig.so" << std::endl;
    gSystem->Exit(1);
    return;
  }

  gInterpreter->AddIncludePath(pwd.Data());
  gInterpreter->AddIncludePath((pwd + "/include").Data());
  gInterpreter->AddIncludePath((pwd + "/StMaker/StCandidateFemtoMaker").Data());
  gInterpreter->AddIncludePath((pwd + "/common/macro/candidate_femto").Data());
  gSystem->AddLinkedLibs(TString::Format("-L%s/lib -lStarAnaConfig -Wl,-rpath,%s/lib", pwd.Data(), pwd.Data()));

  const TString buildDir = TString::Format("%s/tmp/kf-aclic/%s-%d", pwd.Data(), gSystem->HostName(), gSystem->GetPid());
  gSystem->mkdir(buildDir, kTRUE);
  gSystem->SetBuildDir(buildDir, kTRUE);

  if (!gSystem->CompileMacro((pwd + "/common/macro/candidate_femto/checkHistCandidateFemto.C").Data(), "kf")) {
    std::cerr << "ERROR: compiling checkHistCandidateFemto.C" << std::endl;
    gSystem->Exit(1);
    return;
  }
  const TString call = TString::Format("checkHistCandidateFemto(\"%s\",\"%s\",\"%s\",\"%s\")", rootFile, mainconfPath,
                                       outPdf, centGroups);
  Int_t error = 0;
  const Long_t result = gROOT->ProcessLine(call.Data(), &error);
  // This process's ACLiC build directory is not reused: remove it, or tmp/kf-aclic/ fills up.
  gSystem->Exec(TString::Format("rm -rf %s", buildDir.Data()));
  if (error || result != 0) {
    std::cerr << "ERROR: checkHistCandidateFemto failed (status " << result << ")" << std::endl;
    gSystem->Exit(1);
  }
}
