// Loaded before every tools/treecheck macro: adds the project include paths and links the analysis
// libraries, so a check macro can use FemtoPhiTreeSchema.h and the common helpers directly.
// The repository is taken from the working directory, never from $PWD (docs/ai/AGENT_RULES.md,
// "Batch safety").
void prelude() {
  TString repo = gSystem->WorkingDirectory();
  gSystem->AddIncludePath(Form("-I%s/include", repo.Data()));
  gSystem->AddIncludePath(Form("-I%s/StMaker/common", repo.Data()));
  gSystem->Load(Form("%s/lib/libStarAnaConfig.so", repo.Data()));
  gSystem->Load(Form("%s/lib/libStRefMultCorr.so", repo.Data()));
  gSystem->Load(Form("%s/lib/libStCommon.so", repo.Data()));
  gSystem->AddLinkedLibs(Form("-L%s/lib -lStarAnaConfig -lStRefMultCorr -lStCommon -Wl,-rpath,%s/lib",
                              repo.Data(), repo.Data()));
}
