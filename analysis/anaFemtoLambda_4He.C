// KF Lambda x he4. Use the matching new highpurity mainconf explicitly.
#include "StMaker/StFemtoMaker/StFemtoMaker.h"
#include "StMaker/kfparticle/FemtoLambdaKfProvider.h"
#include "FemtoLambdaRunSupport.h"

Int_t anaFemtoLambda_4He(const Char_t* inputFile, const Char_t* outputFile,
    const Char_t* jobid, Long64_t nEventsMax, const Char_t* configPath) {
  FemtoLambdaRunSupport::Run run(inputFile, outputFile, jobid, nEventsMax, configPath);
  if (!run.Prepare("he4")) return 1;
  StChain* chain = new StChain();
  StPicoDstMaker* pico = new StPicoDstMaker(StPicoDstMaker::IoRead, run.input, "picoDst");
  run.EnableBranches(pico);
  StFemtoMaker* femto = new StFemtoMaker("femto", pico, run.output);
  femto->SetMainConfigPath(run.config.Data());
  femto->SetLambdaProvider(createFemtoLambdaKfProvider());
  // StMaker construction registers both Makers; do not AddMaker a second time.
  return run.Execute(chain, pico, femto);
}
