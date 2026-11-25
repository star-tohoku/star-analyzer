// Reconstruction and reduced-tree pairing use the same compiled femto Maker.
#include "StMaker/StFemtoMaker/StFemtoMaker.h"
#include "StMaker/kfparticle/FemtoLambdaKfProvider.h"
#include "FemtoLambdaRunSupport.h"

namespace FemtoLambda4HeTree {
bool Mode(const TString& config, bool& replay) {
  replay = false;
  std::map<std::string, std::string> main;
  if (!YamlParser::ParseFile(config.Data(), main)) return false;
  if (main.find("candidateStore") != main.end()) {
    std::cerr << "ERROR: reducedTree and candidateStore cannot be combined" << std::endl;
    return false;
  }
  const std::map<std::string, std::string>::const_iterator feature = main.find("reducedTree");
  if (feature == main.end() || feature->second.empty()) {
    std::cerr << "ERROR: this entry requires a nonempty mainconf.reducedTree reference" << std::endl;
    return false;
  }
  const std::string path(config.Data());
  const size_t root = path.rfind("/config/");
  if (root == std::string::npos) {
    std::cerr << "ERROR: mainconf must reside below config/" << std::endl;
    return false;
  }
  const std::string source = feature->second[0] == '/' ? feature->second :
      path.substr(0, root) + "/config/" + feature->second;
  std::map<std::string, std::string> values;
  if (!YamlParser::ParseFile(source.c_str(), values)) return false;
  const std::map<std::string, std::string>::const_iterator mode = values.find("mode");
  if (mode == values.end() || (mode->second != "write" && mode->second != "replay")) {
    std::cerr << "ERROR: reducedTree.mode must explicitly be write or replay" << std::endl;
    return false;
  }
  replay = mode->second == "replay";
  return true;
}

bool PrepareReplay(FemtoLambdaRunSupport::Run& run) {
  if (run.input.IsNull() || run.output.IsNull() || run.config.IsNull() || run.job.IsNull()
      || run.requested == 0 || run.requested < -1) {
    std::cerr << "ERROR: explicit input/output/jobid/mainconf are required; nEvents must be positive or -1"
              << std::endl;
    return false;
  }
  std::cout << "[mainconf] (argument) " << run.config << "\n[input-tree] (argument) " << run.input
            << "\n[output] (argument) " << run.output << "\n[jobid] (argument) " << run.job
            << "\n[stage] reduced-tree replay" << std::endl;
  if (!gSystem->AccessPathName(run.output)) {
    std::cerr << "ERROR: refusing to overwrite existing output " << run.output << std::endl;
    return false;
  }
  if (!ConfigManager::GetInstance().LoadConfig(run.config)) return false;
  const FemtoConfig& cfg = ConfigManager::GetInstance().GetFemtoConfig();
  if (cfg.species.find("lambda") == cfg.species.end() || cfg.species.find("he4") == cfg.species.end()
      || cfg.species.size() < 2 || cfg.species.size() > 3
      || (cfg.species.size() == 3 && cfg.species.find("he3") == cfg.species.end())) {
    std::cerr << "ERROR: He4 tree replay requires lambda/he4 and only the optional he3 control" << std::endl;
    return false;
  }
  const TString directory = gSystem->DirName(run.output.Data());
  if (gSystem->mkdir(directory.Data(), kTRUE) != 0 && gSystem->AccessPathName(directory.Data())) {
    std::cerr << "ERROR: cannot create output directory " << directory << std::endl;
    return false;
  }
  return true;
}

Int_t Replay(FemtoLambdaRunSupport::Run& run) {
  if (!PrepareReplay(run)) return 1;
  TStopwatch timer;
  timer.Start();
  std::cout << "[reduced-tree-input] picoDstReader=false lambdaProvider=false" << std::endl;
  StChain* chain = new StChain();
  StFemtoMaker* femto = new StFemtoMaker("femto", 0, run.output);
  femto->SetMainConfigPath(run.config.Data());
  femto->SetReducedTreeProducerId(run.job.Data());
  femto->SetReducedTreeSourceInput(run.input.Data());
  femto->SetReducedTreeInput(run.input.Data());
  if (chain->Init() != kStOK) {
    std::cerr << "ERROR: reduced-tree replay StChain Init failed" << std::endl;
    delete chain;
    return 2;
  }
  const Long64_t available = femto->ReducedTreeInputEntries();
  Int_t status = 0;
  if (available <= 0 || run.requested > available) {
    std::cerr << "ERROR: requested=" << run.requested << " available-tree-events=" << available
              << "; refusing empty or silently truncated replay" << std::endl;
    status = 3;
  }
  if (run.requested == -1) run.requested = available;
  Long64_t completed = 0;
  for (Long64_t i = 0; !status && i < run.requested; ++i) {
    if (i % 1000 == 0) std::cout << "Replaying stored event " << i << std::endl;
    chain->Clear();
    const Int_t result = chain->Make(i);
    if (result != kStOK || femto->ReducedTreeInputFinished()) {
      std::cerr << "ERROR: reduced-tree replay Make=" << result << " entry=" << i
                << " reachedEOF=" << femto->ReducedTreeInputFinished() << std::endl;
      status = 3;
      break;
    }
    ++completed;
  }
  if (completed != run.requested) status = status ? status : 3;
  femto->SetProcessingSucceeded(status == 0);
  if (chain->Finish() != kStOK) status = 4;
  timer.Stop();
  std::cout << "FemtoLambda tree replay job=" << run.job << " requested=" << run.requested
            << " completed=" << completed << " available=" << available << " status=" << status
            << " RealTime=" << timer.RealTime() << " CpuTime=" << timer.CpuTime() << std::endl;
  delete chain;
  return status;
}
}  // namespace FemtoLambda4HeTree

Int_t anaFemtoLambda_4He_tree(const Char_t* inputFile, const Char_t* outputFile,
    const Char_t* jobid, Long64_t nEventsMax, const Char_t* configPath) {
  FemtoLambdaRunSupport::Run run(inputFile, outputFile, jobid, nEventsMax, configPath);
  bool replay = false;
  if (!FemtoLambda4HeTree::Mode(run.config, replay)) return 1;
  if (replay) return FemtoLambda4HeTree::Replay(run);
  if (!run.Prepare("he4")) return 1;
  std::cout << "[stage] reduced-tree production; pairing=false" << std::endl;
  StChain* chain = new StChain();
  StPicoDstMaker* pico = new StPicoDstMaker(StPicoDstMaker::IoRead, run.input, "picoDst");
  run.EnableBranches(pico);
  StFemtoMaker* femto = new StFemtoMaker("femto", pico, run.output);
  femto->SetMainConfigPath(run.config.Data());
  femto->SetReducedTreeProducerId(run.job.Data());
  femto->SetReducedTreeSourceInput(run.input.Data());
  femto->SetLambdaProvider(createFemtoLambdaKfProvider());
  return run.Execute(chain, pico, femto);
}
