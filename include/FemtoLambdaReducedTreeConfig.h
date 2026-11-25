#ifndef FEMTO_LAMBDA_REDUCED_TREE_CONFIG_H
#define FEMTO_LAMBDA_REDUCED_TREE_CONFIG_H

#include "FemtoCandidateEventStoreConfig.h"
#include "Rtypes.h"

namespace femto_lambda_tree {

static const char* const kPrePidControlOmittedMarker =
    "prePID control bank omitted; nonselected prePID pair QA is outside the stored envelope\n";

// Reconstruction provenance deliberately excludes histogram axes, normalization
// and mixing classes: these are choices made when the saved events are paired.
struct Settings {
  bool enabled,requireOriginalMixing,resetPoolAtSegmentBoundary,storePrePidControl;
  Int_t compressionLevel;
  Long64_t autoFlush;
  std::string mode,path,selectionProvenance,mixingProvenance;
  std::vector<std::string> species;
  Settings() : enabled(false),requireOriginalMixing(true),resetPoolAtSegmentBoundary(true),
      storePrePidControl(true),compressionLevel(0),autoFlush(0) {}

  void Load(const std::string& mainPath) {
    const YAML::Node main=YAML::LoadFile(mainPath);
    femto_candidate_store::UniqueKeys(main,"reduced-tree mainconf");
    if(!main["reducedTree"].IsDefined())return;
    if(main["candidateStore"].IsDefined())throw std::runtime_error("reducedTree and candidateStore are mutually exclusive");
    if(!main["reducedTree"].IsScalar())throw std::runtime_error("mainconf.reducedTree must be a YAML path");
    path=femto_candidate_store::ConfigReference(mainPath,main["reducedTree"].as<std::string>());
    const YAML::Node c=YAML::LoadFile(path);femto_candidate_store::UniqueKeys(c,"reduced-tree feature");
    const char* names[]={"mode","speciesKeys","requireOriginalMixing","resetPoolAtSegmentBoundary","compressionLevel","autoFlush","storePrePidControl"};
    std::set<std::string> allowed;
    for(size_t i=0;i<sizeof(names)/sizeof(names[0]);++i) {
      allowed.insert(names[i]);
      if(!c[names[i]].IsDefined())throw std::runtime_error(std::string("reducedTree missing ")+names[i]);
    }
    for(YAML::const_iterator i=c.begin();i!=c.end();++i)
      if(!allowed.count(i->first.as<std::string>()))throw std::runtime_error("unknown reducedTree setting "+i->first.as<std::string>());
    mode=c["mode"].as<std::string>();
    if(mode!="write"&&mode!="replay")throw std::runtime_error("reducedTree.mode must explicitly be write or replay");
    requireOriginalMixing=c["requireOriginalMixing"].as<bool>();
    resetPoolAtSegmentBoundary=c["resetPoolAtSegmentBoundary"].as<bool>();
    storePrePidControl=c["storePrePidControl"].as<bool>();
    compressionLevel=c["compressionLevel"].as<Int_t>();autoFlush=c["autoFlush"].as<Long64_t>();
    if(compressionLevel<0||compressionLevel>9||autoFlush==0)
      throw std::runtime_error("reducedTree requires compressionLevel 0..9 and nonzero autoFlush");
    if(requireOriginalMixing&&!resetPoolAtSegmentBoundary)
      throw std::runtime_error("original mixing requires producer-segment pool resets");
    if(!c["speciesKeys"].IsSequence())throw std::runtime_error("reducedTree.speciesKeys must be a sequence");
    std::set<std::string> seen;
    for(size_t i=0;i<c["speciesKeys"].size();++i) {
      const std::string key=c["speciesKeys"][i].as<std::string>();
      if((key!="lambda"&&key!="he3"&&key!="he4")||!seen.insert(key).second)
        throw std::runtime_error("reducedTree species must be distinct lambda/he3/he4");
      species.push_back(key);
    }
    if(seen.size()!=3)throw std::runtime_error("reducedTree requires lambda, he3 and he4");

    std::ostringstream selection;selection<<"Lambda-helium reconstruction envelope 1\n";
    const char* reco[]={"event","centrality","kf","nuclearid"};
    for(size_t i=0;i<sizeof(reco)/sizeof(reco[0]);++i) {
      if(!main[reco[i]].IsScalar())throw std::runtime_error(std::string("reducedTree requires mainconf.")+reco[i]);
      selection<<reco[i]<<'\n'<<YAML::Dump(YAML::LoadFile(femto_candidate_store::ConfigReference(mainPath,main[reco[i]].as<std::string>())))<<'\n';
    }
    const YAML::Node nuclear=YAML::LoadFile(femto_candidate_store::ConfigReference(mainPath,main["nuclearid"].as<std::string>()));
    if(nuclear["dedxCalibration"].IsDefined()) {
      const YAML::Node calibration=nuclear["dedxCalibration"];
      if(!calibration["enabled"].IsScalar())throw std::runtime_error("dedxCalibration.enabled is required for reducedTree");
      if(calibration["enabled"].as<bool>()) {
        if(!calibration["calibrationFile"].IsScalar())throw std::runtime_error("dedxCalibration.calibrationFile is required");
        selection<<"effective dEdx calibration\n"<<YAML::Dump(YAML::LoadFile(femto_candidate_store::ConfigReference(mainPath,calibration["calibrationFile"].as<std::string>())))<<'\n';
      }
    }
    if(!main["maker"].IsScalar())throw std::runtime_error("reducedTree requires mainconf.maker");
    const YAML::Node maker=YAML::LoadFile(femto_candidate_store::ConfigReference(mainPath,main["maker"].as<std::string>()));
    YAML::Node canonical;
    const char* fixed[]={"lambdaEventPolicy","lambdaHe3ControlEnabled","nuclearSelectionProfile","lambdaPairMassMode","lambdaPairMass"};
    for(size_t i=0;i<sizeof(fixed)/sizeof(fixed[0]);++i) {
      if(!maker[fixed[i]].IsDefined())throw std::runtime_error(std::string("reducedTree canonical setting missing ")+fixed[i]);
      canonical[fixed[i]]=maker[fixed[i]];
    }
    selection<<"canonical candidate construction\n"<<YAML::Dump(canonical)<<'\n';
    if(main["eventPlane"].IsDefined()) {
      if(!main["eventPlane"].IsScalar())throw std::runtime_error("mainconf.eventPlane must be a path");
      const YAML::Node plane=YAML::LoadFile(femto_candidate_store::ConfigReference(mainPath,main["eventPlane"].as<std::string>()));
      YAML::Node definition;
      for(YAML::const_iterator i=plane.begin();i!=plane.end();++i) {
        const std::string key=i->first.as<std::string>();
        if(key!="bins"&&key!="writeDetailedLedger"&&key!="invalidPolicy")definition[key]=i->second;
      }
      selection<<"stored EPD definition/calibration\n"<<YAML::Dump(definition)<<'\n';
    } else selection<<"stored EPD definition/calibration\nunrecorded\n";
    if(!main["mixing"].IsScalar())throw std::runtime_error("reducedTree requires mainconf.mixing");
    mixingProvenance=YAML::Dump(YAML::LoadFile(femto_candidate_store::ConfigReference(mainPath,main["mixing"].as<std::string>())));
    // A full-control profile keeps the original fingerprint, so completed pilot
    // files remain readable. Compact stores have a distinct selection envelope.
    if(!storePrePidControl)selection<<kPrePidControlOmittedMarker;
    selectionProvenance=selection.str();enabled=true;
  }
};
} // namespace femto_lambda_tree
#endif
