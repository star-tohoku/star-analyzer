#ifndef FEMTO_CANDIDATE_EVENT_STORE_CONFIG_H
#define FEMTO_CANDIDATE_EVENT_STORE_CONFIG_H
#include "yaml-cpp/yaml.h"
#include <cmath>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
namespace femto_candidate_store {
inline std::string ConfigReference(const std::string& mainPath,const std::string& value) {
 if(value.empty())throw std::runtime_error("empty candidate-store config reference");
 if(value[0]=='/')return value;
 const size_t p=mainPath.find("/config/");
 return (p==std::string::npos?"":mainPath.substr(0,p+1))+"config/"+value;
}
inline void UniqueKeys(const YAML::Node& node,const std::string& where) {
 if(!node.IsMap())throw std::runtime_error(where+" must be a mapping");
 std::set<std::string> keys;
 for(YAML::const_iterator it=node.begin();it!=node.end();++it)
  if(!it->first.IsScalar()||!keys.insert(it->first.as<std::string>()).second)throw std::runtime_error(where+" has duplicate/invalid keys");
}
struct Settings {
 Settings():enabled(false),requireOriginalMixing(true),resetPoolAtFileBoundary(true),calculateEventPlane(false),excludeCandidateTracks(false),harmonic(0),minTracks(0),minNHitsFit(0),minFitRatio(0),minPt(0),maxPt(0),minEta(0),maxEta(0),maxDca(0){}
 bool enabled,requireOriginalMixing,resetPoolAtFileBoundary,calculateEventPlane,excludeCandidateTracks;
 std::string mode,eventPlaneSource,path,eventPlaneDefinition,selectionProvenance,mixingProvenance;
 std::vector<std::string> species;
 int harmonic,minTracks,minNHitsFit;double minFitRatio,minPt,maxPt,minEta,maxEta,maxDca;
 void Load(const std::string& mainPath) {
  const YAML::Node main=YAML::LoadFile(mainPath);UniqueKeys(main,"candidate-store mainconf");
  if(!main["candidateStore"].IsDefined())return;
  if(!main["candidateStore"].IsScalar())throw std::runtime_error("mainconf.candidateStore must be a YAML path");
  path=ConfigReference(mainPath,main["candidateStore"].as<std::string>());
  const YAML::Node c=YAML::LoadFile(path);UniqueKeys(c,"candidate-store config");
  const char* names[]={"mode","requireOriginalMixing","resetPoolAtFileBoundary","eventPlaneSource","eventPlaneCalculation","speciesKeys"};std::set<std::string> allowed;
  for(unsigned i=0;i<sizeof(names)/sizeof(names[0]);++i){allowed.insert(names[i]);if(!c[names[i]].IsDefined())throw std::runtime_error(std::string("candidate-store missing ")+names[i]);}
  for(YAML::const_iterator it=c.begin();it!=c.end();++it)if(!allowed.count(it->first.as<std::string>()))throw std::runtime_error("unknown candidate-store setting");
  mode=c["mode"].as<std::string>();eventPlaneSource=c["eventPlaneSource"].as<std::string>();
  requireOriginalMixing=c["requireOriginalMixing"].as<bool>();resetPoolAtFileBoundary=c["resetPoolAtFileBoundary"].as<bool>();
  if(!c["speciesKeys"].IsSequence() || c["speciesKeys"].size()==0)throw std::runtime_error("candidate-store speciesKeys must be a nonempty sequence");
  std::set<std::string> speciesSeen;
  for(size_t i=0;i<c["speciesKeys"].size();++i){const std::string key=c["speciesKeys"][i].as<std::string>();if(key.empty()||!speciesSeen.insert(key).second)throw std::runtime_error("invalid/duplicate candidate-store species key");species.push_back(key);}

  if((mode!="write"&&mode!="replay")||(eventPlaneSource!="original"&&eventPlaneSource!="stored"))throw std::runtime_error("invalid candidate-store mode/eventPlaneSource");
  if(mode=="write"&&eventPlaneSource!="original")throw std::runtime_error("candidate writing must preserve original pairing event plane");
  if(requireOriginalMixing && (!resetPoolAtFileBoundary || eventPlaneSource!="original"))throw std::runtime_error("original mixing requires original event plane and per-file pool reset");
  const YAML::Node ep=c["eventPlaneCalculation"];UniqueKeys(ep,"candidate-store event plane");
  const char* epNames[]={"enabled","harmonic","minTracks","minNHitsFit","minFitRatio","minPt","maxPt","minEta","maxEta","maxDca","weight","excludeCandidateTracks"};allowed.clear();
  for(unsigned i=0;i<sizeof(epNames)/sizeof(epNames[0]);++i){allowed.insert(epNames[i]);if(!ep[epNames[i]].IsScalar())throw std::runtime_error(std::string("candidate-store event plane missing ")+epNames[i]);}
  for(YAML::const_iterator it=ep.begin();it!=ep.end();++it)if(!allowed.count(it->first.as<std::string>()))throw std::runtime_error("unknown candidate-store event-plane setting");
  calculateEventPlane=ep["enabled"].as<bool>();harmonic=ep["harmonic"].as<int>();minTracks=ep["minTracks"].as<int>();minNHitsFit=ep["minNHitsFit"].as<int>();
  minFitRatio=ep["minFitRatio"].as<double>();minPt=ep["minPt"].as<double>();maxPt=ep["maxPt"].as<double>();minEta=ep["minEta"].as<double>();maxEta=ep["maxEta"].as<double>();maxDca=ep["maxDca"].as<double>();excludeCandidateTracks=ep["excludeCandidateTracks"].as<bool>();
  if(harmonic!=2||minTracks<1||minNHitsFit<0||!std::isfinite(minFitRatio)||minFitRatio<0||minFitRatio>1||!std::isfinite(minPt)||!std::isfinite(maxPt)||minPt<0||minPt>=maxPt||!std::isfinite(minEta)||!std::isfinite(maxEta)||minEta>=maxEta||!std::isfinite(maxDca)||maxDca<0||ep["weight"].as<std::string>()!="pt")throw std::runtime_error("invalid candidate-store TPC Q2 definition");
  eventPlaneDefinition="primary charged tracks; pT-weighted unflattened TPC Q2; "+YAML::Dump(ep);
  const char* refs[]={"analysis","event","centrality","kf","nuclearid","maker"};std::ostringstream selection;
  selection<<"femto selected-candidate provenance 1\n";
  for(unsigned i=0;i<sizeof(refs)/sizeof(refs[0]);++i){if(!main[refs[i]].IsScalar())throw std::runtime_error(std::string("candidate store requires mainconf.")+refs[i]);const std::string f=ConfigReference(mainPath,main[refs[i]].as<std::string>());selection<<refs[i]<<'\n'<<YAML::Dump(YAML::LoadFile(f))<<'\n';}
  const YAML::Node nuclear=YAML::LoadFile(ConfigReference(mainPath,main["nuclearid"].as<std::string>()));
  if(nuclear["dedxCalibration"].IsDefined()) {
    const YAML::Node cal=nuclear["dedxCalibration"];
    if(!cal["enabled"].IsScalar())throw std::runtime_error("dedxCalibration.enabled is required for provenance");
    if(cal["enabled"].as<bool>()) {
      if(!cal["calibrationFile"].IsScalar())throw std::runtime_error("dedx calibration file is required for provenance");
      selection<<"effective nuclear dEdx calibration\n"<<YAML::Dump(YAML::LoadFile(ConfigReference(mainPath,cal["calibrationFile"].as<std::string>())))<<'\n';
    }
  }
  selection<<"auxiliary event plane\n"<<eventPlaneDefinition<<'\n';selectionProvenance=selection.str();
  if(!main["mixing"].IsScalar())throw std::runtime_error("candidate store requires mainconf.mixing");
  mixingProvenance=YAML::Dump(YAML::LoadFile(ConfigReference(mainPath,main["mixing"].as<std::string>())));
  enabled=true;
 }
};
}
#endif
