#include "cuts/PhiCutConfig.h"
#include "YamlParser.h"
#include <map>
#include <string>
#include <iostream>

PhiCutConfig& PhiCutConfig::GetInstance() {
  static PhiCutConfig instance;
  return instance;
}

PhiCutConfig::PhiCutConfig() {
  SetDefaults();
}

PhiCutConfig::~PhiCutConfig() {
}

void PhiCutConfig::SetDefaults() {
  nSigmaKaon = 2.0;
  minMass2Kaon = 0.15;
  maxMass2Kaon = 0.35;
  maxDCAKaon = 2.0;
  maxDCAKK = 999.0;
  useInvMassCut = kFALSE;
  minInvMass = 0.99;
  maxInvMass = 1.05;
  minOpeningAngle = 0.0;
  maxOpeningAngle = 0.5;
  minPairRapidity = -0.8;
  maxPairRapidity = 0.8;
  minPtEp = 0.15;
  maxPtEp = 2.0;
  maxEtaEp = 1.0;
  maxNTr = 0;  // no limit
}

Bool_t PhiCutConfig::LoadFromFile(const Char_t* filename) {
  return ParseYamlFile(filename);
}

Bool_t PhiCutConfig::ParseYamlFile(const Char_t* filename) {
  std::map<std::string, std::string> values;
  if (!YamlParser::ParseFile(filename, values)) {
    std::cerr << "WARNING: Failed to parse " << filename 
              << ", using default values" << std::endl;
    return kFALSE;
  }
  
  if (values.find("nSigmaKaon") != values.end()) {
    nSigmaKaon = YamlParser::ToDouble(values["nSigmaKaon"], nSigmaKaon);
  }
  if (values.find("minMass2Kaon") != values.end()) {
    minMass2Kaon = YamlParser::ToDouble(values["minMass2Kaon"], minMass2Kaon);
  }
  if (values.find("maxMass2Kaon") != values.end()) {
    maxMass2Kaon = YamlParser::ToDouble(values["maxMass2Kaon"], maxMass2Kaon);
  }
  if (values.find("maxDCAKaon") != values.end()) {
    maxDCAKaon = YamlParser::ToDouble(values["maxDCAKaon"], maxDCAKaon);
  }
  if (values.find("maxDCAKK") != values.end()) {
    maxDCAKK = YamlParser::ToDouble(values["maxDCAKK"], maxDCAKK);
  }
  if (values.find("useInvMassCut") != values.end()) {
    useInvMassCut = YamlParser::ToBool(values["useInvMassCut"], useInvMassCut);
  }
  if (values.find("minInvMass") != values.end()) {
    minInvMass = YamlParser::ToDouble(values["minInvMass"], minInvMass);
  }
  if (values.find("maxInvMass") != values.end()) {
    maxInvMass = YamlParser::ToDouble(values["maxInvMass"], maxInvMass);
  }
  if (values.find("minOpeningAngle") != values.end()) {
    minOpeningAngle = YamlParser::ToDouble(values["minOpeningAngle"], minOpeningAngle);
  }
  if (values.find("maxOpeningAngle") != values.end()) {
    maxOpeningAngle = YamlParser::ToDouble(values["maxOpeningAngle"], maxOpeningAngle);
  }
  if (values.find("minPairRapidity") != values.end()) {
    minPairRapidity = YamlParser::ToDouble(values["minPairRapidity"], minPairRapidity);
  }
  if (values.find("maxPairRapidity") != values.end()) {
    maxPairRapidity = YamlParser::ToDouble(values["maxPairRapidity"], maxPairRapidity);
  }
  if (values.find("minPtEp") != values.end()) {
    minPtEp = YamlParser::ToDouble(values["minPtEp"], minPtEp);
  }
  if (values.find("maxPtEp") != values.end()) {
    maxPtEp = YamlParser::ToDouble(values["maxPtEp"], maxPtEp);
  }
  if (values.find("maxEtaEp") != values.end()) {
    maxEtaEp = YamlParser::ToDouble(values["maxEtaEp"], maxEtaEp);
  }
  if (values.find("maxNTr") != values.end()) {
    maxNTr = YamlParser::ToInt(values["maxNTr"], maxNTr);
  }

  return kTRUE;
}

