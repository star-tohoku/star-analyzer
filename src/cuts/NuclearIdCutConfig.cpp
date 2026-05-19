#include "cuts/NuclearIdCutConfig.h"
#include "YamlParser.h"
#include <map>
#include <string>
#include <iostream>

NuclearIdCutConfig& NuclearIdCutConfig::GetInstance() {
  static NuclearIdCutConfig instance;
  return instance;
}

NuclearIdCutConfig::NuclearIdCutConfig() {
  SetDefaults();
}

NuclearIdCutConfig::~NuclearIdCutConfig() {
}

void NuclearIdCutConfig::SetDefaults() {
  nSigmaFill1D = 1.5;
  nSigmaExclude = 2.0;
  maxNSigmaNuclear = 2.0;
  m2SigmaCut = 1.5;
  minP_M2cut = 0.2;
  minM2_M2cut = 0.01;
  maxM2_M2cut = 10.0;
}

Bool_t NuclearIdCutConfig::LoadFromFile(const Char_t* filename) {
  return ParseYamlFile(filename);
}

Bool_t NuclearIdCutConfig::ParseYamlFile(const Char_t* filename) {
  std::map<std::string, std::string> values;
  if (!YamlParser::ParseFile(filename, values)) {
    std::cerr << "WARNING: Failed to parse " << filename
              << ", using default values" << std::endl;
    return kFALSE;
  }

  if (values.find("nSigmaFill1D") != values.end()) {
    nSigmaFill1D = YamlParser::ToDouble(values["nSigmaFill1D"], nSigmaFill1D);
  }
  if (values.find("nSigmaExclude") != values.end()) {
    nSigmaExclude = YamlParser::ToDouble(values["nSigmaExclude"], nSigmaExclude);
  }
  if (values.find("maxNSigmaNuclear") != values.end()) {
    maxNSigmaNuclear = YamlParser::ToDouble(values["maxNSigmaNuclear"], maxNSigmaNuclear);
  }
  if (values.find("m2SigmaCut") != values.end()) {
    m2SigmaCut = YamlParser::ToDouble(values["m2SigmaCut"], m2SigmaCut);
  }
  if (values.find("minP_M2cut") != values.end()) {
    minP_M2cut = YamlParser::ToDouble(values["minP_M2cut"], minP_M2cut);
  }
  if (values.find("minM2_M2cut") != values.end()) {
    minM2_M2cut = YamlParser::ToDouble(values["minM2_M2cut"], minM2_M2cut);
  }
  if (values.find("maxM2_M2cut") != values.end()) {
    maxM2_M2cut = YamlParser::ToDouble(values["maxM2_M2cut"], maxM2_M2cut);
  }

  return kTRUE;
}
