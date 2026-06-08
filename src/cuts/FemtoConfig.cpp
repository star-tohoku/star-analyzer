#include "cuts/FemtoConfig.h"
#include "YamlParser.h"
#include <iostream>
#include <sstream>

namespace {
std::string Trim(const std::string& s) {
  size_t start = s.find_first_not_of(" \t\r\n");
  if (start == std::string::npos) return "";
  size_t end = s.find_last_not_of(" \t\r\n");
  return s.substr(start, end - start + 1);
}

std::vector<std::string> SplitComma(const std::string& s) {
  std::vector<std::string> out;
  std::stringstream ss(s);
  std::string item;
  while (std::getline(ss, item, ',')) {
    item = Trim(item);
    if (!item.empty()) out.push_back(item);
  }
  return out;
}
}  // namespace

FemtoConfig& FemtoConfig::GetInstance() {
  static FemtoConfig instance;
  return instance;
}

FemtoConfig::FemtoConfig() { SetDefaults(); }

FemtoConfig::~FemtoConfig() {}

void FemtoConfig::SetDefaults() {
  species.clear();
  channels.clear();

  protonChargeMode = "positive";
  protonMaxDca = 1.0;
  protonMinPtPre = 0.2;
  protonMinPtPair = 0.4;
  protonMaxPtPair = 2.0;
  protonMaxAbsEta = 2.0;
  protonMaxAbsNSigma = 2.0;
  protonMinNHitsFit = 15;
  protonMinNHitsRatio = 0.52;
  protonTofMomentumThreshold = 2.0;
  protonMinMass2 = 0.6;
  protonMaxMass2 = 1.2;
  protonMinRapidityCm = -0.5;
  protonMaxRapidityCm = 0.0;

  rotationEnabled = kFALSE;
  rotationSpeciesKey = "phi_rot";
  rotationParticleKey = "phi_rotation";
  rotationN = 10;
  rotationMinAngle = 0.0;
  rotationMaxAngle = 6.283185307179586;
  rotationSeed = 0;
  cfRebinFactor = 1;
  cfCent9Min = 2;
  cfCent9Max = 8;

  SpeciesDef proton;
  proton.key = "proton";
  proton.builderType = "track";
  proton.particleKey = "proton";
  species[proton.key] = proton;

  SpeciesDef phi;
  phi.key = "phi";
  phi.builderType = "resonance";
  phi.particleKey = "phi";
  species[phi.key] = phi;

  ChannelDef ch;
  ch.name = "phi_proton";
  ch.partA = "phi";
  ch.partB = "proton";
  ch.enabled = kTRUE;
  ch.doMixing = kTRUE;
  ch.signalMin = 1.01;
  ch.signalMax = 1.03;
  ch.normQMin = 0.2;
  ch.normQMax = 0.3;
  channels.push_back(ch);
}

void FemtoConfig::ApplyYamlValues(const std::map<std::string, std::string>& values) {
  if (values.find("protonChargeMode") != values.end()) protonChargeMode = values.at("protonChargeMode");
  if (values.find("protonMaxDca") != values.end()) protonMaxDca = YamlParser::ToDouble(values.at("protonMaxDca"), protonMaxDca);
  if (values.find("protonMinPtPre") != values.end()) protonMinPtPre = YamlParser::ToDouble(values.at("protonMinPtPre"), protonMinPtPre);
  if (values.find("protonMinPtPair") != values.end()) protonMinPtPair = YamlParser::ToDouble(values.at("protonMinPtPair"), protonMinPtPair);
  if (values.find("protonMaxPtPair") != values.end()) protonMaxPtPair = YamlParser::ToDouble(values.at("protonMaxPtPair"), protonMaxPtPair);
  if (values.find("protonMaxAbsEta") != values.end()) protonMaxAbsEta = YamlParser::ToDouble(values.at("protonMaxAbsEta"), protonMaxAbsEta);
  if (values.find("protonMaxAbsNSigma") != values.end()) protonMaxAbsNSigma = YamlParser::ToDouble(values.at("protonMaxAbsNSigma"), protonMaxAbsNSigma);
  if (values.find("protonMinNHitsFit") != values.end()) protonMinNHitsFit = YamlParser::ToInt(values.at("protonMinNHitsFit"), protonMinNHitsFit);
  if (values.find("protonMinNHitsRatio") != values.end()) protonMinNHitsRatio = YamlParser::ToDouble(values.at("protonMinNHitsRatio"), protonMinNHitsRatio);
  if (values.find("protonTofMomentumThreshold") != values.end()) {
    protonTofMomentumThreshold = YamlParser::ToDouble(values.at("protonTofMomentumThreshold"), protonTofMomentumThreshold);
  }
  if (values.find("protonMinMass2") != values.end()) protonMinMass2 = YamlParser::ToDouble(values.at("protonMinMass2"), protonMinMass2);
  if (values.find("protonMaxMass2") != values.end()) protonMaxMass2 = YamlParser::ToDouble(values.at("protonMaxMass2"), protonMaxMass2);
  if (values.find("protonMinRapidityCm") != values.end()) protonMinRapidityCm = YamlParser::ToDouble(values.at("protonMinRapidityCm"), protonMinRapidityCm);
  if (values.find("protonMaxRapidityCm") != values.end()) protonMaxRapidityCm = YamlParser::ToDouble(values.at("protonMaxRapidityCm"), protonMaxRapidityCm);

  if (values.find("rotationEnabled") != values.end()) rotationEnabled = YamlParser::ToBool(values.at("rotationEnabled"), rotationEnabled);
  if (values.find("rotationSpeciesKey") != values.end()) rotationSpeciesKey = values.at("rotationSpeciesKey");
  if (values.find("rotationParticleKey") != values.end()) rotationParticleKey = values.at("rotationParticleKey");
  if (values.find("rotationN") != values.end()) rotationN = YamlParser::ToInt(values.at("rotationN"), rotationN);
  if (values.find("rotationMinAngle") != values.end()) rotationMinAngle = YamlParser::ToDouble(values.at("rotationMinAngle"), rotationMinAngle);
  if (values.find("rotationMaxAngle") != values.end()) rotationMaxAngle = YamlParser::ToDouble(values.at("rotationMaxAngle"), rotationMaxAngle);
  if (values.find("rotationSeed") != values.end()) rotationSeed = YamlParser::ToInt(values.at("rotationSeed"), rotationSeed);
  if (values.find("cfRebinFactor") != values.end()) cfRebinFactor = YamlParser::ToInt(values.at("cfRebinFactor"), cfRebinFactor);
  if (values.find("cfCent9Min") != values.end()) cfCent9Min = YamlParser::ToInt(values.at("cfCent9Min"), cfCent9Min);
  if (values.find("cfCent9Max") != values.end()) cfCent9Max = YamlParser::ToInt(values.at("cfCent9Max"), cfCent9Max);
}

Bool_t FemtoConfig::LoadFromFile(const Char_t* filename) {
  return ParseYamlFile(filename);
}

Bool_t FemtoConfig::ParseYamlFile(const Char_t* filename) {
  std::map<std::string, std::string> values;
  if (!YamlParser::ParseFile(filename, values)) {
    std::cerr << "WARNING: Failed to parse " << filename << ", using femto defaults" << std::endl;
    return kFALSE;
  }

  species.clear();
  channels.clear();
  SetDefaults();

  if (values.find("speciesKeys") != values.end()) {
    species.clear();
    std::vector<std::string> keys = SplitComma(values["speciesKeys"]);
    for (size_t i = 0; i < keys.size(); i++) {
      const std::string& key = keys[i];
      SpeciesDef def;
      def.key = key;
      std::string btKey = "species_" + key + "_builderType";
      std::string pkKey = "species_" + key + "_particleKey";
      std::string crKey = "species_" + key + "_cutsRef";
      if (values.find(btKey) != values.end()) def.builderType = values[btKey];
      if (values.find(pkKey) != values.end()) def.particleKey = values[pkKey];
      if (values.find(crKey) != values.end()) def.cutsRef = values[crKey];
      species[key] = def;
    }
  }

  ApplyYamlValues(values);

  Int_t nChannels = 1;
  if (values.find("nChannels") != values.end()) {
    nChannels = YamlParser::ToInt(values["nChannels"], nChannels);
  }
  if (nChannels < 1) nChannels = 1;

  channels.clear();
  for (Int_t ic = 0; ic < nChannels; ic++) {
    std::string prefix = (nChannels == 1) ? "" : ("channel_" + std::to_string(ic) + "_");
    auto getVal = [&](const std::string& k) -> std::string {
      std::string flat = prefix + k;
      if (values.find(flat) != values.end()) return values[flat];
      if (nChannels == 1 && values.find(k) != values.end()) return values[k];
      return "";
    };

    ChannelDef ch;
    ch.name = getVal("channelName");
    if (ch.name.empty()) ch.name = getVal("name");
    ch.partA = getVal("partA");
    ch.partB = getVal("partB");
    ch.enabled = YamlParser::ToBool(getVal("enabled"), kTRUE);
    ch.doMixing = YamlParser::ToBool(getVal("doMixing"), kTRUE);
    ch.signalMin = YamlParser::ToDouble(getVal("signalMin"), 1.01);
    ch.signalMax = YamlParser::ToDouble(getVal("signalMax"), 1.03);
    ch.normQMin = YamlParser::ToDouble(getVal("normQMin"), 0.2);
    ch.normQMax = YamlParser::ToDouble(getVal("normQMax"), 0.3);

    if (ch.name.empty() && !ch.partA.empty() && !ch.partB.empty()) {
      ch.name = ch.partA + "_" + ch.partB;
    }
    if (!ch.partA.empty() && !ch.partB.empty()) {
      channels.push_back(ch);
    }
  }

  if (species.empty()) {
    SetDefaults();
    species.clear();
    channels.clear();
    SpeciesDef proton;
    proton.key = "proton";
    proton.builderType = "track";
    proton.particleKey = "proton";
    species[proton.key] = proton;
    SpeciesDef phi;
    phi.key = "phi";
    phi.builderType = "resonance";
    phi.particleKey = "phi";
    species[phi.key] = phi;
  }
  if (channels.empty()) {
    ChannelDef ch;
    ch.name = "phi_proton";
    ch.partA = "phi";
    ch.partB = "proton";
    ch.enabled = kTRUE;
    ch.doMixing = kTRUE;
    ch.signalMin = 1.01;
    ch.signalMax = 1.03;
    ch.normQMin = 0.2;
    ch.normQMax = 0.3;
    channels.push_back(ch);
  }

  return Validate();
}

Bool_t FemtoConfig::Validate() const {
  Bool_t ok = kTRUE;
  for (size_t i = 0; i < channels.size(); i++) {
    const ChannelDef& ch = channels[i];
    if (species.find(ch.partA) == species.end()) {
      std::cerr << "ERROR: FemtoConfig channel '" << ch.name << "' partA '" << ch.partA
                << "' not found in species map" << std::endl;
      ok = kFALSE;
    }
    if (species.find(ch.partB) == species.end()) {
      std::cerr << "ERROR: FemtoConfig channel '" << ch.name << "' partB '" << ch.partB
                << "' not found in species map" << std::endl;
      ok = kFALSE;
    }
  }
  if (rotationEnabled && species.find(rotationSpeciesKey) == species.end()) {
    std::cerr << "ERROR: rotationEnabled but species '" << rotationSpeciesKey << "' not defined" << std::endl;
    ok = kFALSE;
  }
  if (cfRebinFactor < 1) {
    std::cerr << "ERROR: FemtoConfig cfRebinFactor must be >= 1 (got " << cfRebinFactor << ")" << std::endl;
    ok = kFALSE;
  }
  if (cfCent9Min < 0 || cfCent9Max > 8 || cfCent9Min > cfCent9Max) {
    std::cerr << "ERROR: FemtoConfig cfCent9Min/Max out of range (got [" << cfCent9Min << ", " << cfCent9Max
              << "], valid cent9 is 0-8)" << std::endl;
    ok = kFALSE;
  }
  return ok;
}

const FemtoConfig::SpeciesDef* FemtoConfig::FindSpecies(const std::string& key) const {
  std::map<std::string, SpeciesDef>::const_iterator it = species.find(key);
  if (it == species.end()) return 0;
  return &it->second;
}

const FemtoConfig::ChannelDef* FemtoConfig::FindChannel(const std::string& name) const {
  for (size_t i = 0; i < channels.size(); i++) {
    if (channels[i].name == name) return &channels[i];
  }
  return 0;
}
