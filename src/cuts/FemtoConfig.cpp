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

  if (values.find("speciesKeys") != values.end()) {
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

  Int_t nChannels = 1;
  if (values.find("nChannels") != values.end()) {
    nChannels = YamlParser::ToInt(values["nChannels"], nChannels);
  }
  if (nChannels < 1) nChannels = 1;

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

  if (species.empty()) SetDefaults();
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
