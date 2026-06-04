#ifndef FEMTO_CONFIG_H
#define FEMTO_CONFIG_H

// Species/channel key naming rules: StMaker/StFemtoMaker/README.md

#include "Rtypes.h"
#include <map>
#include <string>
#include <vector>

class FemtoConfig {
public:
  static FemtoConfig& GetInstance();
  Bool_t LoadFromFile(const Char_t* filename);

  struct SpeciesDef {
    std::string key;
    std::string builderType;  // track | resonance
    std::string particleKey;  // proton | phi | ...
    std::string cutsRef;
  };

  struct ChannelDef {
    std::string name;
    std::string partA;
    std::string partB;
    Bool_t enabled;
    Bool_t doMixing;
    Double_t signalMin;
    Double_t signalMax;
    Double_t normQMin;
    Double_t normQMax;
  };

  std::map<std::string, SpeciesDef> species;
  std::vector<ChannelDef> channels;

  Bool_t Validate() const;
  const SpeciesDef* FindSpecies(const std::string& key) const;
  const ChannelDef* FindChannel(const std::string& name) const;

  void SetDefaults();

private:
  FemtoConfig();
  ~FemtoConfig();
  FemtoConfig(const FemtoConfig&);
  FemtoConfig& operator=(const FemtoConfig&);

  Bool_t ParseYamlFile(const Char_t* filename);
};

#endif
