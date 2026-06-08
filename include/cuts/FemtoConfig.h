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
    std::string particleKey;  // proton | phi | phi_rotation | ...
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

  // Zhangwei-like proton bachelor cuts for femto pairing (see 4ReadTree/analysis.cxx).
  std::string protonChargeMode;
  Double_t protonMaxDca;
  Double_t protonMinPtPre;
  Double_t protonMinPtPair;
  Double_t protonMaxPtPair;
  Double_t protonMaxAbsEta;
  Double_t protonMaxAbsNSigma;
  Short_t protonMinNHitsFit;
  Double_t protonMinNHitsRatio;
  Double_t protonTofMomentumThreshold;
  Double_t protonMinMass2;
  Double_t protonMaxMass2;
  Double_t protonMinRapidityCm;
  Double_t protonMaxRapidityCm;

  // Rotation background (phi_rot species).
  Bool_t rotationEnabled;
  std::string rotationSpeciesKey;
  std::string rotationParticleKey;
  Int_t rotationN;
  Double_t rotationMinAngle;
  Double_t rotationMaxAngle;
  Int_t rotationSeed;

  // checkHist CF: merge this many adjacent k* bins after merge (1 = no rebin).
  Int_t cfRebinFactor;

  // checkHist CF cent slice: project hKstar*VsCent over cent9 in [cfCent9Min, cfCent9Max].
  Int_t cfCent9Min;
  Int_t cfCent9Max;

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
  void ApplyYamlValues(const std::map<std::string, std::string>& values);
};

#endif
