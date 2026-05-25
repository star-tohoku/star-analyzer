#ifndef PHI_CUT_CONFIG_H
#define PHI_CUT_CONFIG_H

#include "Rtypes.h"

class PhiCutConfig {
public:
  static PhiCutConfig& GetInstance();
  Bool_t LoadFromFile(const Char_t* filename);
  
  // Cut values (public member variables)
  Double_t nSigmaKaon;
  Double_t minMass2Kaon;
  Double_t maxMass2Kaon;
  Double_t maxDCAKaon;
  Double_t maxDCAKK;
  Bool_t useInvMassCut;
  Double_t minInvMass;
  Double_t maxInvMass;
  Double_t minOpeningAngle;
  Double_t maxOpeningAngle;
  Double_t minPairRapidity;
  Double_t maxPairRapidity;
  // Event-plane (Q-vector) track selection
  Double_t minPtEp;
  Double_t maxPtEp;
  Double_t maxEtaEp;
  // Event-level: skip event if nTracks > maxNTr (0 = no limit)
  Int_t maxNTr;

  // Set default values
  void SetDefaults();
  
private:
  PhiCutConfig();
  ~PhiCutConfig();
  PhiCutConfig(const PhiCutConfig&);
  PhiCutConfig& operator=(const PhiCutConfig&);
  
  Bool_t ParseYamlFile(const Char_t* filename);
};

#endif

