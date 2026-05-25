#ifndef CENTRALITY_CUT_CONFIG_H
#define CENTRALITY_CUT_CONFIG_H

#include "Rtypes.h"
#include <vector>
#include <string>

class CentralityCutConfig {
public:
  static CentralityCutConfig& GetInstance();
  Bool_t LoadFromFile(const Char_t* filename);

  void SetDefaults();

  Bool_t enabled;
  std::string mode;  // "refmult" or "fxtmult"
  Bool_t rejectBadRun;
  Bool_t rejectPileup;
  Bool_t useWeight;
  Bool_t requireValidCentrality;
  Bool_t verbose;
  Bool_t fillCentralityQA;
  std::vector<Int_t> acceptedCentBins;

  Bool_t IsCentBinAccepted(Int_t cent9) const;

private:
  CentralityCutConfig();
  ~CentralityCutConfig();
  CentralityCutConfig(const CentralityCutConfig&);
  CentralityCutConfig& operator=(const CentralityCutConfig&);

  Bool_t ParseYamlFile(const Char_t* filename);
  static void ParseAcceptedBins(const std::string& value, std::vector<Int_t>& out);
};

#endif
