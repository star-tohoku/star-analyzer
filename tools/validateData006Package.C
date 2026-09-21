// Machine validation for a finalized DATA-006 package (ROOT 5 / C++98 compatible).
#include "TFile.h"
#include "TH1.h"
#include "TMath.h"

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace data006_validation {

std::vector<std::string> Split(const std::string& line) {
  std::vector<std::string> out;
  std::stringstream ss(line);
  std::string item;
  while (std::getline(ss, item, ',')) out.push_back(item);
  return out;
}

Bool_t Same(Double_t a, Double_t b) {
  if (TMath::IsNaN(a) && TMath::IsNaN(b)) return kTRUE;
  if (!TMath::Finite(a) || !TMath::Finite(b)) return a == b;
  return std::fabs(a - b) <= 1e-12 + 2e-12 * std::fabs(b);
}

}  // namespace data006_validation

Bool_t validateData006Package(const Char_t* packageDirectory) {
  using namespace data006_validation;
  const std::string dir(packageDirectory ? packageDirectory : "");
  if (dir.empty()) {
    std::cerr << "ERROR: empty DATA-006 package directory" << std::endl;
    return kFALSE;
  }

  std::ifstream status((dir + "/CLOSURE_STATUS.txt").c_str());
  std::string line;
  if (!status || !std::getline(status, line) || line != "FIT_READY") {
    std::cerr << "ERROR: DATA-006 closure status is not FIT_READY" << std::endl;
    return kFALSE;
  }

  TFile f((dir + "/data006_correlations.root").c_str(), "READ");
  if (f.IsZombie()) {
    std::cerr << "ERROR: cannot open finalized DATA-006 ROOT file" << std::endl;
    return kFALSE;
  }

  Long64_t bad = 0, manifestRows = 0, correlationRows = 0;
  std::ifstream manifest((dir + "/object_manifest.csv").c_str());
  if (!manifest || !std::getline(manifest, line)) {
    std::cerr << "ERROR: cannot read DATA-006 object manifest" << std::endl;
    return kFALSE;
  }
  while (std::getline(manifest, line)) {
    const std::vector<std::string> v = Split(line);
    ++manifestRows;
    if (v.size() != 9 || !f.Get(v[4].c_str())) {
      std::cerr << "ERROR: invalid/missing manifest object at row " << manifestRows + 1
                << std::endl;
      ++bad;
    }
  }

  std::ifstream corr((dir + "/correlations.csv").c_str());
  if (!corr || !std::getline(corr, line)) {
    std::cerr << "ERROR: cannot read DATA-006 correlations CSV" << std::endl;
    return kFALSE;
  }
  while (std::getline(corr, line)) {
    const std::vector<std::string> v = Split(line);
    ++correlationRows;
    if (v.size() != 16) {
      std::cerr << "ERROR: correlations CSV has wrong column count at row "
                << correlationRows + 1 << std::endl;
      ++bad;
      continue;
    }
    const std::string base = v[0] + "/" + v[1] + "/";
    TH1* se = (TH1*)f.Get((base + "hSE").c_str());
    TH1* me = (TH1*)f.Get((base + "hME").c_str());
    TH1* cf = (TH1*)f.Get((base + "hCF").c_str());
    if (!se || !me || !cf) {
      std::cerr << "ERROR: missing matched SE/ME/CF for " << base << std::endl;
      ++bad;
      continue;
    }
    Int_t bin = se->GetXaxis()->FindFixBin(std::strtod(v[3].c_str(), 0) / 1000.0 + 1e-12);
    if (v[2] == "underflow") bin = 0;
    if (v[2] == "overflow") bin = se->GetNbinsX() + 1;
    const Double_t csvSE = std::strtod(v[5].c_str(), 0);
    const Double_t csvME = std::strtod(v[7].c_str(), 0);
    const Double_t csvCF = std::strtod(v[13].c_str(), 0);
    if (!Same(csvSE, se->GetBinContent(bin)) ||
        !Same(csvME, me->GetBinContent(bin)) ||
        !Same(csvCF, cf->GetBinContent(bin))) {
      std::cerr << "ERROR: ROOT/CSV mismatch at correlations row " << correlationRows + 1
                << std::endl;
      ++bad;
    }
  }
  f.Close();

  std::cout << "[validateData006Package] manifestObjects=" << manifestRows
            << " correlationRows=" << correlationRows << " mismatches=" << bad << std::endl;
  return bad == 0 && manifestRows > 0 && correlationRows > 0;
}
