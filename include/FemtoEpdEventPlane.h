#ifndef FEMTO_EPD_EVENT_PLANE_H
#define FEMTO_EPD_EVENT_PLANE_H

#include "Rtypes.h"
#include "TMath.h"
#include <cmath>
#include <string>
#include <vector>

class StPicoDst;
class TFile;
class TProfile2D;
class TProfile3D;

namespace femto_epd {

// Reproduce the downstream Float_t conversions, including the rounded bin
// width. The raw physical angle is the corrected first-order EPD plane.
inline Int_t XialeiBin(Double_t physicalPsi1) {
  if (!std::isfinite(physicalPsi1)) return -1;
  const Float_t psi = physicalPsi1;
  const Float_t half = psi / 2;
  const Double_t low = -TMath::Pi()/2.0, high = TMath::Pi()/2.0;
  const Float_t width = (high-low)/Float_t(8);
  if (!(half > low && half < high)) return -1;
  const Int_t bin = Int_t((half-low)/width);
  return bin >= 0 && bin < 8 ? bin : -1;
}

// Reclassify the stored physical first-order plane without changing its
// canonical eight-bin validity. Power-of-two refinements retain the same
// Float_t conversion and rounded-width boundaries as the original binning.
inline Int_t MixingBin(Double_t physicalPsi1,Int_t bins) {
  if(bins<1)return -1;
  const Int_t canonical=XialeiBin(physicalPsi1);
  if(canonical<0)return -1;
  if(bins==1)return 0;
  if(bins==8)return canonical;
  const Float_t psi=physicalPsi1;
  const Float_t half=psi/2;
  const Double_t low=-TMath::Pi()/2.0,high=TMath::Pi()/2.0;
  const Float_t width=(high-low)/Float_t(bins);
  const Int_t bin=Int_t((half-low)/width);
  return bin>=0&&bin<bins?bin:-1;
}

struct Result {
  Bool_t valid;
  Int_t status, bin, runIndex, calibrationCent, hits;
  Double_t psi1, rawQx, rawQy, recenteredQx, recenteredQy, weight;
  Result() : valid(kFALSE), status(0), bin(-1), runIndex(-1), calibrationCent(-1), hits(0),
      psi1(-999.), rawQx(0), rawQy(0), recenteredQx(0), recenteredQy(0), weight(0) {}
};

struct Settings {
  bool enabled, writeDetailedLedger;
  std::string path, recenterFile, shiftFile, runListFile;
  std::string recenterSha256, shiftSha256, runListSha256, definition, pairingDefinition;
  UInt_t randomSeed;
  int bins, harmonic, shiftTerms, tpcMinNHitsFit;
  double minNMip, maxTileWeight, minVz, maxVz, vertexCenterX, vertexCenterY,
      maxVertexRadius, tpcMinPt, tpcMaxPt, tpcMinFitRatio, tpcMaxDca;
  std::vector<int> centralityThresholds, triggerIds;
  std::vector<std::vector<int> > epdRowGroups;
  std::vector<std::vector<double> > tpcEtaRanges;
  std::vector<double> pileupMinCoefficients, pileupMaxCoefficients;
  Settings() : enabled(false), writeDetailedLedger(false), randomSeed(0), bins(0), harmonic(0), shiftTerms(0), tpcMinNHitsFit(0),
      minNMip(0), maxTileWeight(0), minVz(0), maxVz(0), vertexCenterX(0), vertexCenterY(0),
      maxVertexRadius(0), tpcMinPt(0), tpcMaxPt(0), tpcMinFitRatio(0), tpcMaxDca(0) {}
  void Load(const std::string& mainconf);
};

class Provider {
 public:
  Provider();
  ~Provider();
  void Init(const Settings& settings, bool loadCalibration);
  Result Calculate(StPicoDst* pico) const;
  const Settings& Configuration() const { return mSettings; }
 private:
  Settings mSettings;
  std::vector<int> mRuns;
  TFile* mRecenter;
  TFile* mShift;
  TProfile2D* mTpc[2][2][2];
  TProfile2D* mEpd[2][7][2];
  TProfile3D* mShiftSin;
  TProfile3D* mShiftCos;
  Provider(const Provider&);
  Provider& operator=(const Provider&);
};

std::string Sha256File(const std::string& path);
} // namespace femto_epd
#endif
