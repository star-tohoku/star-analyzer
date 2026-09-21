#ifndef DATA006_CONFIG_H
#define DATA006_CONFIG_H

#include "Rtypes.h"
#include "YamlParser.h"
#include <cstdio>

#include <iostream>
#include <map>
#include <string>
#include <vector>

// Configuration used only by the reduced-tree DATA-006 reader/finalizer.  These values do not
// participate in FemtoFlagConfig: they change histogram construction and pair bookkeeping, not
// any producer-side selFlags decision.  The producer mainconf is still loaded and checked first.
struct Data006ChannelConfig {
  std::string name;
  std::string partA;
  std::string partB;
  Bool_t enabled;
  Bool_t doMixing;
  Double_t normQMin;
  Double_t normQMax;
  Double_t closePairDEta;
  Double_t closePairDPhiStar;
};

struct Data006CentralityConfig {
  std::string id;
  Int_t lowPercent;
  Int_t highPercent;
  Int_t cent9Min;
  Int_t cent9Max;
  Int_t refmultLow;
  Int_t refmultHigh;
  Bool_t nativeClass;
};

class Data006Config {
 public:
  Int_t kstarBins;
  Double_t kstarMin;
  Double_t kstarMax;
  Double_t deliveryKstarMax;

  Int_t pairMtBins;
  Double_t pairMtMin;
  Double_t pairMtMax;
  Int_t pairRapidityBins;
  Double_t pairRapidityMin;
  Double_t pairRapidityMax;

  Bool_t closePairEnabled;
  std::string closePairShape;
  Double_t closePairRadiusMin;
  Double_t closePairRadiusMax;
  Double_t closePairRadiusStep;

  std::vector<Data006ChannelConfig> channels;
  std::vector<Data006CentralityConfig> centralities;

  Bool_t Load(const Char_t* filename) {
    std::map<std::string, std::string> v;
    if (!filename || !filename[0] || !YamlParser::ParseFile(filename, v)) return kFALSE;

    Bool_t ok = kTRUE;
    kstarBins = RequiredInt(v, "kstarBins", ok);
    kstarMin = RequiredDouble(v, "kstarMin", ok);
    kstarMax = RequiredDouble(v, "kstarMax", ok);
    deliveryKstarMax = RequiredDouble(v, "deliveryKstarMax", ok);
    pairMtBins = RequiredInt(v, "pairMtBins", ok);
    pairMtMin = RequiredDouble(v, "pairMtMin", ok);
    pairMtMax = RequiredDouble(v, "pairMtMax", ok);
    pairRapidityBins = RequiredInt(v, "pairRapidityBins", ok);
    pairRapidityMin = RequiredDouble(v, "pairRapidityMin", ok);
    pairRapidityMax = RequiredDouble(v, "pairRapidityMax", ok);
    closePairEnabled = RequiredBool(v, "closePairEnabled", ok);
    closePairShape = Required(v, "closePairShape", ok);
    closePairRadiusMin = RequiredDouble(v, "closePairRadiusMin", ok);
    closePairRadiusMax = RequiredDouble(v, "closePairRadiusMax", ok);
    closePairRadiusStep = RequiredDouble(v, "closePairRadiusStep", ok);

    const Int_t nChannels = RequiredInt(v, "nChannels", ok);
    channels.clear();
    for (Int_t i = 0; i < nChannels; ++i) {
      const std::string p = "channel_" + Number(i) + "_";
      Data006ChannelConfig c;
      c.name = Required(v, p + "channelName", ok);
      c.partA = Required(v, p + "partA", ok);
      c.partB = Required(v, p + "partB", ok);
      c.enabled = RequiredBool(v, p + "enabled", ok);
      c.doMixing = RequiredBool(v, p + "doMixing", ok);
      c.normQMin = RequiredDouble(v, p + "normQMin", ok);
      c.normQMax = RequiredDouble(v, p + "normQMax", ok);
      c.closePairDEta = RequiredDouble(v, p + "closePairDEta", ok);
      c.closePairDPhiStar = RequiredDouble(v, p + "closePairDPhiStar", ok);
      channels.push_back(c);
    }

    const Int_t nCentralities = RequiredInt(v, "nCentralityClasses", ok);
    centralities.clear();
    for (Int_t i = 0; i < nCentralities; ++i) {
      const std::string p = "centrality_" + Number(i) + "_";
      Data006CentralityConfig c;
      c.id = Required(v, p + "id", ok);
      c.lowPercent = RequiredInt(v, p + "lowPercent", ok);
      c.highPercent = RequiredInt(v, p + "highPercent", ok);
      c.cent9Min = RequiredInt(v, p + "cent9Min", ok);
      c.cent9Max = RequiredInt(v, p + "cent9Max", ok);
      c.refmultLow = RequiredInt(v, p + "refmultLow", ok);
      c.refmultHigh = RequiredInt(v, p + "refmultHigh", ok);
      c.nativeClass = RequiredBool(v, p + "nativeClass", ok);
      centralities.push_back(c);
    }

    if (!ok) return kFALSE;
    return Validate(filename);
  }

 private:
  static std::string Number(Int_t x) {
    char buf[32];
    std::sprintf(buf, "%d", (int)x);
    return std::string(buf);
  }

  static std::string Required(const std::map<std::string, std::string>& v,
                              const std::string& key, Bool_t& ok) {
    std::map<std::string, std::string>::const_iterator it = v.find(key);
    if (it == v.end() || it->second.empty()) {
      std::cerr << "ERROR: DATA-006 config is missing required key '" << key << "'" << std::endl;
      ok = kFALSE;
      return "";
    }
    return it->second;
  }

  static Int_t RequiredInt(const std::map<std::string, std::string>& v,
                           const std::string& key, Bool_t& ok) {
    const std::string s = Required(v, key, ok);
    return s.empty() ? 0 : YamlParser::ToInt(s, 0);
  }

  static Double_t RequiredDouble(const std::map<std::string, std::string>& v,
                                 const std::string& key, Bool_t& ok) {
    const std::string s = Required(v, key, ok);
    return s.empty() ? 0.0 : YamlParser::ToDouble(s, 0.0);
  }

  static Bool_t RequiredBool(const std::map<std::string, std::string>& v,
                             const std::string& key, Bool_t& ok) {
    const std::string s = Required(v, key, ok);
    return s.empty() ? kFALSE : YamlParser::ToBool(s, kFALSE);
  }

  Bool_t Validate(const Char_t* filename) const {
    Bool_t ok = kTRUE;
    if (kstarBins <= 0 || !(kstarMax > kstarMin) || deliveryKstarMax <= kstarMin ||
        deliveryKstarMax > kstarMax) {
      std::cerr << "ERROR: invalid DATA-006 k* axes in " << filename << std::endl;
      ok = kFALSE;
    }
    if (pairMtBins <= 0 || !(pairMtMax > pairMtMin) || pairRapidityBins <= 0 ||
        !(pairRapidityMax > pairRapidityMin)) {
      std::cerr << "ERROR: invalid DATA-006 pair-kinematic axes in " << filename << std::endl;
      ok = kFALSE;
    }
    if (closePairShape != "box" && closePairShape != "ellipse") {
      std::cerr << "ERROR: DATA-006 closePairShape must be box or ellipse" << std::endl;
      ok = kFALSE;
    }
    if (closePairRadiusMin <= 0 || closePairRadiusMax < closePairRadiusMin ||
        closePairRadiusStep <= 0) {
      std::cerr << "ERROR: invalid DATA-006 close-pair radius scan" << std::endl;
      ok = kFALSE;
    }
    if (channels.empty()) {
      std::cerr << "ERROR: DATA-006 config defines no channels" << std::endl;
      ok = kFALSE;
    }
    std::map<std::string, Bool_t> names;
    for (size_t i = 0; i < channels.size(); ++i) {
      const Data006ChannelConfig& c = channels[i];
      if ((c.partA != "proton" && c.partA != "deuteron") ||
          (c.partB != "proton" && c.partB != "deuteron")) {
        std::cerr << "ERROR: DATA-006 channel '" << c.name
                  << "' supports proton/deuteron track species only" << std::endl;
        ok = kFALSE;
      }
      if (c.name != c.partA + "_" + c.partB) {
        std::cerr << "ERROR: DATA-006 channel name '" << c.name
                  << "' must follow {partA}_{partB}" << std::endl;
        ok = kFALSE;
      }
      if (names[c.name]) {
        std::cerr << "ERROR: duplicate DATA-006 channel '" << c.name << "'" << std::endl;
        ok = kFALSE;
      }
      names[c.name] = kTRUE;
      if (!(c.normQMax > c.normQMin) || c.normQMin < kstarMin || c.normQMax > kstarMax) {
        std::cerr << "ERROR: invalid normalization range for DATA-006 channel '" << c.name << "'"
                  << std::endl;
        ok = kFALSE;
      }
      if (c.closePairDEta < 0 || c.closePairDPhiStar < 0) {
        std::cerr << "ERROR: negative close-pair window for DATA-006 channel '" << c.name << "'"
                  << std::endl;
        ok = kFALSE;
      }
    }
    std::map<std::string, Bool_t> centNames;
    for (size_t i = 0; i < centralities.size(); ++i) {
      const Data006CentralityConfig& c = centralities[i];
      if (c.id.empty() || centNames[c.id] || c.lowPercent < 0 ||
          c.highPercent <= c.lowPercent || c.cent9Min < 0 || c.cent9Max > 8 ||
          c.cent9Min > c.cent9Max) {
        std::cerr << "ERROR: invalid DATA-006 centrality class at index " << i << std::endl;
        ok = kFALSE;
      }
      centNames[c.id] = kTRUE;
    }
    return ok;
  }
};

#endif
