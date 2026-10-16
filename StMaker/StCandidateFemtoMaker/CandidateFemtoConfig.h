#ifndef CANDIDATE_FEMTO_CONFIG_H
#define CANDIDATE_FEMTO_CONFIG_H

#include "FlatYamlConfig.h"

#include <iostream>
#include <set>
#include <string>
#include <vector>

// Configuration of StCandidateFemtoMaker: the YAML referenced by the mainconf key
// "candidateFemto" (flat keys, StFemtoMaker-style species/channel layout). Every key is
// required unless noted. Example: config/cuts/femto/candidateFemto_auau3_k0sXiKaonXi.yaml.
//
// Species (speciesKeys: k0s, xi_minus, kaon_plus; lower_snake_case):
//   species_<s>_source:          name of the CandidateSource that provides it (kfparticle, track_pid)
//   species_<s>_particleKey:     key of that source (k0s_kf, xi_minus_kf, kaon_plus_tof)
//   species_<s>_pdgMass:         mass used for k*, rapidity and m_T [GeV/c^2]
//   species_<s>_massAxis:        nBins, min, max of the mass-like variable (invariant mass, or TOF m^2
//                                for tracks). Candidates outside [min, max) are not used at all.
//   species_<s>_signal:          min, max of the signal window (inclusive)
//   species_<s>_leftSideband:    min, max (optional; no left sideband class without it)
//   species_<s>_rightSideband:   min, max (optional)
//   species_<s>_daughterMasses:  optional; mass hypotheses of the daughters (source order) for the
//                                daughter pT-y acceptance maps
// Keys used only by the analysis macros (common/macro/candidate_femto/); the Maker reads and checks
// them too (same file), so a malformed value stops the Maker job:
//   species_<s>_purityModel:     peak shape of the purity fits: invariantMassGausPol2 (Gaussian + polynomial in the
//                                mass variable), invariantMassStudentTPol2 (Student-t + polynomial; the K0S and Xi-
//                                peaks have non-Gaussian tails), invariantMassDoubleGausPol2 (two
//                                Gaussians, common mean + polynomial) or tofMassSquaredKaonPion (K Gaussian + pion
//                                Gaussian + pol1 in TOF m^2)
//   species_<s>_purityByDefinition: true = purity 1 by definition (e.g. a clean K+ selection), false = fitted
//   species_<s>_purityFitRange:  min, max of the purity fits of the purity-correction macro (inside the mass axis)
//   species_<s>_purityFitSeed:   mean, sigma start values of the purity fits
//   species_<s>_qaPurityFitSigmaRange / species_<s>_qaPurityBackgroundOrder: QA macro (checkHist) purity fit,
//                                local range mean +- N sigma and polynomial order of the background; required
//                                when purityByDefinition is false
//   species_<s>_scoreMin:        optional BDT score threshold of the CFs (absent = no threshold)
// Channels (nChannels: N; channel_<i>_name / _partA / _partB; name = partA_partB by convention).
// Pair histograms: THnSparseF hPairSE_<channel> / hPairME_<channel> with axes
//   (k*, cent9, mass class A, mass class B, mass A, mass B, BDT score A, BDT score B),
//   mass class 0 = none, 1 = signal, 2 = left sideband, 3 = right sideband.
// Shared tracks (always on, no key), per channel: a candidate of species A and a candidate of species B of the
//   same event that share a daughter track are both removed from that channel (its pairs, pools and hChannel*
//   histograms) when at least one of them is in its signal window; other channels are not affected
//   (StCandidateFemtoMaker.cxx, ApplySharedTrackRule). The candidate trees keep every usable candidate and carry
//   one flag per channel of the species (sharedTrackRemoved_<channel>).
// Others: kstarAxis (nBins, min, max); yMidLab; closePairRadiusMin/Max/Step, closePairTpcHalfLength [m];
//   trackQaMinPt / trackQaMinNHitsFit / trackQaMaxDca (all-track PID QA); writeCandidateTrees (BDT training trees);
//   candidateTreeSpecies (required if writeCandidateTrees is true, forbidden otherwise: the species that get a tree);
//   writeCandidateStore: true = reconstruction stage (read the PicoDst once, write every candidate of every part
//   with the event records into the job output, CandidateStore.h; no pairs, no BDT, no trees: the pair stage
//   makes them from the store); false = pair stage, or one PicoDst pass doing everything.
// Analysis-macro keys: cfNormMin/Max (CF normalisation range); qaCfKstarBinWidth (k* bin width of the QA macro
//   CFs); qaClosePairMaxDeltaEta (|delta eta| of the close-pair delta phi* QA); purityFitMinEntries (projections
//   with fewer entries are not fitted); slidingWindowHalfWidth (sliding-window purity: fine k* bins on each side);
//   purityCorrectionMaxMergedBins (k* bins merged at most where a purity fit fails).
struct CandidateFemtoAxis {
  CandidateFemtoAxis() : nBins(0), min(0.), max(0.) {}
  Int_t nBins;
  Double_t min, max;
  bool Contains(Double_t v) const { return v >= min && v < max; }
};

struct CandidateFemtoWindow {
  CandidateFemtoWindow() : defined(false), min(0.), max(0.) {}
  bool defined;
  Double_t min, max;
  bool Contains(Double_t v) const { return defined && v >= min && v <= max; }
};

struct CandidateFemtoSpecies {
  CandidateFemtoSpecies()
      : pdgMass(0.), purityByDefinition(false), qaPurityFitSigmaRange(0.), qaPurityBackgroundOrder(-1), hasScoreMin(false),
        scoreMin(-1e9) {}
  std::string key, source, particleKey;
  Double_t pdgMass;
  CandidateFemtoAxis massAxis;
  CandidateFemtoWindow signal, leftSideband, rightSideband;
  std::vector<Double_t> daughterMasses;
  std::string purityModel;
  bool purityByDefinition;
  CandidateFemtoWindow purityFitRange, purityFitSeed;  // purityFitSeed: min = mean, max = sigma
  Double_t qaPurityFitSigmaRange;  // 0 if purityByDefinition
  Int_t qaPurityBackgroundOrder;   // -1 if purityByDefinition
  bool hasScoreMin;
  Double_t scoreMin;
  // 0 = none, 1 = signal, 2 = left sideband, 3 = right sideband.
  Int_t MassClass(Double_t m) const {
    if (signal.Contains(m)) return 1;
    if (leftSideband.Contains(m)) return 2;
    if (rightSideband.Contains(m)) return 3;
    return 0;
  }
};

struct CandidateFemtoChannel {
  std::string name, partA, partB;
  Int_t indexA, indexB;  // into CandidateFemtoConfig::species
};

class CandidateFemtoConfig {
public:
  CandidateFemtoConfig()
      : cfNormMin(0.), cfNormMax(0.), yMidLab(0.), closePairRadiusMin(0.), closePairRadiusMax(0.),
        closePairRadiusStep(0.), closePairTpcHalfLength(0.), trackQaMinPt(0.), trackQaMinNHitsFit(0.),
        trackQaMaxDca(0.), writeCandidateTrees(false), writeCandidateStore(false), qaCfKstarBinWidth(0.), qaClosePairMaxDeltaEta(0.),
        purityFitMinEntries(0), slidingWindowHalfWidth(-1), purityCorrectionMaxMergedBins(0) {}

  std::vector<CandidateFemtoSpecies> species;
  std::vector<CandidateFemtoChannel> channels;
  CandidateFemtoAxis kstarAxis;
  Double_t cfNormMin, cfNormMax;
  Double_t yMidLab;
  Double_t closePairRadiusMin, closePairRadiusMax, closePairRadiusStep, closePairTpcHalfLength;
  Double_t trackQaMinPt, trackQaMinNHitsFit, trackQaMaxDca;
  bool writeCandidateTrees;
  std::vector<std::string> candidateTreeSpecies;
  bool writeCandidateStore;
  Double_t qaCfKstarBinWidth, qaClosePairMaxDeltaEta;
  Int_t purityFitMinEntries, slidingWindowHalfWidth, purityCorrectionMaxMergedBins;
  std::string path;

  Int_t SpeciesIndex(const std::string& key) const {
    for (size_t i = 0; i < species.size(); ++i) if (species[i].key == key) return static_cast<Int_t>(i);
    return -1;
  }

  bool LoadFromMainConfig(const std::string& mainConfigPath) {
    if (!FlatYamlConfig::ResolveMainConfigReference(mainConfigPath, "candidateFemto", path)) return false;
    return Load(path);
  }

  bool Load(const std::string& file) {
    path = file;
    FlatYamlConfig c;
    if (!c.Load(file)) return false;
    const std::vector<std::string> keys = c.RequireList("speciesKeys");
    std::set<std::string> seen;
    for (size_t i = 0; i < keys.size(); ++i) {
      CandidateFemtoSpecies s;
      s.key = keys[i];
      if (!seen.insert(s.key).second) c.Fail("duplicate species key '" + s.key + "'");
      const std::string p = "species_" + s.key + "_";
      s.source = c.RequireString(p + "source");
      s.particleKey = c.RequireString(p + "particleKey");
      s.pdgMass = c.RequireDouble(p + "pdgMass");
      ReadAxis(c, p + "massAxis", s.massAxis);
      ReadWindow(c, p + "signal", true, s.signal);
      ReadWindow(c, p + "leftSideband", false, s.leftSideband);
      ReadWindow(c, p + "rightSideband", false, s.rightSideband);
      const std::vector<std::string> dm = c.OptionalList(p + "daughterMasses");
      for (size_t j = 0; j < dm.size(); ++j) {
        bool ok = false;
        s.daughterMasses.push_back(FlatYamlConfig::ParseDouble(dm[j], ok));
        if (!ok || s.daughterMasses.back() <= 0.) c.Fail(p + "daughterMasses: invalid value '" + dm[j] + "'");
      }
      s.purityModel = c.RequireString(p + "purityModel");
      if (s.purityModel != "invariantMassGausPol2" && s.purityModel != "invariantMassStudentTPol2" &&
          s.purityModel != "invariantMassDoubleGausPol2" && s.purityModel != "tofMassSquaredKaonPion")
        c.Fail(p + "purityModel must be invariantMassGausPol2, invariantMassStudentTPol2, invariantMassDoubleGausPol2 "
                   "or tofMassSquaredKaonPion");
      s.purityByDefinition = c.RequireBool(p + "purityByDefinition");
      ReadWindow(c, p + "purityFitRange", true, s.purityFitRange);
      if (s.purityFitRange.defined && (s.purityFitRange.min < s.massAxis.min || s.purityFitRange.max > s.massAxis.max))
        c.Fail(p + "purityFitRange must lie inside massAxis");
      // the purity macros project only the fit range: the signal window must lie inside it
      if (s.signal.min < s.purityFitRange.min || s.signal.max > s.purityFitRange.max)
        c.Fail(p + "signal must lie inside purityFitRange");
      {
        const std::vector<std::string> v = c.RequireList(p + "purityFitSeed");
        bool ok0 = false, ok1 = false;
        if (v.size() == 2) {
          s.purityFitSeed.min = FlatYamlConfig::ParseDouble(v[0], ok0);
          s.purityFitSeed.max = FlatYamlConfig::ParseDouble(v[1], ok1);
        }
        s.purityFitSeed.defined = ok0 && ok1 && s.purityFitSeed.max > 0.;
        if (!s.purityFitSeed.defined) c.Fail(p + "purityFitSeed must be: mean, sigma (sigma > 0)");
      }
      if (!s.purityByDefinition) {
        s.qaPurityFitSigmaRange = c.RequireDouble(p + "qaPurityFitSigmaRange");
        s.qaPurityBackgroundOrder = c.RequireInt(p + "qaPurityBackgroundOrder");
        if (!(s.qaPurityFitSigmaRange > 0.) || s.qaPurityBackgroundOrder < 0 || s.qaPurityBackgroundOrder > 9)
          c.Fail(p + "qaPurityFitSigmaRange must be > 0 and qaPurityBackgroundOrder in 0-9");
      }
      s.hasScoreMin = c.Has(p + "scoreMin");
      s.scoreMin = s.hasScoreMin ? c.RequireDouble(p + "scoreMin") : -1e9;
      if (!(s.pdgMass > 0.)) c.Fail(p + "pdgMass must be > 0");
      const CandidateFemtoWindow* w[3] = {&s.signal, &s.leftSideband, &s.rightSideband};
      for (Int_t k = 0; k < 3; ++k)
        if (w[k]->defined && (w[k]->min < s.massAxis.min || w[k]->max >= s.massAxis.max))
          c.Fail(p + ": every window must lie inside massAxis [min, max)");
      if ((s.leftSideband.defined && s.leftSideband.max > s.signal.min) ||
          (s.rightSideband.defined && s.rightSideband.min < s.signal.max))
        c.Fail(p + ": a sideband overlaps the signal window");
      species.push_back(s);
    }
    const Int_t nChannels = c.RequireInt("nChannels");
    if (nChannels < 1) c.Fail("nChannels must be >= 1");
    std::set<std::string> names;
    for (Int_t i = 0; i < nChannels; ++i) {
      CandidateFemtoChannel ch;
      const std::string p = "channel_" + ToString(i) + "_";
      ch.name = c.RequireString(p + "name");
      ch.partA = c.RequireString(p + "partA");
      ch.partB = c.RequireString(p + "partB");
      ch.indexA = SpeciesIndex(ch.partA);
      ch.indexB = SpeciesIndex(ch.partB);
      if (ch.indexA < 0 || ch.indexB < 0) c.Fail(p + ": partA/partB must be listed in speciesKeys");
      if (ch.indexA == ch.indexB) c.Fail(p + ": identical-species channels are not supported");
      if (!names.insert(ch.name).second) c.Fail("duplicate channel name '" + ch.name + "'");
      channels.push_back(ch);
    }
    ReadAxis(c, "kstarAxis", kstarAxis);
    cfNormMin = c.RequireDouble("cfNormMin");
    cfNormMax = c.RequireDouble("cfNormMax");
    yMidLab = c.RequireDouble("yMidLab");
    closePairRadiusMin = c.RequireDouble("closePairRadiusMin");
    closePairRadiusMax = c.RequireDouble("closePairRadiusMax");
    closePairRadiusStep = c.RequireDouble("closePairRadiusStep");
    closePairTpcHalfLength = c.RequireDouble("closePairTpcHalfLength");
    trackQaMinPt = c.RequireDouble("trackQaMinPt");
    trackQaMinNHitsFit = c.RequireDouble("trackQaMinNHitsFit");
    trackQaMaxDca = c.RequireDouble("trackQaMaxDca");
    writeCandidateTrees = c.RequireBool("writeCandidateTrees");
    candidateTreeSpecies.clear();
    if (c.Ok() && writeCandidateTrees) {
      candidateTreeSpecies = c.RequireList("candidateTreeSpecies");
      for (size_t i = 0; i < candidateTreeSpecies.size() && c.Ok(); ++i) {
        if (SpeciesIndex(candidateTreeSpecies[i]) < 0) c.Fail("candidateTreeSpecies: '" + candidateTreeSpecies[i] + "' is not in speciesKeys");
        for (size_t j = 0; j < i; ++j)
          if (candidateTreeSpecies[j] == candidateTreeSpecies[i]) c.Fail("candidateTreeSpecies: '" + candidateTreeSpecies[i] + "' twice");
      }
      if (c.Ok() && candidateTreeSpecies.empty()) c.Fail("candidateTreeSpecies must list at least one species");
    }
    writeCandidateStore = c.RequireBool("writeCandidateStore");
    if (c.Ok() && writeCandidateStore && writeCandidateTrees)
      c.Fail("writeCandidateStore and writeCandidateTrees are exclusive (the pair stage writes the trees from the store)");
    qaCfKstarBinWidth = c.RequireDouble("qaCfKstarBinWidth");
    qaClosePairMaxDeltaEta = c.RequireDouble("qaClosePairMaxDeltaEta");
    purityFitMinEntries = c.RequireInt("purityFitMinEntries");
    slidingWindowHalfWidth = c.RequireInt("slidingWindowHalfWidth");
    purityCorrectionMaxMergedBins = c.RequireInt("purityCorrectionMaxMergedBins");
    if (!(qaCfKstarBinWidth > 0.) || !(qaClosePairMaxDeltaEta > 0.) || purityFitMinEntries < 1 ||
        slidingWindowHalfWidth < 0 || purityCorrectionMaxMergedBins < 1)
      c.Fail("qaCfKstarBinWidth and qaClosePairMaxDeltaEta must be > 0, purityFitMinEntries >= 1, "
             "slidingWindowHalfWidth >= 0, purityCorrectionMaxMergedBins >= 1");
    if (!(cfNormMin >= 0. && cfNormMin < cfNormMax)) c.Fail("invalid cfNormMin/cfNormMax");
    if (!(closePairRadiusMin > 0. && closePairRadiusMin < closePairRadiusMax && closePairRadiusStep > 0. &&
          closePairTpcHalfLength > 0.))
      c.Fail("invalid close-pair radius scan");
    c.RejectUnused();
    return c.Ok();
  }

  void Print(std::ostream& out) const {
    out << "[CandidateFemtoConfig] " << path << "\n";
    for (size_t i = 0; i < species.size(); ++i) {
      const CandidateFemtoSpecies& s = species[i];
      out << "  species " << s.key << ": " << s.source << "/" << s.particleKey << " pdgMass " << s.pdgMass
          << " massAxis " << s.massAxis.nBins << " [" << s.massAxis.min << ", " << s.massAxis.max << ")"
          << " signal [" << s.signal.min << ", " << s.signal.max << "]";
      if (s.leftSideband.defined) out << " leftSB [" << s.leftSideband.min << ", " << s.leftSideband.max << "]";
      if (s.rightSideband.defined) out << " rightSB [" << s.rightSideband.min << ", " << s.rightSideband.max << "]";
      out << "\n";
    }
    for (size_t i = 0; i < channels.size(); ++i)
      out << "  channel " << channels[i].name << " = " << channels[i].partA << " x " << channels[i].partB << "\n";
    out << "  kstarAxis " << kstarAxis.nBins << " [" << kstarAxis.min << ", " << kstarAxis.max << ")"
        << " writeCandidateTrees " << writeCandidateTrees;
    for (size_t i = 0; i < candidateTreeSpecies.size(); ++i) out << (i ? ", " : " (") << candidateTreeSpecies[i];
    out << (candidateTreeSpecies.empty() ? "" : ")") << " writeCandidateStore " << writeCandidateStore << "\n";
  }

private:
  static std::string ToString(Int_t i) {
    std::ostringstream s;
    s << i;
    return s.str();
  }
  static void ReadAxis(FlatYamlConfig& c, const std::string& key, CandidateFemtoAxis& axis) {
    const std::vector<std::string> v = c.RequireList(key);
    bool ok0 = false, ok1 = false, ok2 = false;
    if (v.size() == 3) {
      axis.nBins = static_cast<Int_t>(FlatYamlConfig::ParseDouble(v[0], ok0));
      axis.min = FlatYamlConfig::ParseDouble(v[1], ok1);
      axis.max = FlatYamlConfig::ParseDouble(v[2], ok2);
    }
    if (!(ok0 && ok1 && ok2 && axis.nBins > 0 && axis.min < axis.max)) c.Fail(key + " must be: nBins, min, max");
  }
  static void ReadWindow(FlatYamlConfig& c, const std::string& key, bool required, CandidateFemtoWindow& w) {
    if (!required && !c.Has(key)) return;
    const std::vector<std::string> v = c.RequireList(key);
    bool ok0 = false, ok1 = false;
    if (v.size() == 2) {
      w.min = FlatYamlConfig::ParseDouble(v[0], ok0);
      w.max = FlatYamlConfig::ParseDouble(v[1], ok1);
    }
    w.defined = ok0 && ok1 && w.min < w.max;
    if (!w.defined) c.Fail(key + " must be: min, max (min < max)");
  }
};

#endif
