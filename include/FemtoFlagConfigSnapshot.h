#ifndef FEMTO_FLAG_CONFIG_SNAPSHOT_H
#define FEMTO_FLAG_CONFIG_SNAPSHOT_H

// ---------------------------------------------------------------------------
// Which config values the reduced tree's selFlags depend on, and how to compare them.
//
// The tree stores pre-computed decisions (kSelTrackQualityNom, kSelNominalFemto, ...) so that a
// nominal selection reproduces the maker exactly. The cost of that design is that every flag
// silently BAKES IN the config the maker ran with. A downstream configured differently reads the
// flag and gets the producer's answer with no warning.
//
// This was not hypothetical. The tree-study configs inherited deuteronTofMomentumThreshold: 99.0
// from the dE/dx-only study (TOF effectively not required), while production uses 0.0 (TOF and a
// mass2 window required). Running the downstream with the production value against a tree built
// with the study value returned 554,931 deuterons instead of ~240,912 -- a factor 2.3, silently.
//
// So the maker writes this snapshot into the tree and the downstream refuses to run when its own
// configuration disagrees. Values are formatted as text with %.10g so the comparison is exact and
// covers strings and booleans as well as numbers.
//
// Adding a cut that feeds a flag means adding it here. A key that is absent from the snapshot is
// a key nobody is checking.
// ---------------------------------------------------------------------------

#include "TString.h"
#include <string>
#include <utility>
#include <vector>

#include "ConfigManager.h"
#include "cuts/FemtoConfig.h"
#include "cuts/NuclearIdCutConfig.h"
#include "cuts/PIDCutConfig.h"
#include "cuts/PhiCutConfig.h"
#include "cuts/TrackCutConfig.h"
#include "cuts/MixingConfig.h"
#include "cuts/EventCutConfig.h"

namespace femto_flag_config {

typedef std::pair<TString, TString> Entry;

inline TString Num(Double_t v) { return TString::Format("%.10g", v); }
inline TString Boolean(Bool_t v) { return v ? "true" : "false"; }

// Every value that feeds a selFlags bit, in a fixed order.
inline std::vector<Entry> Collect() {
  ConfigManager& cm = ConfigManager::GetInstance();
  const TrackCutConfig& tr = cm.GetTrackCuts();
  const PIDCutConfig& pid = cm.GetPIDCuts();
  const PhiCutConfig& phi = cm.GetPhiCuts();
  const NuclearIdCutConfig& nuc = cm.GetNuclearIdCuts();
  const FemtoConfig& fc = cm.GetFemtoConfig();
  const MixingConfig& mix = cm.GetMixingConfig();
  const EventCutConfig& evt = cm.GetEventCuts();
  std::vector<Entry> v;

  // The event row stores mixBin, which is a maker decision exactly like a selFlags bit, so the
  // bin geometry that produced it has to match too. bufferSize, mixingMode,
  // maxMixedPairsPerEvent and mixBothDirections are applied downstream and are deliberately NOT
  // here -- they are free to differ.
  v.push_back(Entry("mixing.nVzBins", Num(mix.nVzBins)));
  v.push_back(Entry("mixing.nCentralityBins", Num(mix.nCentralityBins)));
  v.push_back(Entry("mixing.nEventPlaneBins", Num(mix.nEventPlaneBins)));
  v.push_back(Entry("event.minVz", Num(evt.minVz)));
  v.push_back(Entry("event.maxVz", Num(evt.maxVz)));

  // kSelTrackQualityNom, kSelTrackChi2Nom
  v.push_back(Entry("track.minNHitsFit", Num(tr.minNHitsFit)));
  v.push_back(Entry("track.minNHitsRatio", Num(tr.minNHitsRatio)));
  v.push_back(Entry("track.minNHitsDedx", Num(tr.minNHitsDedx)));
  v.push_back(Entry("track.maxDCA", Num(tr.maxDCA)));
  v.push_back(Entry("track.minEta", Num(tr.minEta)));
  v.push_back(Entry("track.maxEta", Num(tr.maxEta)));
  v.push_back(Entry("track.minPt", Num(tr.minPt)));
  v.push_back(Entry("track.maxPt", Num(tr.maxPt)));
  v.push_back(Entry("track.maxChi2", Num(tr.maxChi2)));
  v.push_back(Entry("track.requirePrimaryTrack", Boolean(tr.requirePrimaryTrack)));

  // kSelKaonCutsNom
  v.push_back(Entry("phi.nSigmaKaon", Num(phi.nSigmaKaon)));
  v.push_back(Entry("phi.maxDCAKaon", Num(phi.maxDCAKaon)));

  // kSelLoosePid, kSelNominalPid (kaons)
  v.push_back(Entry("pid.requireTOF", Boolean(pid.requireTOF)));
  v.push_back(Entry("pid.tofFallbackMode", TString(pid.tofFallbackMode.c_str())));
  v.push_back(Entry("pid.pTofFallbackMax", Num(pid.pTofFallbackMax)));
  v.push_back(Entry("pid.tofUseMass2Cut", Boolean(pid.tofUseMass2Cut)));
  v.push_back(Entry("pid.minMass2Kaon", Num(pid.minMass2Kaon)));
  v.push_back(Entry("pid.maxMass2Kaon", Num(pid.maxMass2Kaon)));
  v.push_back(Entry("pid.tofUseDeltaInvBetaCut", Boolean(pid.tofUseDeltaInvBetaCut)));
  v.push_back(Entry("pid.maxAbsDeltaOneOverBetaKaon", Num(pid.maxAbsDeltaOneOverBetaKaon)));
  v.push_back(Entry("pid.pMomKaonPID", Num(pid.pMomKaonPID)));
  v.push_back(Entry("pid.phiDaughterKaonMinusRequireTof", Boolean(pid.phiDaughterKaonMinusRequireTof)));
  v.push_back(Entry("pid.nSigmaProton", Num(pid.nSigmaProton)));

  // kSelNominalPid (nuclei)
  v.push_back(Entry("nuclearId.minNHitsDedxNuclear", Num(nuc.minNHitsDedxNuclear)));
  v.push_back(Entry("nuclearId.maxNSigmaNuclear", Num(nuc.maxNSigmaNuclear)));
  v.push_back(Entry("nuclearId.nSigmaExclude", Num(nuc.nSigmaExclude)));
  v.push_back(Entry("nuclearId.m2_selection", Boolean(nuc.m2_selection)));
  v.push_back(Entry("nuclearId.m2SigmaCut", Num(nuc.m2SigmaCut)));
  v.push_back(Entry("nuclearId.minP_M2cut", Num(nuc.minP_M2cut)));
  v.push_back(Entry("nuclearId.minM2_M2cut", Num(nuc.minM2_M2cut)));
  v.push_back(Entry("nuclearId.maxM2_M2cut", Num(nuc.maxM2_M2cut)));
  v.push_back(Entry("nuclearId.maxPOverQ", Num(nuc.maxPOverQ)));
  v.push_back(Entry("nuclearId.requireBestSpecies", Boolean(nuc.requireBestSpecies)));

  // kSelNominalFemto (deuteron)
  v.push_back(Entry("femto.deuteronMaxDca", Num(fc.deuteronMaxDca)));
  v.push_back(Entry("femto.deuteronMinPMom", Num(fc.deuteronMinPMom)));
  v.push_back(Entry("femto.deuteronMaxPMom", Num(fc.deuteronMaxPMom)));
  v.push_back(Entry("femto.deuteronMinPtPre", Num(fc.deuteronMinPtPre)));
  v.push_back(Entry("femto.deuteronMaxPtPre", Num(fc.deuteronMaxPtPre)));
  v.push_back(Entry("femto.deuteronMinPtPair", Num(fc.deuteronMinPtPair)));
  v.push_back(Entry("femto.deuteronMaxPtPair", Num(fc.deuteronMaxPtPair)));
  v.push_back(Entry("femto.deuteronMaxAbsEta", Num(fc.deuteronMaxAbsEta)));
  v.push_back(Entry("femto.deuteronMaxAbsNSigma", Num(fc.deuteronMaxAbsNSigma)));
  v.push_back(Entry("femto.deuteronMinNHitsFit", Num(fc.deuteronMinNHitsFit)));
  v.push_back(Entry("femto.deuteronMinNHitsRatio", Num(fc.deuteronMinNHitsRatio)));
  v.push_back(Entry("femto.deuteronTofMomentumThreshold", Num(fc.deuteronTofMomentumThreshold)));
  v.push_back(Entry("femto.deuteronMinMass2", Num(fc.deuteronMinMass2)));
  v.push_back(Entry("femto.deuteronMaxMass2", Num(fc.deuteronMaxMass2)));
  v.push_back(Entry("femto.deuteronMinRapidityCm", Num(fc.deuteronMinRapidityCm)));
  v.push_back(Entry("femto.deuteronMaxRapidityCm", Num(fc.deuteronMaxRapidityCm)));

  // kSelNominalFemto (proton)
  v.push_back(Entry("femto.protonChargeMode", TString(fc.protonChargeMode.c_str())));
  v.push_back(Entry("femto.protonMaxDca", Num(fc.protonMaxDca)));
  v.push_back(Entry("femto.protonMinPtPre", Num(fc.protonMinPtPre)));
  v.push_back(Entry("femto.protonMinPtPair", Num(fc.protonMinPtPair)));
  v.push_back(Entry("femto.protonMaxPtPair", Num(fc.protonMaxPtPair)));
  v.push_back(Entry("femto.protonMaxAbsEta", Num(fc.protonMaxAbsEta)));
  v.push_back(Entry("femto.protonMaxAbsNSigma", Num(fc.protonMaxAbsNSigma)));
  v.push_back(Entry("femto.protonMinNHitsFit", Num(fc.protonMinNHitsFit)));
  v.push_back(Entry("femto.protonMinNHitsRatio", Num(fc.protonMinNHitsRatio)));
  v.push_back(Entry("femto.protonTofMomentumThreshold", Num(fc.protonTofMomentumThreshold)));
  v.push_back(Entry("femto.protonMinMass2", Num(fc.protonMinMass2)));
  v.push_back(Entry("femto.protonMaxMass2", Num(fc.protonMaxMass2)));
  v.push_back(Entry("femto.protonMinRapidityCm", Num(fc.protonMinRapidityCm)));
  v.push_back(Entry("femto.protonMaxRapidityCm", Num(fc.protonMaxRapidityCm)));


  // kSelNominalFemto (triton) -- StFemtoPhiTreeMaker::PassFemtoNuclearCuts
  v.push_back(Entry("femto.tritonMaxDca", Num(fc.tritonMaxDca)));
  v.push_back(Entry("femto.tritonMinPMom", Num(fc.tritonMinPMom)));
  v.push_back(Entry("femto.tritonMaxPMom", Num(fc.tritonMaxPMom)));
  v.push_back(Entry("femto.tritonMinPtPre", Num(fc.tritonMinPtPre)));
  v.push_back(Entry("femto.tritonMaxPtPre", Num(fc.tritonMaxPtPre)));
  v.push_back(Entry("femto.tritonMinPtPair", Num(fc.tritonMinPtPair)));
  v.push_back(Entry("femto.tritonMaxPtPair", Num(fc.tritonMaxPtPair)));
  v.push_back(Entry("femto.tritonMaxAbsEta", Num(fc.tritonMaxAbsEta)));
  v.push_back(Entry("femto.tritonMaxAbsNSigma", Num(fc.tritonMaxAbsNSigma)));
  v.push_back(Entry("femto.tritonMinNHitsFit", Num(fc.tritonMinNHitsFit)));
  v.push_back(Entry("femto.tritonMinNHitsRatio", Num(fc.tritonMinNHitsRatio)));
  v.push_back(Entry("femto.tritonTofMomentumThreshold", Num(fc.tritonTofMomentumThreshold)));
  v.push_back(Entry("femto.tritonMinMass2", Num(fc.tritonMinMass2)));
  v.push_back(Entry("femto.tritonMaxMass2", Num(fc.tritonMaxMass2)));
  v.push_back(Entry("femto.tritonMinRapidityCm", Num(fc.tritonMinRapidityCm)));
  v.push_back(Entry("femto.tritonMaxRapidityCm", Num(fc.tritonMaxRapidityCm)));

  // kSelNominalFemto (he3) -- StFemtoPhiTreeMaker::PassFemtoNuclearCuts
  v.push_back(Entry("femto.he3MaxDca", Num(fc.he3MaxDca)));
  v.push_back(Entry("femto.he3MinPMom", Num(fc.he3MinPMom)));
  v.push_back(Entry("femto.he3MaxPMom", Num(fc.he3MaxPMom)));
  v.push_back(Entry("femto.he3MinPtPre", Num(fc.he3MinPtPre)));
  v.push_back(Entry("femto.he3MaxPtPre", Num(fc.he3MaxPtPre)));
  v.push_back(Entry("femto.he3MinPtPair", Num(fc.he3MinPtPair)));
  v.push_back(Entry("femto.he3MaxPtPair", Num(fc.he3MaxPtPair)));
  v.push_back(Entry("femto.he3MaxAbsEta", Num(fc.he3MaxAbsEta)));
  v.push_back(Entry("femto.he3MaxAbsNSigma", Num(fc.he3MaxAbsNSigma)));
  v.push_back(Entry("femto.he3MinNHitsFit", Num(fc.he3MinNHitsFit)));
  v.push_back(Entry("femto.he3MinNHitsRatio", Num(fc.he3MinNHitsRatio)));
  v.push_back(Entry("femto.he3TofMomentumThreshold", Num(fc.he3TofMomentumThreshold)));
  v.push_back(Entry("femto.he3MinMass2", Num(fc.he3MinMass2)));
  v.push_back(Entry("femto.he3MaxMass2", Num(fc.he3MaxMass2)));
  v.push_back(Entry("femto.he3MinRapidityCm", Num(fc.he3MinRapidityCm)));
  v.push_back(Entry("femto.he3MaxRapidityCm", Num(fc.he3MaxRapidityCm)));

  // kSelNominalFemto (he4) -- StFemtoPhiTreeMaker::PassFemtoNuclearCuts
  v.push_back(Entry("femto.he4MaxDca", Num(fc.he4MaxDca)));
  v.push_back(Entry("femto.he4MinPMom", Num(fc.he4MinPMom)));
  v.push_back(Entry("femto.he4MaxPMom", Num(fc.he4MaxPMom)));
  v.push_back(Entry("femto.he4MinPtPre", Num(fc.he4MinPtPre)));
  v.push_back(Entry("femto.he4MaxPtPre", Num(fc.he4MaxPtPre)));
  v.push_back(Entry("femto.he4MinPtPair", Num(fc.he4MinPtPair)));
  v.push_back(Entry("femto.he4MaxPtPair", Num(fc.he4MaxPtPair)));
  v.push_back(Entry("femto.he4MaxAbsEta", Num(fc.he4MaxAbsEta)));
  v.push_back(Entry("femto.he4MaxAbsNSigma", Num(fc.he4MaxAbsNSigma)));
  v.push_back(Entry("femto.he4MinNHitsFit", Num(fc.he4MinNHitsFit)));
  v.push_back(Entry("femto.he4MinNHitsRatio", Num(fc.he4MinNHitsRatio)));
  v.push_back(Entry("femto.he4TofMomentumThreshold", Num(fc.he4TofMomentumThreshold)));
  v.push_back(Entry("femto.he4MinMass2", Num(fc.he4MinMass2)));
  v.push_back(Entry("femto.he4MaxMass2", Num(fc.he4MaxMass2)));
  v.push_back(Entry("femto.he4MinRapidityCm", Num(fc.he4MinRapidityCm)));
  v.push_back(Entry("femto.he4MaxRapidityCm", Num(fc.he4MaxRapidityCm)));

  // The rapidity windows above are evaluated in a frame the maker fixes at Init() from
  // phi.rapidityFrame, phi.sqrtSNNGeV and the centrality config. The window VALUES being equal is
  // not enough: the same deuteronMinRapidityCm means a different selection in a different frame,
  // so what is recorded here is the OUTCOME that ApplyAnalysisRapidity actually uses.
  v.push_back(Entry("phi.rapidityFrameEffective", TString(phi.rapidityFrameEffective.c_str())));
  v.push_back(Entry("phi.rapidityShiftEffective", Num(phi.rapidityShiftEffective)));
  v.push_back(Entry("phi.sqrtSNNGeV", Num(phi.sqrtSNNGeV)));

  return v;
}

}  // namespace femto_flag_config

#endif  // FEMTO_FLAG_CONFIG_SNAPSHOT_H
