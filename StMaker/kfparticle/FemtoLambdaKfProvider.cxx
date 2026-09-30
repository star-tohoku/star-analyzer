#include "FemtoLambdaKfProvider.h"
#include "KfEventSelection.h"
#include "KfLambdaSelector.h"
#include "StPicoKFParticleInterface.h"
#include "ConfigManager.h"
#include "cuts/CentralityCutConfig.h"
#include "cuts/EventCutConfig.h"
#include "cuts/KfParticleCutConfig.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "TSystem.h"
#include "yaml-cpp/yaml.h"

#include <cmath>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>

namespace {
void Require(bool condition, const std::string& message) {
  if (!condition) throw std::runtime_error(message);
}
void UniqueMap(const YAML::Node& node, const std::string& location) {
  Require(node.IsMap(), location + " must be a mapping");
  std::set<std::string> keys;
  for (YAML::const_iterator it = node.begin(); it != node.end(); ++it) {
    Require(it->first.IsScalar(), location + " has a non-scalar key");
    const std::string key = it->first.as<std::string>();
    Require(keys.insert(key).second, location + " has duplicate key: " + key);
  }
}
std::string RequiredString(const YAML::Node& node, const std::string& location) {
  Require(node.IsScalar(), location + " must be an explicit string");
  const std::string result = node.as<std::string>();
  Require(!result.empty(), location + " must not be empty");
  return result;
}
std::string Reference(const std::string& main, const std::string& path) {
  if (path[0] == '/') return path;
  const size_t position = main.find("/config/");
  Require(position != std::string::npos,
          "relative config references require mainconf under project/config/");
  return main.substr(0, position + 1) + "config/" + path;
}
std::string AbsoluteMain(const char* path) {
  Require(path && path[0], "mainconf is required; no implicit default is allowed");
  if (path[0] == '/') return path;
  Require(gSystem && gSystem->WorkingDirectory(), "working directory is unavailable");
  return std::string(gSystem->WorkingDirectory()) + "/" + path;
}

class KfProvider : public FemtoLambdaProvider {
public:
  KfProvider() : mInterface(0), mCuts(0), mPairMass(0.), mInitialized(false) {}
  virtual ~KfProvider() { delete mInterface; }

  virtual bool Init(const char* mainconf, std::ostream& errors) {
    Clear();
    delete mInterface;
    mInterface = 0;
    mCuts = 0;
    mInitialized = false;
    try {
      const std::string mainPath = AbsoluteMain(mainconf);
      const YAML::Node main = YAML::LoadFile(mainPath);
      UniqueMap(main, "mainconf");
      const std::string makerPath = Reference(mainPath, RequiredString(main["maker"], "mainconf.maker"));
      const std::string kfPath = Reference(mainPath, RequiredString(main["kf"], "mainconf.kf"));
      const std::string eventPath = Reference(mainPath, RequiredString(main["event"], "mainconf.event"));
      const std::string centralityPath = Reference(mainPath, RequiredString(main["centrality"], "mainconf.centrality"));
      RequiredString(main["analysis"], "mainconf.analysis");
      const YAML::Node maker = YAML::LoadFile(makerPath);
      UniqueMap(maker, "maker");
      const std::string policy = RequiredString(maker["lambdaEventPolicy"], "maker.lambdaEventPolicy");
      Require(policy == "lambda_imp5_compat" || policy == "mode_vertex",
              "lambdaEventPolicy must be lambda_imp5_compat or mode_vertex");
      Require(RequiredString(maker["lambdaPairMassMode"], "maker.lambdaPairMassMode") == "fixed",
              "Lambda pair mass mode must be fixed (raw KF mass is saved separately)");
      Require(maker["lambdaPairMass"].IsScalar(), "maker.lambdaPairMass is required");
      mPairMass = maker["lambdaPairMass"].as<double>();
      Require(std::isfinite(mPairMass) && mPairMass > 0., "lambdaPairMass must be finite and positive");
      Require(RequiredString(maker["species_lambda_builderType"], "species_lambda_builderType") == "resonance" &&
              RequiredString(maker["species_lambda_particleKey"], "species_lambda_particleKey") == "lambda_kf",
              "Lambda species requires resonance / lambda_kf");

      KfParticleCutConfig& cuts = ConfigManager::GetInstance().GetKfParticleCuts();
      Require(cuts.LoadFromFile(kfPath.c_str()), "KF configuration load/validation failed");
      // Require every effective key, including explicitly disabled optional cuts.
      // KfParticleCutConfig::Dump omits only genuinely profile-inactive fields.
      const YAML::Node kf = YAML::LoadFile(kfPath);
      UniqueMap(kf, "KF");
      std::ostringstream effective;
      cuts.Dump(effective);
      const YAML::Node required = YAML::Load(effective.str());
      for (YAML::const_iterator it = required.begin(); it != required.end(); ++it) {
        const std::string key = it->first.as<std::string>();
        Require(kf[key].IsDefined(), "KF preset is missing effective key: " + key);
      }
      Require(!cuts.reconstructAntiLambda,
              "lambda_kf species is Lambda only; use an explicit future anti_lambda species for anti-Lambda");
      EventCutConfig& eventCuts = ConfigManager::GetInstance().GetEventCuts();
      Require(eventCuts.LoadFromFile(eventPath.c_str()), "event configuration cannot be read");
      const CentralityCutConfig& centrality = ConfigManager::GetInstance().GetCentralityCuts();
      const YAML::Node cent = YAML::LoadFile(centralityPath);
      UniqueMap(cent, "centrality");
      Require(cent["enabled"].IsScalar() && cent["mode"].IsScalar(),
              "centrality.enabled and centrality.mode must be explicit");
      Require(cent["enabled"].as<bool>() == static_cast<bool>(centrality.enabled) &&
              cent["mode"].as<std::string>() == centrality.mode,
              "ConfigManager centrality differs from supplied mainconf; load that mainconf first");
      std::ostringstream eventErrors;
      const bool eventLoaded = mEventSelection.Load(mainPath.c_str(), eventCuts, centrality.mode,
                                                    centrality.enabled, eventErrors,
                                                    policy == "lambda_imp5_compat");
      Require(eventLoaded, eventErrors.str());
      mCuts = &cuts;
      mInterface = new StPicoKFParticleInterface(cuts);
      mInitialized = true;
      std::cout << "FemtoLambdaProvider.mainconf (argument): " << mainPath << '\n'
                << "FemtoLambdaProvider.maker: " << makerPath << '\n'
                << "FemtoLambdaProvider.kf: " << kfPath << '\n'
                << "FemtoLambdaProvider.lambdaPairMass: " << mPairMass << '\n';
      mEventSelection.Dump(std::cout);
      cuts.Dump(std::cout);
      std::cout << StPicoKFParticleInterface::BackendDescription() << std::endl;
      return true;
    } catch (const std::exception& error) {
      mError = error.what();
      errors << "ERROR: Femto Lambda provider: " << mError << '\n';
      return false;
    }
  }
  virtual bool AcceptEvent(const StPicoEvent& event, int nTracks) const {
    return mInitialized && mEventSelection.Check(event, nTracks) == KfEventSelection::kAccepted;
  }
  virtual double ComputeVr(double x, double y) const {
    return mEventSelection.ComputeVr(x, y);
  }
  virtual bool Process(StPicoDst* picoDst, int eventIndex) {
    Clear();
    if (!mInitialized || !mInterface || !mCuts) {
      mError = "provider was not successfully initialized"; return false;
    }
    if (!picoDst || !picoDst->event() || eventIndex < 0) {
      mError = "missing PicoDst/event or invalid event index"; return false;
    }
    if (!mInterface->ProcessEvent(picoDst)) {
      mError = mInterface->LastError(); return false;
    }
    const KfParticleEventStats& s = mInterface->Stats();
#define COPY_STAT(field) mStats.field = s.field
    COPY_STAT(rawTracks); COPY_STAT(qualityTracks); COPY_STAT(covarianceTracks); COPY_STAT(pidTracks);
    COPY_STAT(pidHypotheses); COPY_STAT(protonHypotheses); COPY_STAT(pionHypotheses); COPY_STAT(unknownHypotheses);
    COPY_STAT(primaryHypotheses); COPY_STAT(finderParticles); COPY_STAT(lambdaParticles); COPY_STAT(validCandidates);
    COPY_STAT(invalidCandidates); COPY_STAT(invalidCovariance); COPY_STAT(invalidPid); COPY_STAT(invalidTrackIds);
    COPY_STAT(tofTracks);
#undef COPY_STAT
    const std::vector<KfLambdaCandidate>& raw = mInterface->Candidates();
    for (size_t i = 0; i < raw.size(); ++i) {
      const KfLambdaCandidate& c = raw[i];
      if (!PassKfLambdaCandidateCuts(c, *mCuts)) continue;
      // Validate the explicit array indices, never reinterpret persistent IDs.
      if (c.protonIndex < 0 || c.pionIndex < 0 ||
          static_cast<unsigned>(c.protonIndex) >= picoDst->numberOfTracks() ||
          static_cast<unsigned>(c.pionIndex) >= picoDst->numberOfTracks() ||
          !picoDst->track(c.protonIndex) || !picoDst->track(c.pionIndex) ||
          picoDst->track(c.protonIndex)->id() != c.protonId ||
          picoDst->track(c.pionIndex)->id() != c.pionId) {
        mError = "KF candidate daughter ID/index mapping is inconsistent";
        mCandidates.clear(); mStats.selectedCandidates = 0; return false;
      }
      FemtoLambdaCandidate output;
      if (!ConvertKfLambdaForFemto(c, eventIndex, mPairMass, output, mError)) {
        mCandidates.clear(); mStats.selectedCandidates = 0; return false;
      }
      mCandidates.push_back(output);
    }
    mStats.selectedCandidates = mCandidates.size();
    return true;
  }
  virtual void Clear() {
    mCandidates.clear();
    mStats = FemtoLambdaEventStats();
    mError.clear();
    if (mInterface) mInterface->Clear();
  }
  virtual const std::vector<FemtoLambdaCandidate>& Candidates() const { return mCandidates; }
  virtual const FemtoLambdaEventStats& Stats() const { return mStats; }
  virtual const std::string& LastError() const { return mError; }
private:
  KfProvider(const KfProvider&);
  KfProvider& operator=(const KfProvider&);
  StPicoKFParticleInterface* mInterface;
  const KfParticleCutConfig* mCuts;
  KfEventSelection mEventSelection;
  double mPairMass;
  bool mInitialized;
  std::vector<FemtoLambdaCandidate> mCandidates;
  FemtoLambdaEventStats mStats;
  std::string mError;
};
}

bool ConvertKfLambdaForFemto(const KfLambdaCandidate& source, int eventIndex,
                            double pairMass, FemtoLambdaCandidate& output,
                            std::string& error) {
  output = FemtoLambdaCandidate();
  error.clear();
  if (eventIndex < 0 || source.protonIndex < 0 || source.pionIndex < 0 ||
      source.protonIndex == source.pionIndex || source.pdg != 3122 ||
      !std::isfinite(pairMass) || pairMass <= 0. ||
      !std::isfinite(source.mass) || source.mass <= 0. ||
      !std::isfinite(source.px) || !std::isfinite(source.py) || !std::isfinite(source.pz)) {
    error = "invalid Lambda value/index/pair mass during Femto conversion"; return false;
  }
#define COPY_DIAGNOSTIC(field) output.field = source.field
  COPY_DIAGNOSTIC(x); COPY_DIAGNOSTIC(y); COPY_DIAGNOSTIC(z);
  COPY_DIAGNOSTIC(mass); COPY_DIAGNOSTIC(massError); COPY_DIAGNOSTIC(chi2Ndf); COPY_DIAGNOSTIC(topoChi2Ndf);
  COPY_DIAGNOSTIC(daughterDistance); COPY_DIAGNOSTIC(distanceToPv);
  COPY_DIAGNOSTIC(decayLength); COPY_DIAGNOSTIC(decayLengthError); COPY_DIAGNOSTIC(decayLengthSignificance);
  COPY_DIAGNOSTIC(vertexLineLength); COPY_DIAGNOSTIC(vertexLineLengthError); COPY_DIAGNOSTIC(vertexLineLengthSignificance);
  COPY_DIAGNOSTIC(cosPointing); COPY_DIAGNOSTIC(protonPidPull); COPY_DIAGNOSTIC(pionPidPull);
  COPY_DIAGNOSTIC(protonTofPull); COPY_DIAGNOSTIC(pionTofPull); COPY_DIAGNOSTIC(protonTofM2); COPY_DIAGNOSTIC(pionTofM2);
  COPY_DIAGNOSTIC(protonDcaToPv); COPY_DIAGNOSTIC(pionDcaToPv);
  COPY_DIAGNOSTIC(protonHelixPathLength); COPY_DIAGNOSTIC(pionHelixPathLength);
  COPY_DIAGNOSTIC(imp5PathValid); COPY_DIAGNOSTIC(protonHasTof); COPY_DIAGNOSTIC(pionHasTof);
  COPY_DIAGNOSTIC(protonId); COPY_DIAGNOSTIC(pionId); COPY_DIAGNOSTIC(protonIndex); COPY_DIAGNOSTIC(pionIndex); COPY_DIAGNOSTIC(pdg);
#undef COPY_DIAGNOSTIC
  TLorentzVector p4;
  p4.SetXYZM(source.px, source.py, source.pz, pairMass);
  output.candidate.SetP4(p4);
  output.candidate.eventIndex = eventIndex;
  output.candidate.source = kFemtoCandResonance;
  output.candidate.speciesKey = "lambda";
  output.candidate.charge = 0;
  output.candidate.reso.invMass = source.mass;
  output.candidate.reso.dcaDaughters = source.daughterDistance;
  output.candidate.reso.dau1EventIndex = eventIndex;
  output.candidate.reso.dau2EventIndex = eventIndex;
  output.candidate.reso.dau1Index = source.protonIndex;
  output.candidate.reso.dau2Index = source.pionIndex;
  return true;
}

FemtoLambdaProvider* createFemtoLambdaKfProvider() { return new KfProvider(); }
