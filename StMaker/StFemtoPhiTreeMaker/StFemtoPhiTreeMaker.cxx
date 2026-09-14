#include "StFemtoPhiTreeMaker.h"
#include "ConfigManager.h"
#include "YamlParser.h"
#include "kinematics.h"
#include "cuts/EventCutConfig.h"
#include "cuts/TrackCutConfig.h"
#include "cuts/PIDCutConfig.h"
#include "cuts/PhiCutConfig.h"
#include "cuts/MixingConfig.h"
#include "cuts/CentralityCutConfig.h"
#include "cuts/FemtoConfig.h"
#include "CentralityHelper.h"
#include "StNuclearIdHelper.h"
#include "StPhiKKReconstruction.h"
#include "FemtoFlagConfigSnapshot.h"
#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"

#include "TChain.h"
#include "TFile.h"
#include "TH1D.h"
#include "TNamed.h"
#include "TObjString.h"
#include "TROOT.h"
#include "TString.h"
#include "TSystem.h"
#include "TTree.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "RVersion.h"

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace {
const Double_t kProtonMass = 0.938272;

TString EnvOrEmpty(const char* key) {
  const char* v = gSystem->Getenv(key);
  return v ? TString(v) : TString("");
}

std::string JoinConfigPath(const std::string& mainconf, const std::string& relative) {
  std::string pathStr = mainconf;
  std::string base;
  size_t pos = pathStr.find("/config/");
  if (pos != std::string::npos) {
    base = pathStr.substr(0, pos + 1);
  }
  if (!base.empty() && base[base.size() - 1] != '/') base += "/";
  return base + "config/" + relative;
}
}  // namespace

StFemtoPhiTreeMaker::StFemtoPhiTreeMaker(const char* name, StPicoDstMaker* picoMaker, const char* outName)
    : StMaker(name),
      mPicoDstMaker(picoMaker),
      mPicoDst(0),
      mOutName(outName),
      mOutFile(0),
      mEventTree(0),
      mTrackTree(0),
      mPairTree(0),
      m_centrality(0),
      mWriteEvent(kTRUE),
      mWriteTrack(kTRUE),
      mWritePhiPair(kFALSE),
      mStoreProtons(kTRUE),
      mStoreDeuterons(kTRUE),
      mStoreKaons(kTRUE),
      mCompressLevel(1),
      mAutoFlush(10000),
      mSchemaVersion(femto_phi_tree::kSchemaVersionV2),
      mSubjobId(0),
      mPidCorrectionState("none"),
      mSkipRemaining(0),
      mEnvMinNHitsRatio(0.50),
      mEnvMinNHitsFit(15),
      mEnvMinNHitsDedx(10),
      mEnvMaxDca(3.0),
      mEnvMinPt(0.15),
      mEnvMaxPt(10.0),
      mEnvMinEta(-2.0),
      mEnvMaxEta(0.2),
      mEnvMaxChi2(5.0),
      mEnvMaxAbsNSigmaKaon(5.0),
      mEnvKaonRequireTofOrLowP(kFALSE),
      mEnvKaonLowPMax(0.5),
      mEnvMaxAbsNSigmaDeuteron(5.0),
      mEnvMaxAbsNSigmaProton(4.0),
      mEnvMaxDcaKaon(3.0),
      mEnvDeuteronMaxDca(2.0),
      mEnvDeuteronMinPMom(0.15),
      mEnvDeuteronMaxPMom(4.0),
      mEnvProtonMaxDca(2.0),
      mEnvProtonMinPt(0.15),
      mEnvMinNHitsDedxNuclear(10),
      mEnvKaonRequireDaughterPidReach(kFALSE),
      mWriteKaonOrigin(kTRUE),
      mKaonOriginTree(0),
      mEnvKaonMass2Lo(0.05),
      mEnvKaonMass2Hi(0.50),
      mNInput(0),
      mNNoPico(0),
      mNBadRun(0),
      mNFailEventCuts(0),
      mNPileup(0),
      mNFailCent(0),
      mNFailMaxNTr(0),
      mNAccepted(0),
      mNKp(0),
      mNKm(0),
      mNDeuteron(0),
      mNProton(0),
      mNPhiPair(0) {}

StFemtoPhiTreeMaker::~StFemtoPhiTreeMaker() {
  if (m_centrality) {
    delete m_centrality;
    m_centrality = 0;
  }
  mEventTree = 0;
  mTrackTree = 0;
  mPairTree = 0;
  if (mOutFile) {
    mOutFile->Close();
    delete mOutFile;
    mOutFile = 0;
  }
}

StFemtoPhiTreeMaker* createStFemtoPhiTreeMaker(const char* name, StPicoDstMaker* picoMaker,
                                               const char* outName) {
  return new StFemtoPhiTreeMaker(name, picoMaker, outName);
}

extern "C" void* createStFemtoPhiTreeMakerC(const char* name, void* picoMaker, const char* outName) {
  return (void*)createStFemtoPhiTreeMaker(name, (StPicoDstMaker*)picoMaker, outName);
}

Bool_t StFemtoPhiTreeMaker::LoadTreeConfig() {
  mMainconfPath = EnvOrEmpty("STAR_ANA_MAINCONF").Data();
  if (mMainconfPath.empty()) return kTRUE;

  std::map<std::string, std::string> mainValues;
  if (!YamlParser::ParseFile(mMainconfPath.c_str(), mainValues)) {
    std::cerr << "[StFemtoPhiTreeMaker] failed to parse mainconf for tree keys: " << mMainconfPath << std::endl;
    return kFALSE;
  }
  if (mainValues.find("maker") == mainValues.end()) {
    std::cerr << "[StFemtoPhiTreeMaker] mainconf has no maker: key" << std::endl;
    return kFALSE;
  }
  mMakerYamlPath = JoinConfigPath(mMainconfPath, mainValues["maker"]);

  std::map<std::string, std::string> values;
  if (!YamlParser::ParseFile(mMakerYamlPath.c_str(), values)) {
    std::cerr << "[StFemtoPhiTreeMaker] failed to parse maker YAML: " << mMakerYamlPath << std::endl;
    return kFALSE;
  }

  if (values.find("treeWriteEvent") != values.end()) mWriteEvent = YamlParser::ToBool(values["treeWriteEvent"], mWriteEvent);
  if (values.find("treeWriteTrack") != values.end()) mWriteTrack = YamlParser::ToBool(values["treeWriteTrack"], mWriteTrack);
  if (values.find("treeWritePhiPair") != values.end()) mWritePhiPair = YamlParser::ToBool(values["treeWritePhiPair"], mWritePhiPair);
  if (values.find("treeStoreProtons") != values.end()) mStoreProtons = YamlParser::ToBool(values["treeStoreProtons"], mStoreProtons);
  if (values.find("treeStoreDeuterons") != values.end()) mStoreDeuterons = YamlParser::ToBool(values["treeStoreDeuterons"], mStoreDeuterons);
  if (values.find("treeStoreKaons") != values.end()) mStoreKaons = YamlParser::ToBool(values["treeStoreKaons"], mStoreKaons);
  if (values.find("treeCompressionLevel") != values.end()) mCompressLevel = YamlParser::ToInt(values["treeCompressionLevel"], mCompressLevel);
  if (values.find("treeAutoFlush") != values.end()) mAutoFlush = YamlParser::ToInt(values["treeAutoFlush"], mAutoFlush);
  if (values.find("treeSchemaVersion") != values.end())
    mSchemaVersion = (UInt_t)YamlParser::ToInt(values["treeSchemaVersion"], (Int_t)mSchemaVersion);
  if (values.find("pidCorrectionState") != values.end())
    mPidCorrectionState = values["pidCorrectionState"];
  // Configured trigger ids, comma separated. Empty means "record only, select nothing".
  if (values.find("triggerIds") != values.end()) {
    mTriggerIds.clear();
    std::string raw = values["triggerIds"];
    for (size_t i = 0; i < raw.size(); ++i)
      if (raw[i] == '[' || raw[i] == ']' || raw[i] == ',') raw[i] = ' ';
    std::istringstream iss(raw);
    UInt_t t = 0;
    while (iss >> t) mTriggerIds.push_back(t);
  }

  if (values.find("envMinNHitsFit") != values.end()) mEnvMinNHitsFit = YamlParser::ToInt(values["envMinNHitsFit"], mEnvMinNHitsFit);
  if (values.find("envMinNHitsRatio") != values.end()) mEnvMinNHitsRatio = YamlParser::ToDouble(values["envMinNHitsRatio"], mEnvMinNHitsRatio);
  if (values.find("envMinNHitsDedx") != values.end()) mEnvMinNHitsDedx = YamlParser::ToInt(values["envMinNHitsDedx"], mEnvMinNHitsDedx);
  if (values.find("envMaxDca") != values.end()) mEnvMaxDca = YamlParser::ToDouble(values["envMaxDca"], mEnvMaxDca);
  if (values.find("envMinPt") != values.end()) mEnvMinPt = YamlParser::ToDouble(values["envMinPt"], mEnvMinPt);
  if (values.find("envMaxPt") != values.end()) mEnvMaxPt = YamlParser::ToDouble(values["envMaxPt"], mEnvMaxPt);
  if (values.find("envMinEta") != values.end()) mEnvMinEta = YamlParser::ToDouble(values["envMinEta"], mEnvMinEta);
  if (values.find("envMaxEta") != values.end()) mEnvMaxEta = YamlParser::ToDouble(values["envMaxEta"], mEnvMaxEta);
  if (values.find("envMaxChi2") != values.end()) mEnvMaxChi2 = YamlParser::ToDouble(values["envMaxChi2"], mEnvMaxChi2);
  if (values.find("envMaxAbsNSigmaKaon") != values.end()) mEnvMaxAbsNSigmaKaon = YamlParser::ToDouble(values["envMaxAbsNSigmaKaon"], mEnvMaxAbsNSigmaKaon);
  if (values.find("envKaonRequireTofOrLowP") != values.end())
    mEnvKaonRequireTofOrLowP = YamlParser::ToBool(values["envKaonRequireTofOrLowP"], mEnvKaonRequireTofOrLowP);
  if (values.find("envKaonLowPMax") != values.end())
    mEnvKaonLowPMax = YamlParser::ToDouble(values["envKaonLowPMax"], mEnvKaonLowPMax);
  if (values.find("envMaxAbsNSigmaDeuteron") != values.end()) mEnvMaxAbsNSigmaDeuteron = YamlParser::ToDouble(values["envMaxAbsNSigmaDeuteron"], mEnvMaxAbsNSigmaDeuteron);
  if (values.find("envMaxAbsNSigmaProton") != values.end()) mEnvMaxAbsNSigmaProton = YamlParser::ToDouble(values["envMaxAbsNSigmaProton"], mEnvMaxAbsNSigmaProton);
  if (values.find("envMaxDcaKaon") != values.end()) mEnvMaxDcaKaon = YamlParser::ToDouble(values["envMaxDcaKaon"], mEnvMaxDcaKaon);
  if (values.find("envDeuteronMaxDca") != values.end()) mEnvDeuteronMaxDca = YamlParser::ToDouble(values["envDeuteronMaxDca"], mEnvDeuteronMaxDca);
  if (values.find("envDeuteronMinPMom") != values.end()) mEnvDeuteronMinPMom = YamlParser::ToDouble(values["envDeuteronMinPMom"], mEnvDeuteronMinPMom);
  if (values.find("envDeuteronMaxPMom") != values.end()) mEnvDeuteronMaxPMom = YamlParser::ToDouble(values["envDeuteronMaxPMom"], mEnvDeuteronMaxPMom);
  if (values.find("envProtonMaxDca") != values.end()) mEnvProtonMaxDca = YamlParser::ToDouble(values["envProtonMaxDca"], mEnvProtonMaxDca);
  if (values.find("envProtonMinPt") != values.end()) mEnvProtonMinPt = YamlParser::ToDouble(values["envProtonMinPt"], mEnvProtonMinPt);
  if (values.find("envMinNHitsDedxNuclear") != values.end()) mEnvMinNHitsDedxNuclear = YamlParser::ToInt(values["envMinNHitsDedxNuclear"], mEnvMinNHitsDedxNuclear);
  if (values.find("treeWriteKaonOrigin") != values.end()) mWriteKaonOrigin = YamlParser::ToBool(values["treeWriteKaonOrigin"], mWriteKaonOrigin);
  if (values.find("envKaonRequireDaughterPidReach") != values.end()) mEnvKaonRequireDaughterPidReach = YamlParser::ToBool(values["envKaonRequireDaughterPidReach"], mEnvKaonRequireDaughterPidReach);
  if (values.find("envKaonMass2Lo") != values.end()) mEnvKaonMass2Lo = YamlParser::ToDouble(values["envKaonMass2Lo"], mEnvKaonMass2Lo);
  if (values.find("envKaonMass2Hi") != values.end()) mEnvKaonMass2Hi = YamlParser::ToDouble(values["envKaonMass2Hi"], mEnvKaonMass2Hi);

  TString skipEnv = EnvOrEmpty("STAR_ANA_NSKIP");
  if (skipEnv.Length() > 0) mSkipRemaining = std::atoll(skipEnv.Data());
  if (mSkipRemaining > 0) {
    std::cout << "[StFemtoPhiTreeMaker] skipping first " << mSkipRemaining << " events before fill" << std::endl;
  }


  std::cout << "[StFemtoPhiTreeMaker] tree config from " << mMakerYamlPath
            << " writeEvent=" << mWriteEvent << " writeTrack=" << mWriteTrack
            << " writePhiPair=" << mWritePhiPair << " storeP=" << mStoreProtons
            << " storeD=" << mStoreDeuterons << " compress=" << mCompressLevel
            << " schemaVersion=" << mSchemaVersion
            << " kaonRequireTofOrLowP=" << mEnvKaonRequireTofOrLowP
            << " kaonLowPMax=" << mEnvKaonLowPMax << " nTriggerIds=" << mTriggerIds.size()
            << " pidCorrectionState=" << mPidCorrectionState << std::endl;
  return kTRUE;
}

void StFemtoPhiTreeMaker::BookTrees() {
  if (mSchemaVersion >= femto_phi_tree::kSchemaVersionV2) {
    BookTreesV2();
  } else {
    BookTreesV1();
  }
}

void StFemtoPhiTreeMaker::BookTreesV2() {
  if (mSchemaVersion >= femto_phi_tree::kSchemaVersionV3) {
    BookEventTreeV3();
  } else if (mWriteEvent) {
    mEventTree = new TTree("FemtoEventTree", "accepted femto events (schema 2)");
    mEventTree->SetAutoFlush(mAutoFlush);
    mEventTree->Branch("schemaVersion", &mEvt2.schemaVersion, "schemaVersion/i");
    mEventTree->Branch("eventUID", &mEvt2.eventUID, "eventUID/l");
    mEventTree->Branch("runId", &mEvt2.runId, "runId/I");
    mEventTree->Branch("eventId", &mEvt2.eventId, "eventId/I");
    mEventTree->Branch("sourceFileHash", &mEvt2.sourceFileHash, "sourceFileHash/i");
    mEventTree->Branch("sourceFileIndex", &mEvt2.sourceFileIndex, "sourceFileIndex/s");
    mEventTree->Branch("sourceEntry", &mEvt2.sourceEntry, "sourceEntry/L");
    mEventTree->Branch("subjobId", &mEvt2.subjobId, "subjobId/i");
    mEventTree->Branch("triggerId", &mEvt2.triggerId, "triggerId/i");
    mEventTree->Branch("triggerBits", &mEvt2.triggerBits, "triggerBits/s");
    mEventTree->Branch("vx", &mEvt2.vx, "vx/F");
    mEventTree->Branch("vy", &mEvt2.vy, "vy/F");
    mEventTree->Branch("vz", &mEvt2.vz, "vz/F");
    mEventTree->Branch("vr", &mEvt2.vr, "vr/F");
    mEventTree->Branch("vzVpd", &mEvt2.vzVpd, "vzVpd/F");
    mEventTree->Branch("bField", &mEvt2.bField, "bField/F");
    mEventTree->Branch("refMult", &mEvt2.refMult, "refMult/I");
    mEventTree->Branch("rawMult", &mEvt2.rawMult, "rawMult/I");
    mEventTree->Branch("nBTOFMatch", &mEvt2.nBTOFMatch, "nBTOFMatch/I");
    mEventTree->Branch("nTracks", &mEvt2.nTracks, "nTracks/I");
    mEventTree->Branch("refMultCorr", &mEvt2.refMultCorr, "refMultCorr/F");
    mEventTree->Branch("centWeight", &mEvt2.centWeight, "centWeight/F");
    mEventTree->Branch("centralityPercent", &mEvt2.centralityPercent, "centralityPercent/F");
    mEventTree->Branch("cent9", &mEvt2.cent9, "cent9/I");
    mEventTree->Branch("cent16", &mEvt2.cent16, "cent16/I");
    mEventTree->Branch("qx", &mEvt2.qx, "qx/F");
    mEventTree->Branch("qy", &mEvt2.qy, "qy/F");
    mEventTree->Branch("psi2", &mEvt2.psi2, "psi2/F");
    mEventTree->Branch("mixVzBin", &mEvt2.mixVzBin, "mixVzBin/I");
    mEventTree->Branch("mixCentBin", &mEvt2.mixCentBin, "mixCentBin/I");
    mEventTree->Branch("mixEpBin", &mEvt2.mixEpBin, "mixEpBin/I");
    mEventTree->Branch("mixBin", &mEvt2.mixBin, "mixBin/I");
    mEventTree->Branch("eventFlags", &mEvt2.eventFlags, "eventFlags/i");
    mEventTree->Branch("nKp", &mEvt2.nKp, "nKp/I");
    mEventTree->Branch("nKm", &mEvt2.nKm, "nKm/I");
    mEventTree->Branch("nDeuteron", &mEvt2.nDeuteron, "nDeuteron/I");
    mEventTree->Branch("nProton", &mEvt2.nProton, "nProton/I");
  }
  if (mWriteTrack) {
    mTrackTree = new TTree("FemtoTrackTree", "envelope candidate tracks (packed)");
    mTrackTree->SetAutoFlush(mAutoFlush);
    mTrackTree->Branch("eventUID", &mTrk2.eventUID, "eventUID/l");
    mTrackTree->Branch("trackIndex", &mTrk2.trackIndex, "trackIndex/s");
    mTrackTree->Branch("speciesCode", &mTrk2.speciesCode, "speciesCode/b");
    mTrackTree->Branch("pT", &mTrk2.pT, "pT/s");
    mTrackTree->Branch("eta", &mTrk2.eta, "eta/S");
    mTrackTree->Branch("phi", &mTrk2.phi, "phi/S");
    mTrackTree->Branch("dEdx", &mTrk2.dEdx, "dEdx/s");
    mTrackTree->Branch("nSigmaKaon", &mTrk2.nSigmaKaon, "nSigmaKaon/S");
    mTrackTree->Branch("nSigmaPion", &mTrk2.nSigmaPion, "nSigmaPion/S");
    mTrackTree->Branch("nSigmaProton", &mTrk2.nSigmaProton, "nSigmaProton/S");
    mTrackTree->Branch("nSigmaDeuteron", &mTrk2.nSigmaDeuteron, "nSigmaDeuteron/S");
    mTrackTree->Branch("tofBeta", &mTrk2.tofBeta, "tofBeta/S");
    mTrackTree->Branch("dca", &mTrk2.dca, "dca/s");
    mTrackTree->Branch("nHitsFit", &mTrk2.nHitsFit, "nHitsFit/B");
    mTrackTree->Branch("nHitsMax", &mTrk2.nHitsMax, "nHitsMax/b");
    mTrackTree->Branch("nHitsDedx", &mTrk2.nHitsDedx, "nHitsDedx/b");
    mTrackTree->Branch("selFlags", &mTrk2.selFlags, "selFlags/s");
  }
  if (mWritePhiPair) {
    std::cerr << "[StFemtoPhiTreeMaker] treeWritePhiPair is not supported in schema 2; "
              << "phi pairs are rebuilt downstream (plan sec 2 item 3)" << std::endl;
  }
}

void StFemtoPhiTreeMaker::BookEventTreeV3() {
  if (!mWriteEvent) return;
  mEventTree = new TTree("FemtoEventTree", "accepted femto events (schema 3, packed)");
  mEventTree->SetAutoFlush(mAutoFlush);
  mEventTree->Branch("schemaVersion", &mEvt3.schemaVersion, "schemaVersion/i");
  mEventTree->Branch("eventUID", &mEvt3.eventUID, "eventUID/l");
  mEventTree->Branch("sourceFileHash", &mEvt3.sourceFileHash, "sourceFileHash/i");
  mEventTree->Branch("sourceFileIndex", &mEvt3.sourceFileIndex, "sourceFileIndex/s");
  mEventTree->Branch("sourceEntry", &mEvt3.sourceEntry, "sourceEntry/I");
  mEventTree->Branch("subjobId", &mEvt3.subjobId, "subjobId/i");
  mEventTree->Branch("triggerId", &mEvt3.triggerId, "triggerId/i");
  mEventTree->Branch("triggerBits", &mEvt3.triggerBits, "triggerBits/s");
  mEventTree->Branch("vx", &mEvt3.vx, "vx/S");
  mEventTree->Branch("vy", &mEvt3.vy, "vy/S");
  mEventTree->Branch("vz", &mEvt3.vz, "vz/S");
  mEventTree->Branch("vzVpd", &mEvt3.vzVpd, "vzVpd/S");
  mEventTree->Branch("vr", &mEvt3.vr, "vr/s");
  mEventTree->Branch("bField", &mEvt3.bField, "bField/F");
  mEventTree->Branch("refMult", &mEvt3.refMult, "refMult/s");
  mEventTree->Branch("rawMult", &mEvt3.rawMult, "rawMult/s");
  mEventTree->Branch("nBTOFMatch", &mEvt3.nBTOFMatch, "nBTOFMatch/s");
  mEventTree->Branch("nTracks", &mEvt3.nTracks, "nTracks/s");
  mEventTree->Branch("refMultCorr", &mEvt3.refMultCorr, "refMultCorr/s");
  mEventTree->Branch("centWeight", &mEvt3.centWeight, "centWeight/s");
  mEventTree->Branch("centralityPercent", &mEvt3.centralityPercent, "centralityPercent/s");
  mEventTree->Branch("cent9", &mEvt3.cent9, "cent9/B");
  mEventTree->Branch("cent16", &mEvt3.cent16, "cent16/B");
  mEventTree->Branch("qx", &mEvt3.qx, "qx/S");
  mEventTree->Branch("qy", &mEvt3.qy, "qy/S");
  mEventTree->Branch("psi2", &mEvt3.psi2, "psi2/S");
  mEventTree->Branch("mixBin", &mEvt3.mixBin, "mixBin/s");
  mEventTree->Branch("eventFlags", &mEvt3.eventFlags, "eventFlags/b");
  mEventTree->Branch("nKp", &mEvt3.nKp, "nKp/s");
  mEventTree->Branch("nKm", &mEvt3.nKm, "nKm/s");
  mEventTree->Branch("nDeuteron", &mEvt3.nDeuteron, "nDeuteron/s");
  mEventTree->Branch("nProton", &mEvt3.nProton, "nProton/s");
  BookKaonOriginTree();
}

// Companion tree: one row per stored kaon, carrying the helix origin the KK decay DCA needs.
// Kept out of the track row because only kaons need it and they are 0.66% of stored rows --
// on every row the same three components would cost +32% of the tree instead of +0.2%.
// Record every config value the stored selFlags depend on, so a downstream configured
// differently fails loudly instead of silently returning the producer's selection. One row per
// key per job: after hadd the rows repeat, and a key that then carries two different values means
// trees from two different configurations were merged -- which the reader also catches.
void StFemtoPhiTreeMaker::WriteFlagConfigSnapshot() {
  if (!mOutFile) return;
  mOutFile->cd();
  TTree* t = new TTree("FemtoFlagConfig", "config values the selFlags bits depend on");
  TString key, value;
  TString* pk = &key;
  TString* pv = &value;
  t->Branch("key", &pk);
  t->Branch("value", &pv);
  std::vector<femto_flag_config::Entry> e = femto_flag_config::Collect();
  for (size_t i = 0; i < e.size(); ++i) {
    key = e[i].first;
    value = e[i].second;
    t->Fill();
  }
  t->Write();
  LOG_INFO << "[StFemtoPhiTreeMaker] flag-config snapshot: " << e.size() << " keys" << endm;
}

void StFemtoPhiTreeMaker::BookKaonOriginTree() {
  if (!mWriteKaonOrigin || !mStoreKaons) return;
  mKaonOriginTree = new TTree("FemtoKaonOriginTree", "phi-daughter kaon helix origin");
  mKaonOriginTree->SetAutoFlush(mAutoFlush);
  mKaonOriginTree->Branch("eventUID", &mKaonOrigin.eventUID, "eventUID/l");
  mKaonOriginTree->Branch("trackIndex", &mKaonOrigin.trackIndex, "trackIndex/s");
  mKaonOriginTree->Branch("originDx", &mKaonOrigin.originDx, "originDx/S");
  mKaonOriginTree->Branch("originDy", &mKaonOrigin.originDy, "originDy/S");
  mKaonOriginTree->Branch("originDz", &mKaonOrigin.originDz, "originDz/S");
}

// Origin is stored relative to the primary vertex: two helices translated by the same vector have
// the same distance of closest approach, so a reader works in the vertex frame and the 0.01 cm
// quantisation of the stored vz never reaches the KK DCA.
void StFemtoPhiTreeMaker::FillKaonOriginRow(const TrackState& trk, ULong64_t eventUID,
                                            const TVector3& pVtx) {
  if (!mKaonOriginTree) return;
  using namespace femto_phi_tree;
  mKaonOrigin.Reset();
  mKaonOrigin.eventUID = eventUID;
  mKaonOrigin.trackIndex = (UShort_t)trk.trackIndex;
  mKaonOrigin.originDx = PackI16(trk.originX - pVtx.X(), escale::kOrigin, mPackStats.kaonOrigin);
  mKaonOrigin.originDy = PackI16(trk.originY - pVtx.Y(), escale::kOrigin, mPackStats.kaonOrigin);
  mKaonOrigin.originDz = PackI16(trk.originZ - pVtx.Z(), escale::kOrigin, mPackStats.kaonOrigin);
  mKaonOriginTree->Fill();
}

// Pack the schema 2 event row (already filled from the event) into schema 3.
void StFemtoPhiTreeMaker::FillEventRowV3() {
  using namespace femto_phi_tree;
  mEvt3.Reset();
  mEvt3.eventUID = mEvt2.eventUID;
  mEvt3.sourceFileHash = mEvt2.sourceFileHash;
  mEvt3.sourceFileIndex = mEvt2.sourceFileIndex;
  mEvt3.sourceEntry = (Int_t)mEvt2.sourceEntry;
  mEvt3.subjobId = mEvt2.subjobId;
  mEvt3.triggerId = mEvt2.triggerId;
  mEvt3.triggerBits = mEvt2.triggerBits;
  mEvt3.vx = PackI16(mEvt2.vx, escale::kVxy, mPackStats.evVertex);
  mEvt3.vy = PackI16(mEvt2.vy, escale::kVxy, mPackStats.evVertex);
  mEvt3.vz = PackI16(mEvt2.vz, escale::kVz, mPackStats.evVertex);
  mEvt3.vzVpd = PackI16(mEvt2.vzVpd, escale::kVzVpd, mPackStats.evVertex);
  mEvt3.vr = PackU16(mEvt2.vr, escale::kVr, mPackStats.evVertex);
  mEvt3.bField = mEvt2.bField;
  mEvt3.refMult = PackU16(mEvt2.refMult, 1.0, mPackStats.evMult);
  mEvt3.rawMult = PackU16(mEvt2.rawMult, 1.0, mPackStats.evMult);
  mEvt3.nBTOFMatch = PackU16(mEvt2.nBTOFMatch, 1.0, mPackStats.evMult);
  mEvt3.nTracks = PackU16(mEvt2.nTracks, 1.0, mPackStats.evMult);
  mEvt3.refMultCorr = PackU16(mEvt2.refMultCorr, escale::kRefMultCorr, mPackStats.evCent);
  mEvt3.centWeight = PackU16(mEvt2.centWeight, escale::kCentWeight, mPackStats.evCent);
  mEvt3.centralityPercent =
      PackU16(mEvt2.centralityPercent, escale::kCentPercent, mPackStats.evCent);
  mEvt3.cent9 = (Char_t)mEvt2.cent9;
  mEvt3.cent16 = (Char_t)mEvt2.cent16;
  mEvt3.qx = PackI16(mEvt2.qx, escale::kQ, mPackStats.evQ);
  mEvt3.qy = PackI16(mEvt2.qy, escale::kQ, mPackStats.evQ);
  mEvt3.psi2 = PackI16(mEvt2.psi2, escale::kPsi2, mPackStats.evQ);
  mEvt3.mixBin = (mEvt2.mixBin >= 0 && mEvt2.mixBin <= 65535) ? (UShort_t)mEvt2.mixBin : 0;
  mEvt3.eventFlags = (UChar_t)(mEvt2.eventFlags & 0xFFu);
  mEvt3.nKp = PackU16(mEvt2.nKp, 1.0, mPackStats.evMult);
  mEvt3.nKm = PackU16(mEvt2.nKm, 1.0, mPackStats.evMult);
  mEvt3.nDeuteron = PackU16(mEvt2.nDeuteron, 1.0, mPackStats.evMult);
  mEvt3.nProton = PackU16(mEvt2.nProton, 1.0, mPackStats.evMult);
}

void StFemtoPhiTreeMaker::BookTreesV1() {
  if (mWriteEvent) {
    mEventTree = new TTree("FemtoEventTree", "accepted femto events");
    mEventTree->SetAutoFlush(mAutoFlush);
    mEventTree->Branch("schemaVersion", &mEvt.schemaVersion, "schemaVersion/i");
    mEventTree->Branch("eventUID", &mEvt.eventUID, "eventUID/l");
    mEventTree->Branch("runId", &mEvt.runId, "runId/I");
    mEventTree->Branch("eventId", &mEvt.eventId, "eventId/I");
    mEventTree->Branch("sourceFileHash", &mEvt.sourceFileHash, "sourceFileHash/i");
    mEventTree->Branch("sourceEntry", &mEvt.sourceEntry, "sourceEntry/L");
    mEventTree->Branch("subjobId", &mEvt.subjobId, "subjobId/i");
    mEventTree->Branch("vx", &mEvt.vx, "vx/F");
    mEventTree->Branch("vy", &mEvt.vy, "vy/F");
    mEventTree->Branch("vz", &mEvt.vz, "vz/F");
    mEventTree->Branch("vr", &mEvt.vr, "vr/F");
    mEventTree->Branch("vzVpd", &mEvt.vzVpd, "vzVpd/F");
    mEventTree->Branch("bField", &mEvt.bField, "bField/F");
    mEventTree->Branch("refMult", &mEvt.refMult, "refMult/I");
    mEventTree->Branch("rawMult", &mEvt.rawMult, "rawMult/I");
    mEventTree->Branch("nBTOFMatch", &mEvt.nBTOFMatch, "nBTOFMatch/I");
    mEventTree->Branch("nTracks", &mEvt.nTracks, "nTracks/I");
    mEventTree->Branch("refMultCorr", &mEvt.refMultCorr, "refMultCorr/F");
    mEventTree->Branch("centWeight", &mEvt.centWeight, "centWeight/F");
    mEventTree->Branch("centralityPercent", &mEvt.centralityPercent, "centralityPercent/F");
    mEventTree->Branch("cent9", &mEvt.cent9, "cent9/I");
    mEventTree->Branch("cent16", &mEvt.cent16, "cent16/I");
    mEventTree->Branch("qx", &mEvt.qx, "qx/F");
    mEventTree->Branch("qy", &mEvt.qy, "qy/F");
    mEventTree->Branch("psi2", &mEvt.psi2, "psi2/F");
    mEventTree->Branch("mixVzBin", &mEvt.mixVzBin, "mixVzBin/I");
    mEventTree->Branch("mixCentBin", &mEvt.mixCentBin, "mixCentBin/I");
    mEventTree->Branch("mixEpBin", &mEvt.mixEpBin, "mixEpBin/I");
    mEventTree->Branch("mixBin", &mEvt.mixBin, "mixBin/I");
    mEventTree->Branch("eventFlags", &mEvt.eventFlags, "eventFlags/i");
    mEventTree->Branch("nKp", &mEvt.nKp, "nKp/I");
    mEventTree->Branch("nKm", &mEvt.nKm, "nKm/I");
    mEventTree->Branch("nDeuteron", &mEvt.nDeuteron, "nDeuteron/I");
    mEventTree->Branch("nProton", &mEvt.nProton, "nProton/I");
  }
  if (mWriteTrack) {
    mTrackTree = new TTree("FemtoTrackTree", "envelope candidate tracks");
    mTrackTree->SetAutoFlush(mAutoFlush);
    mTrackTree->Branch("eventUID", &mTrk.eventUID, "eventUID/l");
    mTrackTree->Branch("trackIndex", &mTrk.trackIndex, "trackIndex/I");
    mTrackTree->Branch("speciesCode", &mTrk.speciesCode, "speciesCode/b");
    mTrackTree->Branch("charge", &mTrk.charge, "charge/S");
    mTrackTree->Branch("px", &mTrk.px, "px/F");
    mTrackTree->Branch("py", &mTrk.py, "py/F");
    mTrackTree->Branch("pz", &mTrk.pz, "pz/F");
    mTrackTree->Branch("dEdx", &mTrk.dEdx, "dEdx/F");
    mTrackTree->Branch("nSigmaKaon", &mTrk.nSigmaKaon, "nSigmaKaon/F");
    mTrackTree->Branch("nSigmaDeuteron", &mTrk.nSigmaDeuteron, "nSigmaDeuteron/F");
    mTrackTree->Branch("nSigmaProton", &mTrk.nSigmaProton, "nSigmaProton/F");
    mTrackTree->Branch("tofMatch", &mTrk.tofMatch, "tofMatch/B");
    mTrackTree->Branch("tofBeta", &mTrk.tofBeta, "tofBeta/F");
    mTrackTree->Branch("mass2", &mTrk.mass2, "mass2/F");
    mTrackTree->Branch("deltaOneOverBeta", &mTrk.deltaOneOverBeta, "deltaOneOverBeta/F");
    mTrackTree->Branch("dca", &mTrk.dca, "dca/F");
    mTrackTree->Branch("nHitsFit", &mTrk.nHitsFit, "nHitsFit/S");
    mTrackTree->Branch("nHitsMax", &mTrk.nHitsMax, "nHitsMax/S");
    mTrackTree->Branch("nHitsDedx", &mTrk.nHitsDedx, "nHitsDedx/S");
    mTrackTree->Branch("chi2", &mTrk.chi2, "chi2/F");
    mTrackTree->Branch("originX", &mTrk.originX, "originX/F");
    mTrackTree->Branch("originY", &mTrk.originY, "originY/F");
    mTrackTree->Branch("originZ", &mTrk.originZ, "originZ/F");
    mTrackTree->Branch("bField", &mTrk.bField, "bField/F");
    mTrackTree->Branch("selFlags", &mTrk.selFlags, "selFlags/i");
  }
  if (mWritePhiPair) {
    mPairTree = new TTree("FemtoPhiPairTree", "optional real KK cache");
    mPairTree->SetAutoFlush(mAutoFlush);
    mPairTree->Branch("eventUID", &mPair.eventUID, "eventUID/l");
    mPairTree->Branch("dauPlusIndex", &mPair.dauPlusIndex, "dauPlusIndex/I");
    mPairTree->Branch("dauMinusIndex", &mPair.dauMinusIndex, "dauMinusIndex/I");
    mPairTree->Branch("px", &mPair.px, "px/F");
    mPairTree->Branch("py", &mPair.py, "py/F");
    mPairTree->Branch("pz", &mPair.pz, "pz/F");
    mPairTree->Branch("mKK", &mPair.mKK, "mKK/F");
    mPairTree->Branch("dcaDaughters", &mPair.dcaDaughters, "dcaDaughters/F");
    mPairTree->Branch("openingAngle", &mPair.openingAngle, "openingAngle/F");
    mPairTree->Branch("rapidity", &mPair.rapidity, "rapidity/F");
    mPairTree->Branch("pairFlags", &mPair.pairFlags, "pairFlags/i");
  }
}

Int_t StFemtoPhiTreeMaker::Init() {
  if (!LoadTreeConfig()) return kStErr;

  m_centrality = new CentralityHelper();
  if (!m_centrality->Init(ConfigManager::GetInstance().GetCentralityCuts())) {
    std::cerr << "[StFemtoPhiTreeMaker] CentralityHelper init failed" << std::endl;
  }
  if (!ConfigManager::GetInstance().GetPhiCuts().FinalizeRapidityFrame(
          ConfigManager::GetInstance().GetCentralityCuts())) {
    std::cerr << "[StFemtoPhiTreeMaker] FinalizeRapidityFrame failed" << std::endl;
    return kStErr;
  }

  if (mOutName == "") {
    std::cerr << "[StFemtoPhiTreeMaker] empty output name" << std::endl;
    return kStErr;
  }
  mOutFile = new TFile(mOutName.Data(), "RECREATE");
  if (!mOutFile || mOutFile->IsZombie()) {
    std::cerr << "[StFemtoPhiTreeMaker] cannot create " << mOutName << std::endl;
    return kStErr;
  }
  mOutFile->SetCompressionLevel(mCompressLevel);
#if ROOT_VERSION_CODE >= ROOT_VERSION(6, 0, 0)
  mOutFile->SetCompressionAlgorithm(ROOT::kZLIB);
#endif
  // subjobId identifies the producing job inside a merged tree. STAR jobids are
  // 32-char hex in batch and short integers locally, so hash to a stable UInt.
  TString jobId = EnvOrEmpty("STAR_ANA_JOBID");
  mSubjobId = jobId.Length() ? femto_phi_tree::HashStringFnv(jobId.Data()) : 0u;
  std::cout << "[StFemtoPhiTreeMaker] jobid='" << jobId << "' subjobId=" << mSubjobId << std::endl;
  BookTrees();
  return kStOK;
}

void StFemtoPhiTreeMaker::Clear(Option_t* opt) { (void)opt; }

Bool_t StFemtoPhiTreeMaker::PassEventCuts(Float_t vz, Float_t vr, Int_t refMult, Float_t vzVpd) const {
  EventCutConfig& ev = ConfigManager::GetInstance().GetEventCuts();
  if (vz < ev.minVz || vz > ev.maxVz) return kFALSE;
  if (vr > ev.maxVr) return kFALSE;
  if (refMult < ev.minRefMult) return kFALSE;
  if (refMult > ev.maxRefMult) return kFALSE;
  if (TMath::Abs(vz - vzVpd) > ev.maxVzDiff && TMath::Abs(vzVpd) < ev.maxAbsVzVpd) return kFALSE;
  return kTRUE;
}

Bool_t StFemtoPhiTreeMaker::PassEnvTrackCuts(StPicoTrack* trk, TVector3& pVtx) const {
  if (trk->nHitsFit() < mEnvMinNHitsFit) return kFALSE;
  if (trk->nHitsMax() <= 0) return kFALSE;
  if ((Float_t)trk->nHitsFit() / (Float_t)trk->nHitsMax() < mEnvMinNHitsRatio) return kFALSE;
  if (trk->nHitsDedx() < mEnvMinNHitsDedx) return kFALSE;
  TVector3 pMom = trk->pMom();
  if (pMom.Mag() < 1e-4) return kFALSE;
  Float_t pt = pMom.Perp();
  Float_t eta = pMom.PseudoRapidity();
  if (pt < mEnvMinPt || pt > mEnvMaxPt) return kFALSE;
  if (eta < mEnvMinEta || eta > mEnvMaxEta) return kFALSE;
  if (trk->gDCA(pVtx).Mag() > mEnvMaxDca) return kFALSE;
  TrackCutConfig& tr = ConfigManager::GetInstance().GetTrackCuts();
  if (tr.requirePrimaryTrack && !trk->isPrimary()) return kFALSE;
  if (trk->chi2() > mEnvMaxChi2) return kFALSE;
  return kTRUE;
}

Bool_t StFemtoPhiTreeMaker::PassNominalTrackCuts(StPicoTrack* trk, TVector3& pVtx) const {
  TrackCutConfig& tr = ConfigManager::GetInstance().GetTrackCuts();
  if (trk->nHitsFit() < tr.minNHitsFit) return kFALSE;
  if ((Float_t)trk->nHitsFit() / (Float_t)trk->nHitsMax() < tr.minNHitsRatio) return kFALSE;
  if (trk->nHitsDedx() < tr.minNHitsDedx) return kFALSE;
  TVector3 pMom = trk->pMom();
  if (pMom.Mag() < 1e-4) return kFALSE;
  Float_t pt = pMom.Perp();
  Float_t eta = pMom.PseudoRapidity();
  if (pt < tr.minPt || pt > tr.maxPt) return kFALSE;
  if (eta < tr.minEta || eta > tr.maxEta) return kFALSE;
  if (trk->gDCA(pVtx).Mag() > tr.maxDCA) return kFALSE;
  if (tr.requirePrimaryTrack && !trk->isPrimary()) return kFALSE;
  if (trk->chi2() > tr.maxChi2) return kFALSE;
  return kTRUE;
}

Bool_t StFemtoPhiTreeMaker::PassNominalKaonCuts(StPicoTrack* trk, TVector3& pVtx) const {
  if (!PassNominalTrackCuts(trk, pVtx)) return kFALSE;
  PhiCutConfig& phi = ConfigManager::GetInstance().GetPhiCuts();
  if (trk->gDCA(pVtx).Mag() > phi.maxDCAKaon) return kFALSE;
  if (TMath::Abs(trk->nSigmaKaon()) > phi.nSigmaKaon) return kFALSE;
  return kTRUE;
}

void StFemtoPhiTreeMaker::BuildTrackState(TrackState& track, StPicoTrack* pico, StPicoEvent* event,
                                         TVector3& pVtx, Int_t index) {
  TVector3 pmom = pico->pMom();
  TVector3 org = pico->origin();
  track.originX = org.X();
  track.originY = org.Y();
  track.originZ = org.Z();
  track.momentumX = pmom.X();
  track.momentumY = pmom.Y();
  track.momentumZ = pmom.Z();
  track.BField = event->bField();
  track.pT = pmom.Perp();
  track.eta = pmom.PseudoRapidity();
  track.phi = pmom.Phi();
  track.charge = pico->charge();
  track.nHitsFit = pico->nHitsFit();
  track.nHitsMax = pico->nHitsMax();
  track.nHitsDedx = pico->nHitsDedx();
  track.DCA = pico->gDCA(pVtx).Mag();
  track.chi2 = pico->chi2();
  track.nSigmaKaon = pico->nSigmaKaon();
  track.nSigmaPion = pico->nSigmaPion();
  track.nSigmaProton = pico->nSigmaProton();
  track.nSigmaDeuteron = 0.0f;
  track.dEdx = pico->dEdx();
  track.trackIndex = index;
  track.tofMatch = 0;
  track.mass2 = -999.0f;
  track.deltaOneOverBeta = 999.0f;
  track.tofBeta = -1.0f;
  track.selFlags = 0;
}

PhiKkTrackState StFemtoPhiTreeMaker::ToPhiKkTrack(const TrackState& trk) {
  PhiKkTrackState s;
  s.pT = trk.pT;
  s.eta = trk.eta;
  s.phi = trk.phi;
  s.charge = trk.charge;
  s.originX = trk.originX;
  s.originY = trk.originY;
  s.originZ = trk.originZ;
  s.momentumX = trk.momentumX;
  s.momentumY = trk.momentumY;
  s.momentumZ = trk.momentumZ;
  s.BField = trk.BField;
  s.tofMatch = trk.tofMatch ? kTRUE : kFALSE;
  s.mass2 = trk.mass2;
  s.deltaOneOverBeta = trk.deltaOneOverBeta;
  return s;
}

void StFemtoPhiTreeMaker::FillTofInfo(TrackState& track, StPicoTrack* trk, const TVector3& pMom, Int_t btofIndex) {
  (void)trk;
  track.tofBeta = -1.0f;
  StPicoBTofPidTraits* tof = 0;
  if (btofIndex >= 0) tof = mPicoDst->btofPidTraits(btofIndex);
  Bool_t match = kFALSE;
  StPhiKKReconstruction::FillTofInfo(track.mass2, track.deltaOneOverBeta, match, tof, pMom);
  track.tofMatch = match ? 1 : 0;
  if (tof && match) {
    const Double_t beta = tof->btofBeta();
    if (beta > 1e-4) track.tofBeta = (Float_t)beta;
  }
}

Bool_t StFemtoPhiTreeMaker::PassTofKaonPid(const TrackState& trk) const {
  return StPhiKKReconstruction::PassTofKaonPid(ToPhiKkTrack(trk));
}

Bool_t StFemtoPhiTreeMaker::PassPhiDaughterTofPid(const TrackState& trk) const {
  const Float_t pMag = (Float_t)TrackMomentum(trk).Mag();
  return StPhiKKReconstruction::PassPhiDaughterTofPid(pMag, trk.tofMatch ? kTRUE : kFALSE, trk.mass2,
                                                      trk.deltaOneOverBeta, trk.charge);
}

Bool_t StFemtoPhiTreeMaker::PassTofProtonPid(const TrackState& trk) const {
  const FemtoConfig& fc = ConfigManager::GetInstance().GetFemtoConfig();
  const Double_t pmom = TrackMomentum(trk).Mag();
  if (pmom < fc.protonTofMomentumThreshold) return kTRUE;
  if (!trk.tofMatch) return kFALSE;
  return (trk.mass2 >= fc.protonMinMass2 && trk.mass2 <= fc.protonMaxMass2);
}

TVector3 StFemtoPhiTreeMaker::TrackMomentum(const TrackState& trk) const {
  return TVector3(trk.momentumX, trk.momentumY, trk.momentumZ);
}

Double_t StFemtoPhiTreeMaker::ApplyRapidityFrame(Double_t yLab) const {
  return ConfigManager::GetInstance().GetPhiCuts().ApplyAnalysisRapidity(yLab);
}

Double_t StFemtoPhiTreeMaker::DeuteronRapidityCm(const TrackState& trk) const {
  TLorentzVector lv = StNuclearIdHelper::NuclearP4(TrackMomentum(trk), kNucDeuteron);
  return ApplyRapidityFrame(lv.Rapidity());
}

Double_t StFemtoPhiTreeMaker::ProtonRapidityCm(const TrackState& trk) const {
  TVector3 p = TrackMomentum(trk);
  TLorentzVector lv(p.X(), p.Y(), p.Z(), TMath::Sqrt(kProtonMass * kProtonMass + p.Mag2()));
  return ApplyRapidityFrame(lv.Rapidity());
}

Bool_t StFemtoPhiTreeMaker::PassFemtoDeuteronCuts(const TrackState& trk) const {
  const FemtoConfig& fc = ConfigManager::GetInstance().GetFemtoConfig();
  if (trk.charge <= 0) return kFALSE;
  if (trk.DCA >= fc.deuteronMaxDca) return kFALSE;
  TVector3 p = TrackMomentum(trk);
  Double_t pmom = p.Mag();
  if (pmom < fc.deuteronMinPMom || pmom > fc.deuteronMaxPMom) return kFALSE;
  if (trk.pT < fc.deuteronMinPtPre || trk.pT > fc.deuteronMaxPtPre) return kFALSE;
  if (TMath::Abs(trk.eta) >= fc.deuteronMaxAbsEta) return kFALSE;
  if (TMath::Abs(trk.nSigmaDeuteron) >= fc.deuteronMaxAbsNSigma) return kFALSE;
  if (trk.nHitsFit < fc.deuteronMinNHitsFit) return kFALSE;
  if (trk.nHitsMax <= 0) return kFALSE;
  if ((Float_t)trk.nHitsFit / (Float_t)trk.nHitsMax < fc.deuteronMinNHitsRatio) return kFALSE;
  const Bool_t passTofRule =
      (pmom < fc.deuteronTofMomentumThreshold) ||
      (pmom > fc.deuteronTofMomentumThreshold && trk.tofMatch && trk.mass2 >= fc.deuteronMinMass2 &&
       trk.mass2 <= fc.deuteronMaxMass2);
  if (!passTofRule) return kFALSE;
  if (trk.pT < fc.deuteronMinPtPair || trk.pT > fc.deuteronMaxPtPair) return kFALSE;
  Double_t yCm = DeuteronRapidityCm(trk);
  if (yCm < fc.deuteronMinRapidityCm || yCm > fc.deuteronMaxRapidityCm) return kFALSE;
  return kTRUE;
}

Bool_t StFemtoPhiTreeMaker::PassFemtoProtonCuts(const TrackState& trk) const {
  const FemtoConfig& fc = ConfigManager::GetInstance().GetFemtoConfig();
  if (fc.protonChargeMode == "positive" && trk.charge <= 0) return kFALSE;
  if (fc.protonChargeMode == "negative" && trk.charge >= 0) return kFALSE;
  if (trk.DCA >= fc.protonMaxDca) return kFALSE;
  if (trk.pT < fc.protonMinPtPre) return kFALSE;
  if (TMath::Abs(trk.eta) >= fc.protonMaxAbsEta) return kFALSE;
  if (TMath::Abs(trk.nSigmaProton) >= fc.protonMaxAbsNSigma) return kFALSE;
  if (trk.nHitsFit < fc.protonMinNHitsFit) return kFALSE;
  if (trk.nHitsMax <= 0) return kFALSE;
  if ((Float_t)trk.nHitsFit / (Float_t)trk.nHitsMax < fc.protonMinNHitsRatio) return kFALSE;
  if (!PassTofProtonPid(trk)) return kFALSE;
  if (trk.pT < fc.protonMinPtPair || trk.pT > fc.protonMaxPtPair) return kFALSE;
  Double_t yCm = ProtonRapidityCm(trk);
  if (yCm < fc.protonMinRapidityCm || yCm > fc.protonMaxRapidityCm) return kFALSE;
  return kTRUE;
}

// Packed (schema 2) track row. Every narrowing conversion goes through Pack*, which
// counts saturations so that Finish() can assert none happened.
void StFemtoPhiTreeMaker::FillTrackTreeV2(const TrackState& src, UChar_t species,
                                          ULong64_t eventUID) {
  if (!mTrackTree) return;
  using namespace femto_phi_tree;
  mTrk2.Reset();
  mTrk2.eventUID = eventUID;
  if (src.trackIndex < 0 || src.trackIndex > 65535) {
    ++mPackStats.trackIndex;
    mTrk2.trackIndex = 0;
  } else {
    mTrk2.trackIndex = (UShort_t)src.trackIndex;
  }
  mTrk2.speciesCode = species;

  TVector3 pmom(src.momentumX, src.momentumY, src.momentumZ);
  mTrk2.pT = PackU16(pmom.Perp(), scale::kPt, mPackStats.pT);
  mTrk2.eta = PackI16(pmom.PseudoRapidity(), scale::kEta, mPackStats.eta);
  mTrk2.phi = PackI16(pmom.Phi(), scale::kPhi, mPackStats.phi);
  mTrk2.dEdx = PackU16(src.dEdx, scale::kDedx, mPackStats.dEdx);
  mTrk2.nSigmaKaon = PackNSigma(src.nSigmaKaon, mPackStats.nSigmaKaon, mPackStats.nSigmaSentinel);
  mTrk2.nSigmaPion = PackNSigma(src.nSigmaPion, mPackStats.nSigmaPion, mPackStats.nSigmaSentinel);
  mTrk2.nSigmaProton =
      PackNSigma(src.nSigmaProton, mPackStats.nSigmaProton, mPackStats.nSigmaSentinel);
  mTrk2.nSigmaDeuteron =
      PackNSigma(src.nSigmaDeuteron, mPackStats.nSigmaDeuteron, mPackStats.nSigmaSentinel);
  mTrk2.tofBeta = PackI16(src.tofBeta, scale::kTofBeta, mPackStats.tofBeta);
  mTrk2.dca = PackU16(src.DCA, scale::kDca, mPackStats.dca);

  Int_t nhf = src.nHitsFit < 0 ? -src.nHitsFit : src.nHitsFit;
  if (nhf > 127) { nhf = 127; ++mPackStats.nHits; }
  mTrk2.nHitsFit = (Char_t)(src.charge >= 0 ? nhf : -nhf);
  Int_t nhm = src.nHitsMax, nhd = src.nHitsDedx;
  if (nhm < 0 || nhm > 255) { nhm = (nhm < 0) ? 0 : 255; ++mPackStats.nHits; }
  if (nhd < 0 || nhd > 255) { nhd = (nhd < 0) ? 0 : 255; ++mPackStats.nHits; }
  mTrk2.nHitsMax = (UChar_t)nhm;
  mTrk2.nHitsDedx = (UChar_t)nhd;

  UInt_t flags = src.selFlags;
  if (src.tofMatch) flags |= kSelTofMatch;
  mTrk2.selFlags = (UShort_t)(flags & 0xFFFFu);
  mTrackTree->Fill();
}

void StFemtoPhiTreeMaker::FillTrackTree(const TrackState& src, UChar_t species, ULong64_t eventUID) {
  if (mSchemaVersion >= femto_phi_tree::kSchemaVersionV2) {
    FillTrackTreeV2(src, species, eventUID);
    return;
  }
  if (!mTrackTree) return;
  mTrk.Reset();
  mTrk.eventUID = eventUID;
  mTrk.trackIndex = src.trackIndex;
  mTrk.speciesCode = species;
  mTrk.charge = src.charge;
  mTrk.px = src.momentumX;
  mTrk.py = src.momentumY;
  mTrk.pz = src.momentumZ;
  mTrk.dEdx = src.dEdx;
  mTrk.nSigmaKaon = src.nSigmaKaon;
  mTrk.nSigmaDeuteron = src.nSigmaDeuteron;
  mTrk.nSigmaProton = src.nSigmaProton;
  mTrk.tofMatch = src.tofMatch;
  mTrk.tofBeta = src.tofBeta;
  mTrk.mass2 = src.mass2;
  mTrk.deltaOneOverBeta = src.deltaOneOverBeta;
  mTrk.dca = src.DCA;
  mTrk.nHitsFit = src.nHitsFit;
  mTrk.nHitsMax = src.nHitsMax;
  mTrk.nHitsDedx = src.nHitsDedx;
  mTrk.chi2 = src.chi2;
  mTrk.originX = src.originX;
  mTrk.originY = src.originY;
  mTrk.originZ = src.originZ;
  mTrk.bField = src.BField;
  mTrk.selFlags = src.selFlags;
  mTrackTree->Fill();
}

void StFemtoPhiTreeMaker::FillOptionalPhiPairs(const std::vector<TrackState>& kp,
                                               const std::vector<TrackState>& km, ULong64_t eventUID) {
  if (!mPairTree) return;
  PhiCutConfig& phiCfg = ConfigManager::GetInstance().GetPhiCuts();
  for (size_t i = 0; i < kp.size(); ++i) {
    for (size_t j = 0; j < km.size(); ++j) {
      TVector3 dcaPosPlus, dcaPosMinus;
      Double_t invMass = 0.0;
      TVector3 phiMom;
      if (!StPhiKKReconstruction::ReconstructPhi(ToPhiKkTrack(kp[i]), ToPhiKkTrack(km[j]), invMass, phiMom,
                                                 dcaPosPlus, dcaPosMinus)) {
        continue;
      }
      Double_t dcaKK = StPhiKKReconstruction::CalculateDCA(ToPhiKkTrack(kp[i]), ToPhiKkTrack(km[j]), dcaPosPlus,
                                                          dcaPosMinus);
      Double_t opening = StPhiKKReconstruction::CalculateOpeningAngle(ToPhiKkTrack(kp[i]), ToPhiKkTrack(km[j]));
      Double_t yLab = StPhiKKReconstruction::CalculatePairRapidity(invMass, phiMom);
      Double_t yPair = ApplyRapidityFrame(yLab);
      UInt_t flags = 0;
      if (StPhiKKReconstruction::PassPairTofCut(ToPhiKkTrack(kp[i]), ToPhiKkTrack(km[j]))) {
        flags |= femto_phi_tree::kPairPassTof;
      }
      if (opening >= phiCfg.minOpeningAngle && opening <= phiCfg.maxOpeningAngle &&
          yPair >= phiCfg.minPairRapidity && yPair <= phiCfg.maxPairRapidity) {
        flags |= femto_phi_tree::kPairPassKinematics;
      }
      if (!(flags & femto_phi_tree::kPairPassKinematics)) continue;
      mPair.Reset();
      mPair.eventUID = eventUID;
      mPair.dauPlusIndex = kp[i].trackIndex;
      mPair.dauMinusIndex = km[j].trackIndex;
      mPair.px = (Float_t)phiMom.X();
      mPair.py = (Float_t)phiMom.Y();
      mPair.pz = (Float_t)phiMom.Z();
      mPair.mKK = (Float_t)invMass;
      mPair.dcaDaughters = (Float_t)dcaKK;
      mPair.openingAngle = (Float_t)opening;
      mPair.rapidity = (Float_t)yPair;
      mPair.pairFlags = flags;
      mPairTree->Fill();
      mNPhiPair++;
    }
  }
}

void StFemtoPhiTreeMaker::MixBins(Float_t vz, Int_t cent9, Double_t psi2, Int_t& vzBin, Int_t& centBin,
                                  Int_t& epBin, Int_t& mixBin) const {
  const EventCutConfig& ev = ConfigManager::GetInstance().GetEventCuts();
  const MixingConfig& mix = ConfigManager::GetInstance().GetMixingConfig();
  vzBin = 0;
  if (mix.nVzBins > 0) {
    Double_t vzSpan = ev.maxVz - ev.minVz;
    if (vzSpan > 0) {
      vzBin = (Int_t)((vz - ev.minVz) / vzSpan * mix.nVzBins);
      if (vzBin < 0) vzBin = 0;
      if (vzBin >= mix.nVzBins) vzBin = mix.nVzBins - 1;
    }
  }
  centBin = 0;
  if (mix.nCentralityBins > 0 && cent9 >= 0) {
    centBin = cent9;
    if (centBin >= mix.nCentralityBins) centBin = mix.nCentralityBins - 1;
  }
  epBin = 0;
  if (mix.nEventPlaneBins > 1 && psi2 >= 0) {
    epBin = (Int_t)(psi2 / TMath::Pi() * mix.nEventPlaneBins);
    if (epBin < 0) epBin = 0;
    if (epBin >= mix.nEventPlaneBins) epBin = mix.nEventPlaneBins - 1;
  }
  mixBin = vzBin + mix.nVzBins * (centBin + mix.nCentralityBins * epBin);
}

TString StFemtoPhiTreeMaker::CurrentPicoPath() const {
  if (!mPicoDstMaker || !mPicoDstMaker->chain() || !mPicoDstMaker->chain()->GetFile()) return "";
  return mPicoDstMaker->chain()->GetFile()->GetName();
}

Long64_t StFemtoPhiTreeMaker::CurrentEntry() const {
  if (!mPicoDstMaker || !mPicoDstMaker->chain()) return -1;
  return mPicoDstMaker->chain()->GetReadEntry();
}

Int_t StFemtoPhiTreeMaker::Make() {
  if (mSkipRemaining > 0) {
    mSkipRemaining--;
    return kStOK;
  }
  mNInput++;
  if (!mPicoDstMaker) {
    mNNoPico++;
    return kStWarn;
  }
  mPicoDst = mPicoDstMaker->picoDst();
  if (!mPicoDst) {
    mNNoPico++;
    return kStWarn;
  }
  StPicoEvent* event = mPicoDst->event();
  if (!event) {
    mNNoPico++;
    return kStWarn;
  }

  TVector3 pVtx = event->primaryVertex();
  Float_t vzVpd = event->vzVpd();
  Int_t refMult = event->refMult();
  EventCutConfig& evCuts = ConfigManager::GetInstance().GetEventCuts();
  Float_t vr = evCuts.ComputeVr(pVtx.X(), pVtx.Y());
  const Double_t vz = pVtx.Z();
  const Int_t runId = event->runId();
  const Int_t eventId = event->eventId();
  const Int_t nBTOFMatch = event->nBTOFMatch();

  CentralityCutConfig& centCfg = ConfigManager::GetInstance().GetCentralityCuts();
  Int_t rawMult = refMult;
  if (centCfg.enabled) {
    TString mode(centCfg.mode.c_str());
    mode.ToLower();
    if (mode == "fxtmult") rawMult = event->fxtMult();
  }

  Int_t cent9 = -1, cent16 = -1;
  Double_t refMultCorr = -1.0, centWeight = 1.0, centralityPercent = -1.0, psi2 = -1.0;

  CentralityRejectReason centReason = kCentralityOk;
  if (m_centrality && m_centrality->IsEnabled()) {
    if (!m_centrality->CheckBadRun(runId, centReason)) {
      mNBadRun++;
      return kStOK;
    }
  }
  if (!PassEventCuts(pVtx.Z(), vr, refMult, vzVpd)) {
    mNFailEventCuts++;
    return kStOK;
  }
  if (m_centrality && m_centrality->IsEnabled()) {
    if (!m_centrality->CheckPileup(rawMult, nBTOFMatch, vz, centReason)) {
      mNPileup++;
      return kStOK;
    }
    if (!m_centrality->ComputeBins(event, rawMult, vz, cent9, cent16, refMultCorr, centWeight, centReason)) {
      mNFailCent++;
      return kStOK;
    }
    if (!m_centrality->AcceptCentBin(cent9, refMultCorr, centReason)) {
      mNFailCent++;
      return kStOK;
    }
    centralityPercent = CentralityHelper::Cent9ToPercentile(cent9);
  }

  PhiCutConfig& phiCfg = ConfigManager::GetInstance().GetPhiCuts();
  const Int_t nTracks = mPicoDst->numberOfTracks();
  if (phiCfg.maxNTr > 0 && nTracks > phiCfg.maxNTr) {
    mNFailMaxNTr++;
    return kStOK;
  }

  std::vector<TrackState> kaonsPlus;
  std::vector<TrackState> kaonsMinus;
  Double_t Qx = 0.0, Qy = 0.0;
  const FemtoConfig& femtoCfg = ConfigManager::GetInstance().GetFemtoConfig();
  const PIDCutConfig& pidCfg = ConfigManager::GetInstance().GetPIDCuts();

  Int_t nKp = 0, nKm = 0, nD = 0, nP = 0;
  const ULong64_t eventUID = femto_phi_tree::MakeEventUID(runId, eventId);

  for (Int_t itrk = 0; itrk < nTracks; itrk++) {
    StPicoTrack* trk = mPicoDst->track(itrk);
    if (!trk) continue;
    TVector3 pMom = trk->pMom();
    Float_t pt = pMom.Perp();
    Float_t eta = pMom.PseudoRapidity();
    Float_t phi = pMom.Phi();
    if (pt >= phiCfg.minPtEp && pt <= phiCfg.maxPtEp && TMath::Abs(eta) < phiCfg.maxEtaEp) {
      Qx += TMath::Cos(2.0 * phi);
      Qy += TMath::Sin(2.0 * phi);
    }
    if (!PassEnvTrackCuts(trk, pVtx)) continue;

    TrackState ts;
    BuildTrackState(ts, trk, event, pVtx, itrk);
    FillTofInfo(ts, trk, pMom, trk->bTofPidTraitsIndex());
    if (PassNominalTrackCuts(trk, pVtx)) ts.selFlags |= femto_phi_tree::kSelTrackQualityNom;
    if (trk->chi2() <= ConfigManager::GetInstance().GetTrackCuts().maxChi2) {
      ts.selFlags |= femto_phi_tree::kSelTrackChi2Nom;
    }
    if (ts.tofMatch) ts.selFlags |= femto_phi_tree::kSelTofMatch;

    // Storage rule for kaons. With envKaonRequireTofOrLowP the store is restricted to what
    // the production daughter PID can ever use: K+ is TPC-only up to pMomKaonPID and needs
    // TOF above it, K- always needs TOF. Signal, ROT, MIX and the K-p kaon species all run
    // through PassPhiDaughterTofPid, so nothing downstream sees a kaon outside this set.
    // envKaonRequireDaughterPidReach narrows it further to the rows PassPhiDaughterTofPid can
    // ever accept: TOF-matched inside a widened m2 window, or TPC-only below pMomKaonPID. The
    // TPC-only branch deliberately keeps BOTH charges, so the K- TPC-only variation
    // (phiDaughterKaonMinusRequireTof: false) stays reproducible. Rows outside the reach only
    // ever fed phi-daughter PID QA histograms, never an observable.
    // The m2 window is a storage envelope and must stay wider than any planned m2 systematic:
    // tightening it forecloses that systematic and needs a full re-production.
    Bool_t kaonStoreOk =
        !mEnvKaonRequireTofOrLowP || ts.tofMatch || (pMom.Mag() < mEnvKaonLowPMax);
    if (kaonStoreOk && mEnvKaonRequireDaughterPidReach) {
      kaonStoreOk = ts.tofMatch
                        ? (ts.mass2 >= mEnvKaonMass2Lo && ts.mass2 <= mEnvKaonMass2Hi)
                        : (pMom.Mag() <= mEnvKaonLowPMax);
    }
    if (mStoreKaons && kaonStoreOk && TMath::Abs(ts.nSigmaKaon) <= mEnvMaxAbsNSigmaKaon &&
        ts.DCA <= mEnvMaxDcaKaon) {
      if (PassNominalKaonCuts(trk, pVtx)) ts.selFlags |= femto_phi_tree::kSelKaonCutsNom;
      if (PassTofKaonPid(ts)) ts.selFlags |= femto_phi_tree::kSelLoosePid;
      if (PassPhiDaughterTofPid(ts)) ts.selFlags |= femto_phi_tree::kSelNominalPid;
      if (ts.charge > 0) {
        kaonsPlus.push_back(ts);
        FillTrackTree(ts, femto_phi_tree::kSpeciesKp, eventUID);
        FillKaonOriginRow(ts, eventUID, pVtx);
        nKp++;
        mNKp++;
      } else if (ts.charge < 0) {
        kaonsMinus.push_back(ts);
        FillTrackTree(ts, femto_phi_tree::kSpeciesKm, eventUID);
        FillKaonOriginRow(ts, eventUID, pVtx);
        nKm++;
        mNKm++;
      }
    }

    if (mStoreDeuterons && trk->nHitsDedx() >= mEnvMinNHitsDedxNuclear && ts.dEdx > 0 && ts.charge > 0) {
      NuclearTrackState nucState;
      StNuclearIdHelper::FillFromPico(nucState, trk, mPicoDst);
      ts.nSigmaDeuteron = (Float_t)StNuclearIdHelper::GetNSigma(kNucDeuteron, nucState.pMag, nucState.dedx);
      const Double_t pmom = pMom.Mag();
      if (TMath::Abs(ts.nSigmaDeuteron) <= mEnvMaxAbsNSigmaDeuteron && ts.DCA <= mEnvDeuteronMaxDca &&
          pmom >= mEnvDeuteronMinPMom && pmom <= mEnvDeuteronMaxPMom) {
        TrackState dTrk = ts;
        if (StNuclearIdHelper::IsDeuteron(nucState)) dTrk.selFlags |= femto_phi_tree::kSelNominalPid;
        dTrk.selFlags |= femto_phi_tree::kSelLoosePid;
        if (PassFemtoDeuteronCuts(dTrk)) dTrk.selFlags |= femto_phi_tree::kSelNominalFemto;
        FillTrackTree(dTrk, femto_phi_tree::kSpeciesDeuteron, eventUID);
        nD++;
        mNDeuteron++;
      }
    }

    if (mStoreProtons && TMath::Abs(ts.nSigmaProton) <= mEnvMaxAbsNSigmaProton && ts.DCA <= mEnvProtonMaxDca &&
        ts.pT >= mEnvProtonMinPt) {
      Bool_t passCharge = kTRUE;
      if (femtoCfg.protonChargeMode == "positive") passCharge = (ts.charge > 0);
      else if (femtoCfg.protonChargeMode == "negative") passCharge = (ts.charge < 0);
      if (passCharge && TMath::Abs(ts.nSigmaProton) <= pidCfg.nSigmaProton * 2.0) {
        TrackState pTrk = ts;
        if (PassTofProtonPid(pTrk)) pTrk.selFlags |= femto_phi_tree::kSelLoosePid;
        if (PassTofProtonPid(pTrk) && TMath::Abs(pTrk.nSigmaProton) <= pidCfg.nSigmaProton) {
          pTrk.selFlags |= femto_phi_tree::kSelNominalPid;
        }
        if (PassFemtoProtonCuts(pTrk)) pTrk.selFlags |= femto_phi_tree::kSelNominalFemto;
        FillTrackTree(pTrk, femto_phi_tree::kSpeciesProton, eventUID);
        nP++;
        mNProton++;
      }
    }
  }

  TVector2 Q(Qx, Qy);
  if (Q.Mod() > 0) {
    psi2 = 0.5 * TMath::ATan2(Qy, Qx);
    if (psi2 < 0) psi2 += TMath::Pi();
  }

  if (mWritePhiPair) FillOptionalPhiPairs(kaonsPlus, kaonsMinus, eventUID);

  Int_t mixVz = 0, mixCent = 0, mixEp = 0, mixBin = 0;
  MixBins((Float_t)vz, cent9, psi2, mixVz, mixCent, mixEp, mixBin);

  mEvt.Reset();
  mEvt.eventUID = eventUID;
  mEvt.runId = runId;
  mEvt.eventId = eventId;
  mEvt.sourceFileHash = femto_phi_tree::HashStringFnv(CurrentPicoPath().Data());
  mEvt.sourceEntry = CurrentEntry();
  mEvt.subjobId = mSubjobId;
  mEvt.vx = (Float_t)pVtx.X();
  mEvt.vy = (Float_t)pVtx.Y();
  mEvt.vz = (Float_t)vz;
  mEvt.vr = vr;
  mEvt.vzVpd = vzVpd;
  mEvt.bField = event->bField();
  mEvt.refMult = refMult;
  mEvt.rawMult = rawMult;
  mEvt.nBTOFMatch = nBTOFMatch;
  mEvt.nTracks = nTracks;
  mEvt.refMultCorr = (Float_t)refMultCorr;
  mEvt.centWeight = (Float_t)centWeight;
  mEvt.centralityPercent = (Float_t)centralityPercent;
  mEvt.cent9 = cent9;
  mEvt.cent16 = cent16;
  mEvt.qx = (Float_t)Qx;
  mEvt.qy = (Float_t)Qy;
  mEvt.psi2 = (Float_t)psi2;
  mEvt.mixVzBin = mixVz;
  mEvt.mixCentBin = mixCent;
  mEvt.mixEpBin = mixEp;
  mEvt.mixBin = mixBin;
  mEvt.eventFlags = femto_phi_tree::kEvtPassEventCuts | femto_phi_tree::kEvtPassPileup |
                    femto_phi_tree::kEvtPassCent | femto_phi_tree::kEvtPassMaxNTr;
  mEvt.nKp = nKp;
  mEvt.nKm = nKm;
  mEvt.nDeuteron = nD;
  mEvt.nProton = nP;

  if (mSchemaVersion >= femto_phi_tree::kSchemaVersionV2) {
    mEvt2.Reset();
    mEvt2.eventUID = mEvt.eventUID;
    mEvt2.runId = mEvt.runId;
    mEvt2.eventId = mEvt.eventId;
    mEvt2.sourceFileHash = mEvt.sourceFileHash;
    mEvt2.sourceFileIndex = SourceFileIndex(CurrentPicoPath());
    mEvt2.sourceEntry = mEvt.sourceEntry;
    mEvt2.subjobId = mEvt.subjobId;
    FillTriggerInfo(event);
    mEvt2.vx = mEvt.vx;  mEvt2.vy = mEvt.vy;  mEvt2.vz = mEvt.vz;
    mEvt2.vr = mEvt.vr;  mEvt2.vzVpd = mEvt.vzVpd;  mEvt2.bField = mEvt.bField;
    mEvt2.refMult = mEvt.refMult;  mEvt2.rawMult = mEvt.rawMult;
    mEvt2.nBTOFMatch = mEvt.nBTOFMatch;  mEvt2.nTracks = mEvt.nTracks;
    mEvt2.refMultCorr = mEvt.refMultCorr;  mEvt2.centWeight = mEvt.centWeight;
    mEvt2.centralityPercent = mEvt.centralityPercent;
    mEvt2.cent9 = mEvt.cent9;  mEvt2.cent16 = mEvt.cent16;
    mEvt2.qx = mEvt.qx;  mEvt2.qy = mEvt.qy;  mEvt2.psi2 = mEvt.psi2;
    mEvt2.mixVzBin = mEvt.mixVzBin;  mEvt2.mixCentBin = mEvt.mixCentBin;
    mEvt2.mixEpBin = mEvt.mixEpBin;  mEvt2.mixBin = mEvt.mixBin;
    mEvt2.eventFlags = mEvt.eventFlags;
    mEvt2.nKp = nKp;  mEvt2.nKm = nKm;  mEvt2.nDeuteron = nD;  mEvt2.nProton = nP;
  }

  if (mSchemaVersion >= femto_phi_tree::kSchemaVersionV3) FillEventRowV3();
  if (mEventTree) mEventTree->Fill();
  mNAccepted++;
  return kStOK;
}

// sourceFileHash is a 32-bit FNV of the PicoDst path and cannot be inverted; with
// 99,246 files in the dataset the expected number of colliding pairs is 1.15.
// The index into this per-output-file table is the invertible provenance key
// (plan sec 10.7).
UShort_t StFemtoPhiTreeMaker::SourceFileIndex(const TString& path) {
  const std::string key = path.Data();
  for (size_t i = 0; i < mSourceFiles.size(); ++i)
    if (mSourceFiles[i] == key) return (UShort_t)i;
  if (mSourceFiles.size() >= 65535) {
    std::cerr << "[StFemtoPhiTreeMaker] more than 65535 source files in one job" << std::endl;
    return 65535;
  }
  mSourceFiles.push_back(key);
  return (UShort_t)(mSourceFiles.size() - 1);
}

// Record which of the configured triggers this event fired. No selection is applied
// here: the tree stores the information so that a trigger variation can be run
// downstream (plan sec 3.1).
void StFemtoPhiTreeMaker::FillTriggerInfo(StPicoEvent* event) {
  mEvt2.triggerId = 0;
  mEvt2.triggerBits = 0;
  if (!event) return;
  const std::vector<unsigned int>& ids = event->triggerIds();
  for (size_t k = 0; k < mTriggerIds.size() && k < 16; ++k)
    for (size_t j = 0; j < ids.size(); ++j)
      if (ids[j] == mTriggerIds[k]) {
        mEvt2.triggerBits |= (UShort_t)(1u << k);
        if (mEvt2.triggerId == 0) mEvt2.triggerId = mTriggerIds[k];
      }
}

// The hash -> path table that makes sourceFileHash / sourceFileIndex invertible.
void StFemtoPhiTreeMaker::WriteSourceFileTable() {
  if (!mOutFile) return;
  mOutFile->cd();
  UShort_t idx = 0;
  UInt_t hash = 0;
  TString path;
  TString* pathPtr = &path;
  TTree* t = new TTree("SourceFileTable", "sourceFileIndex -> PicoDst path");
  t->Branch("sourceFileIndex", &idx, "sourceFileIndex/s");
  t->Branch("sourceFileHash", &hash, "sourceFileHash/i");
  t->Branch("path", &pathPtr);
  for (size_t i = 0; i < mSourceFiles.size(); ++i) {
    idx = (UShort_t)i;
    path = mSourceFiles[i].c_str();
    hash = femto_phi_tree::HashStringFnv(mSourceFiles[i].c_str());
    t->Fill();
  }
  t->Write();
}

static void WriteNamed(const char* name, const TString& title) {
  TNamed obj(name, title.Data());
  obj.Write();
}

void StFemtoPhiTreeMaker::WriteMetadata() {
  if (!mOutFile) return;
  mOutFile->cd();
  WriteNamed("schemaVersion", TString::Format("%u", mSchemaVersion));
  WriteNamed("pidCorrectionState", mPidCorrectionState.c_str());
  {
    TString trg;
    for (size_t i = 0; i < mTriggerIds.size(); ++i)
      trg += TString::Format("%s%u", i ? "," : "", mTriggerIds[i]);
    WriteNamed("triggerIds", trg);
  }
  WriteNamed("nSourceFiles", TString::Format("%d", (int)mSourceFiles.size()));
  WriteNamed("packSaturationTotal", TString::Format("%lld", (long long)mPackStats.Total()));
  WriteNamed("packNSigmaSentinel", TString::Format("%lld", (long long)mPackStats.nSigmaSentinel));
  WriteNamed("packSaturationDetail",
             TString::Format("pT=%lld eta=%lld phi=%lld dEdx=%lld nSigma=%lld tofBeta=%lld "
                             "dca=%lld trackIndex=%lld nHits=%lld nSigmaK=%lld nSigmaPi=%lld "
                             "nSigmaP=%lld nSigmaD=%lld sentinel=%lld evVertex=%lld "
                             "evQ=%lld evMult=%lld evCent=%lld",
                             (long long)mPackStats.pT, (long long)mPackStats.eta,
                             (long long)mPackStats.phi, (long long)mPackStats.dEdx,
                             (long long)0, (long long)mPackStats.tofBeta,
                             (long long)mPackStats.dca, (long long)mPackStats.trackIndex,
                             (long long)mPackStats.nHits, (long long)mPackStats.nSigmaKaon,
                             (long long)mPackStats.nSigmaPion, (long long)mPackStats.nSigmaProton,
                             (long long)mPackStats.nSigmaDeuteron,
                             (long long)mPackStats.nSigmaSentinel, (long long)mPackStats.evVertex,
                             (long long)mPackStats.evQ, (long long)mPackStats.evMult,
                             (long long)mPackStats.evCent));
  WriteNamed("gitRevision", EnvOrEmpty("STAR_ANA_GIT_REV"));
  WriteNamed("gitDirty", EnvOrEmpty("STAR_ANA_GIT_DIRTY"));
  WriteNamed("mainconfPath", mMainconfPath.c_str());
  WriteNamed("makerYamlPath", mMakerYamlPath.c_str());
  WriteNamed("starLibraryTag", EnvOrEmpty("STAR_ANA_LIBRARY_TAG"));
  WriteNamed("rootVersion", gROOT->GetVersion());
  WriteNamed("compressionLevel", TString::Format("%d", mCompressLevel));
  WriteNamed("nInputEvents", TString::Format("%lld", (long long)mNInput));
  WriteNamed("nAcceptedEvents", TString::Format("%lld", (long long)mNAccepted));
  WriteNamed("nBadRun", TString::Format("%lld", (long long)mNBadRun));
  WriteNamed("nFailEventCuts", TString::Format("%lld", (long long)mNFailEventCuts));
  WriteNamed("nPileup", TString::Format("%lld", (long long)mNPileup));
  WriteNamed("nFailCent", TString::Format("%lld", (long long)mNFailCent));
  WriteNamed("nFailMaxNTr", TString::Format("%lld", (long long)mNFailMaxNTr));
  WriteNamed("nKp", TString::Format("%lld", (long long)mNKp));
  WriteNamed("nKm", TString::Format("%lld", (long long)mNKm));
  WriteNamed("nDeuteron", TString::Format("%lld", (long long)mNDeuteron));
  WriteNamed("nProton", TString::Format("%lld", (long long)mNProton));
  WriteNamed("nPhiPair", TString::Format("%lld", (long long)mNPhiPair));
  WriteNamed("derivedCandidatesStored", "none: no phi_rot, phi_mix, SE/ME pair rows");
  WriteCounterHistogram();
  if (mSchemaVersion >= femto_phi_tree::kSchemaVersionV2) WriteSourceFileTable();
}

// TNamed metadata is overwritten, not summed, by `hadd`: a merged file reports the
// counters of one input only. The same counters therefore also go into a labelled
// TH1D, which hadd adds bin by bin, so closure can be checked after merging.
void StFemtoPhiTreeMaker::WriteCounterHistogram() {
  if (!mOutFile) return;
  mOutFile->cd();
  const char* label[13] = {"nInputEvents", "nAcceptedEvents", "nBadRun",    "nFailEventCuts",
                           "nPileup",      "nFailCent",       "nFailMaxNTr", "nKp",
                           "nKm",          "nDeuteron",       "nProton",     "nPhiPair",
                           "nFiles"};
  const Long64_t value[13] = {mNInput,      mNAccepted, mNBadRun, mNFailEventCuts,
                              mNPileup,     mNFailCent, mNFailMaxNTr, mNKp,
                              mNKm,         mNDeuteron, mNProton, mNPhiPair,
                              1};
  TH1D h("hCounters", "summable production counters;;count", 13, 0.5, 13.5);
  h.SetDirectory(0);
  for (Int_t i = 0; i < 13; ++i) {
    h.GetXaxis()->SetBinLabel(i + 1, label[i]);
    h.SetBinContent(i + 1, (Double_t)value[i]);
    h.SetBinError(i + 1, 0.0);
  }
  h.SetEntries(mNInput > 0 ? (Double_t)mNInput : 1.0);
  h.Write();
}

Int_t StFemtoPhiTreeMaker::Finish() {
  if (m_centrality) m_centrality->Finish();
  std::cout << "[StFemtoPhiTreeMaker] nInput=" << mNInput << " nAccepted=" << mNAccepted
            << " nBadRun=" << mNBadRun << " nFailEventCuts=" << mNFailEventCuts
            << " nPileup=" << mNPileup << " nFailCent=" << mNFailCent
            << " nKp=" << mNKp << " nKm=" << mNKm << " nD=" << mNDeuteron << " nP=" << mNProton
            << " nPhiPair=" << mNPhiPair << std::endl;
  if (mSchemaVersion >= femto_phi_tree::kSchemaVersionV2) {
    std::cout << "[StFemtoPhiTreeMaker] schema=" << mSchemaVersion
              << " sourceFiles=" << mSourceFiles.size()
              << " packSaturation total=" << mPackStats.Total()
              << " (pT=" << mPackStats.pT << " eta=" << mPackStats.eta
              << " phi=" << mPackStats.phi << " dEdx=" << mPackStats.dEdx
              << " tofBeta=" << mPackStats.tofBeta << " dca=" << mPackStats.dca
              << " trackIndex=" << mPackStats.trackIndex << " nHits=" << mPackStats.nHits
              << " nSigmaK=" << mPackStats.nSigmaKaon << " nSigmaPi=" << mPackStats.nSigmaPion
              << " nSigmaP=" << mPackStats.nSigmaProton
              << " nSigmaD=" << mPackStats.nSigmaDeuteron
              << " evVertex=" << mPackStats.evVertex << " evQ=" << mPackStats.evQ
              << " evMult=" << mPackStats.evMult << " evCent=" << mPackStats.evCent << " kaonOrigin=" << mPackStats.kaonOrigin << ")"
              << "  nSigmaSentinel=" << mPackStats.nSigmaSentinel << " (intended)" << std::endl;
    if (mPackStats.Total() > 0)
      std::cerr << "[StFemtoPhiTreeMaker] WARNING: packed values were saturated; a stored "
                << "quantity left its representable range (plan Step 2 requires zero)"
                << std::endl;
  }
  if (mOutFile) {
    mOutFile->cd();
    if (mEventTree) mEventTree->Write();
    if (mTrackTree) mTrackTree->Write();
    if (mPairTree) mPairTree->Write();
    // Without this the companion loses everything after its last autoflush -- 7,877 of 37,877
    // rows in the first run, which the reader saw only as a suspiciously round 30,000.
    if (mKaonOriginTree) mKaonOriginTree->Write();
    WriteFlagConfigSnapshot();
    WriteMetadata();
    mOutFile->Close();
    delete mOutFile;
    mOutFile = 0;
    mEventTree = 0;
    mTrackTree = 0;
    mPairTree = 0;
  }
  return kStOK;
}
