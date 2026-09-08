// ROOT5/CINT-compatible, read-only Imp5 validation layered on structural QA.
// root4star -b -q 'tests/check_kfparticle_imp5_output.C("result.root")'
// Thresholds come exclusively from the output's effective configuration.
// This checks stored candidate quantities, not an independent PicoDst rerun.
#include "check_kfparticle_output.C"

Bool_t KfImp5MetadataBool(const TString& configuration, const char* key, Bool_t& value)
{
  const TString token = KfOutputMetadataValue(configuration, key);
  if (token == "true" || token == "1") { value = kTRUE; return kTRUE; }
  if (token == "false" || token == "0") { value = kFALSE; return kTRUE; }
  return kFALSE;
}

Int_t KfImp5OutputQaInspect(const char* filename, Bool_t requireCandidates)
{
  const Int_t structural = KfOutputQaInspect(filename, requireCandidates);
  if (structural != 0) return structural;
  TFile* file = TFile::Open(filename, "READ");
  if (!file || file->IsZombie())
    return KfOutputQaFailure(file, "cannot reopen Imp5 output read-only", 7);
  TObjString* configuration = dynamic_cast<TObjString*>(file->Get("KFParticleEffectiveConfiguration"));
  TObjString* eventConfiguration = dynamic_cast<TObjString*>(file->Get("KFEventSelectionConfiguration"));
  if (!configuration || !eventConfiguration)
    return KfOutputQaFailure(file, "Imp5 effective configuration is missing", 7);
  const TString config = configuration->GetString();
  Bool_t legacy = kFALSE, useTof = kTRUE, antiEnabled = kFALSE;
  Bool_t strictTof = kTRUE, cleanKaons = kTRUE, hftOnly = kTRUE;
  if (!KfOutputEventPolicy(eventConfiguration->GetString(), config, legacy) || !legacy ||
      KfOutputMetadataValue(config, "selectionProfile") != "lambda_imp5" ||
      KfOutputMetadataValue(config, "pidProfile") != "pico_nsigma" ||
      !KfImp5MetadataBool(config, "useTof", useTof) || useTof ||
      !KfImp5MetadataBool(config, "strictTofPid", strictTof) || strictTof ||
      !KfImp5MetadataBool(config, "cleanKaonsWithTof", cleanKaons) || cleanKaons ||
      !KfImp5MetadataBool(config, "useHftTracksOnly", hftOnly) || hftOnly ||
      !KfImp5MetadataBool(config, "reconstructAntiLambda", antiEnabled))
    return KfOutputQaFailure(file, "output is not a consistent completed TPC-only Imp5 profile", 7);

  enum { nThresholds = 16, nLeaves = 21 };
  const char* thresholdNames[nThresholds] = {
    "nSigmaProton", "nSigmaPion", "imp5MinDCAProton", "imp5MinDCAPion",
    "imp5MaxPathLength", "maxDaughterDistance", "maxDistanceToPv", "minCosPointing",
    "minMass", "maxMass", "maxMassError", "maxChi2Ndf", "maxTopoChi2Ndf",
    "minDecayLength", "minDecayLengthSignificance", "minVertexLineSignificance"
  };
  Double_t cut[nThresholds];
  Int_t i = 0;
  for (i = 0; i < nThresholds; ++i) {
    if (!KfOutputMetadataDouble(config, thresholdNames[i], cut[i])) {
      std::cerr << "Missing/nonfinite threshold: " << thresholdNames[i] << std::endl;
      return KfOutputQaFailure(file, "Imp5 threshold metadata are incomplete", 7);
    }
  }
  if (cut[0] <= 0. || cut[1] <= 0. || cut[2] < 0. || cut[3] < 0. || cut[4] <= 0. ||
      cut[5] < 0. || cut[6] < 0. || cut[7] <= -1. || cut[7] > 1. ||
      cut[8] <= 0. || cut[9] <= cut[8])
    return KfOutputQaFailure(file, "Imp5 common-cut metadata contain disabled/invalid thresholds", 7);

  TTree* tree = dynamic_cast<TTree*>(file->Get("KfLambdaCandidates"));
  if (!tree) return KfOutputQaFailure(file, "Imp5 candidate tree is missing", 7);
  const char* leafNames[nLeaves] = {
    "protonPidPull", "pionPidPull", "protonDcaToPv", "pionDcaToPv",
    "protonHelixPathLength", "pionHelixPathLength", "daughterDistance",
    "distanceToPv", "cosPointing", "mass", "massError", "chi2Ndf", "topoChi2Ndf",
    "decayLength", "decayLengthSignificance", "vertexLineLengthSignificance",
    "pdg", "selected", "imp5PathValid", "protonHasTof", "pionHasTof"
  };
  TLeaf* leaf[nLeaves];
  for (i = 0; i < nLeaves; ++i) {
    leaf[i] = tree->GetLeaf(leafNames[i]);
    if (!leaf[i]) {
      std::cerr << "Missing leaf: " << leafNames[i] << std::endl;
      return KfOutputQaFailure(file, "Imp5 candidate audit leaf is missing", 7);
    }
  }
  // The original gDCA/path diagnostics must keep Double_t precision.
  for (i = 2; i <= 5; ++i)
    if (TString(leaf[i]->GetTypeName()) != "Double_t")
      return KfOutputQaFailure(file, "Imp5 DCA/path diagnostic is not Double_t", 7);
  for (i = 17; i < nLeaves; ++i)
    if (TString(leaf[i]->GetTypeName()) != "Bool_t")
      return KfOutputQaFailure(file, "Imp5 selected/path/TOF flag is not Bool_t", 7);

  Long64_t survived[7] = {0,0,0,0,0,0,0};
  Long64_t selectedCount = 0, invalidPaths = 0;
  const Double_t floatTolerance = 1.e-6; // Float_t serialization/rounding QA only.
  for (Long64_t entry = 0; entry < tree->GetEntries(); ++entry) {
    if (tree->GetEntry(entry) <= 0)
      return KfOutputQaFailure(file, "cannot read Imp5 candidate entry", 7);
    Double_t v[nLeaves];
    for (i = 0; i < nLeaves; ++i) v[i] = leaf[i]->GetValue();
    const Bool_t selected = v[17] == 1.;
    const Bool_t validPath = v[18] == 1.;
    for (i = 17; i < nLeaves; ++i)
      if (v[i] != 0. && v[i] != 1.)
        return KfOutputQaFailure(file, "non-boolean Imp5 selection/path/TOF flag", 7);
    if (v[19] != 0. || v[20] != 0.)
      return KfOutputQaFailure(file, "TOF PID appears in a TPC-only Imp5 candidate", 7);
    if (!antiEnabled && v[16] != 3122.)
      return KfOutputQaFailure(file, "anti-Lambda appears although reconstruction is disabled", 7);
    // PID and original gDCA are input-hypothesis cuts, so every raw candidate
    // must already pass them. These tests do not measure track-hit efficiency.
    if (!TMath::Finite(v[0]) || !TMath::Finite(v[1]) ||
        TMath::Abs(v[0]) > cut[0] + floatTolerance ||
        TMath::Abs(v[1]) > cut[1] + floatTolerance ||
        !TMath::Finite(v[2]) || !TMath::Finite(v[3]) || v[2] < cut[2] || v[3] < cut[3]) {
      std::cerr << "Imp5 PID/gDCA violation at raw entry " << entry << std::endl;
      return KfOutputQaFailure(file, "raw candidate violates Imp5 daughter PID/gDCA cuts", 7);
    }
    if (validPath && (!TMath::Finite(v[4]) || !TMath::Finite(v[5])))
      return KfOutputQaFailure(file, "valid Imp5 path flag has a nonfinite path length", 7);
    if (!validPath) ++invalidPaths;
    const Bool_t passPath = validPath && TMath::Finite(v[4]) && TMath::Finite(v[5]) &&
        TMath::Abs(v[4]) <= cut[4] && TMath::Abs(v[5]) <= cut[4];
    const Bool_t passDaughter = v[6] <= cut[5];
    const Bool_t passParent = v[7] <= cut[6];
    const Bool_t passPointing = v[8] >= cut[7];
    const Bool_t passMass = v[9] >= cut[8] && v[9] <= cut[9];
    Bool_t passAdditional = kTRUE;
    if (cut[10] >= 0. && v[10] > cut[10]) passAdditional = kFALSE;
    if (cut[11] >= 0. && v[11] > cut[11]) passAdditional = kFALSE;
    if (cut[12] >= 0. && v[12] > cut[12]) passAdditional = kFALSE;
    if (cut[13] >= 0. && v[13] < cut[13]) passAdditional = kFALSE;
    if (cut[14] >= 0. && v[14] < cut[14]) passAdditional = kFALSE;
    if (cut[15] >= 0. && v[15] < cut[15]) passAdditional = kFALSE;
    const Bool_t decisions[6] = {passPath, passDaughter, passParent,
                                 passPointing, passMass, passAdditional};
    Bool_t cumulative = kTRUE;
    ++survived[0];
    for (i = 0; i < 6; ++i) {
      cumulative = cumulative && decisions[i];
      if (cumulative) ++survived[i + 1];
    }
    if (selected) {
      ++selectedCount;
      if (!passPath || v[6] > cut[5] + floatTolerance ||
          v[7] > cut[6] + floatTolerance || v[8] < cut[7] - floatTolerance ||
          v[9] < cut[8] - floatTolerance || v[9] > cut[9] + floatTolerance ||
          !passAdditional) {
        std::cerr << "Imp5 final-selection violation at entry " << entry << std::endl;
        return KfOutputQaFailure(file, "selected candidate violates stored Imp5 final cuts", 7);
      }
    }
  }
  std::cout << "Imp5 cumulative candidate survival (reordered diagnostic cuts, not Maker stages):"
            << "\n  raw after Finder/Topo and input PID/gDCA = " << survived[0]
            << "\n  + original helix path guard = " << survived[1]
            << "\n  + KF daughter distance = " << survived[2]
            << "\n  + KF parent distance to PV = " << survived[3]
            << "\n  + KF pointing cosine = " << survived[4]
            << "\n  + mass range = " << survived[5]
            << "\n  + any additional configured KF final cuts = " << survived[6]
            << "\n  stored selected = " << selectedCount
            << "; raw candidates with invalid path flag = " << invalidPaths << std::endl;
  // Replayed decisions use the exact stored Float_t values and Double_t cuts.
  // No mass constraint, fit, background subtraction, or S/B estimate is made.
  if (survived[6] != selectedCount)
    return KfOutputQaFailure(file, "replayed Imp5 final selection disagrees with stored selected count", 7);
  std::cout << "PASS Imp5 metadata, TPC-only PID, original gDCA/path and final topology/mass cuts."
            << "\nThis is stored-output consistency QA, not a PID calibration or S/B measurement."
            << std::endl;
  file->Close();
  delete file;
  return 0;
}

void check_kfparticle_imp5_output(const char* filename, Bool_t requireCandidates = kTRUE)
{
  const Int_t result = KfImp5OutputQaInspect(filename, requireCandidates);
  gSystem->Exit(result);
}
