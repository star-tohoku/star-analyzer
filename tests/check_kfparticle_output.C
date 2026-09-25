// ROOT5/CINT-compatible, read-only structural/output QA. No fit or cut tuning.
// Example: root4star -b -q 'tests/check_kfparticle_output.C("result.root")'
// requireCandidates=false allows a COMPLETED zero-candidate sample only;
// it never treats KFRunStatus=incomplete or zero reconstructed events as success.
// Current-format event selection metadata/QA are required. Older outputs without
// them fail explicitly; they must not be reported as mode-validated runs.
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TLeaf.h"
#include "TMath.h"
#include "TNamed.h"
#include "TObjString.h"
#include "TString.h"
#include "TSystem.h"
#include "TTree.h"

#include <cfloat>
#include <iostream>

Int_t KfOutputQaFailure(TFile* file, const char* reason, Int_t code)
{
  std::cerr << "FAIL KF output QA: " << reason << std::endl;
  if (file) { file->Close(); delete file; }
  return code;
}

TString KfOutputMetadataValue(const TString& configuration, const char* key)
{
  const TString content = TString("\n") + configuration + "\n";
  const TString prefix = TString("\n") + key + ":";
  Ssiz_t begin = content.Index(prefix);
  if (begin == kNPOS) return TString("");
  begin += prefix.Length();
  const Ssiz_t end = content.Index("\n", begin);
  return TString(content(begin, end - begin)).Strip(TString::kBoth);
}

Bool_t KfOutputMetadataDouble(const TString& configuration, const char* key, Double_t& value)
{
  const TString token = KfOutputMetadataValue(configuration, key);
  if (!token.IsFloat()) return kFALSE;
  value = token.Atof();
  return TMath::Finite(value);
}

// Absence of the policy/flags denotes the original mode-aware output format.
// Only an explicit, self-consistent Imp5 policy may bypass vertex-cut QA.
Bool_t KfOutputEventPolicy(const TString& eventConfiguration,
                          const TString& kfConfiguration, Bool_t& legacy)
{
  legacy = kFALSE;
  const TString policy = KfOutputMetadataValue(eventConfiguration, "KFEventSelection.policy");
  const TString vertex = KfOutputMetadataValue(eventConfiguration, "KFEventSelection.vertexCutsApplied");
  const TString refMult = KfOutputMetadataValue(eventConfiguration, "KFEventSelection.refMultCutsApplied");
  const TString vpd = KfOutputMetadataValue(eventConfiguration, "KFEventSelection.vpdCutsApplied");
  const TString profile = KfOutputMetadataValue(kfConfiguration, "selectionProfile");
  if (profile != "" && profile != "kf_reference" && profile != "lambda_imp5") return kFALSE;
  if (policy == "") {
    return vertex == "" && refMult == "" && vpd == "" && profile != "lambda_imp5";
  }
  if (policy != "mode_vertex" && policy != "legacy_lambda_imp5") return kFALSE;
  legacy = policy == "legacy_lambda_imp5";
  const TString expected = legacy ? "false" : "true";
  if (vertex != expected || refMult != expected || vpd != expected) return kFALSE;
  if (legacy) return profile == "lambda_imp5";
  return profile != "lambda_imp5";
}

// Vertex histograms are unweighted, one fill per eligible event. Include flow
// cells: axis range choices must not silently hide entries from this check.
Bool_t KfOutputVertexHistogram(TH1* histogram, Int_t dimensions, Double_t expected)
{
  if (!histogram || histogram->GetDimension() != dimensions ||
      !TMath::Finite(histogram->GetEntries()) || histogram->GetEntries() != expected)
    return kFALSE;
  Double_t total = 0.;
  // ROOT 5 does not expose TH1::GetNcells(); use the public axis sizes.
  const Int_t cells = (histogram->GetNbinsX() + 2) *
      (dimensions == 2 ? histogram->GetNbinsY() + 2 : 1);
  for (Int_t cell = 0; cell < cells; ++cell) {
    const Double_t count = histogram->GetBinContent(cell);
    if (!TMath::Finite(count) || count < 0. || count != TMath::Floor(count)) return kFALSE;
    total += count;
  }
  if (total != expected) return kFALSE;
  for (Int_t axis = 1; axis <= dimensions; ++axis)
    if (!TMath::Finite(histogram->GetMean(axis)) ||
        !TMath::Finite(histogram->GetRMS(axis))) return kFALSE;
  return kTRUE;
}

// A binned histogram cannot reconstruct the exact extrema. Check that no
// populated cell lies wholly outside the inclusive configured interval.
Bool_t KfOutputHistogramWithin(TH1* histogram, Double_t low, Double_t high)
{
  if (!histogram || histogram->GetDimension() != 1 ||
      !TMath::Finite(low) || !TMath::Finite(high) || high < low) return kFALSE;
  const Double_t tolerance = 64. * DBL_EPSILON *
      (1. + TMath::Max(TMath::Abs(low), TMath::Abs(high)));
  const Int_t bins = histogram->GetNbinsX();
  const TAxis* axis = histogram->GetXaxis();
  if (histogram->GetBinContent(0) > 0. && axis->GetXmin() <= low) return kFALSE;
  if (histogram->GetBinContent(bins + 1) > 0. && axis->GetXmax() > high + tolerance) return kFALSE;
  for (Int_t bin = 1; bin <= bins; ++bin) {
    if (histogram->GetBinContent(bin) <= 0.) continue;
    if (axis->GetBinUpEdge(bin) < low - tolerance ||
        axis->GetBinLowEdge(bin) > high + tolerance) return kFALSE;
  }
  if (histogram->GetEntries() > 0. &&
      (histogram->GetMean() < low - tolerance || histogram->GetMean() > high + tolerance))
    return kFALSE;
  return kTRUE;
}

Int_t KfOutputQaInspect(const char* filename, Bool_t requireCandidates)
{
  if (!filename || !*filename)
    return KfOutputQaFailure(0, "an input filename is required", 1);
  TFile* file = TFile::Open(filename, "READ");
  if (!file || file->IsZombie())
    return KfOutputQaFailure(file, "cannot open the output file read-only", 1);
  TNamed* status = dynamic_cast<TNamed*>(file->Get("KFRunStatus"));
  if (!status || TString(status->GetTitle()) != "completed")
    return KfOutputQaFailure(file, "KFRunStatus is missing or is not completed", 2);
  TNamed* backend = dynamic_cast<TNamed*>(file->Get("KFParticleBackend"));
  TObjString* configuration = dynamic_cast<TObjString*>(file->Get("KFParticleEffectiveConfiguration"));
  if (!backend || !configuration)
    return KfOutputQaFailure(file, "backend/effective configuration provenance is missing", 2);
  std::cout << "File: " << filename << "\nKFRunStatus: " << status->GetTitle()
            << "\nBackend: " << backend->GetTitle() << std::endl;

  TH1* stages = dynamic_cast<TH1*>(file->Get("hKfStages"));
  TTree* tree = dynamic_cast<TTree*>(file->Get("KfLambdaCandidates"));
  if (!stages || !tree || stages->GetNbinsX() != 14)
    return KfOutputQaFailure(file, "candidate tree or 14-bin stage histogram is missing", 3);
  const char* stageLabels[14] = {"read events", "reconstructed events", "raw tracks", "quality tracks",
      "covariance tracks", "PID tracks", "PID hypotheses", "primary hypotheses",
      "Topo particle slots", "Topo Lambdas", "invalid candidates", "valid raw Lambdas",
      "selected Lambdas", "selected anti-Lambdas"};
  Double_t stageCounts[14];
  std::cout << "Processing stages:" << std::endl;
  Int_t i = 0;
  for (i = 0; i < 14; ++i) {
    stageCounts[i] = stages->GetBinContent(i + 1);
    const TString label(stages->GetXaxis()->GetBinLabel(i + 1));
    if (label != stageLabels[i] || !TMath::Finite(stageCounts[i]) || stageCounts[i] < 0.)
      return KfOutputQaFailure(file, "stage labels/counts are inconsistent", 3);
    std::cout << "  " << i << " " << label << ": " << stageCounts[i] << std::endl;
  }
  if (stageCounts[0] <= 0. || stageCounts[1] <= 0. || stageCounts[1] > stageCounts[0])
    return KfOutputQaFailure(file, "no successfully reconstructed events or invalid event counts", 3);

  TObjString* eventConfiguration = dynamic_cast<TObjString*>(file->Get("KFEventSelectionConfiguration"));
  if (!eventConfiguration)
    return KfOutputQaFailure(file, "current-format KFEventSelectionConfiguration is missing; rerun the mode-aware Maker", 6);
  const TString eventText = eventConfiguration->GetString();
  const TString mode = KfOutputMetadataValue(eventText, "KFEventSelection.mode");
  if (KfOutputMetadataValue(eventText, "KFEventSelection.loaded") != "true" ||
      (mode != "refmult" && mode != "fxtmult"))
    return KfOutputQaFailure(file, "event selection provenance has no valid loaded mode", 6);
  Bool_t legacyLambdaCuts = kFALSE;
  if (!KfOutputEventPolicy(eventText, configuration->GetString(), legacyLambdaCuts))
    return KfOutputQaFailure(file, "unknown or inconsistent KF event-selection policy/flags/profile", 6);
  const TString vertexSource = KfOutputMetadataValue(eventText, "KFEventSelection.vertexSource");
  const Bool_t fullVertexSource = vertexSource == TString("event.vertexByMode.") + mode;
  const Bool_t qaVertexSource = vertexSource == TString("event.qaVertexByMode.") + mode;
  if (!fullVertexSource && !(legacyLambdaCuts && qaVertexSource))
    return KfOutputQaFailure(file, "event vertex provenance is inconsistent with its selection policy", 6);
  Double_t centerX = 0., centerY = 0., radius = 0., minimumVz = 0., maximumVz = 0.;
  if (!KfOutputMetadataDouble(eventText, "KFEventSelection.vtxCenterX", centerX) ||
      !KfOutputMetadataDouble(eventText, "KFEventSelection.vtxCenterY", centerY))
    return KfOutputQaFailure(file, "effective vertex QA center metadata are invalid", 6);
  // Imp5 has only a diagnostic center. Old Imp5 files can still contain unused
  // radius/Vz fields, but neither old nor sparse Imp5 metadata claim those cuts.
  if (!legacyLambdaCuts &&
      (!KfOutputMetadataDouble(eventText, "KFEventSelection.maxVr", radius) ||
       !KfOutputMetadataDouble(eventText, "KFEventSelection.minVz", minimumVz) ||
       !KfOutputMetadataDouble(eventText, "KFEventSelection.maxVz", maximumVz) ||
       radius <= 0. || minimumVz >= maximumVz))
    return KfOutputQaFailure(file, "effective vertex radius/Vz cut metadata are invalid", 6);
  std::cout << "Effective event selection (stored by the Maker):\n" << eventText << std::endl;

  TH1* eventSelection = dynamic_cast<TH1*>(file->Get("hKfEventSelection"));
  const char* eventLabels[13] = {"read events", "bad run", "invalid vertex", "Vz", "Vr",
      "refMult", "VPD difference", "track count", "pileup", "centrality invalid",
      "centrality bin", "KF error", "reconstructed"};
  if (!eventSelection || eventSelection->GetDimension() != 1 || eventSelection->GetNbinsX() != 13)
    return KfOutputQaFailure(file, "current-format 13-bin hKfEventSelection is missing", 6);
  Double_t eventCounts[13];
  Double_t terminalTotal = 0., afterEvent = 0.;
  std::cout << "Event selection (first bin is total; other bins are mutually exclusive outcomes):" << std::endl;
  for (i = 0; i < 13; ++i) {
    eventCounts[i] = eventSelection->GetBinContent(i + 1);
    if (TString(eventSelection->GetXaxis()->GetBinLabel(i + 1)) != eventLabels[i] ||
        !TMath::Finite(eventCounts[i]) || eventCounts[i] < 0. ||
        eventCounts[i] != TMath::Floor(eventCounts[i]))
      return KfOutputQaFailure(file, "event-selection labels/counts are inconsistent", 6);
    if (i > 0) terminalTotal += eventCounts[i];
    if (i >= 8) afterEvent += eventCounts[i];
    std::cout << "  " << eventLabels[i] << ": " << eventCounts[i] << std::endl;
  }
  if (eventSelection->GetBinContent(0) != 0. || eventSelection->GetBinContent(14) != 0. ||
      terminalTotal != eventCounts[0] || eventCounts[0] != stageCounts[0] ||
      eventCounts[12] != stageCounts[1] || eventCounts[11] != 0.)
    return KfOutputQaFailure(file, "event outcomes do not partition read events, disagree with stages, or contain KF errors in a completed run", 6);

  if (legacyLambdaCuts && (eventCounts[3] != 0. || eventCounts[4] != 0. ||
      eventCounts[5] != 0. || eventCounts[6] != 0.))
    return KfOutputQaFailure(file, "Imp5 policy contains rejections from disabled Vz/Vr/RefMult/VPD cuts", 6);

  TH2* vertexBefore = dynamic_cast<TH2*>(file->Get("hKfVertexXYBefore"));
  TH2* vertexAfter = dynamic_cast<TH2*>(file->Get("hKfVertexXYAfterEvent"));
  TH1* radiusBefore = dynamic_cast<TH1*>(file->Get("hKfVertexRadiusBefore"));
  TH1* radiusAfter = dynamic_cast<TH1*>(file->Get("hKfVertexRadiusAfterEvent"));
  TH1* vzAfter = dynamic_cast<TH1*>(file->Get("hKfVertexZAfterEvent"));
  if (!vertexBefore || !vertexAfter || !radiusBefore || !radiusAfter || !vzAfter)
    return KfOutputQaFailure(file, "required before/after-event vertex QA histograms are missing", 6);
  const Double_t beforeEntries = vertexBefore->GetEntries();
  if (!TMath::Finite(beforeEntries) || beforeEntries < afterEvent || beforeEntries > eventCounts[0] ||
      !KfOutputVertexHistogram(vertexBefore, 2, beforeEntries) ||
      !KfOutputVertexHistogram(radiusBefore, 1, beforeEntries) ||
      !KfOutputVertexHistogram(vertexAfter, 2, afterEvent) ||
      !KfOutputVertexHistogram(radiusAfter, 1, afterEvent) ||
      !KfOutputVertexHistogram(vzAfter, 1, afterEvent))
    return KfOutputQaFailure(file, "vertex histogram entries/cells/statistics disagree with event selection", 6);
  if (radiusBefore->GetMean() < 0. || radiusAfter->GetMean() < 0. ||
      (!legacyLambdaCuts && (!KfOutputHistogramWithin(radiusAfter, 0., radius) ||
                            !KfOutputHistogramWithin(vzAfter, minimumVz, maximumVz))))
    return KfOutputQaFailure(file, "vertex radius/Vz QA is outside the stored effective cuts or has a negative radius", 6);
  std::cout << "Vertex QA: policy=" << (legacyLambdaCuts ? "legacy_lambda_imp5" : "mode_vertex")
            << "; vertex cuts=" << (legacyLambdaCuts ? "disabled (diagnostic center only)" : "applied")
            << "; mode=" << mode << "; center=(" << centerX << "," << centerY << ") cm";
  if (!legacyLambdaCuts)
    std::cout << "; maxVr=" << radius << " cm; Vz=[" << minimumVz << "," << maximumVz << "] cm";
  std::cout << "\n  finite-XY entries before bad-run/event cuts=" << beforeEntries
            << "; after bad-run/event cuts, before centrality=" << afterEvent
            << "\n  after-event mean (x,y)=(" << vertexAfter->GetMean(1) << ","
            << vertexAfter->GetMean(2) << ") cm; mean effective radius=" << radiusAfter->GetMean()
            << " cm; mean Vz=" << vzAfter->GetMean() << " cm" << std::endl;

  const char* integerNames[5] = {"pdg", "protonId", "pionId", "protonIndex", "pionIndex"};
  TLeaf* integerLeaves[5];
  for (i = 0; i < 5; ++i) {
    integerLeaves[i] = tree->GetLeaf(integerNames[i]);
    if (!integerLeaves[i])
      return KfOutputQaFailure(file, "required PDG/daughter identity leaf is missing", 4);
  }
  TLeaf* selectedLeaf = tree->GetLeaf("selected");
  if (!selectedLeaf)
    return KfOutputQaFailure(file, "selected leaf is missing", 4);
  const char* floatNames[26] = {"mass", "massError", "chi2Ndf", "topoChi2Ndf",
      "x", "y", "z", "px", "py", "pz", "magneticField",
      "daughterDistance", "distanceToPv", "decayLength", "decayLengthError",
      "decayLengthSignificance", "vertexLineLength", "vertexLineLengthError",
      "vertexLineLengthSignificance", "cosPointing", "protonPidPull", "pionPidPull",
      "protonTofM2", "pionTofM2", "protonTofPull", "pionTofPull"};
  TLeaf* floatLeaves[26];
  for (i = 0; i < 26; ++i) {
    floatLeaves[i] = tree->GetLeaf(floatNames[i]);
    if (!floatLeaves[i]) {
      std::cerr << "Missing leaf: " << floatNames[i] << std::endl;
      return KfOutputQaFailure(file, "required fit/PID leaf is missing", 4);
    }
  }

  const Long64_t entries = tree->GetEntries();
  if (entries < 0 || stageCounts[11] != Double_t(entries))
    return KfOutputQaFailure(file, "tree count differs from valid raw Lambda stage", 4);
  Long64_t rawCounts[2] = {0, 0}, selectedCounts[2] = {0, 0};
  Long64_t offPdgCounts[2] = {0, 0};
  Long64_t windowCounts[2][3] = {{0, 0, 0}, {0, 0, 0}};
  Double_t minimumMass[2] = {DBL_MAX, DBL_MAX};
  Double_t maximumMass[2] = {-DBL_MAX, -DBL_MAX};
  Double_t sumMass[2] = {0., 0.}, sumMass2[2] = {0., 0.};
  Double_t minimumError[2] = {DBL_MAX, DBL_MAX}, maximumError[2] = {0., 0.};
  const Double_t pdgMass = 1.115683;
  // Equal-width diagnostic windows, NOT a signal definition or a peak test.
  const Double_t windowLow[3] = {1.086, 1.110, 1.134};
  const Double_t windowHigh[3] = {1.098, 1.122, 1.146};
  for (Long64_t entry = 0; entry < entries; ++entry) {
    if (tree->GetEntry(entry) <= 0)
      return KfOutputQaFailure(file, "failed reading a candidate tree entry", 4);
    const Int_t pdg = Int_t(integerLeaves[0]->GetValue());
    if (pdg != 3122 && pdg != -3122)
      return KfOutputQaFailure(file, "candidate has a non-Lambda PDG", 4);
    const Int_t signIndex = pdg > 0 ? 0 : 1;
    const Int_t protonId = Int_t(integerLeaves[1]->GetValue());
    const Int_t pionId = Int_t(integerLeaves[2]->GetValue());
    const Int_t protonIndex = Int_t(integerLeaves[3]->GetValue());
    const Int_t pionIndex = Int_t(integerLeaves[4]->GetValue());
    if (protonId < 0 || pionId < 0 || protonId == pionId ||
        protonIndex < 0 || pionIndex < 0 || protonIndex == pionIndex)
      return KfOutputQaFailure(file, "invalid or repeated daughter IDs/indices", 4);
    for (i = 0; i < 26; ++i) {
      if (!TMath::Finite(floatLeaves[i]->GetValue())) {
        std::cerr << "Non-finite entry=" << entry << " field=" << floatNames[i] << std::endl;
        return KfOutputQaFailure(file, "candidate contains a non-finite fit/PID value", 4);
      }
    }
    const Double_t mass = floatLeaves[0]->GetValue();
    const Double_t massError = floatLeaves[1]->GetValue();
    if (mass <= 0. || massError <= 0. || floatLeaves[2]->GetValue() < 0. ||
        floatLeaves[3]->GetValue() < 0. || floatLeaves[14]->GetValue() <= 0. ||
        floatLeaves[17]->GetValue() <= 0.)
      return KfOutputQaFailure(file, "invalid mass/uncertainty or fit quality", 4);
    const Double_t selected = selectedLeaf->GetValue();
    if (selected != 0. && selected != 1.)
      return KfOutputQaFailure(file, "selected flag is not boolean", 4);
    ++rawCounts[signIndex];
    if (selected != 0.) ++selectedCounts[signIndex];
    if (TMath::Abs(mass - pdgMass) > 0.0001) ++offPdgCounts[signIndex];
    if (mass < minimumMass[signIndex]) minimumMass[signIndex] = mass;
    if (mass > maximumMass[signIndex]) maximumMass[signIndex] = mass;
    if (massError < minimumError[signIndex]) minimumError[signIndex] = massError;
    if (massError > maximumError[signIndex]) maximumError[signIndex] = massError;
    sumMass[signIndex] += mass;
    sumMass2[signIndex] += mass * mass;
    for (i = 0; i < 3; ++i)
      if (mass >= windowLow[i] && mass < windowHigh[i]) ++windowCounts[signIndex][i];
  }
  if (stageCounts[12] != Double_t(selectedCounts[0]) ||
      stageCounts[13] != Double_t(selectedCounts[1]))
    return KfOutputQaFailure(file, "selected tree counts differ from stage counters", 5);
  if (requireCandidates && (entries == 0 || selectedCounts[0] + selectedCounts[1] == 0))
    return KfOutputQaFailure(file, "at least one selected Lambda/anti-Lambda was required", 5);

  const char* signLabels[2] = {"Lambda (+3122)", "anti-Lambda (-3122)"};
  const char* rawHistogramNames[2] = {"hKfLambdaMassRaw", "hKfAntiLambdaMassRaw"};
  const char* selectedHistogramNames[2] = {"hKfLambdaMassSelected", "hKfAntiLambdaMassSelected"};
  std::cout << "KfLambdaCandidates entries: " << entries << std::endl;
  for (Int_t species = 0; species < 2; ++species) {
    TH1* rawHistogram = dynamic_cast<TH1*>(file->Get(rawHistogramNames[species]));
    TH1* selectedHistogram = dynamic_cast<TH1*>(file->Get(selectedHistogramNames[species]));
    if (!rawHistogram || !selectedHistogram)
      return KfOutputQaFailure(file, "signed mass histograms are missing", 5);
    if (rawHistogram->GetEntries() != Double_t(rawCounts[species]) ||
        selectedHistogram->GetEntries() != Double_t(selectedCounts[species]))
      return KfOutputQaFailure(file, "signed histogram and tree counts disagree", 5);
    std::cout << signLabels[species] << ": raw=" << rawCounts[species]
              << ", selected=" << selectedCounts[species] << std::endl;
    if (rawCounts[species] > 0) {
      const Double_t mean = sumMass[species] / Double_t(rawCounts[species]);
      const Double_t variance = sumMass2[species] / Double_t(rawCounts[species]) - mean * mean;
      std::cout << "  raw mass range [GeV/c2]: " << minimumMass[species] << " .. "
                << maximumMass[species] << "; mean=" << mean
                << ", RMS=" << TMath::Sqrt(TMath::Max(0., variance))
                << "; mass-error range=" << minimumError[species] << " .. "
                << maximumError[species] << std::endl;
      const Int_t maximumBin = rawHistogram->GetMaximumBin();
      std::cout << "  highest-populated raw histogram bin: ["
                << rawHistogram->GetXaxis()->GetBinLowEdge(maximumBin) << ", "
                << rawHistogram->GetXaxis()->GetBinUpEdge(maximumBin) << "), count="
                << rawHistogram->GetBinContent(maximumBin) << std::endl;
    }
    std::cout << "  raw diagnostic counts: left [1.086,1.098)=" << windowCounts[species][0]
              << ", central [1.110,1.122)=" << windowCounts[species][1]
              << ", right [1.134,1.146)=" << windowCounts[species][2]
              << "; |mass-PDG|>0.0001 GeV/c2: " << offPdgCounts[species] << std::endl;
  }
  std::cout << "Diagnostic windows/highest bins are descriptive only: no peak significance, "
               "background fit, calibrated efficiency, or physics-yield validation is asserted.\n"
               "Raw mass means no parent mass constraint but AFTER upstream Topo selection. "
               "An off-PDG spread is a diagnostic; code review and the synthetic copy-constraint "
               "test separately check parent-mass handling.\n"
               "PASS KF output structure, mode/event/vertex QA, counters/finite values and signed mass summaries."
            << std::endl;
  file->Close();
  delete file;
  return 0;
}

void check_kfparticle_output(const char* filename, Bool_t requireCandidates = kTRUE)
{
  const Int_t result = KfOutputQaInspect(filename, requireCandidates);
  gSystem->Exit(result);
}
