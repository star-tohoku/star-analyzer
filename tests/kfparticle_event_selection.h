#ifndef STAR_ANALYZER_KF_EVENT_SELECTION_TEST_H
#define STAR_ANALYZER_KF_EVENT_SELECTION_TEST_H

// Compiled Pico-adapter fixture: temporary synthetic data only.
#include "KfEventSelection.h"
#include "cuts/EventCutConfig.h"
#include "StPicoEvent/StPicoEvent.h"
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <sys/stat.h>
#include <unistd.h>

namespace kf_event_selection_test {
inline void Require(bool condition, const std::string& message) {
  if (!condition) throw std::runtime_error("KF event fixture: " + message);
}

class RestoreEventCuts {
public:
  explicit RestoreEventCuts(EventCutConfig& cuts) : c(cuts),
      minVz(c.minVz), maxVz(c.maxVz), maxVr(c.maxVr), x(c.vtxCenterX), y(c.vtxCenterY),
      minRef(c.minRefMult), maxRef(c.maxRefMult), diff(c.maxVzDiff), vpd(c.maxAbsVzVpd), n(c.maxNTr) {}
  ~RestoreEventCuts() {
    c.minVz = minVz; c.maxVz = maxVz; c.maxVr = maxVr;
    c.vtxCenterX = x; c.vtxCenterY = y;
    c.minRefMult = minRef; c.maxRefMult = maxRef;
    c.maxVzDiff = diff; c.maxAbsVzVpd = vpd; c.maxNTr = n;
  }
private:
  EventCutConfig& c;
  double minVz, maxVz, maxVr, x, y, minRef, maxRef, diff, vpd;
  int n;
};

class ConfigTree {
public:
  ConfigTree() {
    char pattern[] = "/tmp/star_analyzer_kf_event_XXXXXX";
    char* result = mkdtemp(pattern);
    Require(result != 0, "cannot create temporary configuration directory");
    root = result;
    Require(mkdir((root + "/config").c_str(), 0700) == 0, "cannot create config directory");
    Require(mkdir((root + "/config/mainconf").c_str(), 0700) == 0, "cannot create mainconf directory");
    main = root + "/config/mainconf/main.yaml";
    analysis = root + "/config/analysis.yaml";
    event = root + "/config/event.yaml";
  }
  ~ConfigTree() {
    unlink(main.c_str()); unlink(analysis.c_str()); unlink(event.c_str());
    rmdir((root + "/config/mainconf").c_str());
    rmdir((root + "/config").c_str()); rmdir(root.c_str());
  }
  void Write(const std::string& path, const std::string& contents) {
    std::ofstream output(path.c_str()); output << contents; output.close();
    Require(static_cast<bool>(output), "cannot write temporary configuration: " + path);
  }
  void Set(const std::string& mode, const std::string& profiles, bool absolute = false) {
    Write(main, "analysis: " + (absolute ? analysis : "analysis.yaml") + "\nevent: " +
          (absolute ? event : "event.yaml") + "\n");
    Write(analysis, "analysis:\n  mode: " + mode + "\n"); Write(event, profiles);
  }
  std::string root, main, analysis, event;
};

inline std::string Profiles() {
  return "minVz: -100\nmaxVz: 100\nmaxVr: 99\nvtxCenterX: 77\nvtxCenterY: 88\n"
         "vertexByMode:\n"
         "  refmult:\n    vzRange: [-100, 100]\n    center: [0, 0]\n    radius: 2\n"
         "  fxtmult:\n    vzRange: [198, 202]\n    center: [-0.4, -2.0]\n    radius: 2\n";
}

inline std::string QaProfiles() {
  return "maxNTr: 20\nqaVertexByMode:\n"
         "  refmult:\n    center: [0, 0]\n"
         "  fxtmult:\n    center: [-0.4, -2.0]\n";
}

inline std::string Replace(const std::string& text, const std::string& from, const std::string& to) {
  std::string result(text);
  const size_t position = result.find(from);
  Require(position != std::string::npos, "missing replacement fixture token");
  result.replace(position, from.size(), to); return result;
}

inline void Invalid(ConfigTree& tree, EventCutConfig& base, const std::string& mode,
                    const std::string& event, const std::string& message) {
  tree.Set(mode, event);
  KfEventSelection selection;
  std::ostringstream errors;
  Require(!selection.Load(tree.main.c_str(), base, "fxtmult", false, errors), message);
  Require(!errors.str().empty(), "configuration rejection had no explanation");
}

inline void InvalidSparse(ConfigTree& tree, EventCutConfig& base, const std::string& mode,
                          const std::string& event, const std::string& message) {
  tree.Set(mode, event);
  KfEventSelection selection;
  std::ostringstream errors;
  Require(!selection.Load(tree.main.c_str(), base, "fxtmult", false, errors, true), message);
  Require(!errors.str().empty(), "sparse configuration rejection had no explanation");
}

inline void CheckLegacyDump(const std::string& dump, const std::string& vertexSource) {
  Require(dump.find("KFEventSelection.policy: legacy_lambda_imp5") != std::string::npos &&
          dump.find("KFEventSelection.vertexCutsApplied: false") != std::string::npos &&
          dump.find("KFEventSelection.refMultCutsApplied: false") != std::string::npos &&
          dump.find("KFEventSelection.vpdCutsApplied: false") != std::string::npos &&
          dump.find("KFEventSelection.vertexSource: " + vertexSource) != std::string::npos &&
          dump.find("KFEventSelection.vtxCenterX:") != std::string::npos &&
          dump.find("KFEventSelection.vtxCenterY:") != std::string::npos &&
          dump.find("KFEventSelection.maxNTr:") != std::string::npos,
          "legacy dump lost its policy, diagnostic center, source or active maxNTr");
  const char* inactive[] = {"minVz", "maxVz", "maxVr", "minRefMult", "maxRefMult",
                           "maxVzDiff", "maxAbsVzVpd"};
  for (size_t i = 0; i < sizeof(inactive) / sizeof(inactive[0]); ++i)
    Require(dump.find(std::string("KFEventSelection.") + inactive[i] + ":") == std::string::npos,
            std::string("legacy dump advertises an unused event cut: ") + inactive[i]);
}

inline void Run() {
  EventCutConfig& base = EventCutConfig::GetInstance();
  RestoreEventCuts restore(base);
  base.SetDefaults();
  // Conflicting legacy vertex values must neither win nor mutate.
  base.minVz = -10.; base.maxVz = 10.; base.maxVr = 0.1;
  base.vtxCenterX = 77.; base.vtxCenterY = 88.;
  base.minRefMult = 2.; base.maxRefMult = 10.; base.maxNTr = 20;
  ConfigTree tree;
  tree.Set("FXTMULT", Profiles());
  KfEventSelection selection;
  std::ostringstream errors;
  Require(selection.Load(tree.main.c_str(), base, "fxtmult", true, errors), errors.str());
  Require(selection.Mode() == "fxtmult", "mode is not canonicalized");
  Require(selection.ComputeVr(-0.4, -2.) == 0., "fixed-target center was not applied");
  Require(base.minVz == -10. && base.maxVr == 0.1 && base.vtxCenterX == 77. && base.vtxCenterY == 88.,
          "profile selection mutated shared EventCutConfig");
  StPicoEvent event;
  event.setPrimaryVertexPosition(-0.23f, -2.28f, 200.f);
  event.setRefMultPos(3); event.setRefMultNeg(2); event.setVzVpd(-999.f);
  Require(selection.Check(event, 20) == KfEventSelection::kAccepted, "fixed-target vertex rejected");
  event.setPrimaryVertexPosition(0.f, 0.f, 0.f);
  Require(selection.Check(event, 20) == KfEventSelection::kVz, "collider vertex accepted in fixed-target mode");
  event.setPrimaryVertexPosition(0.f, 0.f, 200.f);
  Require(selection.Check(event, 20) == KfEventSelection::kVr, "fixed-target XY center ignored");
  event.setPrimaryVertexPosition(-0.4f, -2.f, 198.f);
  Require(selection.Check(event, 20) == KfEventSelection::kAccepted, "inclusive lower Vz boundary rejected");
  event.setPrimaryVertexPosition(-0.4f, -2.f, 202.f);
  Require(selection.Check(event, 20) == KfEventSelection::kAccepted, "inclusive upper Vz boundary rejected");
  base.maxNTr = 1;
  Require(selection.Check(event, 20) == KfEventSelection::kAccepted, "base mutation changed loaded cuts");
  base.maxNTr = 20;

  tree.Set("ReFmUlT", Profiles(), true);
  Require(selection.Load(tree.main.c_str(), base, "REFMULT", true, errors), errors.str());
  Require(selection.Mode() == "refmult" && selection.ComputeVr(0., 0.) == 0., "collider/absolute refs failed");
  event.setPrimaryVertexPosition(0.f, 0.f, 0.f);
  Require(selection.Check(event, 20) == KfEventSelection::kAccepted, "collider vertex rejected");
  event.setPrimaryVertexPosition(-0.23f, -2.28f, 200.f);
  Require(selection.Check(event, 20) == KfEventSelection::kVz, "fixed-target vertex accepted in collider mode");
  event.setPrimaryVertexPosition(2.f, 0.f, 100.f);
  Require(selection.Check(event, 20) == KfEventSelection::kAccepted, "inclusive collider Vz/Vr boundaries changed");
  event.setPrimaryVertexPosition(2.01f, 0.f, 0.f);
  Require(selection.Check(event, 20) == KfEventSelection::kVr, "outside-radius vertex accepted");
  event.setPrimaryVertexPosition(std::numeric_limits<float>::quiet_NaN(), 0.f, 0.f);
  Require(selection.Check(event, 20) == KfEventSelection::kInvalidVertex, "nonfinite vertex accepted");
  event.setPrimaryVertexPosition(0.f, 0.f, 0.f);
  event.setRefMultPos(0); event.setRefMultNeg(0);
  Require(selection.Check(event, 20) == KfEventSelection::kRefMult, "low RefMult accepted");
  event.setRefMultPos(9); event.setRefMultNeg(2);
  Require(selection.Check(event, 20) == KfEventSelection::kRefMult, "high RefMult accepted");
  event.setRefMultPos(3); event.setRefMultNeg(2); event.setVzVpd(4.f);
  Require(selection.Check(event, 20) == KfEventSelection::kVpd, "valid-VPD difference cut ignored");
  event.setVzVpd(3.f);
  Require(selection.Check(event, 20) == KfEventSelection::kAccepted, "inclusive VPD boundary changed");
  event.setVzVpd(200.f);
  Require(selection.Check(event, 20) == KfEventSelection::kAccepted, "invalid-VPD threshold changed");
  event.setVzVpd(std::numeric_limits<float>::quiet_NaN());
  Require(selection.Check(event, 21) == KfEventSelection::kTrackCount, "track-count cut ignored");
  Require(selection.Check(event, 20) == KfEventSelection::kAccepted, "missing VPD stopped being optional");
  std::ostringstream effective; selection.Dump(effective);
  Require(effective.str().find("KFEventSelection.mode: refmult") != std::string::npos &&
          effective.str().find(tree.analysis) != std::string::npos &&
          effective.str().find("KFEventSelection.maxNTr: 20") != std::string::npos,
          "effective config lacks mode, paths or values");
  Require(!selection.Load(tree.main.c_str(), base, "fxtmult", true, errors), "centrality mismatch accepted");
  Require(selection.Check(event, 20) != KfEventSelection::kAccepted && selection.Mode().empty() &&
          !std::isfinite(selection.ComputeVr(0., 0.)), "failed reload retained valid state");
  Require(selection.Load(tree.main.c_str(), base, "fxtmult", false, errors), "disabled centrality blocked mode");

  // Explicit Imp5 event compatibility must not silently change the default.
  // Keep the selected mode/profile for QA, while the opt-in selection itself
  // follows the old Lambda maxNTr-only cut plus required finite-PV safety.
  tree.Set("fxtmult", Profiles());
  Require(selection.Load(tree.main.c_str(), base, "fxtmult", true, errors, true), errors.str());
  Require(selection.Mode() == "fxtmult" && selection.ComputeVr(-0.4, -2.) == 0.,
          "Imp5 discarded mode-aware diagnostic coordinates");
  event.setPrimaryVertexPosition(8.f, 9.f, 250.f);
  event.setRefMultPos(0); event.setRefMultNeg(0); event.setVzVpd(0.f);
  Require(selection.Check(event, 20) == KfEventSelection::kAccepted,
          "Imp5 unexpectedly applied Vz/Vr/RefMult/VPD cuts");
  Require(selection.Check(event, 21) == KfEventSelection::kTrackCount,
          "Imp5 failed to enforce the legacy maxNTr cut");
  event.setPrimaryVertexPosition(std::numeric_limits<float>::quiet_NaN(), 0.f, 200.f);
  Require(selection.Check(event, 20) == KfEventSelection::kInvalidVertex,
          "Imp5 accepted nonfinite transverse vertex");
  event.setPrimaryVertexPosition(0.f, 0.f, std::numeric_limits<float>::quiet_NaN());
  Require(selection.Check(event, 20) == KfEventSelection::kInvalidVertex,
          "Imp5 accepted nonfinite longitudinal vertex");
  event.setPrimaryVertexPosition(8.f, 9.f, 250.f);
  base.maxNTr = 0;
  Require(selection.Load(tree.main.c_str(), base, "fxtmult", true, errors, true), errors.str());
  Require(selection.Check(event, 100000) == KfEventSelection::kAccepted,
          "Imp5 maxNTr=0 did not disable the legacy track-count cut");
  base.maxNTr = 20;
  Require(selection.Load(tree.main.c_str(), base, "fxtmult", true, errors, true), errors.str());
  std::ostringstream imp5Effective; selection.Dump(imp5Effective);
  Require(imp5Effective.str().find("KFEventSelection.policy: legacy_lambda_imp5") != std::string::npos &&
          imp5Effective.str().find("KFEventSelection.vertexCutsApplied: false") != std::string::npos &&
          imp5Effective.str().find("KFEventSelection.refMultCutsApplied: false") != std::string::npos &&
          imp5Effective.str().find("KFEventSelection.vpdCutsApplied: false") != std::string::npos,
          "Imp5 provenance incorrectly advertises mode vertex/refMult/VPD cuts");
  CheckLegacyDump(imp5Effective.str(), "event.vertexByMode.fxtmult");
  Require(base.minVz == -10. && base.maxVr == 0.1 && base.vtxCenterX == 77. && base.vtxCenterY == 88.,
          "Imp5 modified the shared EventCutConfig");
  Require(!selection.Load(tree.main.c_str(), base, "refmult", true, errors, true),
          "Imp5 bypassed centrality-mode consistency validation");
  Require(selection.Check(event, 20) != KfEventSelection::kAccepted,
          "failed Imp5 reload retained accepting state");
  Require(selection.Load(tree.main.c_str(), base, "fxtmult", true, errors), errors.str());
  Require(selection.Check(event, 20) == KfEventSelection::kVz,
          "omitting the optional Imp5 flag retained its old event selection");
  std::ostringstream restoredEffective; selection.Dump(restoredEffective);
  Require(restoredEffective.str().find("KFEventSelection.policy: mode_vertex") != std::string::npos &&
          restoredEffective.str().find("KFEventSelection.vertexCutsApplied: true") != std::string::npos,
          "default mode event selection was not restored in provenance");
  event.setPrimaryVertexPosition(-0.23f, -2.28f, 200.f);
  event.setRefMultPos(3); event.setRefMultNeg(2); event.setVzVpd(-999.f);
  Require(selection.Check(event, 20) == KfEventSelection::kAccepted,
          "valid fixed-target event failed after restoring default cuts");

  // New sparse schema stores only active maxNTr plus a diagnostic mode center.
  // It must not acquire numeric cuts from omitted fields or mutable base state.
  tree.Set("FXTMULT", QaProfiles());
  base.maxNTr = 7;
  Require(selection.Load(tree.main.c_str(), base, "fxtmult", true, errors, true), errors.str());
  Require(selection.Mode() == "fxtmult" && selection.ComputeVr(-0.4, -2.) == 0.,
          "sparse fixed-target QA center was not selected");
  Require(base.maxNTr == 7, "sparse config mutated the shared base track count");
  event.setPrimaryVertexPosition(8.f, 9.f, 250.f);
  event.setRefMultPos(0); event.setRefMultNeg(0); event.setVzVpd(0.f);
  Require(selection.Check(event, 20) == KfEventSelection::kAccepted &&
          selection.Check(event, 21) == KfEventSelection::kTrackCount,
          "sparse maxNTr was not read from YAML or unused vertex/refMult/VPD cuts became active");
  event.setPrimaryVertexPosition(0.f, std::numeric_limits<float>::infinity(), 200.f);
  Require(selection.Check(event, 20) == KfEventSelection::kInvalidVertex,
          "sparse event policy bypassed finite vertex safety");
  std::ostringstream sparseEffective; selection.Dump(sparseEffective);
  CheckLegacyDump(sparseEffective.str(), "event.qaVertexByMode.fxtmult");
  Require(sparseEffective.str().find("KFEventSelection.maxNTr: 20") != std::string::npos,
          "sparse dump reported stale maxNTr from the base singleton");
  Require(!selection.Load(tree.main.c_str(), base, "fxtmult", true, errors),
          "standard mode accepted a diagnostic-only QA vertex profile");
  Require(selection.Mode().empty() && !std::isfinite(selection.ComputeVr(0., 0.)),
          "failed standard reload retained the preceding sparse profile");

  tree.Set("ReFmUlT", QaProfiles(), true);
  Require(selection.Load(tree.main.c_str(), base, "refmult", true, errors, true), errors.str());
  Require(selection.Mode() == "refmult" && selection.ComputeVr(0., 0.) == 0.,
          "sparse collider QA center or absolute paths failed");
  tree.Set("fxtmult", Replace(QaProfiles(), "maxNTr: 20", "maxNTr: 0"));
  Require(selection.Load(tree.main.c_str(), base, "fxtmult", true, errors, true), errors.str());
  event.setPrimaryVertexPosition(8.f, 9.f, 250.f);
  Require(selection.Check(event, 100000) == KfEventSelection::kAccepted,
          "sparse maxNTr=0 did not disable track-count selection");
  base.maxNTr = 20;

  InvalidSparse(tree, base, "fxmult", QaProfiles(), "sparse unknown analysis mode accepted");
  InvalidSparse(tree, base, "fxtmult", Replace(QaProfiles(), "maxNTr: 20\n", ""),
                "sparse missing maxNTr accepted");
  const char* badSparseNTr[] = {"-1", "20.5", ".inf", "null", "true"};
  for (size_t i = 0; i < sizeof(badSparseNTr) / sizeof(badSparseNTr[0]); ++i)
    InvalidSparse(tree, base, "fxtmult",
                  Replace(QaProfiles(), "maxNTr: 20", std::string("maxNTr: ") + badSparseNTr[i]),
                  "sparse invalid/noninteger maxNTr accepted");
  InvalidSparse(tree, base, "fxtmult", QaProfiles() + "maxNTr: 10\n",
                "sparse duplicate root maxNTr accepted");
  InvalidSparse(tree, base, "fxtmult", QaProfiles() + "minVz: 198\n",
                "sparse schema accepted an inactive top-level vertex cut");
  InvalidSparse(tree, base, "fxtmult",
                QaProfiles() + "vertexByMode:\n  fxtmult:\n    center: [-0.4, -2.0]\n    vzRange: [198, 202]\n    radius: 2\n",
                "ambiguous qaVertexByMode and vertexByMode were accepted");
  InvalidSparse(tree, base, "fxtmult", "maxNTr: 20\nqaVertexByMode:\n  refmult:\n    center: [0, 0]\n",
                "sparse selected mode profile missing");
  InvalidSparse(tree, base, "fxtmult", Replace(QaProfiles(), "  refmult:", "  refmull:"),
                "sparse unrecognized unselected mode accepted");
  InvalidSparse(tree, base, "fxtmult", QaProfiles() + "  fxtmult:\n    center: [0, 0]\n",
                "sparse duplicate selected mode accepted");
  InvalidSparse(tree, base, "fxtmult", Replace(QaProfiles(), "center: [-0.4, -2.0]",
                "center: [-0.4, -2.0]\n    center: [0, 0]"), "sparse duplicate center accepted");
  InvalidSparse(tree, base, "fxtmult", Replace(QaProfiles(), "center: [-0.4, -2.0]",
                "center: [-0.4, -2.0]\n    radius: 2"), "sparse inactive radius accepted");
  InvalidSparse(tree, base, "fxtmult", Replace(QaProfiles(), "center: [-0.4, -2.0]",
                "center: [-0.4, -2.0]\n    typo: 2"), "sparse unrecognized profile key accepted");
  InvalidSparse(tree, base, "fxtmult", Replace(QaProfiles(), "center: [-0.4, -2.0]", "center: [-0.4]"),
                "sparse wrong center dimension accepted");
  InvalidSparse(tree, base, "fxtmult", Replace(QaProfiles(), "center: [-0.4, -2.0]", "center: [.nan, -2.0]"),
                "sparse nonfinite selected center accepted");
  InvalidSparse(tree, base, "fxtmult", Replace(QaProfiles(), "center: [0, 0]", "center: [0, .inf]"),
                "sparse invalid unselected profile was not validated");

  // Archived full profiles and standard selection remain usable after sparse
  // loads and failures, with no changes to their original geometry semantics.
  tree.Set("fxtmult", Profiles());
  Require(selection.Load(tree.main.c_str(), base, "fxtmult", true, errors, true), errors.str());
  Require(selection.Check(event, 20) == KfEventSelection::kAccepted,
          "archived full-profile Imp5 compatibility was lost after sparse reloads");
  Require(selection.Load(tree.main.c_str(), base, "fxtmult", true, errors), errors.str());
  Require(selection.Check(event, 20) == KfEventSelection::kVz,
          "standard vertex cuts were not restored after sparse Imp5 selection");

  Invalid(tree, base, "fxmult", Profiles(), "unknown mode accepted");
  Invalid(tree, base, "null", Profiles(), "null mode accepted");
  Invalid(tree, base, "[]", Profiles(), "sequence mode accepted");
  Invalid(tree, base, "fxtmult", "maxVr: 2\n", "missing profiles accepted");
  Invalid(tree, base, "fxtmult", "vertexByMode:\n  refmult:\n    vzRange: [-1, 1]\n    center: [0, 0]\n    radius: 2\n",
          "missing selected profile accepted");
  Invalid(tree, base, "fxtmult", Replace(Profiles(), "  fxtmult:", "  fxmult:"), "unknown profile mode accepted");
  Invalid(tree, base, "fxtmult", Replace(Profiles(), "center: [-0.4, -2.0]", "center: [-0.4, -2.0]\n    typo: 1"),
          "unknown profile key accepted");
  Invalid(tree, base, "fxtmult", Replace(Profiles(), "center: [-0.4, -2.0]", "center: [-0.4, -2.0]\n    center: [0, 0]"),
          "duplicate profile key accepted");
  Invalid(tree, base, "fxtmult", Replace(Profiles(), "center: [-0.4, -2.0]", "center: [.nan, -2.0]"),
          "nonfinite center accepted");
  Invalid(tree, base, "fxtmult", Replace(Profiles(), "vzRange: [198, 202]", "vzRange: [198, .inf]"),
          "nonfinite Vz accepted");
  Invalid(tree, base, "fxtmult", Replace(Profiles(), "center: [-0.4, -2.0]", "center: [-0.4]"),
          "wrong center dimension accepted");
  Invalid(tree, base, "fxtmult", Replace(Profiles(), "vzRange: [198, 202]", "vzRange: [198, 202, 203]"),
          "wrong Vz dimension accepted");
  Invalid(tree, base, "fxtmult", Replace(Profiles(), "vzRange: [198, 202]", "vzRange: [202, 198]"),
          "reversed Vz interval accepted");
  const char* badRadius[] = {"0", "-2", ".inf", "nope"};
  for (unsigned int i = 0; i < sizeof(badRadius) / sizeof(badRadius[0]); ++i)
    Invalid(tree, base, "fxtmult", Replace(Profiles(), "    radius: 2", std::string("    radius: ") + badRadius[i]),
            "invalid radius accepted");
  Invalid(tree, base, "fxtmult", Profiles() + "maxVr: 2\n", "duplicate event key accepted");
  tree.Set("fxtmult", Profiles());
  tree.Write(tree.analysis, "analysis:\n  mode: fxtmult\n  mode: refmult\n");
  Require(!selection.Load(tree.main.c_str(), base, "fxtmult", false, errors), "duplicate analysis mode accepted");
  tree.Write(tree.analysis, "analysis:\n  anaName: test\n");
  Require(!selection.Load(tree.main.c_str(), base, "fxtmult", false, errors), "missing analysis mode accepted");
  tree.Set("fxtmult", Profiles());
  tree.Write(tree.main, "analysis: analysis.yaml\nevent: event.yaml\nevent: event.yaml\n");
  Require(!selection.Load(tree.main.c_str(), base, "fxtmult", false, errors), "duplicate mainconf key accepted");
  Require(!selection.Load("", base, "fxtmult", false, errors), "empty mainconf accepted");
  std::cout << "PASS KF event selection: mode profiles, immutable precedence, boundaries, strict YAML, "
               "centrality consistency, sparse/archived Imp5 schemas and inactive dump omission, "
               "default restoration and failed-reload safety" << std::endl;
}
}

inline void TestKfEventSelection() { kf_event_selection_test::Run(); }
#endif


