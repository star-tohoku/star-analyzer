#include "OfficialEventSelection.h"

#include "ConfigManager.h"
#include "HistManager.h"
#include "KfEventSelection.h"
#include "cuts/CentralityCutConfig.h"
#include "cuts/KfParticleCutConfig.h"

#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"

#include "TH1.h"

#include <fstream>
#include <iostream>
#include <set>
#include <sstream>

namespace {
class OfficialEventSelection : public FemtoEventSelector {
public:
  OfficialEventSelection() : mHist(0), mEventSelection(0) {
    for (Int_t i = 0; i < 13; ++i) mEventCounts[i] = 0;
    for (Int_t i = 0; i < 8; ++i) mOfficialCounts[i] = 0;
  }
  virtual ~OfficialEventSelection() { delete mEventSelection; }

  virtual const char* Name() const { return "official_run21_3gev"; }

  virtual Bool_t Init(const char* mainConfigPath, HistManager* hist) {
    mHist = hist;
    const KfParticleCutConfig& k = ConfigManager::GetInstance().GetKfParticleCuts();
    if (!k.officialEventSelection) {
      std::cerr << "ERROR: OfficialEventSelection requires officialEventSelection: true in the kf YAML" << std::endl;
      return kFALSE;
    }
    const CentralityCutConfig& centralityCuts = ConfigManager::GetInstance().GetCentralityCuts();
    mEventSelection = new KfEventSelection();
    if (!mEventSelection->Load(mainConfigPath, ConfigManager::GetInstance().GetEventCuts(),
                               centralityCuts.mode, centralityCuts.enabled, std::cerr, false))
      return kFALSE;
    mEventSelection->Dump(std::cout);
    if (!k.officialBadRunFile.empty()) {
      std::ifstream in(k.officialBadRunFile.c_str());
      if (!in) { std::cerr << "ERROR: cannot open officialBadRunFile " << k.officialBadRunFile << std::endl; return kFALSE; }
      // One run number per line; '#' starts a comment. Anything else is an error (a partly read list
      // would silently accept bad runs).
      std::string line;
      Int_t lineNumber = 0;
      while (std::getline(in, line)) {
        ++lineNumber;
        const std::string::size_type hash = line.find('#');
        if (hash != std::string::npos) line = line.substr(0, hash);
        std::istringstream words(line);
        Int_t run = 0;
        std::string extra;
        if (!(words >> run)) {
          if (line.find_first_not_of(" \t\r") == std::string::npos) continue;
          std::cerr << "ERROR: officialBadRunFile " << k.officialBadRunFile << " line " << lineNumber
                    << " is not a run number" << std::endl;
          return kFALSE;
        }
        if (words >> extra) {
          std::cerr << "ERROR: officialBadRunFile " << k.officialBadRunFile << " line " << lineNumber
                    << " has more than one value" << std::endl;
          return kFALSE;
        }
        mBadRuns.insert(run);
      }
      if (mBadRuns.empty()) { std::cerr << "ERROR: empty officialBadRunFile" << std::endl; return kFALSE; }
      std::cout << "[OfficialEventSelection] " << mBadRuns.size() << " bad runs from " << k.officialBadRunFile << std::endl;
    }
    if (mHist) {
      TH1* ev = mHist->Get("hKfEventSelection");
      const char* eventLabels[] = {"read events", "bad run", "invalid vertex", "Vz", "Vr",
          "refMult", "VPD difference", "track count", "pileup", "centrality invalid",
          "centrality bin", "KF error", "reconstructed"};
      if (ev && ev->GetNbinsX() == 13)
        for (Int_t i = 0; i < 13; ++i) ev->GetXaxis()->SetBinLabel(i + 1, eventLabels[i]);
      TH1* off = mHist->Get("hKfOfficialEventSelection");
      const char* offLabels[] = {"read events", "bad run", "trigger", "event cuts", "reconstructed (trees)",
          "pileup", "below 80%", "accepted"};
      if (off && off->GetNbinsX() == 8)
        for (Int_t i = 0; i < 8; ++i) off->GetXaxis()->SetBinLabel(i + 1, offLabels[i]);
    }
    std::cout << "[OfficialEventSelection] " << mBadRuns.size() << " bad runs, " << k.officialTriggerIds.size()
              << " trigger ids, pileup " << (k.officialPileupLow.empty() ? "off" : "on") << ", "
              << k.officialBatchRunEdges.size() + 1 << " run batch(es)" << std::endl;
    return kTRUE;
  }

  virtual Bool_t Accept(StPicoDst* picoDst, Int_t& cent9) {
    cent9 = -1;
    const StPicoEvent* event = picoDst ? picoDst->event() : 0;
    if (!event || !mEventSelection) return kFALSE;
    const KfParticleCutConfig& k = ConfigManager::GetInstance().GetKfParticleCuts();
    ++mEventCounts[0];
    ++mOfficialCounts[0];
    const Int_t runId = event->runId();
    if (mBadRuns.count(runId)) { ++mOfficialCounts[1]; ++mEventCounts[1]; return kFALSE; }
    Int_t triggerMask = 0;
    for (size_t i = 0; i < k.officialTriggerIds.size(); ++i)
      if (event->isTrigger(static_cast<unsigned int>(k.officialTriggerIds[i]))) triggerMask |= (1 << i);
    if (!k.officialTriggerIds.empty() && triggerMask == 0) { ++mOfficialCounts[2]; return kFALSE; }
    const KfEventSelection::Result result = mEventSelection->Check(*event, picoDst->numberOfTracks());
    if (result != KfEventSelection::kAccepted) {
      ++mOfficialCounts[3];
      ++mEventCounts[static_cast<Int_t>(result) + 1];
      return kFALSE;
    }
    const Int_t fxt = event->fxtMult();
    const Int_t nTof = event->nBTOFMatch();
    if (k.officialRequirePositiveMult && (nTof <= 0 || fxt <= 0)) { ++mOfficialCounts[3]; return kFALSE; }
    const Int_t batch = k.OfficialBatch(runId);
    const Bool_t pileupPass = PileupPass(k, fxt, nTof, batch);
    const Int_t c9 = Cent9(k, fxt, batch);
    ++mOfficialCounts[4];
    if (!pileupPass) { ++mOfficialCounts[5]; ++mEventCounts[8]; return kFALSE; }
    if (c9 < 0) { ++mOfficialCounts[6]; ++mEventCounts[10]; return kFALSE; }
    ++mOfficialCounts[7];
    ++mEventCounts[12];
    cent9 = c9;
    return kTRUE;
  }

  virtual void Finish() {
    if (!mHist) return;
    TH1* ev = mHist->Get("hKfEventSelection");
    if (ev) for (Int_t i = 0; i < 13 && i < ev->GetNbinsX(); ++i) ev->SetBinContent(i + 1, mEventCounts[i]);
    TH1* off = mHist->Get("hKfOfficialEventSelection");
    if (off) for (Int_t i = 0; i < 8 && i < off->GetNbinsX(); ++i) off->SetBinContent(i + 1, mOfficialCounts[i]);
    std::cout << "[OfficialEventSelection]";
    for (Int_t i = 0; i < 8; ++i) std::cout << " " << mOfficialCounts[i];
    std::cout << " (read, bad run, trigger, event cuts, passed cuts, pileup, below 80%, accepted)" << std::endl;
  }

  virtual std::string Describe() const {
    std::ostringstream s;
    if (mEventSelection) mEventSelection->Dump(s);
    s << "official bad runs: " << mBadRuns.size() << "\n";
    return s.str();
  }

private:
  // Pile-up band: low(n) < fxtMult < high(n), n = nBTOFMatch, polynomials per run batch.
  static Bool_t PileupPass(const KfParticleCutConfig& k, Int_t fxtMult, Int_t nTof, Int_t batch) {
    if (k.officialPileupLow.empty()) return kTRUE;
    const Double_t x = nTof;
    Double_t lo = 0., hi = 0., xp = 1.;
    const std::vector<double>& cl = k.officialPileupLow[batch];
    const std::vector<double>& ch = k.officialPileupHigh[batch];
    for (size_t i = 0; i < cl.size(); ++i) { lo += cl[i] * xp; xp *= x; }
    xp = 1.;
    for (size_t i = 0; i < ch.size(); ++i) { hi += ch[i] * xp; xp *= x; }
    return fxtMult > lo && fxtMult < hi;
  }
  // cent9 from fxtMult with strict lower bounds (fxtMult > edge): 0 = 70-80% ... 8 = 0-5%; -1 below.
  static Int_t Cent9(const KfParticleCutConfig& k, Int_t fxtMult, Int_t batch) {
    const std::vector<double>& e = k.officialCentralityEdges[batch];
    for (Int_t c = static_cast<Int_t>(e.size()) - 1; c >= 0; --c)
      if (fxtMult > e[c]) return c;
    return -1;
  }

  HistManager* mHist;
  KfEventSelection* mEventSelection;
  std::set<Int_t> mBadRuns;
  Long64_t mEventCounts[13];
  Long64_t mOfficialCounts[8];
};
}  // namespace

FemtoEventSelector* createOfficialEventSelection() { return new OfficialEventSelection(); }
