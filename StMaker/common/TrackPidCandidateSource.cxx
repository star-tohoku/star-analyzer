#include "TrackPidCandidateSource.h"

#include "FlatYamlConfig.h"
#include "HistManager.h"

#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"

#include "TMath.h"
#include "TVector3.h"

#include <iostream>
#include <map>

namespace {
struct TrackSelection {
  Int_t charge;
  Double_t minPt, maxP, maxDca, minNHitsFit, maxAbsNSigma;
  std::string nSigmaHypothesis;
  Bool_t requireTof;
};

class TrackPidCandidateSource : public CandidateSource {
public:
  TrackPidCandidateSource() : mHist(0) {}

  virtual const char* Name() const { return "track_pid"; }

  virtual Bool_t Init(const char* mainConfigPath, HistManager* hist) {
    mHist = hist;
    std::string path;
    if (!FlatYamlConfig::ResolveMainConfigReference(mainConfigPath ? mainConfigPath : "", "trackPid", path))
      return kFALSE;
    FlatYamlConfig cfg;
    if (!cfg.Load(path)) return kFALSE;
    mKeys = cfg.RequireList("particleKeys");
    for (size_t i = 0; i < mKeys.size(); ++i) {
      const std::string& k = mKeys[i];
      TrackSelection s;
      s.charge = cfg.RequireInt(k + "_charge");
      s.minPt = cfg.RequireDouble(k + "_minPt");
      s.maxP = cfg.RequireDouble(k + "_maxP");
      s.maxDca = cfg.RequireDouble(k + "_maxDca");
      s.minNHitsFit = cfg.RequireDouble(k + "_minNHitsFit");
      s.nSigmaHypothesis = cfg.RequireString(k + "_nSigmaHypothesis");
      s.maxAbsNSigma = cfg.RequireDouble(k + "_maxAbsNSigma");
      s.requireTof = cfg.RequireBool(k + "_requireTof");
      if (s.charge != 1 && s.charge != -1) cfg.Fail(k + "_charge must be +1 or -1");
      if (s.nSigmaHypothesis != "pion" && s.nSigmaHypothesis != "kaon" && s.nSigmaHypothesis != "proton" &&
          s.nSigmaHypothesis != "electron")
        cfg.Fail(k + "_nSigmaHypothesis must be pion, kaon, proton or electron");
      if (!(s.minPt >= 0. && s.maxP > s.minPt && s.maxDca > 0. && s.maxAbsNSigma > 0.))
        cfg.Fail(k + ": invalid kinematic/PID cuts");
      mSelections[k] = s;
      mCandidates[k];
      std::cout << "[TrackPidCandidateSource] " << k << ": charge " << s.charge << ", pT >= " << s.minPt
                << ", |p| < " << s.maxP << ", DCA < " << s.maxDca << ", nHitsFit >= " << s.minNHitsFit
                << ", |nSigma(" << s.nSigmaHypothesis << ")| < " << s.maxAbsNSigma
                << (s.requireTof ? ", TOF required" : "") << std::endl;
    }
    cfg.RejectUnused();
    return cfg.Ok();
  }

  virtual std::vector<std::string> SpeciesKeys() const { return mKeys; }
  virtual Bool_t ProvidesSpecies(const std::string& key) const { return mSelections.count(key) > 0; }
  virtual std::vector<std::string> DaughterNames(const std::string&) const {
    return std::vector<std::string>(1, "track");
  }
  virtual std::vector<std::string> VariableNames(const std::string&) const {
    const char* v[] = {"pt", "p", "dca", "nHitsFit", "nSigmaPion", "nSigmaKaon", "nSigmaProton", "beta",
                       "isPrimary"};
    return std::vector<std::string>(v, v + sizeof(v) / sizeof(v[0]));
  }

  virtual Bool_t ProcessEvent(StPicoDst* picoDst, Int_t /*cent9*/) {
    if (!picoDst || !picoDst->event()) { mError = "no event"; return kFALSE; }
    const TVector3 pv = picoDst->event()->primaryVertex();
    const Int_t nTracks = picoDst->numberOfTracks();
    for (std::map<std::string, TrackSelection>::const_iterator sel = mSelections.begin();
         sel != mSelections.end(); ++sel) {
      const TrackSelection& s = sel->second;
      std::vector<FemtoSourceCandidate>& out = mCandidates[sel->first];
      for (Int_t i = 0; i < nTracks; ++i) {
        const StPicoTrack* t = picoDst->track(i);
        if (!t || t->charge() != s.charge) continue;
        const TVector3 g = t->gMom();
        if (g.Pt() < s.minPt || g.Mag() >= s.maxP) continue;
        if (t->nHitsFit() < s.minNHitsFit) continue;
        if (TMath::Abs(NSigma(t, s.nSigmaHypothesis)) >= s.maxAbsNSigma) continue;
        const Double_t dca = t->gDCA(pv).Mag();
        if (dca >= s.maxDca) continue;
        Double_t beta = -1., m2 = -999.;
        if (t->isTofTrack()) {
          const StPicoBTofPidTraits* tof = picoDst->btofPidTraits(t->bTofPidTraitsIndex());
          if (tof && tof->btofBeta() > 0.) {
            beta = tof->btofBeta();
            m2 = g.Mag2() * (1. / (beta * beta) - 1.);
          }
        }
        if (s.requireTof && beta <= 0.) continue;
        FemtoSourceCandidate c;
        const TVector3 p = t->isPrimary() ? t->pMom() : g;
        c.px = p.X(); c.py = p.Y(); c.pz = p.Z();
        c.massValue = m2;
        c.x = pv.X(); c.y = pv.Y(); c.z = pv.Z();
        c.decayLength = 0.f;
        c.daughters.push_back(i);
        const Float_t v[] = {static_cast<Float_t>(g.Pt()), static_cast<Float_t>(g.Mag()), static_cast<Float_t>(dca),
                             static_cast<Float_t>(t->nHitsFit()), t->nSigmaPion(), t->nSigmaKaon(),
                             t->nSigmaProton(), static_cast<Float_t>(beta), t->isPrimary() ? 1.f : 0.f};
        c.variables.assign(v, v + sizeof(v) / sizeof(v[0]));
        out.push_back(c);
      }
    }
    return kTRUE;
  }

  virtual const std::vector<FemtoSourceCandidate>& Candidates(const std::string& key) const {
    std::map<std::string, std::vector<FemtoSourceCandidate> >::const_iterator it = mCandidates.find(key);
    return it != mCandidates.end() ? it->second : mEmpty;
  }

  virtual void Clear() {
    for (std::map<std::string, std::vector<FemtoSourceCandidate> >::iterator it = mCandidates.begin();
         it != mCandidates.end(); ++it)
      it->second.clear();
  }

  virtual const std::string& LastError() const { return mError; }

private:
  static Float_t NSigma(const StPicoTrack* t, const std::string& hypothesis) {
    if (hypothesis == "pion") return t->nSigmaPion();
    if (hypothesis == "kaon") return t->nSigmaKaon();
    if (hypothesis == "proton") return t->nSigmaProton();
    return t->nSigmaElectron();
  }

  HistManager* mHist;
  std::vector<std::string> mKeys;
  std::map<std::string, TrackSelection> mSelections;
  std::map<std::string, std::vector<FemtoSourceCandidate> > mCandidates;
  std::vector<FemtoSourceCandidate> mEmpty;
  std::string mError;
};
}  // namespace

CandidateSource* createTrackPidCandidateSource() { return new TrackPidCandidateSource(); }
