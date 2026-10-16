#include "StCandidateFemtoMaker.h"

#include "CandidateBdtScorer.h"
#include "CandidateFemtoConfig.h"
#include "CandidateSource.h"
#include "CandidateStore.h"
#include "ConfigManager.h"
#include "FlatYamlConfig.h"
#include "HistManager.h"
#include "cuts/EventCutConfig.h"
#include "cuts/MixingConfig.h"

#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"

#include "THnSparse.h"
#include "TFile.h"
#include "TH1.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TNamed.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TTree.h"
#include "TVector2.h"
#include "TVector3.h"

#include <algorithm>
#include <deque>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <utility>
#include <vector>

namespace {
// Daughter-track kinematics kept with each candidate so that close-pair QA also works for
// mixed-event pairs, whose PicoDst is gone. Global momentum at the DCA to the PV; rMax =
// transverse radius [m] at which the track leaves the TPC through an endcap (straight r-z line).
struct DauKin {
  Float_t eta, phi, pt, rMax;
  Short_t charge;
};

// One usable candidate of a species in one event (mass inside the species mass axis).
struct Entry {
  Float_t px, py, pz, mass, score;
  Int_t massClass;
  Float_t rel[3];  // decay vertex minus the event primary vertex [cm]
  Float_t decayLength;
  std::vector<Int_t> tracks;
  std::vector<DauKin> dau;
  Float_t bT;      // magnetic field [T]
  std::vector<char> channelRemoved;  // per channel: removed by that channel's shared-track rule (ApplySharedTrackRule)
  std::vector<Float_t> treeValues;  // FillCandidateModelValues output, kept until the tree is filled
};

// One past event of a channel's mixing pool: its A and B candidates of this channel (both non-empty).
struct PoolEvent {
  std::vector<Entry> a, b;
};

struct ChannelState {
  ChannelState() : se(0), me(0), seCount(0), meCount(0), shareRejectSE(0), eventsWithBoth(0) {
    for (Int_t i = 0; i < 4; ++i) removedA[i] = removedB[i] = 0;
  }
  std::map<Int_t, std::deque<PoolEvent> > pool;  // per mixing class, at most bufferSize events
  std::vector<Entry> a, b;  // this event's candidates of the channel after the channel's shared-track rule
  THnSparseF* se;
  THnSparseF* me;
  Long64_t seCount, meCount, shareRejectSE, eventsWithBoth;
  Long64_t removedA[4], removedB[4];  // removed by this channel's shared-track rule, per mass class
};

struct SpeciesState {
  SpeciesState() : source(0), bdt(0), tree(0), usable(0), signal(0) {}
  CandidateSource* source;
  CandidateBdtScorer* bdt;
  std::vector<Double_t> scoreEdges;
  std::vector<std::string> daughterNames;
  std::vector<std::string> variableNames;
  std::vector<std::string> modelVariableNames;  // CandidateModelVariableNames(variableNames)
  std::vector<Float_t> modelFloats;             // FillCandidateModelValues: tree and BDT see the same values
  std::vector<double> modelValues;
  std::vector<Entry> entries;
  TTree* tree;
  std::vector<Float_t> treeVariables;
  Float_t treeMass, treePx, treePy, treePz, treeScore;
  std::vector<Int_t> treeChannels;       // channels with this species: one flag branch sharedTrackRemoved_<channel> each
  std::vector<Char_t> treeChannelFlags;  // one byte each, written as /O (std::vector<bool> is packed: no addresses)
  Long64_t usable, signal;   // every usable candidate (the shared-track rule acts per channel)
};

Bool_t SharesTrack(const Entry& a, const Entry& b) {
  for (size_t i = 0; i < a.tracks.size(); ++i)
    for (size_t j = 0; j < b.tracks.size(); ++j)
      if (a.tracks[i] >= 0 && a.tracks[i] == b.tracks[j]) return kTRUE;
  return kFALSE;
}

// Index i * nB + j of the first shared (A daughter i, B daughter j) combination, -1 if none.
Int_t SharedCombo(const Entry& a, const Entry& b) {
  for (size_t i = 0; i < a.tracks.size(); ++i)
    for (size_t j = 0; j < b.tracks.size(); ++j)
      if (a.tracks[i] >= 0 && a.tracks[i] == b.tracks[j]) return static_cast<Int_t>(i * b.tracks.size() + j);
  return -1;
}

Double_t ComputeKStar(const Entry& a, Double_t massA, const Entry& b, Double_t massB) {
  TLorentzVector va; va.SetXYZM(a.px, a.py, a.pz, massA);
  TLorentzVector vb; vb.SetXYZM(b.px, b.py, b.pz, massB);
  const TLorentzVector pair = va + vb;
  TLorentzVector q = va - vb;
  q.Boost(-pair.BoostVector());
  return 0.5 * q.Vect().Mag();
}
}  // namespace

class StCandidateFemtoMakerImpl {
public:
  StCandidateFemtoMakerImpl()
      : selector(0), hist(0), outputFile(0), treeFile(0), picoDst(0), storeWriter(0), storeReader(0), storeFinished(kFALSE),
        readEvents(0), acceptedEvents(0), cent9(-1), runId(0), eventId(0), trackQaInit(kFALSE) {
    for (Int_t i = 0; i < 9; ++i) trackQaDedx[i] = trackQaM2[i] = 0;
  }
  ~StCandidateFemtoMakerImpl() {
    for (size_t i = 0; i < species.size(); ++i) delete species[i].bdt;
    for (size_t i = 0; i < sources.size(); ++i) delete sources[i];
    delete selector;
    delete storeWriter;
    delete storeReader;
    delete hist;
  }

  CandidateFemtoConfig cfg;
  FemtoEventSelector* selector;
  std::vector<CandidateSource*> sources;
  std::vector<SpeciesState> species;
  std::vector<ChannelState> channels;
  HistManager* hist;
  TFile* outputFile;
  TFile* treeFile;
  StPicoDst* picoDst;
  CandidateStoreWriter* storeWriter;  // reconstruction stage (writeCandidateStore: true)
  CandidateStoreReader* storeReader;  // pair stage (SetCandidateStoreInput)
  Bool_t storeFinished;
  std::string bdtConfigPath;
  Long64_t readEvents, acceptedEvents;
  Int_t cent9, runId, eventId;
  Bool_t trackQaInit;
  TH1* trackQaDedx[9];
  TH1* trackQaM2[9];

  Bool_t InitBdt(const std::string& mainConfigPath);
  THnSparseF* BookPairHistogram(const CandidateFemtoChannel& ch, const char* tag);
  FemtoDaughterTrack PicoDaughterTrack(Int_t trackIndex) const;
  DauKin MakeDauKin(const FemtoDaughterTrack& t, Double_t vzCm) const;
  Bool_t DaughterTracks(const FemtoSourceCandidate& c, std::vector<FemtoDaughterTrack>& out) const;
  Double_t DPhiStarMin(const DauKin& a, Float_t bTA, const DauKin& b, Float_t bTB) const;
  Bool_t BuildEntries(Int_t s, const TVector3& pv, Float_t bT);
  void ApplySharedTrackRule(Int_t c);
  void FillSpecies(Int_t s);
  void FillTrackPidQa(const TVector3& pv);
  void FillEventQa(Int_t c, const std::vector<Entry>& a, const std::vector<Entry>& b, Int_t fxtMult, Int_t nTof,
                   Double_t vr);
  void ProcessChannel(Int_t c, Int_t mixBin);
  void FillPair(Int_t c, const Entry& a, const Entry& b, Bool_t sameEvent);
  void FillMixedEventPairs(Int_t c, const std::deque<PoolEvent>& pool);
  Int_t GetMixingBin(Double_t vz) const;
};

THnSparseF* StCandidateFemtoMakerImpl::BookPairHistogram(const CandidateFemtoChannel& ch, const char* tag) {
  const CandidateFemtoSpecies& a = cfg.species[ch.indexA];
  const CandidateFemtoSpecies& b = cfg.species[ch.indexB];
  const SpeciesState& sa = species[ch.indexA];
  const SpeciesState& sb = species[ch.indexB];
  Int_t nBins[8] = {cfg.kstarAxis.nBins, 9, 4, 4, a.massAxis.nBins, b.massAxis.nBins,
                    static_cast<Int_t>(sa.scoreEdges.size()) - 1, static_cast<Int_t>(sb.scoreEdges.size()) - 1};
  Double_t lo[8] = {cfg.kstarAxis.min, -0.5, -0.5, -0.5, a.massAxis.min, b.massAxis.min,
                    sa.scoreEdges.front(), sb.scoreEdges.front()};
  Double_t hi[8] = {cfg.kstarAxis.max, 8.5, 3.5, 3.5, a.massAxis.max, b.massAxis.max,
                    sa.scoreEdges.back(), sb.scoreEdges.back()};
  const TString name = TString::Format("hPair%s_%s", tag, ch.name.c_str());
  const TString title = TString::Format("%s pairs %s (%s x %s)", tag, ch.name.c_str(), a.key.c_str(), b.key.c_str());
  THnSparseF* h = new THnSparseF(name, title, 8, nBins, lo, hi);
  h->GetAxis(6)->Set(nBins[6], &sa.scoreEdges[0]);
  h->GetAxis(7)->Set(nBins[7], &sb.scoreEdges[0]);
  const char* axisNames[8] = {"kstar", "cent9", "massClassA", "massClassB", "massA", "massB", "scoreA", "scoreB"};
  const TString axisTitles[8] = {"k* [GeV/c]", "cent9 (0 = 70-80%, 8 = 0-5%)",
      "mass class " + TString(a.key.c_str()) + " (0 none, 1 signal, 2 left SB, 3 right SB)",
      "mass class " + TString(b.key.c_str()) + " (0 none, 1 signal, 2 left SB, 3 right SB)",
      "mass variable " + TString(a.key.c_str()), "mass variable " + TString(b.key.c_str()),
      "BDT score " + TString(a.key.c_str()), "BDT score " + TString(b.key.c_str())};
  for (Int_t i = 0; i < 8; ++i) {
    h->GetAxis(i)->SetName(axisNames[i]);
    h->GetAxis(i)->SetTitle(axisTitles[i]);
  }
  h->Sumw2();
  return h;
}

Bool_t StCandidateFemtoMakerImpl::InitBdt(const std::string& mainConfigPath) {
  for (size_t i = 0; i < species.size(); ++i) {
    species[i].modelVariableNames = CandidateModelVariableNames(species[i].variableNames);
    const std::vector<std::string>& n = species[i].modelVariableNames;
    for (size_t a = 0; a < n.size(); ++a)
      for (size_t b = a + 1; b < n.size(); ++b)
        if (n[a] == n[b]) {
          std::cerr << "ERROR: species " << cfg.species[i].key << ": variable name '" << n[a]
                    << "' is used twice (event and source variables of the candidate trees must be unique)" << std::endl;
          return kFALSE;
        }
    // No BDT: one score bin around the default score (-2).
    species[i].scoreEdges.clear();
    species[i].scoreEdges.push_back(-3.);
    species[i].scoreEdges.push_back(3.);
  }
  std::map<std::string, std::string> main;
  if (!YamlParser::ParseFile(mainConfigPath.c_str(), main)) return kFALSE;
  if (main.find("candidateBdt") == main.end()) {
    std::cout << "[StCandidateFemtoMaker] no candidateBdt key in the mainconf: no BDT scores" << std::endl;
    return kTRUE;
  }
  if (!FlatYamlConfig::ResolveMainConfigReference(mainConfigPath, "candidateBdt", bdtConfigPath)) return kFALSE;
  FlatYamlConfig c;
  if (!c.Load(bdtConfigPath)) return kFALSE;
  const std::string configDir = mainConfigPath.substr(0, mainConfigPath.rfind("config/mainconf/")) + "config/";
  const std::vector<std::string> keys = c.RequireList("speciesKeys");
  const Double_t outside = c.RequireDouble("outsideTrainingScore");
  if (c.Ok() && outside >= -1.) c.Fail("outsideTrainingScore must be below -1 (BDT scores are in [-1, 1])");
  for (size_t k = 0; k < keys.size() && c.Ok(); ++k) {
    const Int_t s = cfg.SpeciesIndex(keys[k]);
    if (s < 0) { c.Fail("species '" + keys[k] + "' is not in the candidateFemto speciesKeys"); break; }
    if (species[s].bdt) { c.Fail("species '" + keys[k] + "' is listed twice in speciesKeys"); break; }
    const std::string model = configDir + c.RequireString(keys[k] + "_model");
    const std::vector<std::string> edges = c.RequireList(keys[k] + "_scoreEdges");
    std::vector<Double_t> e;
    for (size_t i = 0; i < edges.size(); ++i) {
      bool ok = false;
      e.push_back(FlatYamlConfig::ParseDouble(edges[i], ok));
      if (!ok || (i > 0 && e[i] <= e[i - 1])) c.Fail(keys[k] + "_scoreEdges must be increasing numbers");
    }
    // first bin: the outside-training score only; the other bins cover [-1, 1]
    if (c.Ok() && (e.size() < 3 || !(e[0] < outside && outside < e[1]) || e[1] > -1. || e.back() < 1.))
      c.Fail(keys[k] + "_scoreEdges must start with a bin around outsideTrainingScore, then cover [-1, 1] "
             "(e.g. -3, -1, ..., 1)");
    if (!c.Ok()) break;
    SpeciesState& st = species[s];
    st.scoreEdges = e;
    st.bdt = new CandidateBdtScorer();
    if (!st.bdt->Init(keys[k], st.modelVariableNames, model, static_cast<Float_t>(outside))) return kFALSE;
  }
  c.RejectUnused();
  return c.Ok();
}

// The PicoDst values of a final-state track (all zeros if the index is not a track of this event).
FemtoDaughterTrack StCandidateFemtoMakerImpl::PicoDaughterTrack(Int_t trackIndex) const {
  FemtoDaughterTrack d;
  if (!picoDst || trackIndex < 0 || trackIndex >= static_cast<Int_t>(picoDst->numberOfTracks())) return d;
  const StPicoTrack* t = picoDst->track(trackIndex);
  if (!t) return d;
  const TVector3 g = t->gMom();  // built from the stored Float_t components: the casts below are exact
  d.gx = static_cast<Float_t>(g.X()); d.gy = static_cast<Float_t>(g.Y()); d.gz = static_cast<Float_t>(g.Z());
  d.charge = t->charge();
  return d;
}

// The candidate's own daughter tracks (candidate store) or, if it has none, those of the PicoDst.
Bool_t StCandidateFemtoMakerImpl::DaughterTracks(const FemtoSourceCandidate& c, std::vector<FemtoDaughterTrack>& out) const {
  out.clear();
  if (!c.daughterTracks.empty()) {
    if (c.daughterTracks.size() != c.daughters.size()) return kFALSE;
    out = c.daughterTracks;
    return kTRUE;
  }
  if (!picoDst && !c.daughters.empty()) return kFALSE;
  for (size_t d = 0; d < c.daughters.size(); ++d) out.push_back(PicoDaughterTrack(c.daughters[d]));
  return kTRUE;
}

DauKin StCandidateFemtoMakerImpl::MakeDauKin(const FemtoDaughterTrack& t, Double_t vzCm) const {
  DauKin k;
  k.eta = 0.f; k.phi = 0.f; k.pt = -1.f; k.rMax = 0.f; k.charge = 0;
  const TVector3 g(t.gx, t.gy, t.gz);
  if (g.Pt() <= 0.) return k;
  k.eta = g.PseudoRapidity(); k.phi = g.Phi(); k.pt = g.Pt(); k.charge = static_cast<Short_t>(t.charge);
  const Double_t halfCm = 100. * cfg.closePairTpcHalfLength;
  const Double_t dz = (g.Z() < 0.) ? vzCm + halfCm : halfCm - vzCm;  // cm to the endcap
  k.rMax = (TMath::Abs(g.Z()) > 1e-9) ? 0.01 * TMath::Max(0., dz) * g.Pt() / TMath::Abs(g.Z()) : 1e3;
  return k;
}

// phi* = phi + asin(-0.3 q B r / (2 pT)) [B in T, r in m, pT in GeV/c], minimised over the
// configured radii up to the smaller rMax of the two tracks.
Double_t StCandidateFemtoMakerImpl::DPhiStarMin(const DauKin& a, Float_t bTA, const DauKin& b, Float_t bTB) const {
  const Double_t rEnd = TMath::Min(cfg.closePairRadiusMax, static_cast<Double_t>(TMath::Min(a.rMax, b.rMax)));
  Double_t best = 1e9;
  for (Double_t r = cfg.closePairRadiusMin;; r += cfg.closePairRadiusStep) {
    Double_t argA = -0.3 * a.charge * bTA * r / (2. * a.pt);
    Double_t argB = -0.3 * b.charge * bTB * r / (2. * b.pt);
    argA = TMath::Max(-1., TMath::Min(1., argA));
    argB = TMath::Max(-1., TMath::Min(1., argB));
    Double_t d = (a.phi + TMath::ASin(argA)) - (b.phi + TMath::ASin(argB));
    while (d > TMath::Pi()) d -= TMath::TwoPi();
    while (d < -TMath::Pi()) d += TMath::TwoPi();
    if (TMath::Abs(d) < TMath::Abs(best)) best = d;
    if (r + cfg.closePairRadiusStep > rEnd + 1e-9) break;
  }
  return best;
}

Bool_t StCandidateFemtoMakerImpl::BuildEntries(Int_t s, const TVector3& pv, Float_t bT) {
  const CandidateFemtoSpecies& sp = cfg.species[s];
  SpeciesState& st = species[s];
  st.entries.clear();
  const std::vector<FemtoSourceCandidate>& cands = st.source->Candidates(sp.particleKey);
  std::vector<FemtoDaughterTrack> dauTracks;
  for (size_t i = 0; i < cands.size(); ++i) {
    const FemtoSourceCandidate& c = cands[i];
    if (!sp.massAxis.Contains(c.massValue)) continue;
    Entry e;
    e.px = c.px; e.py = c.py; e.pz = c.pz;
    e.mass = c.massValue;
    e.massClass = sp.MassClass(c.massValue);
    if (c.variables.size() != st.variableNames.size()) {
      std::cerr << "ERROR: source " << st.source->Name() << " gave " << c.variables.size() << " variables for "
                << sp.key << ", VariableNames() has " << st.variableNames.size() << std::endl;
      return kFALSE;
    }
    if (st.bdt || st.tree) FillCandidateModelValues(c, runId, eventId, cent9, st.modelFloats, st.modelValues);
    e.score = st.bdt ? st.bdt->Score(st.modelValues, runId) : -2.f;
    e.rel[0] = c.x - pv.X(); e.rel[1] = c.y - pv.Y(); e.rel[2] = c.z - pv.Z();
    e.decayLength = c.decayLength;
    e.tracks = c.daughters;
    if (!DaughterTracks(c, dauTracks)) {
      std::cerr << "ERROR: " << sp.key << " candidate without daughter tracks (no PicoDst and none stored)" << std::endl;
      return kFALSE;
    }
    for (size_t d = 0; d < dauTracks.size(); ++d) e.dau.push_back(MakeDauKin(dauTracks[d], pv.Z()));
    e.bT = bT;
    e.channelRemoved.assign(cfg.channels.size(), 0);
    if (st.tree) e.treeValues = st.modelFloats;
    st.entries.push_back(e);
  }
  return kTRUE;
}

// Shared-track rule, per channel: the rule concerns the two particles whose
// correlation the channel measures. When a candidate of species A and a candidate of species B of the same
// event share a daughter track and at least one of them is in its signal window, both are removed from this
// channel: from its SE and ME pairs, its mixing pool and its candidate histograms (hChannel*, the purity
// inputs of the CF macros). Sharing between candidates of the same species is not considered (identical-species
// channels are not supported; if they ever are, the rule applies inside such a channel). A candidate removed in
// any channel is NOT removed elsewhere: each correlation function is computed independently.
// The per-species histograms and the candidate trees keep every usable candidate; the trees carry one flag per
// channel (sharedTrackRemoved_<channel>), so a BDT recipe decides itself whether to train without them. Rationale: a candidate that looks good
// only because it borrowed a track of the partner may come from a pathological reconstruction; the rule depends
// only on the mass windows, not on BDT thresholds chosen later. Pairs that still share a track afterwards (both
// candidates outside the signal windows) are rejected by the SE veto in FillPair.
void StCandidateFemtoMakerImpl::ApplySharedTrackRule(Int_t c) {
  const CandidateFemtoChannel& ch = cfg.channels[c];
  ChannelState& cs = channels[c];
  std::vector<Entry>& ea = species[ch.indexA].entries;
  std::vector<Entry>& eb = species[ch.indexB].entries;
  std::vector<char> ra(ea.size(), 0), rb(eb.size(), 0);
  for (size_t i = 0; i < ea.size(); ++i)
    for (size_t j = 0; j < eb.size(); ++j)
      if ((ea[i].massClass == 1 || eb[j].massClass == 1) && SharesTrack(ea[i], eb[j])) ra[i] = rb[j] = 1;
  cs.a.clear();
  cs.b.clear();
  const char* name = ch.name.c_str();
  for (Int_t side = 0; side < 2; ++side) {
    std::vector<Entry>& e = side == 0 ? ea : eb;
    const std::vector<char>& r = side == 0 ? ra : rb;
    std::vector<Entry>& kept = side == 0 ? cs.a : cs.b;
    Long64_t* removed = side == 0 ? cs.removedA : cs.removedB;
    const char* key = cfg.species[side == 0 ? ch.indexA : ch.indexB].key.c_str();
    for (size_t i = 0; i < e.size(); ++i) {
      if (!r[i]) { kept.push_back(e[i]); continue; }
      e[i].channelRemoved[c] = 1;
      ++removed[(e[i].massClass >= 0 && e[i].massClass <= 3) ? e[i].massClass : 0];
      if (hist) hist->Fill(Form("hChannelSharedTrackRemovedMass_%s_%s", name, key), e[i].mass);
    }
  }
}

// Counters, candidate tree and QA histograms of every usable candidate of the species (the shared-track rule
// acts per channel; the tree records each channel's decision in sharedTrackRemoved_<channel>).
void StCandidateFemtoMakerImpl::FillSpecies(Int_t s) {
  const CandidateFemtoSpecies& sp = cfg.species[s];
  SpeciesState& st = species[s];
  for (size_t i = 0; i < st.entries.size(); ++i) {
    const Entry& e = st.entries[i];
    ++st.usable;
    if (e.massClass == 1) ++st.signal;
    if (!st.tree) continue;
    st.treeMass = e.treeValues[3]; st.treePx = e.treeValues[4]; st.treePy = e.treeValues[5];
    st.treePz = e.treeValues[6]; st.treeScore = e.score;
    for (size_t v = 0; v < st.treeVariables.size(); ++v)
      st.treeVariables[v] = e.treeValues[kCandidateModelFirstSourceVariable + v];
    for (size_t k = 0; k < st.treeChannels.size(); ++k) st.treeChannelFlags[k] = e.channelRemoved[st.treeChannels[k]] ? 1 : 0;
    st.tree->Fill();
  }
  if (!hist) return;
  const char* k = sp.key.c_str();
  for (size_t i = 0; i < st.entries.size(); ++i) {
    const Entry& e = st.entries[i];
    const TVector3 p(e.px, e.py, e.pz);
    hist->Fill(Form("hCandidateMass_%s", k), e.mass);
    hist->Fill(Form("hCandidateMassVsP_%s", k), p.Mag(), e.mass);
    hist->Fill(Form("hCandidateMassVsCent9_%s", k), e.mass, cent9);
    if (st.bdt) {
      hist->Fill(Form("hCandidateScore_%s", k), e.score);
      hist->Fill(Form("hCandidateScoreVsMass_%s", k), e.mass, e.score);
      // per-cent9 mass spectra above any score threshold (purity fits of the CF macros)
      hist->Fill(Form("hCandidateMassVsCent9VsScore_%s", k), e.mass, cent9, e.score);
    }
    if (e.massClass != 1) continue;
    TLorentzVector v; v.SetVectM(p, sp.pdgMass);
    if (v.E() > TMath::Abs(v.Pz())) {
      hist->Fill(Form("hCandidatePtVsYLab_%s_signal", k), v.Rapidity(), v.Pt());
      hist->Fill(Form("hCandidatePtVsYCm_%s_signal", k), -(v.Rapidity() - cfg.yMidLab), v.Pt());
    }
    for (size_t d = 0; d < e.dau.size() && d < st.daughterNames.size(); ++d) {
      if (e.dau[d].pt <= 0.f) continue;
      hist->Fill(Form("hDaughterEtaLab_%s_%s_signal", k, st.daughterNames[d].c_str()), e.dau[d].eta);
      if (d < sp.daughterMasses.size()) {
        TLorentzVector dv; dv.SetPtEtaPhiM(e.dau[d].pt, e.dau[d].eta, e.dau[d].phi, sp.daughterMasses[d]);
        if (dv.E() > TMath::Abs(dv.Pz()))
          hist->Fill(Form("hDaughterPtVsYLab_%s_%s_signal", k, st.daughterNames[d].c_str()), dv.Rapidity(), dv.Pt());
      }
    }
  }
}

// dE/dx and TOF m^2 vs q|p| of all charged tracks of accepted events, per cent9 (QA only).
void StCandidateFemtoMakerImpl::FillTrackPidQa(const TVector3& pv) {
  if (!hist || cent9 < 0 || cent9 > 8) return;
  if (!trackQaInit) {
    for (Int_t i = 0; i < 9; ++i) {
      trackQaDedx[i] = hist->Get(TString::Format("hTrackQA_dEdxVsQP_cent%d", i).Data());
      trackQaM2[i] = hist->Get(TString::Format("hTrackQA_m2VsQP_cent%d", i).Data());
    }
    trackQaInit = kTRUE;
  }
  TH1* hD = trackQaDedx[cent9];
  TH1* hM = trackQaM2[cent9];
  if (!hD && !hM) return;
  const Int_t nTracks = picoDst->numberOfTracks();
  for (Int_t i = 0; i < nTracks; ++i) {
    const StPicoTrack* t = picoDst->track(i);
    if (!t || t->charge() == 0) continue;
    const TVector3 g = t->gMom();
    if (g.Pt() < cfg.trackQaMinPt) continue;
    if (t->nHitsFit() < cfg.trackQaMinNHitsFit) continue;
    if (t->gDCA(pv).Mag() >= cfg.trackQaMaxDca) continue;
    const Double_t qp = t->charge() * g.Mag();
    if (hD) hD->Fill(qp, t->dEdx());
    if (!hM || !t->isTofTrack()) continue;
    const StPicoBTofPidTraits* tof = picoDst->btofPidTraits(t->bTofPidTraitsIndex());
    if (!tof || tof->btofBeta() <= 0.) continue;
    const Double_t beta = tof->btofBeta();
    hM->Fill(qp, g.Mag2() * (1. / (beta * beta) - 1.));
  }
}

// Candidate multiplicities and masses vs event variables, per channel (counts clamped).
void StCandidateFemtoMakerImpl::FillEventQa(Int_t c, const std::vector<Entry>& a, const std::vector<Entry>& b,
                                            Int_t fxtMult, Int_t nTof, Double_t vr) {
  if (!hist) return;
  const char* ch = cfg.channels[c].name.c_str();
  const Int_t nA = a.size(), nB = b.size();
  const Bool_t withB = nB > 0;
  hist->Fill(Form("hEvtQA_NAVsNB_%s", ch), std::min(nA, 40), std::min(nB, 4));
  if (withB) hist->Fill(Form("hEvtQA_MultVsNTof_%s_withB", ch), fxtMult, nTof);
  for (Int_t i = 0; i < nA; ++i) {
    const Double_t m = a[i].mass;
    hist->Fill(Form("hEvtQA_MassAVsNB_%s", ch), m, std::min(nB, 3));
    hist->Fill(Form("hEvtQA_MassAVsMult_%s_all", ch), m, fxtMult);
    hist->Fill(Form("hEvtQA_MassAVsNTof_%s_all", ch), m, nTof);
    hist->Fill(Form("hEvtQA_MassAVsVr_%s_all", ch), m, vr);
    if (!withB) continue;
    hist->Fill(Form("hEvtQA_MassAVsMult_%s_withB", ch), m, fxtMult);
    hist->Fill(Form("hEvtQA_MassAVsNTof_%s_withB", ch), m, nTof);
    hist->Fill(Form("hEvtQA_MassAVsVr_%s_withB", ch), m, vr);
  }
  for (Int_t j = 0; j < nB; ++j) {
    hist->Fill(Form("hEvtQA_MassBVsNA_%s", ch), b[j].mass, std::min(nA, 40));
    hist->Fill(Form("hEvtQA_MassBVsMult_%s", ch), b[j].mass, fxtMult);
  }
}

Int_t StCandidateFemtoMakerImpl::GetMixingBin(Double_t vz) const {
  const EventCutConfig& ev = ConfigManager::GetInstance().GetEventCuts();
  const MixingConfig& mix = ConfigManager::GetInstance().GetMixingConfig();
  Int_t vzBin = 0;
  if (mix.nVzBins > 0) {
    const Double_t span = ev.maxVz - ev.minVz;
    if (span > 0) {
      vzBin = static_cast<Int_t>((vz - ev.minVz) / span * mix.nVzBins);
      if (vzBin < 0) vzBin = 0;
      if (vzBin >= mix.nVzBins) vzBin = mix.nVzBins - 1;
    }
  }
  Int_t centBin = 0;
  if (mix.nCentralityBins > 1 && cent9 >= 0) {
    centBin = cent9;
    if (centBin >= mix.nCentralityBins) centBin = mix.nCentralityBins - 1;
  }
  return vzBin + mix.nVzBins * centBin;
}

void StCandidateFemtoMakerImpl::FillPair(Int_t c, const Entry& a, const Entry& b, Bool_t sameEvent) {
  const CandidateFemtoChannel& ch = cfg.channels[c];
  const CandidateFemtoSpecies& sa = cfg.species[ch.indexA];
  const CandidateFemtoSpecies& sb = cfg.species[ch.indexB];
  ChannelState& cs = channels[c];
  const Bool_t windows = a.massClass > 0 && b.massClass > 0;
  const Bool_t signalSignal = a.massClass == 1 && b.massClass == 1;
  const Double_t kstar = ComputeKStar(a, sa.pdgMass, b, sb.pdgMass);
  const char* name = ch.name.c_str();
  // PicoDst track indices are unique only within one event: the shared-track veto applies to SE
  // pairs only (ME pairs always come from different events).
  if (sameEvent && SharesTrack(a, b)) {
    if (!windows) return;
    ++cs.shareRejectSE;
    if (hist) {
      hist->Fill(Form("hShareRejectSE_%s", name), SharedCombo(a, b));
      hist->Fill(Form("hKstarShareRejectSE_%s", name), kstar);
    }
    return;
  }
  const Double_t x[8] = {kstar, static_cast<Double_t>(cent9), static_cast<Double_t>(a.massClass),
                         static_cast<Double_t>(b.massClass), a.mass, b.mass, a.score, b.score};
  (sameEvent ? cs.se : cs.me)->Fill(x);
  if (windows) { if (sameEvent) ++cs.seCount; else ++cs.meCount; }
  if (!hist) return;
  const char* tag = sameEvent ? "SE" : "ME";
  if (b.massClass == 1) {
    // Split/duplicate-track QA: smallest (eta, phi) distance between same-charge daughters of A
    // and B, vs mass A and k*; common-secondary-vertex QA: distance between the two decay
    // vertices, each relative to its own event PV.
    Double_t dMin = 1e9, dPt = 0.;
    for (size_t i = 0; i < a.dau.size(); ++i)
      for (size_t j = 0; j < b.dau.size(); ++j) {
        const DauKin& u = a.dau[i];
        const DauKin& v = b.dau[j];
        if (u.pt <= 0.f || v.pt <= 0.f || u.charge != v.charge) continue;
        const Double_t dEta = u.eta - v.eta;
        const Double_t dPhi = TVector2::Phi_mpi_pi(u.phi - v.phi);
        const Double_t d = TMath::Sqrt(dEta * dEta + dPhi * dPhi);
        if (d < dMin) { dMin = d; dPt = 2. * (u.pt - v.pt) / (u.pt + v.pt + 1e-9); }
      }
    if (dMin < 1e8) {
      hist->Fill(Form("hDaughterMinDist%s_%s", tag, name), dMin, a.mass, kstar);
      if (kstar < 0.5) hist->Fill(Form("hDaughterMinDistVsDPt%s_%s", tag, name), dMin, dPt);
    }
    const Double_t dx = a.rel[0] - b.rel[0], dy = a.rel[1] - b.rel[1], dz = a.rel[2] - b.rel[2];
    hist->Fill(Form("hDecayVertexDistance%s_%s", tag, name), TMath::Sqrt(dx * dx + dy * dy + dz * dz), a.mass, kstar);
    hist->Fill(Form("hDecayLengthA%s_%s", tag, name), a.decayLength, a.mass, kstar);
    if (sameEvent) {
      const Double_t r = TMath::Sqrt(a.rel[0] * a.rel[0] + a.rel[1] * a.rel[1] + a.rel[2] * a.rel[2]);
      hist->Fill(Form("hDecayVertexVsDecayLengthA_%s", name), r, a.decayLength);
    }
  }
  if (!signalSignal) return;
  hist->Fill(Form("hKstar%s_%s", tag, name), kstar);
  if (cent9 >= 0) hist->Fill(Form("hKstar%sVsCent_%s", tag, name), kstar, cent9);
  const SpeciesState& stA = species[ch.indexA];
  const SpeciesState& stB = species[ch.indexB];
  for (size_t i = 0; i < a.dau.size() && i < stA.daughterNames.size(); ++i)
    for (size_t j = 0; j < b.dau.size() && j < stB.daughterNames.size(); ++j) {
      if (a.dau[i].pt <= 0.f || b.dau[j].pt <= 0.f) continue;
      hist->Fill(Form("hDEtaDPhiStar%s_%s_%s_%s", tag, name, stA.daughterNames[i].c_str(), stB.daughterNames[j].c_str()),
                 DPhiStarMin(a.dau[i], a.bT, b.dau[j], b.bT), a.dau[i].eta - b.dau[j].eta, kstar);
    }
  TLorentzVector va; va.SetXYZM(a.px, a.py, a.pz, sa.pdgMass);
  TLorentzVector vb; vb.SetXYZM(b.px, b.py, b.pz, sb.pdgMass);
  const TLorentzVector pair = va + vb;
  const Double_t kT = 0.5 * pair.Pt();
  const Double_t mAvg = 0.5 * (sa.pdgMass + sb.pdgMass);
  hist->Fill(Form("hPairYCmVsKstar%s_%s", tag, name), kstar, -(pair.Rapidity() - cfg.yMidLab));
  hist->Fill(Form("hPairMtVsKstar%s_%s", tag, name), kstar, TMath::Sqrt(kT * kT + mAvg * mAvg));
}

// Mixed pairs of an event that has both species of the channel: current A x every B of the pool events, then
// (mixBothDirections) every A of the pool events x current B. All pairs are used (no subsampling), so every
// event enters ME with the same per-pair weight as in SE.
void StCandidateFemtoMakerImpl::FillMixedEventPairs(Int_t c, const std::deque<PoolEvent>& pool) {
  const MixingConfig& mix = ConfigManager::GetInstance().GetMixingConfig();
  const std::vector<Entry>& a = channels[c].a;
  const std::vector<Entry>& b = channels[c].b;
  for (size_t ie = 0; ie < pool.size(); ++ie)
    for (size_t i = 0; i < a.size(); ++i)
      for (size_t j = 0; j < pool[ie].b.size(); ++j) FillPair(c, a[i], pool[ie].b[j], kFALSE);
  if (!mix.mixBothDirections) return;
  for (size_t ie = 0; ie < pool.size(); ++ie)
    for (size_t i = 0; i < pool[ie].a.size(); ++i)
      for (size_t j = 0; j < b.size(); ++j) FillPair(c, pool[ie].a[i], b[j], kFALSE);
}

// Per channel, with the candidates left by the channel's shared-track rule. SE pairs exist only in events with
// both species; ME is made from the same kind of events: an event takes part in mixing (as the current event and
// in the pool) only if it has at least one A and one B candidate. Order: SE -> ME against the past pool -> store
// the current event (no self-mixing). The channel candidate histograms (hChannel*) are filled from the same
// events, so the purity inputs are exactly the candidates that enter the channel's pairs.
void StCandidateFemtoMakerImpl::ProcessChannel(Int_t c, Int_t mixBin) {
  const CandidateFemtoChannel& ch = cfg.channels[c];
  ChannelState& cs = channels[c];
  const std::vector<Entry>& a = cs.a;
  const std::vector<Entry>& b = cs.b;
  if (a.empty() || b.empty()) return;
  ++cs.eventsWithBoth;
  if (hist) {
    const char* name = ch.name.c_str();
    for (Int_t side = 0; side < 2; ++side) {
      const Int_t s = side == 0 ? ch.indexA : ch.indexB;
      const std::vector<Entry>& e = side == 0 ? a : b;
      const char* key = cfg.species[s].key.c_str();
      for (size_t i = 0; i < e.size(); ++i) {
        hist->Fill(Form("hChannelMassVsCent9_%s_%s", name, key), e[i].mass, cent9);
        if (species[s].bdt) hist->Fill(Form("hChannelMassVsCent9VsScore_%s_%s", name, key), e[i].mass, cent9, e[i].score);
      }
    }
  }
  for (size_t i = 0; i < a.size(); ++i)
    for (size_t j = 0; j < b.size(); ++j) FillPair(c, a[i], b[j], kTRUE);
  std::deque<PoolEvent>& pool = cs.pool[mixBin];
  FillMixedEventPairs(c, pool);
  PoolEvent evt;
  evt.a = a;
  evt.b = b;
  pool.push_back(evt);
  const Int_t depth = ConfigManager::GetInstance().GetMixingConfig().bufferSize;
  while (static_cast<Int_t>(pool.size()) > depth) pool.pop_front();
}

// ---------------------------------------------------------------------------------------------

StCandidateFemtoMaker::StCandidateFemtoMaker(const char* name, StPicoDstMaker* picoMaker, const char* outName)
    : StMaker(name), mPicoDstMaker(picoMaker), mOutName(outName), mMainConfigPath(""),
      mImpl(new StCandidateFemtoMakerImpl()) {}

StCandidateFemtoMaker::~StCandidateFemtoMaker() {
  if (mImpl->outputFile) { mImpl->outputFile->Close(); delete mImpl->outputFile; }
  if (mImpl->treeFile) { mImpl->treeFile->Close(); delete mImpl->treeFile; }
  delete mImpl;
}

StCandidateFemtoMaker* createStCandidateFemtoMaker(const char* name, StPicoDstMaker* picoMaker, const char* outName) {
  return new StCandidateFemtoMaker(name, picoMaker, outName);
}

void StCandidateFemtoMaker::SetEventSelector(FemtoEventSelector* selector) {
  delete mImpl->selector;
  mImpl->selector = selector;
}

void StCandidateFemtoMaker::AddCandidateSource(CandidateSource* source) {
  if (source) mImpl->sources.push_back(source);
}

Long64_t StCandidateFemtoMaker::GetAcceptedEvents() const { return mImpl->acceptedEvents; }

Bool_t StCandidateFemtoMaker::CandidateStoreInputFinished() const { return mImpl->storeFinished; }

Long64_t StCandidateFemtoMaker::GetSameEventPairs() const {
  Long64_t n = 0;
  for (size_t i = 0; i < mImpl->channels.size(); ++i) n += mImpl->channels[i].seCount;
  return n;
}

Long64_t StCandidateFemtoMaker::GetMixedEventPairs() const {
  Long64_t n = 0;
  for (size_t i = 0; i < mImpl->channels.size(); ++i) n += mImpl->channels[i].meCount;
  return n;
}

Int_t StCandidateFemtoMaker::Init() {
  StCandidateFemtoMakerImpl& m = *mImpl;
  const std::string mainConfig(mMainConfigPath.Data());
  std::cout << "[StCandidateFemtoMaker] mainconf " << mainConfig << std::endl;
  const Bool_t storeInput = mStoreList.Length() > 0;
  if (storeInput) {
    if (mPicoDstMaker || m.selector || !m.sources.empty()) {
      std::cerr << "ERROR: StCandidateFemtoMaker pair stage (candidate store input) takes no PicoDst maker, event selector "
                   "or candidate source: the stores provide them" << std::endl;
      return kStErr;
    }
  } else if (!mPicoDstMaker || !m.selector || m.sources.empty()) {
    std::cerr << "ERROR: StCandidateFemtoMaker requires a PicoDst maker, an event selector and >= 1 candidate source"
              << std::endl;
    return kStErr;
  }
  if (!m.cfg.LoadFromMainConfig(mainConfig)) return kStErr;
  if (storeInput && m.cfg.writeCandidateStore) {
    std::cerr << "ERROR: writeCandidateStore must be false in the pair stage (candidate store input)" << std::endl;
    return kStErr;
  }
  // The event-cut YAML fixes the mixing classes and the vertex radius: the pair stage checks it against the store.
  std::string eventCutsPath;
  if ((storeInput || m.cfg.writeCandidateStore) &&
      !FlatYamlConfig::ResolveMainConfigReference(mainConfig, "event", eventCutsPath))
    return kStErr;
  m.cfg.Print(std::cout);
  const std::string histPath = ConfigManager::GetInstance().GetHistConfigPath(GetName());
  if (histPath.empty()) { std::cerr << "ERROR: StCandidateFemtoMaker histogram config missing" << std::endl; return kStErr; }
  m.hist = new HistManager();
  if (!m.hist->LoadFromFile(histPath.c_str())) return kStErr;

  if (storeInput) {
    std::cout << "[StCandidateFemtoMaker] pair stage, candidate stores listed in " << mStoreList << std::endl;
    m.storeReader = new CandidateStoreReader();
    if (!m.storeReader->Open(mStoreList.Data(), eventCutsPath, m.sources)) return kStErr;
  } else if (!m.selector->Init(mainConfig.c_str(), m.hist)) {
    return kStErr;
  }
  for (size_t i = 0; i < m.sources.size(); ++i) {
    if (!m.sources[i]->Init(mainConfig.c_str(), m.hist)) {
      std::cerr << "ERROR: candidate source " << m.sources[i]->Name() << ": " << m.sources[i]->LastError() << std::endl;
      return kStErr;
    }
  }
  m.species.assign(m.cfg.species.size(), SpeciesState());
  for (size_t s = 0; s < m.cfg.species.size(); ++s) {
    const CandidateFemtoSpecies& sp = m.cfg.species[s];
    for (size_t i = 0; i < m.sources.size(); ++i)
      if (sp.source == m.sources[i]->Name()) m.species[s].source = m.sources[i];
    if (!m.species[s].source || !m.species[s].source->ProvidesSpecies(sp.particleKey)) {
      std::cerr << "ERROR: species " << sp.key << ": no source '" << sp.source << "' providing '" << sp.particleKey
                << "' (check the macro and the source YAML)" << std::endl;
      return kStErr;
    }
    m.species[s].daughterNames = m.species[s].source->DaughterNames(sp.particleKey);
    m.species[s].variableNames = m.species[s].source->VariableNames(sp.particleKey);
    if (!sp.daughterMasses.empty() && sp.daughterMasses.size() != m.species[s].daughterNames.size()) {
      std::cerr << "ERROR: species " << sp.key << ": daughterMasses needs " << m.species[s].daughterNames.size()
                << " values" << std::endl;
      return kStErr;
    }
  }
  if (!m.InitBdt(mainConfig)) return kStErr;
  if (m.cfg.writeCandidateStore)
    for (size_t s = 0; s < m.species.size(); ++s)
      if (m.species[s].bdt) {
        std::cerr << "ERROR: writeCandidateStore: no candidateBdt in the reconstruction stage (the pair stage scores)" << std::endl;
        return kStErr;
      }
  {
    // Mixing classes: Vz bin x cent9. Only "no centrality classes" (<= 1) or one class per cent9 bin (9) are
    // defined; any other number would clamp cent9 into the last class instead of grouping bins.
    const MixingConfig& mix = ConfigManager::GetInstance().GetMixingConfig();
    // Mixing classes need at least one Vz bin (with nVzBins 0 every centrality would share one pool).
    if (mix.nVzBins < 1) {
      std::cerr << "ERROR: mixing nVzBins = " << mix.nVzBins << ": StCandidateFemtoMaker needs >= 1" << std::endl;
      return kStErr;
    }
    if (mix.nCentralityBins > 1 && mix.nCentralityBins != 9) {
      std::cerr << "ERROR: mixing nCentralityBins = " << mix.nCentralityBins << ": StCandidateFemtoMaker supports 1 (none) or 9 (cent9)" << std::endl;
      return kStErr;
    }
    // ME uses every pair of the pools (no random subset): the mixing YAML must say so.
    if (!mix.IsBufferAllMode() || mix.bufferSize < 1) {
      std::cerr << "ERROR: mixing YAML: StCandidateFemtoMaker needs mixingMode: bufferAll and bufferSize >= 1 (got "
                << mix.mixingMode << ", " << mix.bufferSize << ")" << std::endl;
      return kStErr;
    }
    std::cout << "[StCandidateFemtoMaker] mixing: " << mix.nVzBins << " vz x " << mix.nCentralityBins
              << " centrality classes, " << mix.bufferSize << " events per pool, all pairs, both directions "
              << (mix.mixBothDirections ? "yes" : "no") << std::endl;
  }

  // Every output path is checked before anything is created, and a failure after creating the output removes it again, so
  // that a failed Init leaves no empty file that would block a rerun with the same path.
  TString treeName = mOutName;
  if (treeName.EndsWith(".root")) treeName.Remove(treeName.Length() - 5);
  treeName += ".tree.root";
  if (m.cfg.writeCandidateTrees && !gSystem->AccessPathName(treeName.Data())) {
    std::cerr << "ERROR: candidate tree output already exists (choose a new path): " << treeName << std::endl;
    return kStErr;
  }
  m.outputFile = TFile::Open(mOutName.Data(), "CREATE");
  if (!m.outputFile || m.outputFile->IsZombie()) {
    std::cerr << "ERROR: cannot create output (choose a new path): " << mOutName << std::endl;
    return kStErr;
  }
  m.channels.assign(m.cfg.channels.size(), ChannelState());
  for (size_t c = 0; c < m.cfg.channels.size() && !m.cfg.writeCandidateStore; ++c) {
    m.channels[c].se = m.BookPairHistogram(m.cfg.channels[c], "SE");
    m.channels[c].me = m.BookPairHistogram(m.cfg.channels[c], "ME");
  }
  if (m.cfg.writeCandidateStore) {
    m.storeWriter = new CandidateStoreWriter();
    if (!m.storeWriter->Init(m.outputFile, m.sources, eventCutsPath)) {
      m.outputFile->Close(); delete m.outputFile; m.outputFile = 0;
      gSystem->Unlink(mOutName.Data());
      return kStErr;
    }
    std::cout << "[StCandidateFemtoMaker] reconstruction stage: candidate store in " << mOutName << " (no pairs)" << std::endl;
  }
  if (m.cfg.writeCandidateTrees) {
    m.treeFile = TFile::Open(treeName.Data(), "CREATE");
    if (!m.treeFile || m.treeFile->IsZombie()) {
      std::cerr << "ERROR: cannot create candidate tree output: " << treeName << std::endl;
      delete m.treeFile; m.treeFile = 0;
      m.outputFile->Close(); delete m.outputFile; m.outputFile = 0;
      gSystem->Unlink(mOutName.Data());
      return kStErr;
    }
    for (size_t s = 0; s < m.species.size(); ++s) {
      SpeciesState& st = m.species[s];
      if (std::find(m.cfg.candidateTreeSpecies.begin(), m.cfg.candidateTreeSpecies.end(), m.cfg.species[s].key) ==
          m.cfg.candidateTreeSpecies.end())
        continue;
      m.treeFile->cd();
      st.tree = new TTree(("candidates_" + m.cfg.species[s].key).c_str(),
                          ("usable " + m.cfg.species[s].key + " candidates of accepted events").c_str());
      st.tree->Branch("runId", &m.runId, "runId/I");
      st.tree->Branch("eventId", &m.eventId, "eventId/I");
      st.tree->Branch("cent9", &m.cent9, "cent9/I");
      st.tree->Branch("mass", &st.treeMass, "mass/F");
      st.tree->Branch("px", &st.treePx, "px/F");
      st.tree->Branch("py", &st.treePy, "py/F");
      st.tree->Branch("pz", &st.treePz, "pz/F");
      st.tree->Branch("score", &st.treeScore, "score/F");
      st.treeVariables.assign(st.variableNames.size(), 0.f);
      for (size_t v = 0; v < st.variableNames.size(); ++v)
        st.tree->Branch(st.variableNames[v].c_str(), &st.treeVariables[v], (st.variableNames[v] + "/F").c_str());
      for (size_t c = 0; c < m.cfg.channels.size(); ++c)
        if (m.cfg.channels[c].indexA == static_cast<Int_t>(s) || m.cfg.channels[c].indexB == static_cast<Int_t>(s))
          st.treeChannels.push_back(c);
      st.treeChannelFlags.assign(st.treeChannels.size(), 0);
      for (size_t k = 0; k < st.treeChannels.size(); ++k) {
        const std::string b = "sharedTrackRemoved_" + m.cfg.channels[st.treeChannels[k]].name;
        st.tree->Branch(b.c_str(), &st.treeChannelFlags[k], (b + "/O").c_str());
      }
    }
    std::cout << "[StCandidateFemtoMaker] candidate trees: " << treeName << std::endl;
  }
  m.outputFile->cd();
  return kStOK;
}

void StCandidateFemtoMaker::Clear(Option_t* option) {
  StMaker::Clear(option);
  for (size_t i = 0; i < mImpl->sources.size(); ++i) mImpl->sources[i]->Clear();
}

Int_t StCandidateFemtoMaker::Make() {
  StCandidateFemtoMakerImpl& m = *mImpl;
  CandidateStoreEvent ev;
  if (m.storeReader) {
    // pair stage: the next stored event (accepted by the reconstruction stage) and its candidates
    m.picoDst = 0;
    if (!m.storeReader->Next(ev)) {
      if (!m.storeReader->Failed().empty()) {
        std::cerr << "ERROR: candidate store: " << m.storeReader->Failed() << std::endl;
        return kStErr;
      }
      m.storeFinished = kTRUE;
      return kStEOF;
    }
    ++m.readEvents;
  } else {
    if (!mPicoDstMaker) return kStErr;
    m.picoDst = mPicoDstMaker->picoDst();
    if (!m.picoDst) return kStWarn;
    StPicoEvent* event = m.picoDst->event();
    if (!event) return kStWarn;
    ++m.readEvents;
    const TVector3 pv = event->primaryVertex();
    // hVz: every read event; everything else: accepted events only (hN = events behind the pairs).
    m.hist->Fill("hVz", pv.Z());
    if (!m.selector->Accept(m.picoDst, m.cent9)) return kStOK;
    // The event record keeps the PicoDst Float_t / integer values (exact in the candidate store).
    ev.runId = event->runId();
    ev.eventId = event->eventId();
    ev.cent9 = m.cent9;
    ev.pvX = static_cast<Float_t>(pv.X()); ev.pvY = static_cast<Float_t>(pv.Y()); ev.pvZ = static_cast<Float_t>(pv.Z());
    ev.bField = event->bField();
    ev.fxtMult = event->fxtMult();
    ev.nBTOFMatch = event->nBTOFMatch();
    for (size_t i = 0; i < m.sources.size(); ++i) {
      if (!m.sources[i]->ProcessEvent(m.picoDst, m.cent9)) {
        std::cerr << "ERROR: candidate source " << m.sources[i]->Name() << ": " << m.sources[i]->LastError() << std::endl;
        return kStErr;
      }
    }
  }
  ++m.acceptedEvents;
  m.runId = ev.runId;
  m.eventId = ev.eventId;
  m.cent9 = ev.cent9;
  const TVector3 pv(ev.pvX, ev.pvY, ev.pvZ);
  if (m.storeWriter) {
    // reconstruction stage: every candidate of every part, with the PicoDst values of its daughter tracks
    m.storeWriter->BeginEvent(ev);
    for (size_t g = 0; g < m.storeWriter->NumberOfGroups(); ++g) {
      const std::vector<FemtoSourceCandidate>& cands = m.storeWriter->SourceOfGroup(g)->Candidates(m.storeWriter->ParticleKeyOfGroup(g));
      for (size_t i = 0; i < cands.size(); ++i) {
        FemtoSourceCandidate c = cands[i];
        std::vector<FemtoDaughterTrack> tracks;
        if (!m.DaughterTracks(c, tracks)) return kStErr;
        c.daughterTracks = tracks;
        if (!m.storeWriter->AddCandidate(g, c)) return kStErr;
      }
    }
    m.storeWriter->EndEvent();
  }
  if (m.picoDst) {
    // PicoDst-only QA (in the pair stage it comes from the stores)
    m.FillTrackPidQa(pv);
    m.hist->Fill("hEvtQA_MultVsNTof_all", ev.fxtMult, ev.nBTOFMatch);
  }
  if (m.storeWriter) return kStOK;
  const Float_t bT = 0.1 * ev.bField;  // kG -> T
  for (size_t s = 0; s < m.species.size(); ++s)
    if (!m.BuildEntries(s, pv, bT)) return kStErr;
  for (size_t c = 0; c < m.channels.size(); ++c) m.ApplySharedTrackRule(c);
  for (size_t s = 0; s < m.species.size(); ++s) m.FillSpecies(s);
  const EventCutConfig& evCuts = ConfigManager::GetInstance().GetEventCuts();
  const Double_t vr = TMath::Sqrt((pv.X() - evCuts.vtxCenterX) * (pv.X() - evCuts.vtxCenterX) +
                                  (pv.Y() - evCuts.vtxCenterY) * (pv.Y() - evCuts.vtxCenterY));
  const Int_t mixBin = m.GetMixingBin(pv.Z());
  for (size_t c = 0; c < m.channels.size(); ++c) {
    m.FillEventQa(c, m.channels[c].a, m.channels[c].b, ev.fxtMult, ev.nBTOFMatch, vr);
    m.ProcessChannel(c, mixBin);
  }
  m.hist->Fill("hN", 0.);
  return kStOK;
}

Int_t StCandidateFemtoMaker::Finish() {
  StCandidateFemtoMakerImpl& m = *mImpl;
  Int_t result = kStOK;
  if (m.selector) m.selector->Finish();
  for (size_t i = 0; i < m.sources.size(); ++i) m.sources[i]->Finish();
  if (m.outputFile && !m.outputFile->IsZombie()) {
    m.outputFile->cd();
    // pair stage: the reconstruction-stage histograms of the stores (event selection, parts, all-track QA, hVz)
    if (m.storeReader && m.hist && !m.storeReader->AddReconstructionHistograms(m.hist)) result = kStErr;
    if (m.hist) m.hist->Write();
    for (size_t c = 0; c < m.channels.size(); ++c) {
      if (!m.channels[c].se) continue;  // reconstruction stage: no pairs
      if (m.channels[c].se->Write() <= 0) result = kStErr;
      if (m.channels[c].me->Write() <= 0) result = kStErr;
    }
    if (m.storeWriter && !m.storeWriter->Write()) result = kStErr;
    std::ostringstream conf;
    m.cfg.Print(conf);
    const char* files[2] = {m.cfg.path.c_str(), m.bdtConfigPath.c_str()};
    for (Int_t f = 0; f < 2; ++f) {
      if (!*files[f]) continue;
      std::ifstream in(files[f]);
      conf << "\n# ---- " << files[f] << "\n" << in.rdbuf();
    }
    for (size_t s = 0; s < m.species.size(); ++s)
      if (m.species[s].bdt) conf << "\n# ---- BDT of " << m.cfg.species[s].key << "\n" << m.species[s].bdt->Model().Describe() << "\n";
    if (m.selector) conf << "\n# ---- event selection (" << m.selector->Name() << ")\n" << m.selector->Describe();
    if (m.storeWriter) conf << "\n# ---- candidate store written (reconstruction stage)\n" << m.storeWriter->LayoutText();
    if (m.storeReader)
      conf << "\n# ---- pair stage: " << m.storeReader->Events() << " events from the candidate stores listed in " << mStoreList
           << "\n" << m.storeReader->FirstLayout() << "\n# ---- reconstruction-stage configuration (first store)\n"
           << m.storeReader->FirstConfiguration();
    // Provenance from the batch joblist (the runtime bundle has no git repository); "unrecorded" if absent.
    const char* envKeys[3] = {"STAR_ANA_GIT_REV", "STAR_ANA_GIT_DIRTY", "STAR_ANA_LIBRARY_TAG"};
    conf << "\n# ---- provenance";
    for (Int_t e = 0; e < 3; ++e) {
      const char* v = gSystem->Getenv(envKeys[e]);
      conf << "\n" << envKeys[e] << ": " << ((v && *v) ? v : "unrecorded");
    }
    conf << "\n";
    TObjString confObj(conf.str().c_str());
    if (confObj.Write("CandidateFemtoConfiguration") <= 0) result = kStErr;
    m.outputFile->Flush();
    if (m.outputFile->TestBit(TFile::kWriteError)) result = kStErr;
    m.outputFile->Close();
    if (m.outputFile->TestBit(TFile::kWriteError)) result = kStErr;
    delete m.outputFile;
    m.outputFile = 0;
  }
  if (m.treeFile && !m.treeFile->IsZombie()) {
    m.treeFile->cd();
    for (size_t s = 0; s < m.species.size(); ++s)
      if (m.species[s].tree && m.species[s].tree->Write("", TObject::kOverwrite) <= 0) result = kStErr;
    m.treeFile->Close();
    delete m.treeFile;
    m.treeFile = 0;
  }
  std::cout << "[StCandidateFemtoMaker] readEvents=" << m.readEvents << " acceptedEvents=" << m.acceptedEvents
            << (m.storeWriter ? TString::Format(" storedEvents=%lld", m.storeWriter->Entries()).Data() : "")
            << (m.storeReader ? " (pair stage: events read from the candidate stores)" : "") << std::endl;
  for (size_t s = 0; s < m.species.size(); ++s) {
    const SpeciesState& st = m.species[s];
    std::cout << "[StCandidateFemtoMaker] species " << m.cfg.species[s].key << ": usable=" << st.usable
              << " signal=" << st.signal
              << (st.bdt ? TString::Format(" outsideTrainingDomain=%lld", st.bdt->OutsideCount()).Data() : "")
              << std::endl;
  }
  for (size_t c = 0; c < m.channels.size(); ++c) {
    const ChannelState& cs = m.channels[c];
    const CandidateFemtoChannel& ch = m.cfg.channels[c];
    std::cout << "[StCandidateFemtoMaker] channel " << ch.name << ": eventsWithBoth=" << cs.eventsWithBoth
              << " sePairs(windows)=" << cs.seCount << " mePairs(windows)=" << cs.meCount
              << " ME/SE=" << (cs.seCount > 0 ? TString::Format("%.1f", Double_t(cs.meCount) / cs.seCount).Data() : "n/a")
              << " shareRejectSE=" << cs.shareRejectSE << " sharedTrackRemoved(signal/leftSB/rightSB/other) "
              << ch.partA << "=" << cs.removedA[1] << "/" << cs.removedA[2] << "/" << cs.removedA[3] << "/" << cs.removedA[0]
              << " " << ch.partB << "=" << cs.removedB[1] << "/" << cs.removedB[2] << "/" << cs.removedB[3] << "/"
              << cs.removedB[0] << std::endl;
  }
  return result;
}
