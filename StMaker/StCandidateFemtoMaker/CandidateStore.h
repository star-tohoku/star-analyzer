#ifndef CANDIDATE_STORE_H
#define CANDIDATE_STORE_H

// Candidate store: the output of the reconstruction stage of StCandidateFemtoMaker (one PicoDst pass) and the
// input of its pair stage (no PicoDst). Header-only, compiled into libStCandidateFemtoMaker.
//
// One ROOT file per reconstruction job holds
//   - TTree "candidateEvents": one entry per accepted event (event record + every candidate of every species
//     of every reconstruction part, before any femto selection), in processing order;
//   - TObjString "CandidateStoreLayout": the parts, species keys, daughter and variable names (the branch layout);
//   - TObjString "CandidateStoreEventCuts": the event-cut YAML of the reconstruction job (the pair stage uses the
//     same values for the mixing classes and the vertex radius and checks that they agree, ignoring comments);
//   - the reconstruction-stage histograms (event selection, reconstruction parts, all-track QA, hVz): the pair
//     stage adds them to its own output, so its output has the same content as a one-pass job.
// Nothing here knows any particle: species, daughters and variables come from the parts (CandidateSource).
// Every stored value is the Float_t / Int_t the parts and the PicoDst give, so the pair stage reproduces a
// one-pass job exactly.

#include "CandidateSource.h"
#include "HistManager.h"

#include "TClass.h"
#include "TFile.h"
#include "TH1.h"
#include "TKey.h"
#include "TObjString.h"
#include "TTree.h"

#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

static const char* const kCandidateStoreTreeName = "candidateEvents";
static const char* const kCandidateStoreLayoutName = "CandidateStoreLayout";
static const char* const kCandidateStoreEventCutsName = "CandidateStoreEventCuts";
static const char* const kCandidateStoreLayoutVersion = "candidate store layout 1";

// Event record: everything the pair stage needs from the PicoDst event.
struct CandidateStoreEvent {
  CandidateStoreEvent() : runId(0), eventId(0), cent9(-1), pvX(0.f), pvY(0.f), pvZ(0.f), bField(0.f), fxtMult(0), nBTOFMatch(0) {}
  Int_t runId, eventId, cent9;
  Float_t pvX, pvY, pvZ;  // primary vertex [cm]
  Float_t bField;         // [kG]
  Int_t fxtMult, nBTOFMatch;
};

// Branches of one species of one part. Prefix "<part>_<particleKey>_".
struct CandidateStoreGroup {
  CandidateStoreGroup()
      : px(new std::vector<Float_t>), py(new std::vector<Float_t>), pz(new std::vector<Float_t>), mass(new std::vector<Float_t>),
        x(new std::vector<Float_t>), y(new std::vector<Float_t>), z(new std::vector<Float_t>), decayLength(new std::vector<Float_t>),
        variables(new std::vector<Float_t>), daughterTrack(new std::vector<Int_t>), daughterGx(new std::vector<Float_t>),
        daughterGy(new std::vector<Float_t>), daughterGz(new std::vector<Float_t>), daughterCharge(new std::vector<Int_t>) {}
  std::string source, particleKey;
  std::vector<std::string> daughterNames, variableNames;
  // per candidate
  std::vector<Float_t>* px; std::vector<Float_t>* py; std::vector<Float_t>* pz; std::vector<Float_t>* mass;
  std::vector<Float_t>* x; std::vector<Float_t>* y; std::vector<Float_t>* z; std::vector<Float_t>* decayLength;
  // per candidate x variable / x daughter, flattened
  std::vector<Float_t>* variables;
  std::vector<Int_t>* daughterTrack;
  std::vector<Float_t>* daughterGx; std::vector<Float_t>* daughterGy; std::vector<Float_t>* daughterGz;
  std::vector<Int_t>* daughterCharge;

  std::string Prefix() const { return source + "_" + particleKey + "_"; }
  void ClearValues() {
    px->clear(); py->clear(); pz->clear(); mass->clear(); x->clear(); y->clear(); z->clear(); decayLength->clear();
    variables->clear(); daughterTrack->clear(); daughterGx->clear(); daughterGy->clear(); daughterGz->clear(); daughterCharge->clear();
  }
  void DeleteBuffers() {
    delete px; delete py; delete pz; delete mass; delete x; delete y; delete z; delete decayLength;
    delete variables; delete daughterTrack; delete daughterGx; delete daughterGy; delete daughterGz; delete daughterCharge;
  }
  // The vector branches of this group, in a fixed order (writer and reader use the same list).
  void FloatBranches(std::vector<std::pair<std::string, std::vector<Float_t>**> >& f) {
    const std::string p = Prefix();
    f.push_back(std::make_pair(p + "px", &px)); f.push_back(std::make_pair(p + "py", &py));
    f.push_back(std::make_pair(p + "pz", &pz)); f.push_back(std::make_pair(p + "mass", &mass));
    f.push_back(std::make_pair(p + "x", &x)); f.push_back(std::make_pair(p + "y", &y)); f.push_back(std::make_pair(p + "z", &z));
    f.push_back(std::make_pair(p + "decayLength", &decayLength)); f.push_back(std::make_pair(p + "variables", &variables));
    f.push_back(std::make_pair(p + "daughterGx", &daughterGx)); f.push_back(std::make_pair(p + "daughterGy", &daughterGy));
    f.push_back(std::make_pair(p + "daughterGz", &daughterGz));
  }
  void IntBranches(std::vector<std::pair<std::string, std::vector<Int_t>**> >& f) {
    const std::string p = Prefix();
    f.push_back(std::make_pair(p + "daughterTrack", &daughterTrack));
    f.push_back(std::make_pair(p + "daughterCharge", &daughterCharge));
  }
};

inline std::string CandidateStoreJoin(const std::vector<std::string>& v) {
  std::string s;
  for (size_t i = 0; i < v.size(); ++i) s += (i ? "," : "") + v[i];
  return s;
}

inline std::vector<std::string> CandidateStoreSplit(const std::string& s) {
  std::vector<std::string> v;
  std::string cur;
  for (size_t i = 0; i < s.size(); ++i) {
    if (s[i] == ',') { v.push_back(cur); cur.clear(); } else cur += s[i];
  }
  if (!s.empty()) v.push_back(cur);
  return v;
}

// Text form of the layout: "candidate store layout 1", then per group
//   group <part> <particleKey>
//   daughters <comma-separated names>
//   variables <comma-separated names>
inline std::string CandidateStoreLayoutText(const std::vector<CandidateStoreGroup>& groups) {
  std::ostringstream o;
  o << kCandidateStoreLayoutVersion << "\n";
  for (size_t g = 0; g < groups.size(); ++g)
    o << "group " << groups[g].source << " " << groups[g].particleKey << "\n"
      << "daughters " << CandidateStoreJoin(groups[g].daughterNames) << "\n"
      << "variables " << CandidateStoreJoin(groups[g].variableNames) << "\n";
  return o.str();
}

inline Bool_t CandidateStoreParseLayout(const std::string& text, std::vector<CandidateStoreGroup>& groups, std::string& error) {
  std::istringstream in(text);
  std::string line;
  if (!std::getline(in, line) || line != kCandidateStoreLayoutVersion) { error = "unknown layout version '" + line + "'"; return kFALSE; }
  while (std::getline(in, line)) {
    if (line.empty()) continue;
    std::istringstream f(line);
    std::string word, a, b;
    f >> word >> a >> b;
    if (word != "group" || a.empty() || b.empty()) { error = "bad layout line '" + line + "'"; return kFALSE; }
    CandidateStoreGroup g;
    g.source = a; g.particleKey = b;
    std::string d, v;
    if (!std::getline(in, d) || d.compare(0, 10, "daughters ") != 0 || !std::getline(in, v) || v.compare(0, 10, "variables ") != 0) {
      error = "group " + a + " " + b + ": daughters/variables lines missing";
      return kFALSE;
    }
    g.daughterNames = CandidateStoreSplit(d.substr(10));
    g.variableNames = CandidateStoreSplit(v.substr(10));
    groups.push_back(g);
  }
  if (groups.empty()) { error = "no species in the layout"; return kFALSE; }
  return kTRUE;
}

inline std::string CandidateStoreReadFile(const std::string& path) {
  std::ifstream in(path.c_str());
  std::ostringstream o;
  o << in.rdbuf();
  return o.str();
}

// The content of a YAML text without comments and blank lines: two event-cut files that differ only in comments or
// layout select the same events, so stores are compared on this, not on the raw text. A '#' starts a comment at the
// beginning of a line or after a space or tab, outside quotes.
inline std::string CandidateStoreYamlContent(const std::string& text) {
  std::istringstream in(text);
  std::string line, out;
  while (std::getline(in, line)) {
    char quote = 0;
    size_t cut = line.size();
    for (size_t i = 0; i < line.size(); ++i) {
      const char c = line[i];
      if (quote) { if (c == quote) quote = 0; continue; }
      if (c == '"' || c == '\'') { quote = c; continue; }
      if (c == '#' && (i == 0 || line[i - 1] == ' ' || line[i - 1] == '\t')) { cut = i; break; }
    }
    line.erase(cut);
    while (!line.empty() && (line[line.size() - 1] == ' ' || line[line.size() - 1] == '\t' || line[line.size() - 1] == '\r'))
      line.erase(line.size() - 1);
    if (line.find_first_not_of(" \t") == std::string::npos) continue;
    out += line;
    out += '\n';
  }
  return out;
}

// Reconstruction stage: one store per job, written into the job output file.
class CandidateStoreWriter {
public:
  CandidateStoreWriter() : mTree(0), mEntries(0) {}
  ~CandidateStoreWriter() { for (size_t g = 0; g < mGroups.size(); ++g) mGroups[g].DeleteBuffers(); }

  // Every species of every part. file: the job output (current directory for the tree).
  Bool_t Init(TFile* file, const std::vector<CandidateSource*>& sources, const std::string& eventCutsPath) {
    mEventCuts = CandidateStoreReadFile(eventCutsPath);
    if (mEventCuts.empty()) { std::cerr << "ERROR: candidate store: cannot read the event-cut YAML " << eventCutsPath << std::endl; return kFALSE; }
    for (size_t s = 0; s < sources.size(); ++s) {
      const std::vector<std::string> keys = sources[s]->SpeciesKeys();
      for (size_t k = 0; k < keys.size(); ++k) {
        CandidateStoreGroup g;
        g.source = sources[s]->Name();
        g.particleKey = keys[k];
        g.daughterNames = sources[s]->DaughterNames(keys[k]);
        g.variableNames = sources[s]->VariableNames(keys[k]);
        if (g.source.find_first_of(" ,\n") != std::string::npos || g.particleKey.find_first_of(" ,\n") != std::string::npos) {
          std::cerr << "ERROR: candidate store: part/species names must not contain blanks or commas" << std::endl;
          return kFALSE;
        }
        mGroups.push_back(g);
        mSourceOfGroup.push_back(sources[s]);
      }
    }
    file->cd();
    mTree = new TTree(kCandidateStoreTreeName, "accepted events with every reconstructed candidate (candidate store)");
    mTree->Branch("runId", &mEvent.runId, "runId/I");
    mTree->Branch("eventId", &mEvent.eventId, "eventId/I");
    mTree->Branch("cent9", &mEvent.cent9, "cent9/I");
    mTree->Branch("pvX", &mEvent.pvX, "pvX/F");
    mTree->Branch("pvY", &mEvent.pvY, "pvY/F");
    mTree->Branch("pvZ", &mEvent.pvZ, "pvZ/F");
    mTree->Branch("bField", &mEvent.bField, "bField/F");
    mTree->Branch("fxtMult", &mEvent.fxtMult, "fxtMult/I");
    mTree->Branch("nBTOFMatch", &mEvent.nBTOFMatch, "nBTOFMatch/I");
    for (size_t g = 0; g < mGroups.size(); ++g) {
      std::vector<std::pair<std::string, std::vector<Float_t>**> > fb;
      std::vector<std::pair<std::string, std::vector<Int_t>**> > ib;
      mGroups[g].FloatBranches(fb);
      mGroups[g].IntBranches(ib);
      for (size_t i = 0; i < fb.size(); ++i) mTree->Branch(fb[i].first.c_str(), fb[i].second);
      for (size_t i = 0; i < ib.size(); ++i) mTree->Branch(ib[i].first.c_str(), ib[i].second);
    }
    std::cout << "[CandidateStore] writing " << mGroups.size() << " species:";
    for (size_t g = 0; g < mGroups.size(); ++g) std::cout << " " << mGroups[g].source << "/" << mGroups[g].particleKey;
    std::cout << std::endl;
    return kTRUE;
  }

  size_t NumberOfGroups() const { return mGroups.size(); }
  CandidateSource* SourceOfGroup(size_t g) const { return mSourceOfGroup[g]; }
  const std::string& ParticleKeyOfGroup(size_t g) const { return mGroups[g].particleKey; }

  void BeginEvent(const CandidateStoreEvent& ev) {
    mEvent = ev;
    for (size_t g = 0; g < mGroups.size(); ++g) mGroups[g].ClearValues();
  }

  // c.daughterTracks must be filled (one per daughter).
  Bool_t AddCandidate(size_t g, const FemtoSourceCandidate& c) {
    CandidateStoreGroup& s = mGroups[g];
    if (c.variables.size() != s.variableNames.size() || c.daughters.size() != s.daughterNames.size() ||
        c.daughterTracks.size() != c.daughters.size()) {
      std::cerr << "ERROR: candidate store: " << s.source << "/" << s.particleKey << " candidate has " << c.variables.size()
                << " variables, " << c.daughters.size() << " daughters, " << c.daughterTracks.size() << " daughter tracks (layout: "
                << s.variableNames.size() << " variables, " << s.daughterNames.size() << " daughters)" << std::endl;
      return kFALSE;
    }
    s.px->push_back(c.px); s.py->push_back(c.py); s.pz->push_back(c.pz); s.mass->push_back(c.massValue);
    s.x->push_back(c.x); s.y->push_back(c.y); s.z->push_back(c.z); s.decayLength->push_back(c.decayLength);
    s.variables->insert(s.variables->end(), c.variables.begin(), c.variables.end());
    for (size_t d = 0; d < c.daughters.size(); ++d) {
      s.daughterTrack->push_back(c.daughters[d]);
      s.daughterGx->push_back(c.daughterTracks[d].gx); s.daughterGy->push_back(c.daughterTracks[d].gy);
      s.daughterGz->push_back(c.daughterTracks[d].gz); s.daughterCharge->push_back(c.daughterTracks[d].charge);
    }
    return kTRUE;
  }

  void EndEvent() { mTree->Fill(); }

  // Into the job output file (current directory).
  // The tree belongs to the output file: Entries() uses the count cached here, since the file (and the tree with
  // it) is closed and deleted right after (GetEntries() on the deleted tree in Finish).
  Bool_t Write() {
    mEntries = mTree->GetEntries();
    TObjString layout(CandidateStoreLayoutText(mGroups).c_str());
    TObjString cuts(mEventCuts.c_str());
    return mTree->Write("", TObject::kOverwrite) > 0 && layout.Write(kCandidateStoreLayoutName) > 0 &&
           cuts.Write(kCandidateStoreEventCutsName) > 0;
  }
  Long64_t Entries() const { return mEntries; }  // valid after Write()
  std::string LayoutText() const { return CandidateStoreLayoutText(mGroups); }

private:
  TTree* mTree;
  Long64_t mEntries;
  CandidateStoreEvent mEvent;
  std::vector<CandidateStoreGroup> mGroups;
  std::vector<CandidateSource*> mSourceOfGroup;
  std::string mEventCuts;
};

// A reconstruction part replayed from the store: same name, species, daughters and variables as the part that
// wrote it. ProcessEvent does nothing: CandidateStoreReader::Next fills the candidates of the current event.
class StoredCandidateSource : public CandidateSource {
public:
  explicit StoredCandidateSource(const std::string& name) : mName(name) {}
  virtual const char* Name() const { return mName.c_str(); }
  virtual Bool_t Init(const char*, HistManager*) { return kTRUE; }
  virtual std::vector<std::string> SpeciesKeys() const { return mKeys; }
  virtual Bool_t ProvidesSpecies(const std::string& key) const { return mDaughters.count(key) > 0; }
  virtual std::vector<std::string> DaughterNames(const std::string& key) const {
    std::map<std::string, std::vector<std::string> >::const_iterator it = mDaughters.find(key);
    return it != mDaughters.end() ? it->second : std::vector<std::string>();
  }
  virtual std::vector<std::string> VariableNames(const std::string& key) const {
    std::map<std::string, std::vector<std::string> >::const_iterator it = mVariables.find(key);
    return it != mVariables.end() ? it->second : std::vector<std::string>();
  }
  virtual Bool_t ProcessEvent(StPicoDst*, Int_t) { return kTRUE; }
  virtual const std::vector<FemtoSourceCandidate>& Candidates(const std::string& key) const {
    std::map<std::string, std::vector<FemtoSourceCandidate> >::const_iterator it = mCandidates.find(key);
    return it != mCandidates.end() ? it->second : mEmpty;
  }
  virtual void Clear() {
    for (std::map<std::string, std::vector<FemtoSourceCandidate> >::iterator it = mCandidates.begin(); it != mCandidates.end(); ++it)
      it->second.clear();
  }
  virtual const std::string& LastError() const { return mError; }

  void AddSpecies(const CandidateStoreGroup& g) {
    mKeys.push_back(g.particleKey);
    mDaughters[g.particleKey] = g.daughterNames;
    mVariables[g.particleKey] = g.variableNames;
    mCandidates[g.particleKey];
  }
  std::vector<FemtoSourceCandidate>& MutableCandidates(const std::string& key) { return mCandidates[key]; }

private:
  std::string mName, mError;
  std::vector<std::string> mKeys;
  std::map<std::string, std::vector<std::string> > mDaughters, mVariables;
  std::map<std::string, std::vector<FemtoSourceCandidate> > mCandidates;
  std::vector<FemtoSourceCandidate> mEmpty;
};

// Pair stage: reads the stores of a list, file after file, entry after entry (the order of the reconstruction
// jobs), so mixing pools continue across files exactly as in one job over the same PicoDst files.
class CandidateStoreReader {
public:
  CandidateStoreReader() : mFile(0), mTree(0), mFileIndex(-1), mEntry(0), mEntries(0), mEvents(0) {}
  ~CandidateStoreReader() {
    CloseFile();
    for (size_t g = 0; g < mGroups.size(); ++g) mGroups[g].DeleteBuffers();
    for (std::map<std::string, TH1*>::iterator it = mHistSums.begin(); it != mHistSums.end(); ++it) delete it->second;
  }

  // list: one store .root file, or a text file with one store path per line ('#' comments). eventCutsPath: event-cut YAML of this job, must
  // equal the one stored. Creates one StoredCandidateSource per stored part (ownership passes to the caller).
  Bool_t Open(const std::string& list, const std::string& eventCutsPath, std::vector<CandidateSource*>& sources) {
    std::ifstream in(list.c_str());
    if (!in) { std::cerr << "ERROR: candidate store: cannot open the list " << list << std::endl; return kFALSE; }
    std::string line;
    if (list.size() > 5 && list.compare(list.size() - 5, 5, ".root") == 0) mFiles.push_back(list);  // one store file
    else while (std::getline(in, line)) {
      const std::string::size_type c = line.find('#');
      if (c != std::string::npos) line.erase(c);
      std::istringstream f(line);
      std::string path;
      if (f >> path) mFiles.push_back(path);
    }
    if (mFiles.empty()) { std::cerr << "ERROR: candidate store: the list " << list << " is empty" << std::endl; return kFALSE; }
    mEventCuts = CandidateStoreReadFile(eventCutsPath);
    if (mEventCuts.empty()) { std::cerr << "ERROR: candidate store: cannot read the event-cut YAML " << eventCutsPath << std::endl; return kFALSE; }
    if (!OpenFile(0)) return kFALSE;
    std::map<std::string, StoredCandidateSource*> byName;
    for (size_t g = 0; g < mGroups.size(); ++g) {
      StoredCandidateSource*& s = byName[mGroups[g].source];
      if (!s) { s = new StoredCandidateSource(mGroups[g].source); sources.push_back(s); }
      s->AddSpecies(mGroups[g]);
      mSourceOfGroup.push_back(s);
    }
    std::cout << "[CandidateStore] reading " << mFiles.size() << " store file(s), species:";
    for (size_t g = 0; g < mGroups.size(); ++g) std::cout << " " << mGroups[g].source << "/" << mGroups[g].particleKey;
    std::cout << std::endl;
    return kTRUE;
  }

  // Next stored event: fills ev and the candidates of the stored parts. kFALSE at the end or on error (see Failed).
  Bool_t Next(CandidateStoreEvent& ev) {
    while (mEntry >= mEntries) {
      if (mFileIndex + 1 >= static_cast<Int_t>(mFiles.size())) { CloseFile(); return kFALSE; }
      if (!OpenFile(mFileIndex + 1)) { mFailed = "cannot use " + mFiles[mFileIndex + 1]; return kFALSE; }
    }
    if (mTree->GetEntry(mEntry++) <= 0) { mFailed = "cannot read an entry of " + mFiles[mFileIndex]; return kFALSE; }
    ev = mEvent;
    for (size_t g = 0; g < mGroups.size(); ++g) {
      const CandidateStoreGroup& s = mGroups[g];
      std::vector<FemtoSourceCandidate>& out = mSourceOfGroup[g]->MutableCandidates(s.particleKey);
      out.clear();
      const size_t n = s.px->size(), nv = s.variableNames.size(), nd = s.daughterNames.size();
      if (s.py->size() != n || s.pz->size() != n || s.mass->size() != n || s.x->size() != n || s.y->size() != n ||
          s.z->size() != n || s.decayLength->size() != n || s.variables->size() != n * nv || s.daughterTrack->size() != n * nd ||
          s.daughterGx->size() != n * nd || s.daughterGy->size() != n * nd || s.daughterGz->size() != n * nd ||
          s.daughterCharge->size() != n * nd) {
        mFailed = "inconsistent vector sizes of " + s.Prefix() + " in " + mFiles[mFileIndex];
        return kFALSE;
      }
      for (size_t i = 0; i < n; ++i) {
        FemtoSourceCandidate c;
        c.px = (*s.px)[i]; c.py = (*s.py)[i]; c.pz = (*s.pz)[i]; c.massValue = (*s.mass)[i];
        c.x = (*s.x)[i]; c.y = (*s.y)[i]; c.z = (*s.z)[i]; c.decayLength = (*s.decayLength)[i];
        c.variables.assign(s.variables->begin() + i * nv, s.variables->begin() + (i + 1) * nv);
        for (size_t d = 0; d < nd; ++d) {
          const size_t k = i * nd + d;
          c.daughters.push_back((*s.daughterTrack)[k]);
          FemtoDaughterTrack t;
          t.gx = (*s.daughterGx)[k]; t.gy = (*s.daughterGy)[k]; t.gz = (*s.daughterGz)[k]; t.charge = (*s.daughterCharge)[k];
          c.daughterTracks.push_back(t);
        }
        out.push_back(c);
      }
    }
    ++mEvents;
    return kTRUE;
  }

  const std::string& Failed() const { return mFailed; }
  Long64_t Events() const { return mEvents; }
  const std::string& FirstLayout() const { return mLayout; }
  const std::string& FirstConfiguration() const { return mFirstConfiguration; }

  // Add the reconstruction-stage histograms of every store read so far into the pair-stage HistManager. The two
  // stages fill disjoint histograms; one filled by both stops the job.
  Bool_t AddReconstructionHistograms(HistManager* hist) {
    for (std::map<std::string, TH1*>::iterator it = mHistSums.begin(); it != mHistSums.end(); ++it) {
      TH1* mine = hist->Get(it->first.c_str());
      if (!mine) { std::cerr << "ERROR: candidate store: histogram " << it->first << " is not in the hist YAML of this job" << std::endl; return kFALSE; }
      if (mine->GetEntries() > 0 && it->second->GetEntries() > 0) {
        std::cerr << "ERROR: candidate store: histogram " << it->first << " was filled in both stages" << std::endl;
        return kFALSE;
      }
      // bin labels set by the reconstruction parts in their Init (e.g. stage counters) come with the histogram
      TAxis* to[3] = {mine->GetXaxis(), mine->GetYaxis(), mine->GetZaxis()};
      TAxis* from[3] = {it->second->GetXaxis(), it->second->GetYaxis(), it->second->GetZaxis()};
      for (Int_t a = 0; a < 3; ++a)
        if (from[a]->GetLabels() && from[a]->GetNbins() == to[a]->GetNbins())
          for (Int_t b = 1; b <= from[a]->GetNbins(); ++b) to[a]->SetBinLabel(b, from[a]->GetBinLabel(b));
      mine->Add(it->second);
    }
    return kTRUE;
  }

private:
  std::vector<std::string> mFiles;
  TFile* mFile;
  TTree* mTree;
  Int_t mFileIndex;
  Long64_t mEntry, mEntries, mEvents;
  CandidateStoreEvent mEvent;
  std::vector<CandidateStoreGroup> mGroups;
  std::vector<StoredCandidateSource*> mSourceOfGroup;
  std::string mLayout, mEventCuts, mFailed, mFirstConfiguration;
  std::map<std::string, TH1*> mHistSums;

  void CloseFile() {
    if (mFile) { mFile->Close(); delete mFile; }
    mFile = 0; mTree = 0;
  }

  Bool_t OpenFile(Int_t index) {
    CloseFile();
    mFileIndex = index;
    mEntry = mEntries = 0;
    const std::string& path = mFiles[index];
    mFile = TFile::Open(path.c_str(), "READ");
    if (!mFile || mFile->IsZombie()) { std::cerr << "ERROR: candidate store: cannot open " << path << std::endl; return kFALSE; }
    TObjString* layout = dynamic_cast<TObjString*>(mFile->Get(kCandidateStoreLayoutName));
    TObjString* cuts = dynamic_cast<TObjString*>(mFile->Get(kCandidateStoreEventCutsName));
    mTree = dynamic_cast<TTree*>(mFile->Get(kCandidateStoreTreeName));
    if (!layout || !cuts || !mTree) { std::cerr << "ERROR: candidate store: " << path << " is not a candidate store" << std::endl; return kFALSE; }
    if (CandidateStoreYamlContent(cuts->GetString().Data()) != CandidateStoreYamlContent(mEventCuts)) {
      std::cerr << "ERROR: candidate store: the event cuts of " << path << " differ from the event-cut YAML of this job"
                << " (comments and blank lines are ignored)" << std::endl;
      return kFALSE;
    }
    if (index == 0) {
      mLayout = layout->GetString().Data();
      TObjString* conf = dynamic_cast<TObjString*>(mFile->Get("CandidateFemtoConfiguration"));
      mFirstConfiguration = conf ? conf->GetString().Data() : "(not recorded)";
      std::string error;
      if (!CandidateStoreParseLayout(mLayout, mGroups, error)) { std::cerr << "ERROR: candidate store: " << path << ": " << error << std::endl; return kFALSE; }
    } else if (mLayout != layout->GetString().Data()) {
      std::cerr << "ERROR: candidate store: the layout of " << path << " differs from the first file" << std::endl;
      return kFALSE;
    }
    mTree->SetBranchAddress("runId", &mEvent.runId);
    mTree->SetBranchAddress("eventId", &mEvent.eventId);
    mTree->SetBranchAddress("cent9", &mEvent.cent9);
    mTree->SetBranchAddress("pvX", &mEvent.pvX);
    mTree->SetBranchAddress("pvY", &mEvent.pvY);
    mTree->SetBranchAddress("pvZ", &mEvent.pvZ);
    mTree->SetBranchAddress("bField", &mEvent.bField);
    mTree->SetBranchAddress("fxtMult", &mEvent.fxtMult);
    mTree->SetBranchAddress("nBTOFMatch", &mEvent.nBTOFMatch);
    for (size_t g = 0; g < mGroups.size(); ++g) {
      std::vector<std::pair<std::string, std::vector<Float_t>**> > fb;
      std::vector<std::pair<std::string, std::vector<Int_t>**> > ib;
      mGroups[g].FloatBranches(fb);
      mGroups[g].IntBranches(ib);
      for (size_t i = 0; i < fb.size(); ++i)
        if (mTree->SetBranchAddress(fb[i].first.c_str(), fb[i].second) < 0) { std::cerr << "ERROR: candidate store: no branch " << fb[i].first << " in " << path << std::endl; return kFALSE; }
      for (size_t i = 0; i < ib.size(); ++i)
        if (mTree->SetBranchAddress(ib[i].first.c_str(), ib[i].second) < 0) { std::cerr << "ERROR: candidate store: no branch " << ib[i].first << " in " << path << std::endl; return kFALSE; }
    }
    mEntries = mTree->GetEntries();
    // reconstruction-stage histograms of this file
    TIter next(mFile->GetListOfKeys());
    while (TKey* key = static_cast<TKey*>(next())) {
      TClass* cl = TClass::GetClass(key->GetClassName());
      if (!cl || !cl->InheritsFrom(TH1::Class())) continue;
      TH1* h = static_cast<TH1*>(key->ReadObj());
      if (!h) continue;
      std::map<std::string, TH1*>::iterator it = mHistSums.find(key->GetName());
      if (it == mHistSums.end()) {
        h->SetDirectory(0);
        mHistSums[key->GetName()] = h;
      } else {
        it->second->Add(h);
        delete h;
      }
    }
    return kTRUE;
  }
};

#endif
