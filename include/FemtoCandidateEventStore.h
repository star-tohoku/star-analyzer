#ifndef FEMTO_CANDIDATE_EVENT_STORE_H
#define FEMTO_CANDIDATE_EVENT_STORE_H

#include "FemtoCandidate.h"
#include "FemtoNuclearPidPairQaCandidate.h"
#include "FemtoEpdEventPlane.h"
#include "TBranch.h"
#include "TFile.h"
#include "TLeaf.h"
#include "TList.h"
#include "TMath.h"
#include "TNamed.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TTree.h"
#include "TUUID.h"
#include <cctype>
#include <fstream>
#include <limits>
#include <set>
#include <sstream>
#include <string>
#include <vector>

// Exact persistence of the existing StFemtoMaker candidate contract. This layout
// keeps the stored energy and derived Float_t values; replay never calls SetP4.
// The file owns its tree. Write must finish before that file is closed.
// Replay original per-job stores: repacking/merging changes the file UUID and
// discards original mixing-pool boundaries, so a reader rejects those files.
namespace femto_candidate_store {
static const char* const kTreeName = "femtoCandidateEvents";
static const char* const kLayoutName = "FemtoCandidateEventStoreLayout";
static const char* const kSelectionName = "FemtoCandidateEventStoreSelection";
static const char* const kMixingName = "FemtoCandidateEventStoreMixing";
static const char* const kFileUuidName = "FemtoCandidateEventStoreFileUUID";
static const char* const kProducerStatusName = "FemtoLambdaRunStatus";
static const char* const kSchemaVersion = "femto candidate event store 1";
static const char* const kEpdSchemaVersion = "femto candidate event store 2 epd";
enum { kFloatFields = 20, kIntFields = 10, kPhysical = 0, kPidQa = 1 };

struct Event {
  Int_t inputIndex, runId, eventId, cent9, cent16, refMult, fxtMult, nTracks, nBTOFMatch;
  Double_t pvX, pvY, pvZ, bField, refMultCorr, centWeight, psi2;
  Bool_t eventPlaneValid;
  Double_t eventPlane2, q2x, q2y;
  Int_t eventPlane2Tracks;
  Bool_t eventPlane2Valid;
  Double_t epdPsi1, epdRawQx, epdRawQy, epdRecenteredQx, epdRecenteredQy, epdWeight;
  Int_t epdStatus, epdBin, epdRunIndex, epdCalibrationCent, epdHits;
  Bool_t epdValid;
  Event() : inputIndex(-1), runId(0), eventId(0), cent9(-1), cent16(-1),
      refMult(0), fxtMult(0), nTracks(0), nBTOFMatch(0), pvX(0), pvY(0), pvZ(0),
      bField(0), refMultCorr(0), centWeight(1), psi2(-1), eventPlaneValid(kFALSE),
      eventPlane2(-1), q2x(0), q2y(0), eventPlane2Tracks(0), eventPlane2Valid(kFALSE), epdPsi1(-999), epdRawQx(0), epdRawQy(0),
      epdRecenteredQx(0), epdRecenteredQy(0), epdWeight(0), epdStatus(0), epdBin(-1),
      epdRunIndex(-1), epdCalibrationCent(-1), epdHits(0), epdValid(kFALSE) {}
};

inline bool ValidSpecies(const std::vector<std::string>& species, std::string& error) {
  if (species.empty()) { error = "candidate store requires explicit species"; return false; }
  std::set<std::string> seen;
  for (size_t i = 0; i < species.size(); ++i) {
    if (species[i].empty() || species[i].find_first_not_of("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-") != std::string::npos || !seen.insert(species[i]).second) {
      error = "candidate store has empty, duplicate or invalid species key"; return false;
    }
  }
  return true;
}

inline std::string Layout(const std::vector<std::string>& species, bool epd=false) {
  std::ostringstream o;
  o << (epd?kEpdSchemaVersion:kSchemaVersion) << '\n'
    << "floatFields px py pz energy y pt eta phi nSigmaProton nSigmaKaon nSigmaHe4 nSigmaDeuteron nSigmaTriton nSigmaHe3 mass2 deltaOneOverBeta dca invMass dcaDaughters betaGamma\n"
    << "intFields eventIndex source charge tofMatch nHitsFit trackIndex dau1EventIndex dau1Index dau2EventIndex dau2Index\n"
    << "uintFields qualityFlags\ncollections physical=0 nuclearPidQa=1\n";
  if(epd)o << "epdFields psi1 rawQx rawQy recenteredQx recenteredQy weight status bin runIndex calibrationCent hits valid\n";
  for (size_t i = 0; i < species.size(); ++i) o << "species " << species[i] << '\n';
  return o.str();
}

inline bool ValidEvent(const Event& e, std::string& error, bool epd=false) {
  if (e.inputIndex < 0 || e.cent9 < 0 || e.cent9 > 8 || e.cent16 < -1 || e.cent16 > 15 ||
      e.refMult < 0 || e.fxtMult < 0 || e.nTracks < 0 || e.nBTOFMatch < 0 || e.eventPlane2Tracks < 0) {
    error = "invalid candidate-store event indices/counts"; return false;
  }
  const Double_t values[] = {e.pvX,e.pvY,e.pvZ,e.bField,e.refMultCorr,e.centWeight,e.psi2,e.eventPlane2,e.q2x,e.q2y};
  for (unsigned i = 0; i < sizeof(values)/sizeof(values[0]); ++i)
    if (!TMath::Finite(values[i])) { error = "nonfinite candidate-store event value"; return false; }
  if(epd) {
    const double epdValues[]={e.epdPsi1,e.epdRawQx,e.epdRawQy,e.epdRecenteredQx,e.epdRecenteredQy,e.epdWeight};
    for(unsigned i=0;i<sizeof(epdValues)/sizeof(epdValues[0]);++i)if(!TMath::Finite(epdValues[i])){error="nonfinite EPD event value";return false;}
    if(e.epdHits<0 || (e.epdValid && (e.epdBin<0 || e.epdBin>=8 || e.epdRunIndex<0 || e.epdCalibrationCent<0 || e.epdCalibrationCent>8 || e.epdWeight<=0 || e.epdHits<=0 || femto_epd::XialeiBin(e.epdPsi1)!=e.epdBin || e.epdStatus!=0))){error="invalid EPD event payload";return false;}
  }
  return true;
}

inline bool ValidReference(Int_t event, Int_t track, Int_t input, Int_t nTracks) {
  return (event == -1 && track == -1) || (event == input && track >= 0 && track < nTracks);
}

inline bool ValidCandidate(const FemtoCandidate& c, const Event& e, std::string& error) {
  if (c.eventIndex != e.inputIndex || (c.source != kFemtoCandTrack && c.source != kFemtoCandResonance) ||
      c.trk.trackIndex < -1 || !ValidReference(c.reso.dau1EventIndex,c.reso.dau1Index,e.inputIndex,e.nTracks) ||
      !ValidReference(c.reso.dau2EventIndex,c.reso.dau2Index,e.inputIndex,e.nTracks)) {
    error = "invalid candidate source/event/track reference"; return false;
  }
  if ((c.source == kFemtoCandTrack && (c.trk.trackIndex < 0 || c.trk.trackIndex >= e.nTracks)) ||
      (c.source == kFemtoCandResonance && (c.reso.dau1Index < 0 || c.reso.dau2Index < 0))) {
    error = "candidate has no required daughter/track identity"; return false;
  }
  if (!TMath::Finite(c.px) || !TMath::Finite(c.py) || !TMath::Finite(c.pz) ||
      !TMath::Finite(c.energy) || c.energy <= 0 || !TMath::Finite(c.reso.invMass)) {
    error = "nonfinite/invalid candidate pairing four-vector or mass"; return false;
  }
  return true;
}

struct Buffers {
  std::vector<Int_t>* species;
  std::vector<Int_t>* collection;
  std::vector<Int_t>* qaSelected;
  std::vector<Int_t>* physicalPresent;
  std::vector<Int_t>* qaPresent;
  std::vector<Float_t>* floats;
  std::vector<Int_t>* ints;
  std::vector<UInt_t>* quality;
  Buffers() : species(new std::vector<Int_t>), collection(new std::vector<Int_t>),
      qaSelected(new std::vector<Int_t>), physicalPresent(new std::vector<Int_t>),
      qaPresent(new std::vector<Int_t>), floats(new std::vector<Float_t>),
      ints(new std::vector<Int_t>), quality(new std::vector<UInt_t>) {}
  ~Buffers() { delete species;delete collection;delete qaSelected;delete physicalPresent;delete qaPresent;delete floats;delete ints;delete quality; }
  void Clear() { species->clear();collection->clear();qaSelected->clear();physicalPresent->clear();qaPresent->clear();floats->clear();ints->clear();quality->clear(); }
private:
  Buffers(const Buffers&);
  Buffers& operator=(const Buffers&);
};

inline void Append(Buffers& b, const FemtoCandidate& c, Int_t species, Int_t collection, bool selected) {
  b.species->push_back(species);b.collection->push_back(collection);b.qaSelected->push_back(selected ? 1 : 0);b.quality->push_back(c.qualityFlags);
  const Float_t f[] = {c.px,c.py,c.pz,c.energy,c.y,c.pt,c.eta,c.phi,
    c.trk.nSigmaProton,c.trk.nSigmaKaon,c.trk.nSigmaHe4,c.trk.nSigmaDeuteron,c.trk.nSigmaTriton,c.trk.nSigmaHe3,
    c.trk.mass2,c.trk.deltaOneOverBeta,c.trk.dca,c.reso.invMass,c.reso.dcaDaughters,c.reso.betaGamma};
  const Int_t n[] = {c.eventIndex,static_cast<Int_t>(c.source),c.charge,c.trk.tofMatch ? 1 : 0,c.trk.nHitsFit,c.trk.trackIndex,
    c.reso.dau1EventIndex,c.reso.dau1Index,c.reso.dau2EventIndex,c.reso.dau2Index};
  b.floats->insert(b.floats->end(),f,f+kFloatFields);b.ints->insert(b.ints->end(),n,n+kIntFields);
}

inline bool Decode(const Buffers& b, size_t row, const std::string& key, FemtoCandidate& c, std::string& error) {
  const Float_t* f = &(*b.floats)[row*kFloatFields];const Int_t* n = &(*b.ints)[row*kIntFields];
  if (n[1] < 0 || n[1] > 1 || n[2] < std::numeric_limits<Short_t>::min() || n[2] > std::numeric_limits<Short_t>::max() ||
      (n[3] != 0 && n[3] != 1) || n[4] < std::numeric_limits<Short_t>::min() || n[4] > std::numeric_limits<Short_t>::max()) {
    error = "invalid candidate enum, boolean or short integer"; return false;
  }
  c.speciesKey=key;c.eventIndex=n[0];c.source=static_cast<FemtoCandSource>(n[1]);c.charge=static_cast<Short_t>(n[2]);
  c.trk.tofMatch=n[3] != 0;c.trk.nHitsFit=static_cast<Short_t>(n[4]);c.trk.trackIndex=n[5];
  c.reso.dau1EventIndex=n[6];c.reso.dau1Index=n[7];c.reso.dau2EventIndex=n[8];c.reso.dau2Index=n[9];c.qualityFlags=(*b.quality)[row];
  c.px=f[0];c.py=f[1];c.pz=f[2];c.energy=f[3];c.y=f[4];c.pt=f[5];c.eta=f[6];c.phi=f[7];
  c.trk.nSigmaProton=f[8];c.trk.nSigmaKaon=f[9];c.trk.nSigmaHe4=f[10];c.trk.nSigmaDeuteron=f[11];
  c.trk.nSigmaTriton=f[12];c.trk.nSigmaHe3=f[13];c.trk.mass2=f[14];c.trk.deltaOneOverBeta=f[15];c.trk.dca=f[16];
  c.reso.invMass=f[17];c.reso.dcaDaughters=f[18];c.reso.betaGamma=f[19];return true;
}

inline bool CheckScalar(TTree* tree, const char* name, const char* type, std::string& error) {
  TBranch* branch=tree->GetBranch(name);TLeaf* leaf=branch ? branch->GetLeaf(name) : 0;
  if (!branch || !leaf || std::string(leaf->GetTypeName())!=type || leaf->GetLeafCount() || leaf->GetLenStatic()!=1) {
    error=std::string("missing/wrong scalar branch ")+name;return false;
  }
  return true;
}

inline std::string CompactClass(const char* name) {
  std::string out;const std::string raw=name ? name : "";
  for(size_t i=0;i<raw.size();++i)if(!std::isspace(static_cast<unsigned char>(raw[i])))out+=raw[i];
  const std::string prefix="std::";size_t at=0;while((at=out.find(prefix))!=std::string::npos)out.erase(at,prefix.size());
  return out;
}

inline bool CheckVector(TTree* tree,const char* name,const char* type,std::string& error) {
  TBranch* branch=tree->GetBranch(name);
  if(!branch || CompactClass(branch->GetClassName())!=type) {error=std::string("missing/wrong vector branch ")+name;return false;}
  return true;
}

#define FCES_EVENT_INT_FIELDS(M) M(inputIndex) M(runId) M(eventId) M(cent9) M(cent16) M(refMult) M(fxtMult) M(nTracks) M(nBTOFMatch) M(eventPlane2Tracks)
#define FCES_EVENT_DOUBLE_FIELDS(M) M(pvX) M(pvY) M(pvZ) M(bField) M(refMultCorr) M(centWeight) M(psi2) M(eventPlane2) M(q2x) M(q2y)
#define FCES_EVENT_BOOL_FIELDS(M) M(eventPlaneValid) M(eventPlane2Valid)

inline bool Bind(TTree* tree,Event& e,Buffers& b,bool writing,std::string& error,bool epd=false) {
#define FCES_BIND_INT(name) if(writing){if(!tree->Branch(#name,&e.name,#name "/I")){error="cannot create " #name;return false;}}else{if(!CheckScalar(tree,#name,"Int_t",error))return false;if(tree->SetBranchAddress(#name,&e.name)<0){error="cannot bind " #name;return false;}}
  FCES_EVENT_INT_FIELDS(FCES_BIND_INT)
  if(epd){FCES_BIND_INT(epdStatus) FCES_BIND_INT(epdBin) FCES_BIND_INT(epdRunIndex) FCES_BIND_INT(epdCalibrationCent) FCES_BIND_INT(epdHits)}
#undef FCES_BIND_INT
#define FCES_BIND_DOUBLE(name) if(writing){if(!tree->Branch(#name,&e.name,#name "/D")){error="cannot create " #name;return false;}}else{if(!CheckScalar(tree,#name,"Double_t",error))return false;if(tree->SetBranchAddress(#name,&e.name)<0){error="cannot bind " #name;return false;}}
  FCES_EVENT_DOUBLE_FIELDS(FCES_BIND_DOUBLE)
  if(epd){FCES_BIND_DOUBLE(epdPsi1) FCES_BIND_DOUBLE(epdRawQx) FCES_BIND_DOUBLE(epdRawQy) FCES_BIND_DOUBLE(epdRecenteredQx) FCES_BIND_DOUBLE(epdRecenteredQy) FCES_BIND_DOUBLE(epdWeight)}
#undef FCES_BIND_DOUBLE
#define FCES_BIND_BOOL(name) if(writing){if(!tree->Branch(#name,&e.name,#name "/O")){error="cannot create " #name;return false;}}else{if(!CheckScalar(tree,#name,"Bool_t",error))return false;if(tree->SetBranchAddress(#name,&e.name)<0){error="cannot bind " #name;return false;}}
  FCES_EVENT_BOOL_FIELDS(FCES_BIND_BOOL)
  if(epd){FCES_BIND_BOOL(epdValid)}
#undef FCES_BIND_BOOL
#define FCES_BIND_VECTOR(name,member,type) if(writing){if(!tree->Branch(name,&b.member)){error=std::string("cannot create ")+name;return false;}}else{if(!CheckVector(tree,name,type,error))return false;if(tree->SetBranchAddress(name,&b.member)<0){error=std::string("cannot bind ")+name;return false;}tree->GetBranch(name)->SetAutoDelete(kFALSE);}
  FCES_BIND_VECTOR("candidateSpecies",species,"vector<int>")
  FCES_BIND_VECTOR("candidateCollection",collection,"vector<int>")
  FCES_BIND_VECTOR("candidateQaSelected",qaSelected,"vector<int>")
  FCES_BIND_VECTOR("physicalPresent",physicalPresent,"vector<int>")
  FCES_BIND_VECTOR("qaPresent",qaPresent,"vector<int>")
  FCES_BIND_VECTOR("candidateFloat",floats,"vector<float>")
  FCES_BIND_VECTOR("candidateInt",ints,"vector<int>")
  FCES_BIND_VECTOR("candidateQuality",quality,"vector<unsignedint>")
#undef FCES_BIND_VECTOR
  return true;
}

class Writer {
public:
  Writer():mFile(0),mTree(0),mEntries(0),mReady(false),mWritten(false),mEpd(false){}
  bool Init(TFile* file,const std::vector<std::string>& species,const std::string& selectionProvenance,const std::string& mixingProvenance,bool epd=false) {
    if(mReady || mTree || !mError.empty())return Fail("writer already initialized or failed");
    if(!file || !file->IsOpen() || file->IsZombie() || !file->IsWritable())return Fail("candidate-store output file is not writable");
    if(!ValidSpecies(species,mError))return false;
    if(selectionProvenance.empty() || mixingProvenance.empty())return Fail("candidate-store provenance is required");
    if(file->GetKey(kTreeName) || file->GetKey(kLayoutName) || file->GetKey(kSelectionName) || file->GetKey(kMixingName) || file->GetKey(kFileUuidName) || file->GetList()->FindObject(kTreeName))return Fail("refuse existing candidate-store output");
    mSpecies=species;mSelection=selectionProvenance;mMixing=mixingProvenance;mEpd=epd;mFile=file;file->cd();
    mTree=new TTree(kTreeName,"accepted events and exact StFemtoMaker candidates in processing order");mTree->SetDirectory(file);
    if(!Bind(mTree,mEvent,mBuffers,true,mError,mEpd))return false;
    mReady=true;return true;
  }
  bool Fill(const Event& event,const FemtoCandidateStore& physical,const FemtoNuclearPidPairQaStore& qa) {
    if(!mReady || mWritten || !mError.empty())return Fail("writer is not ready for Fill");
    if(!ValidEvent(event,mError,mEpd))return false;
    for(FemtoCandidateStore::const_iterator it=physical.begin();it!=physical.end();++it)if(SpeciesIndex(it->first)<0)return Fail("unknown physical species "+it->first);
    for(FemtoNuclearPidPairQaStore::const_iterator it=qa.begin();it!=qa.end();++it)if(SpeciesIndex(it->first)<0)return Fail("unknown QA species "+it->first);
    mEvent=event;mBuffers.Clear();
    for(size_t i=0;i<mSpecies.size();++i){
      FemtoCandidateStore::const_iterator p=physical.find(mSpecies[i]);FemtoNuclearPidPairQaStore::const_iterator q=qa.find(mSpecies[i]);
      mBuffers.physicalPresent->push_back(p!=physical.end()?1:0);mBuffers.qaPresent->push_back(q!=qa.end()?1:0);
      if(p!=physical.end())for(size_t j=0;j<p->second.size();++j){
        const FemtoCandidate& c=p->second[j];if(c.speciesKey!=mSpecies[i])return Fail("candidate species/key mismatch");
        if(!ValidCandidate(c,event,mError))return false;Append(mBuffers,c,static_cast<Int_t>(i),kPhysical,false);
      }
      if(q!=qa.end())for(size_t j=0;j<q->second.size();++j){
        const FemtoCandidate& c=q->second[j].candidate;if(c.speciesKey!=mSpecies[i])return Fail("QA candidate species/key mismatch");
        if(!ValidCandidate(c,event,mError))return false;Append(mBuffers,c,static_cast<Int_t>(i),kPidQa,q->second[j].selected);
      }
    }
    if(mTree->Fill()<=0)return Fail("candidate-store TTree Fill failed");++mEntries;return true;
  }
  bool Write() {
    if(!mReady || mWritten || !mError.empty() || !mFile || !mFile->IsOpen())return Fail("candidate-store writer cannot Write");
    mFile->cd();TObjString layout(Layout(mSpecies,mEpd).c_str()),selection(mSelection.c_str()),mixing(mMixing.c_str()),uuid(mFile->GetUUID().AsString());
    if(mTree->Write("",TObject::kOverwrite)<=0 || layout.Write(kLayoutName)<=0 || selection.Write(kSelectionName)<=0 || mixing.Write(kMixingName)<=0 || uuid.Write(kFileUuidName)<=0)return Fail("candidate-store metadata/tree Write failed");
    mWritten=true;return true;
  }
  Long64_t Entries()const{return mEntries;}
  const std::string& LastError()const{return mError;}
private:
  TFile* mFile;TTree* mTree;Event mEvent;Buffers mBuffers;std::vector<std::string> mSpecies;
  std::string mSelection,mMixing,mError;Long64_t mEntries;bool mReady,mWritten,mEpd;
  bool Fail(const std::string& text){if(mError.empty())mError=text;return false;}
  Int_t SpeciesIndex(const std::string& key)const{for(size_t i=0;i<mSpecies.size();++i)if(mSpecies[i]==key)return static_cast<Int_t>(i);return -1;}
  Writer(const Writer&);Writer& operator=(const Writer&);
};

class Reader {
public:
  Reader():mFile(0),mTree(0),mFileIndex(-1),mEntry(0),mEntries(0),mEntriesTotal(0),mReady(false),mFinished(false),mNewFile(false),mRequireMixing(true),mRequireCompletedProducer(false),mEpd(false),mRequireEpd(false){}
  ~Reader(){Close();}
  bool Init(const char* input,const std::vector<std::string>& species,const std::string& selectionProvenance,const std::string& mixingProvenance,bool requireOriginalMixing,bool requireCompletedProducer=false,bool requireEpd=false) {
    if(mReady || !mFiles.empty() || !mError.empty())return Fail("reader already initialized or failed");
    if(!input || !*input || !ValidSpecies(species,mError))return Fail("candidate-store input/species are required");
    if(selectionProvenance.empty() || mixingProvenance.empty())return Fail("candidate-store provenance is required");
    mSpecies=species;mSelection=selectionProvenance;mMixing=mixingProvenance;mRequireMixing=requireOriginalMixing;mRequireCompletedProducer=requireCompletedProducer;mRequireEpd=requireEpd;
    const std::string path=Absolute(input,gSystem->WorkingDirectory());
    if(path.size()>=5 && path.compare(path.size()-5,5,".root")==0)mFiles.push_back(path);
    else{
      std::ifstream in(path.c_str());if(!in)return Fail("cannot read candidate-store list "+path);
      const size_t slash=path.find_last_of('/');const std::string base=slash==std::string::npos?gSystem->WorkingDirectory():path.substr(0,slash);
      std::string line;while(std::getline(in,line)){
        const size_t comment=line.find('#');if(comment!=std::string::npos)line.erase(comment);
        std::istringstream fields(line);std::string item,extra;if(!(fields>>item))continue;
        if(fields>>extra)return Fail("candidate-store list requires one path per line");mFiles.push_back(Absolute(item,base));
      }
      if(in.bad())return Fail("candidate-store list read failed");
    }
    if(mFiles.empty())return Fail("candidate-store input list is empty");std::set<std::string> seen,seenUuid;
    for(size_t i=0;i<mFiles.size();++i){
      if(!seen.insert(mFiles[i]).second)return Fail("duplicate candidate-store input file "+mFiles[i]);
      TFile* file=TFile::Open(mFiles[i].c_str(),"READ");if(!file || file->IsZombie()){if(file)delete file;return Fail("cannot open candidate store "+mFiles[i]);}
      TTree* tree=0;bool good=CheckFile(file,tree);
      if(good && !seenUuid.insert(file->GetUUID().AsString()).second){good=false;Fail("duplicate candidate-store stream UUID");}
      Long64_t count=good?tree->GetEntries():-1;
      if(good && (count<0 || count>std::numeric_limits<Long64_t>::max()-mEntriesTotal)){good=false;Fail("invalid candidate-store entry count");}
      if(good){mFileEntries.push_back(count);mEntriesTotal+=count;}
      if(tree)tree->ResetBranchAddresses();file->Close();delete file;
      if(!good)return false;
    }
    mReady=true;return true;
  }
  bool Next(Event& event,FemtoCandidateStore& physical,FemtoNuclearPidPairQaStore& qa) {
    physical.clear();qa.clear();mNewFile=false;
    if(!mError.empty())return false;if(!mReady)return Fail("reader is not initialized");if(mFinished)return false;
    while(!mTree || mEntry>=mEntries){
      if(mFileIndex+1>=static_cast<Int_t>(mFiles.size())){Close();mFinished=true;return false;}
      if(!Open(mFileIndex+1))return false;
      if(mEntries==0)continue;
    }
    const bool first=mEntry==0;
    if(mTree->GetEntry(mEntry)<=0)return Fail("candidate-store entry read failed in "+mFiles[mFileIndex]);
    if(!ValidEvent(mEvent,mError,mEpd))return false;
    if(!ValidateBuffers())return false;
    FemtoCandidateStore p;FemtoNuclearPidPairQaStore q;
    for(size_t i=0;i<mSpecies.size();++i){if((*mBuffers.physicalPresent)[i])p[mSpecies[i]];if((*mBuffers.qaPresent)[i])q[mSpecies[i]];}
    for(size_t i=0;i<mBuffers.species->size();++i){
      const Int_t si=(*mBuffers.species)[i],collection=(*mBuffers.collection)[i],selected=(*mBuffers.qaSelected)[i];
      if(si<0 || si>=static_cast<Int_t>(mSpecies.size()) || (collection!=kPhysical && collection!=kPidQa) || (selected!=0 && selected!=1) || (collection==kPhysical && selected!=0))return Fail("invalid candidate collection/species/QA flag");
      if((collection==kPhysical && !(*mBuffers.physicalPresent)[si]) || (collection==kPidQa && !(*mBuffers.qaPresent)[si]))return Fail("candidate belongs to an absent collection");
      FemtoCandidate c;if(!Decode(mBuffers,i,mSpecies[si],c,mError) || !ValidCandidate(c,mEvent,mError))return false;
      if(collection==kPhysical)p[mSpecies[si]].push_back(c);else{FemtoNuclearPidPairQaCandidate item;item.candidate=c;item.selected=selected!=0;q[mSpecies[si]].push_back(item);}
    }
    ++mEntry;event=mEvent;physical.swap(p);qa.swap(q);mNewFile=first;return true;
  }
  bool Finished()const{return mFinished && mError.empty();}
  bool NewFile()const{return mNewFile;}
  bool HasEpdPayload()const{return mEpd;}
  Long64_t EntriesTotal()const{return mEntriesTotal;}
  const std::string& LastError()const{return mError;}
private:
  TFile* mFile;TTree* mTree;Int_t mFileIndex;Long64_t mEntry,mEntries,mEntriesTotal;Event mEvent;Buffers mBuffers;
  std::vector<std::string> mFiles,mSpecies;std::vector<Long64_t> mFileEntries;std::string mSelection,mMixing,mError;
  bool mReady,mFinished,mNewFile,mRequireMixing,mRequireCompletedProducer,mEpd,mRequireEpd;
  bool Fail(const std::string& text){if(mError.empty())mError=text;return false;}
  static std::string Absolute(const std::string& path,const std::string& base){return (!path.empty() && path[0]=='/') || path.find("://")!=std::string::npos?path:base+"/"+path;}
  void Close(){if(mTree)mTree->ResetBranchAddresses();if(mFile){mFile->Close();delete mFile;}mFile=0;mTree=0;mEntry=mEntries=0;}
  bool CheckFile(TFile* file,TTree*& tree){
    TObject* producer=file->Get(kProducerStatusName);TNamed* status=dynamic_cast<TNamed*>(producer);
    if((mRequireCompletedProducer && !status) || (producer && (!status || std::string(status->GetTitle())!="completed")))return Fail("candidate-store producer is missing, malformed or incomplete");
    TObjString* layout=dynamic_cast<TObjString*>(file->Get(kLayoutName));TObjString* selection=dynamic_cast<TObjString*>(file->Get(kSelectionName));TObjString* mixing=dynamic_cast<TObjString*>(file->Get(kMixingName));tree=dynamic_cast<TTree*>(file->Get(kTreeName));
    TObjString* uuid=dynamic_cast<TObjString*>(file->Get(kFileUuidName));
    if(!layout || !selection || !mixing || !uuid || !tree)return Fail("missing candidate-store tree/metadata");
    if(uuid->GetString()!=file->GetUUID().AsString())return Fail("candidate store was merged/repacked; replay original per-job files");
    const std::string storedLayout=layout->GetString().Data();
    if(storedLayout==Layout(mSpecies,true))mEpd=true;
    else if(storedLayout==Layout(mSpecies,false))mEpd=false;
    else return Fail("candidate-store schema/species mismatch");
    if(mRequireEpd&&!mEpd)return Fail("EPD replay requires an EPD schema-2 payload");
    if(selection->GetString().Data()!=mSelection)return Fail("candidate-store selection provenance mismatch");
    if(mixing->GetString().Length()==0 || (mRequireMixing && mixing->GetString().Data()!=mMixing))return Fail("candidate-store mixing provenance mismatch");
    mEvent=Event();
    return Bind(tree,mEvent,mBuffers,false,mError,mEpd);
  }
  bool Open(Int_t index){
    Close();mFileIndex=index;mFile=TFile::Open(mFiles[index].c_str(),"READ");
    if(!mFile || mFile->IsZombie())return Fail("cannot reopen candidate store "+mFiles[index]);
    if(!CheckFile(mFile,mTree))return false;mEntries=mTree->GetEntries();
    if(mEntries!=mFileEntries[index])return Fail("candidate-store entry count changed after Init");return true;
  }
  bool ValidateBuffers(){
    if(!mBuffers.species || !mBuffers.collection || !mBuffers.qaSelected || !mBuffers.quality || !mBuffers.floats || !mBuffers.ints || !mBuffers.physicalPresent || !mBuffers.qaPresent)return Fail("null candidate-store vector branch");
    const size_t n=mBuffers.species->size();
    if(n>std::numeric_limits<size_t>::max()/kFloatFields || n>std::numeric_limits<size_t>::max()/kIntFields ||
       mBuffers.collection->size()!=n || mBuffers.qaSelected->size()!=n || mBuffers.quality->size()!=n ||
       mBuffers.floats->size()!=n*kFloatFields || mBuffers.ints->size()!=n*kIntFields ||
       mBuffers.physicalPresent->size()!=mSpecies.size() || mBuffers.qaPresent->size()!=mSpecies.size())return Fail("inconsistent candidate-store vector lengths");
    for(size_t i=0;i<mSpecies.size();++i)if(((*mBuffers.physicalPresent)[i]!=0 && (*mBuffers.physicalPresent)[i]!=1) || ((*mBuffers.qaPresent)[i]!=0 && (*mBuffers.qaPresent)[i]!=1))return Fail("invalid collection-presence flag");
    return true;
  }
  Reader(const Reader&);Reader& operator=(const Reader&);
};

#undef FCES_EVENT_INT_FIELDS
#undef FCES_EVENT_DOUBLE_FIELDS
#undef FCES_EVENT_BOOL_FIELDS
} // namespace femto_candidate_store
#endif
