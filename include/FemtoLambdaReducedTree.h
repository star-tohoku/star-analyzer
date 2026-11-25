#ifndef FEMTO_LAMBDA_REDUCED_TREE_H
#define FEMTO_LAMBDA_REDUCED_TREE_H

#include "FemtoCandidateEventStore.h"
#include "FemtoLambdaProvider.h"
#include "TVector3.h"
#include <algorithm>
#include <map>
#include <utility>

// Value-only persistence for Lambda--helium reconstruction. Physical four-vectors
// retain the original Float_t fields; measured rigidity and geometry are separate.
namespace femto_lambda_tree {
typedef femto_candidate_store::Event Event;
static const char* const kTreeName = "femtoLambdaReducedEvents";
static const char* const kMetadataName = "FemtoLambdaReducedMetadata";
static const char* const kVersion = "lambda-helium reduced tree 1";

enum TrackRole { kProtonDaughter=1, kPionDaughter=2, kHe3Track=4, kHe4Track=8 };

struct Track {
  Int_t index, id, charge, nHitsFit, nHitsMax, nHitsDedx;
  Int_t tofIndex, tofCellId, tofMatchFlag, helixH;
  UInt_t roleFlags;
  TVector3 pMom, gMom, origin, tofHitPosition;
  Float_t dEdx, nSigmaElectron, nSigmaPion, nSigmaKaon, nSigmaProton;
  Float_t beta, tofTime, tofPathLength, tofLocalY, tofLocalZ, tofDeltaY, tofDeltaZ;
  Bool_t hasTof, tofValid, helixValid;
  Double_t helixCurvature, helixDipAngle, helixPhase;
  Track() : index(-1), id(-1), charge(0), nHitsFit(0), nHitsMax(0), nHitsDedx(0),
      tofIndex(-1), tofCellId(-1), tofMatchFlag(0), helixH(0), roleFlags(0),
      dEdx(-999), nSigmaElectron(-999), nSigmaPion(-999), nSigmaKaon(-999), nSigmaProton(-999),
      beta(-999), tofTime(-999), tofPathLength(-999), tofLocalY(-999), tofLocalZ(-999),
      tofDeltaY(-999), tofDeltaZ(-999), hasTof(kFALSE), tofValid(kFALSE), helixValid(kFALSE),
      helixCurvature(-999), helixDipAngle(-999), helixPhase(-999) {}
};

struct Payload {
  std::vector<Track> tracks;
  std::vector<FemtoLambdaCandidate> lambdas;
  FemtoLambdaEventStats stats;
  Bool_t hasAnyLegacyCandidate, baselinePoolEligible;
  Int_t originalMixBin, originalDeuteronCount, originalTritonCount;
  Int_t sourceFileIndex;
  Long64_t sourceEntry;
  Payload() : hasAnyLegacyCandidate(kFALSE), baselinePoolEligible(kFALSE),
      originalMixBin(-1), originalDeuteronCount(0), originalTritonCount(0),
      sourceFileIndex(-1), sourceEntry(-1) {}
};

inline std::string Layout(const std::vector<std::string>& species) {
  std::ostringstream s;
  s << kVersion << '\n' << femto_candidate_store::Layout(species,true)
    << "measured track TVector3 primary/global rigidity GeV/c; origin cm; B kG\n"
    << "physical lambdaP/he3P/he4P TVector3 GeV/c; canonical Float energy GeV\n"
    << "trackInt index id charge nHitsFit nHitsMax nHitsDedx tofIndex tofCellId tofMatchFlag helixH\n"
    << "trackFloat dEdx_keV_per_cm nSigmaElectron nSigmaPion nSigmaKaon nSigmaProton beta tofTime tofPathLength tofLocalY tofLocalZ tofDeltaY tofDeltaZ\n"
    << "trackDouble helixCurvature helixDipAngle helixPhase\n"
    << "trackFlags roleFlags hasTof tofValid helixValid; optional unavailable values -999\n"
    << "lambdaFloat x y z mass massError chi2Ndf topoChi2Ndf daughterDistance distanceToPv decayLength decayLengthError decayLengthSignificance vertexLineLength vertexLineLengthError vertexLineLengthSignificance cosPointing protonPidPull pionPidPull protonTofPull pionTofPull protonTofM2 pionTofM2\n"
    << "lambdaDouble protonDcaToPv pionDcaToPv protonHelixPathLength pionHelixPathLength\n"
    << "lambdaInt protonId pionId protonIndex pionIndex pdg imp5PathValid protonHasTof pionHasTof\n"
    << "stats rawTracks qualityTracks covarianceTracks pidTracks pidHypotheses protonHypotheses pionHypotheses unknownHypotheses primaryHypotheses finderParticles lambdaParticles validCandidates invalidCandidates invalidCovariance invalidPid invalidTrackIds tofTracks selectedCandidates\n"
    << "identity segmentHigh segmentLow ordinal; original inputIndex is retained\n";
  return s.str();
}

// Sparse bit blocks detect repeated physical run/event keys without a tree-node
// allocation for every accepted event. The key remains independent of segment ID.
class EventIdentities {
  struct Block { ULong64_t words[64]; Block(){for(int i=0;i<64;++i)words[i]=0;} };
  std::map<ULong64_t,Block> mBlocks;
 public:
  bool Insert(Int_t run,Int_t event) {
    if(run<=0 || event<0)return false;
    const UInt_t id=UInt_t(event);
    const ULong64_t key=(ULong64_t(UInt_t(run))<<32)|(id>>12);
    Block& block=mBlocks[key];const UInt_t word=(id&4095u)>>6;
    const ULong64_t mask=ULong64_t(1)<<(id&63u);
    if(block.words[word]&mask)return false;block.words[word]|=mask;return true;
  }
};

enum { kTrackInts=10, kTrackFloats=12, kTrackDoubles=3, kLambdaFloats=22,
       kLambdaDoubles=4, kLambdaInts=8, kStats=18 };

struct ExtraBuffers {
  std::vector<Int_t>* trackInt;
  std::vector<Float_t>* trackFloat;
  std::vector<Double_t>* trackDouble;
  std::vector<UInt_t>* trackFlags;
  std::vector<TVector3>* trackPrimaryMomentum;
  std::vector<TVector3>* trackGlobalMomentum;
  std::vector<TVector3>* trackOrigin;
  std::vector<TVector3>* trackTofHitPosition;
  std::vector<Float_t>* lambdaFloat;
  std::vector<Double_t>* lambdaDouble;
  std::vector<Int_t>* lambdaInt;
  std::vector<UInt_t>* stats;
  std::vector<TVector3>* lambdaP;
  std::vector<TVector3>* he3P;
  std::vector<TVector3>* he4P;
  std::vector<Int_t>* lambdaProtonTrackSlot;
  std::vector<Int_t>* lambdaPionTrackSlot;
  std::vector<Int_t>* he3TrackSlot;
  std::vector<Int_t>* he4TrackSlot;
  ExtraBuffers() : trackInt(new std::vector<Int_t>), trackFloat(new std::vector<Float_t>),
      trackDouble(new std::vector<Double_t>), trackFlags(new std::vector<UInt_t>),
      trackPrimaryMomentum(new std::vector<TVector3>), trackGlobalMomentum(new std::vector<TVector3>),
      trackOrigin(new std::vector<TVector3>), trackTofHitPosition(new std::vector<TVector3>),
      lambdaFloat(new std::vector<Float_t>), lambdaDouble(new std::vector<Double_t>),
      lambdaInt(new std::vector<Int_t>), stats(new std::vector<UInt_t>),
      lambdaP(new std::vector<TVector3>), he3P(new std::vector<TVector3>), he4P(new std::vector<TVector3>),
      lambdaProtonTrackSlot(new std::vector<Int_t>), lambdaPionTrackSlot(new std::vector<Int_t>),
      he3TrackSlot(new std::vector<Int_t>), he4TrackSlot(new std::vector<Int_t>) {}
  ~ExtraBuffers() {
    delete trackInt;delete trackFloat;delete trackDouble;delete trackFlags;
    delete trackPrimaryMomentum;delete trackGlobalMomentum;delete trackOrigin;delete trackTofHitPosition;
    delete lambdaFloat;delete lambdaDouble;delete lambdaInt;delete stats;
    delete lambdaP;delete he3P;delete he4P;delete lambdaProtonTrackSlot;delete lambdaPionTrackSlot;
    delete he3TrackSlot;delete he4TrackSlot;
  }
  void Clear() {
    trackInt->clear();trackFloat->clear();trackDouble->clear();trackFlags->clear();
    trackPrimaryMomentum->clear();trackGlobalMomentum->clear();trackOrigin->clear();trackTofHitPosition->clear();
    lambdaFloat->clear();lambdaDouble->clear();lambdaInt->clear();stats->clear();
    lambdaP->clear();he3P->clear();he4P->clear();lambdaProtonTrackSlot->clear();lambdaPionTrackSlot->clear();
    he3TrackSlot->clear();he4TrackSlot->clear();
  }
 private:
  ExtraBuffers(const ExtraBuffers&);ExtraBuffers& operator=(const ExtraBuffers&);
};

struct Frame {
  ULong64_t segmentHigh,segmentLow;
  Long64_t ordinal,sourceEntry;
  Int_t originalMixBin,originalDeuteronCount,originalTritonCount,sourceFileIndex;
  Bool_t hasAnyLegacyCandidate,baselinePoolEligible;
  Frame() : segmentHigh(0),segmentLow(0),ordinal(0),sourceEntry(-1),originalMixBin(-1),
      originalDeuteronCount(0),originalTritonCount(0),sourceFileIndex(-1),
      hasAnyLegacyCandidate(kFALSE),baselinePoolEligible(kFALSE) {}
};
typedef std::pair<ULong64_t,ULong64_t> SegmentKey;

inline bool BindExtra(TTree* tree,Frame& f,ExtraBuffers& b,bool writing,std::string& error) {
#define FLRT_SCALAR(name,type,leaf) if(writing){if(!tree->Branch(#name,&f.name,#name "/" leaf)){error="cannot create " #name;return false;}}else{if(!femto_candidate_store::CheckScalar(tree,#name,type,error))return false;if(tree->SetBranchAddress(#name,&f.name)<0){error="cannot bind " #name;return false;}}
  FLRT_SCALAR(segmentHigh,"ULong64_t","l") FLRT_SCALAR(segmentLow,"ULong64_t","l")
  FLRT_SCALAR(ordinal,"Long64_t","L") FLRT_SCALAR(sourceEntry,"Long64_t","L")
  FLRT_SCALAR(originalMixBin,"Int_t","I") FLRT_SCALAR(originalDeuteronCount,"Int_t","I")
  FLRT_SCALAR(originalTritonCount,"Int_t","I") FLRT_SCALAR(sourceFileIndex,"Int_t","I")
  FLRT_SCALAR(hasAnyLegacyCandidate,"Bool_t","O") FLRT_SCALAR(baselinePoolEligible,"Bool_t","O")
#undef FLRT_SCALAR
#define FLRT_VECTOR(member,type) if(writing){if(!tree->Branch(#member,&b.member)){error="cannot create " #member;return false;}}else{if(!femto_candidate_store::CheckVector(tree,#member,type,error))return false;if(tree->SetBranchAddress(#member,&b.member)<0){error="cannot bind " #member;return false;}tree->GetBranch(#member)->SetAutoDelete(kFALSE);}
  FLRT_VECTOR(trackInt,"vector<int>") FLRT_VECTOR(trackFloat,"vector<float>")
  FLRT_VECTOR(trackDouble,"vector<double>") FLRT_VECTOR(trackFlags,"vector<unsignedint>")
  FLRT_VECTOR(trackPrimaryMomentum,"vector<TVector3>") FLRT_VECTOR(trackGlobalMomentum,"vector<TVector3>")
  FLRT_VECTOR(trackOrigin,"vector<TVector3>") FLRT_VECTOR(trackTofHitPosition,"vector<TVector3>")
  FLRT_VECTOR(lambdaFloat,"vector<float>") FLRT_VECTOR(lambdaDouble,"vector<double>") FLRT_VECTOR(lambdaInt,"vector<int>")
  FLRT_VECTOR(stats,"vector<unsignedint>")
  FLRT_VECTOR(lambdaP,"vector<TVector3>") FLRT_VECTOR(he3P,"vector<TVector3>") FLRT_VECTOR(he4P,"vector<TVector3>")
  FLRT_VECTOR(lambdaProtonTrackSlot,"vector<int>") FLRT_VECTOR(lambdaPionTrackSlot,"vector<int>")
  FLRT_VECTOR(he3TrackSlot,"vector<int>") FLRT_VECTOR(he4TrackSlot,"vector<int>")
#undef FLRT_VECTOR
  return true;
}

inline bool FiniteVector(const TVector3& p) {
  return TMath::Finite(p.X())&&TMath::Finite(p.Y())&&TMath::Finite(p.Z());
}

inline bool EncodePayload(const Event& event,const FemtoCandidateStore& physical,const Payload& p,
                          Frame& frame,ExtraBuffers& b,std::string& error) {
  b.Clear();frame.hasAnyLegacyCandidate=p.hasAnyLegacyCandidate;frame.baselinePoolEligible=p.baselinePoolEligible;
  frame.originalMixBin=p.originalMixBin;frame.originalDeuteronCount=p.originalDeuteronCount;
  frame.originalTritonCount=p.originalTritonCount;frame.sourceFileIndex=p.sourceFileIndex;frame.sourceEntry=p.sourceEntry;
  if(p.originalDeuteronCount<0 || p.originalTritonCount<0 || p.sourceFileIndex< -1 || p.sourceEntry< -1 ||
      (p.baselinePoolEligible&&!p.hasAnyLegacyCandidate)){error="invalid reduced-tree event admission/provenance";return false;}
  std::map<Int_t,Int_t> slots;
  for(size_t i=0;i<p.tracks.size();++i) {
    const Track& t=p.tracks[i];
    if(t.index<0 || t.index>=event.nTracks || !slots.insert(std::make_pair(t.index,Int_t(i))).second ||
       !FiniteVector(t.pMom)||!FiniteVector(t.gMom)||!FiniteVector(t.origin)||!FiniteVector(t.tofHitPosition)||
       (t.charge!= -1&&t.charge!=1)||t.nHitsFit<0||t.nHitsMax<0||t.nHitsDedx<0||
       (t.tofValid&&(!t.hasTof||!TMath::Finite(t.beta)||t.beta<=0))) {
      error="invalid/duplicate reduced-tree track";return false;
    }
    const Int_t ints[]={t.index,t.id,t.charge,t.nHitsFit,t.nHitsMax,t.nHitsDedx,t.tofIndex,t.tofCellId,t.tofMatchFlag,t.helixH};
    const Float_t floats[]={t.dEdx,t.nSigmaElectron,t.nSigmaPion,t.nSigmaKaon,t.nSigmaProton,
      t.beta,t.tofTime,t.tofPathLength,t.tofLocalY,t.tofLocalZ,t.tofDeltaY,t.tofDeltaZ};
    const Double_t doubles[]={t.helixCurvature,t.helixDipAngle,t.helixPhase};
    b.trackInt->insert(b.trackInt->end(),ints,ints+kTrackInts);
    b.trackFloat->insert(b.trackFloat->end(),floats,floats+kTrackFloats);
    b.trackDouble->insert(b.trackDouble->end(),doubles,doubles+kTrackDoubles);
    b.trackFlags->push_back(t.roleFlags);b.trackFlags->push_back(t.hasTof?1:0);
    b.trackFlags->push_back(t.tofValid?1:0);b.trackFlags->push_back(t.helixValid?1:0);
    b.trackPrimaryMomentum->push_back(t.pMom);b.trackGlobalMomentum->push_back(t.gMom);
    b.trackOrigin->push_back(t.origin);b.trackTofHitPosition->push_back(t.tofHitPosition);
  }
  FemtoCandidateStore::const_iterator li=physical.find("lambda");
  const size_t nl=li==physical.end()?0:li->second.size();
  if(p.lambdas.size()!=nl){error="Lambda detail/canonical collection length mismatch";return false;}
  for(size_t i=0;i<nl;++i) {
    const FemtoLambdaCandidate& l=p.lambdas[i];const FemtoCandidate& c=li->second[i];
    if(l.candidate.px!=c.px||l.candidate.py!=c.py||l.candidate.pz!=c.pz||l.candidate.energy!=c.energy||
       l.protonIndex!=c.reso.dau1Index||l.pionIndex!=c.reso.dau2Index||
       !slots.count(l.protonIndex)||!slots.count(l.pionIndex)) {
      error="Lambda detail/track references disagree with canonical candidate";return false;
    }
    const Float_t f[]={l.x,l.y,l.z,l.mass,l.massError,l.chi2Ndf,l.topoChi2Ndf,l.daughterDistance,l.distanceToPv,
      l.decayLength,l.decayLengthError,l.decayLengthSignificance,l.vertexLineLength,l.vertexLineLengthError,
      l.vertexLineLengthSignificance,l.cosPointing,l.protonPidPull,l.pionPidPull,l.protonTofPull,l.pionTofPull,l.protonTofM2,l.pionTofM2};
    const Double_t d[]={l.protonDcaToPv,l.pionDcaToPv,l.protonHelixPathLength,l.pionHelixPathLength};
    const Int_t n[]={l.protonId,l.pionId,l.protonIndex,l.pionIndex,l.pdg,l.imp5PathValid?1:0,l.protonHasTof?1:0,l.pionHasTof?1:0};
    b.lambdaFloat->insert(b.lambdaFloat->end(),f,f+kLambdaFloats);b.lambdaDouble->insert(b.lambdaDouble->end(),d,d+kLambdaDoubles);
    b.lambdaInt->insert(b.lambdaInt->end(),n,n+kLambdaInts);b.lambdaP->push_back(TVector3(c.px,c.py,c.pz));
    b.lambdaProtonTrackSlot->push_back(slots[l.protonIndex]);b.lambdaPionTrackSlot->push_back(slots[l.pionIndex]);
  }
  const char* names[]={"he3","he4"};
  for(int species=0;species<2;++species) {
    FemtoCandidateStore::const_iterator it=physical.find(names[species]);
    if(it==physical.end())continue;
    for(size_t i=0;i<it->second.size();++i) {
      const FemtoCandidate& c=it->second[i];
      if(!slots.count(c.trk.trackIndex)){error="selected helium has no saved raw track";return false;}
      (species==0?b.he3P:b.he4P)->push_back(TVector3(c.px,c.py,c.pz));
      (species==0?b.he3TrackSlot:b.he4TrackSlot)->push_back(slots[c.trk.trackIndex]);
    }
  }
  const FemtoLambdaEventStats& s=p.stats;
  const UInt_t stats[]={s.rawTracks,s.qualityTracks,s.covarianceTracks,s.pidTracks,s.pidHypotheses,s.protonHypotheses,
    s.pionHypotheses,s.unknownHypotheses,s.primaryHypotheses,s.finderParticles,s.lambdaParticles,s.validCandidates,
    s.invalidCandidates,s.invalidCovariance,s.invalidPid,s.invalidTrackIds,s.tofTracks,s.selectedCandidates};
  b.stats->assign(stats,stats+kStats);return true;
}

inline bool DecodePayload(const Event& event,const FemtoCandidateStore& physical,const Frame& f,
                          const ExtraBuffers& b,Payload& p,std::string& error) {
  if(!b.trackInt||!b.trackFloat||!b.trackDouble||!b.trackFlags||!b.trackPrimaryMomentum||!b.trackGlobalMomentum||
     !b.trackOrigin||!b.trackTofHitPosition||!b.lambdaFloat||!b.lambdaDouble||!b.lambdaInt||!b.stats||
     !b.lambdaP||!b.he3P||!b.he4P||!b.lambdaProtonTrackSlot||!b.lambdaPionTrackSlot||!b.he3TrackSlot||!b.he4TrackSlot) {
    error="null reduced-tree detail vector";return false;
  }
  const size_t nt=b.trackPrimaryMomentum->size(),nl=b.lambdaP->size();
  if(b.trackInt->size()!=nt*kTrackInts||b.trackFloat->size()!=nt*kTrackFloats||b.trackDouble->size()!=nt*kTrackDoubles||
     b.trackFlags->size()!=nt*4||b.trackGlobalMomentum->size()!=nt||b.trackOrigin->size()!=nt||b.trackTofHitPosition->size()!=nt||
     b.lambdaFloat->size()!=nl*kLambdaFloats||b.lambdaDouble->size()!=nl*kLambdaDoubles||b.lambdaInt->size()!=nl*kLambdaInts||
     b.lambdaProtonTrackSlot->size()!=nl||b.lambdaPionTrackSlot->size()!=nl||b.he3TrackSlot->size()!=b.he3P->size()||
     b.he4TrackSlot->size()!=b.he4P->size()||b.stats->size()!=kStats){error="inconsistent reduced-tree detail vector lengths";return false;}
  p=Payload();p.hasAnyLegacyCandidate=f.hasAnyLegacyCandidate;p.baselinePoolEligible=f.baselinePoolEligible;
  p.originalMixBin=f.originalMixBin;p.originalDeuteronCount=f.originalDeuteronCount;p.originalTritonCount=f.originalTritonCount;
  p.sourceFileIndex=f.sourceFileIndex;p.sourceEntry=f.sourceEntry;
  for(size_t i=0;i<nt;++i) {
    Track t;const Int_t* n=&(*b.trackInt)[i*kTrackInts];const Float_t* v=&(*b.trackFloat)[i*kTrackFloats];
    const Double_t* d=&(*b.trackDouble)[i*kTrackDoubles];const UInt_t* flags=&(*b.trackFlags)[i*4];
    if(flags[1]>1||flags[2]>1||flags[3]>1){error="invalid reduced-tree track flags";return false;}
    t.index=n[0];t.id=n[1];t.charge=n[2];t.nHitsFit=n[3];t.nHitsMax=n[4];t.nHitsDedx=n[5];
    t.tofIndex=n[6];t.tofCellId=n[7];t.tofMatchFlag=n[8];t.helixH=n[9];
    t.dEdx=v[0];t.nSigmaElectron=v[1];t.nSigmaPion=v[2];t.nSigmaKaon=v[3];t.nSigmaProton=v[4];
    t.beta=v[5];t.tofTime=v[6];t.tofPathLength=v[7];t.tofLocalY=v[8];t.tofLocalZ=v[9];t.tofDeltaY=v[10];t.tofDeltaZ=v[11];
    t.helixCurvature=d[0];t.helixDipAngle=d[1];t.helixPhase=d[2];
    t.roleFlags=flags[0];t.hasTof=flags[1];t.tofValid=flags[2];t.helixValid=flags[3];
    t.pMom=(*b.trackPrimaryMomentum)[i];t.gMom=(*b.trackGlobalMomentum)[i];t.origin=(*b.trackOrigin)[i];
    t.tofHitPosition=(*b.trackTofHitPosition)[i];p.tracks.push_back(t);
  }
  FemtoCandidateStore::const_iterator it=physical.find("lambda");
  if(nl!=(it==physical.end()?0:it->second.size())){error="Lambda detail/canonical count mismatch";return false;}
  for(size_t i=0;i<nl;++i) {
    FemtoLambdaCandidate l;l.candidate=it->second[i];const Float_t* v=&(*b.lambdaFloat)[i*kLambdaFloats];
    const Double_t* d=&(*b.lambdaDouble)[i*kLambdaDoubles];const Int_t* n=&(*b.lambdaInt)[i*kLambdaInts];
    if(n[5]<0||n[5]>1||n[6]<0||n[6]>1||n[7]<0||n[7]>1){error="invalid Lambda detail flags";return false;}
    l.x=v[0];l.y=v[1];l.z=v[2];l.mass=v[3];l.massError=v[4];l.chi2Ndf=v[5];l.topoChi2Ndf=v[6];
    l.daughterDistance=v[7];l.distanceToPv=v[8];l.decayLength=v[9];l.decayLengthError=v[10];l.decayLengthSignificance=v[11];
    l.vertexLineLength=v[12];l.vertexLineLengthError=v[13];l.vertexLineLengthSignificance=v[14];l.cosPointing=v[15];
    l.protonPidPull=v[16];l.pionPidPull=v[17];l.protonTofPull=v[18];l.pionTofPull=v[19];l.protonTofM2=v[20];l.pionTofM2=v[21];
    l.protonDcaToPv=d[0];l.pionDcaToPv=d[1];l.protonHelixPathLength=d[2];l.pionHelixPathLength=d[3];
    l.protonId=n[0];l.pionId=n[1];l.protonIndex=n[2];l.pionIndex=n[3];l.pdg=n[4];
    l.imp5PathValid=n[5];l.protonHasTof=n[6];l.pionHasTof=n[7];p.lambdas.push_back(l);
  }
  const UInt_t* s=&(*b.stats)[0];
  p.stats.rawTracks=s[0];p.stats.qualityTracks=s[1];p.stats.covarianceTracks=s[2];p.stats.pidTracks=s[3];
  p.stats.pidHypotheses=s[4];p.stats.protonHypotheses=s[5];p.stats.pionHypotheses=s[6];p.stats.unknownHypotheses=s[7];
  p.stats.primaryHypotheses=s[8];p.stats.finderParticles=s[9];p.stats.lambdaParticles=s[10];p.stats.validCandidates=s[11];
  p.stats.invalidCandidates=s[12];p.stats.invalidCovariance=s[13];p.stats.invalidPid=s[14];p.stats.invalidTrackIds=s[15];
  p.stats.tofTracks=s[16];p.stats.selectedCandidates=s[17];
  ExtraBuffers check;Frame frame;
  if(!EncodePayload(event,physical,p,frame,check,error))return false;
  if(*check.lambdaProtonTrackSlot!=*b.lambdaProtonTrackSlot||*check.lambdaPionTrackSlot!=*b.lambdaPionTrackSlot||
      *check.he3TrackSlot!=*b.he3TrackSlot||*check.he4TrackSlot!=*b.he4TrackSlot||
      *check.lambdaP!=*b.lambdaP||*check.he3P!=*b.he3P||*check.he4P!=*b.he4P) {
    error="physical TVector3/track slots disagree with canonical payload";return false;
  }
  return true;
}

struct Metadata {
  ULong64_t segmentHigh,segmentLow;
  Long64_t entries;
  Bool_t completed;
  std::string* producerId;
  std::string* layout;
  std::string* selection;
  std::string* mixing;
  Metadata() : segmentHigh(0),segmentLow(0),entries(0),completed(kFALSE),producerId(new std::string),
      layout(new std::string),selection(new std::string),mixing(new std::string) {}
  ~Metadata(){delete producerId;delete layout;delete selection;delete mixing;}
 private:
  Metadata(const Metadata&);Metadata& operator=(const Metadata&);
};

inline bool BindMetadata(TTree* tree,Metadata& m,bool writing,std::string& error) {
#define FLRT_META_SCALAR(name,type,leaf) if(writing){if(!tree->Branch(#name,&m.name,#name "/" leaf)){error="cannot create metadata " #name;return false;}}else{if(!femto_candidate_store::CheckScalar(tree,#name,type,error))return false;if(tree->SetBranchAddress(#name,&m.name)<0){error="cannot bind metadata " #name;return false;}}
  FLRT_META_SCALAR(segmentHigh,"ULong64_t","l") FLRT_META_SCALAR(segmentLow,"ULong64_t","l")
  FLRT_META_SCALAR(entries,"Long64_t","L") FLRT_META_SCALAR(completed,"Bool_t","O")
#undef FLRT_META_SCALAR
#define FLRT_META_STRING(name) if(writing){if(!tree->Branch(#name,&m.name)){error="cannot create metadata " #name;return false;}}else{TBranch* branch=tree->GetBranch(#name);if(!branch||femto_candidate_store::CompactClass(branch->GetClassName())!="string"){error="missing/wrong metadata string " #name;return false;}if(tree->SetBranchAddress(#name,&m.name)<0){error="cannot bind metadata " #name;return false;}branch->SetAutoDelete(kFALSE);}
  FLRT_META_STRING(producerId) FLRT_META_STRING(layout) FLRT_META_STRING(selection) FLRT_META_STRING(mixing)
#undef FLRT_META_STRING
  return true;
}

class Writer {
 public:
  Writer() : mFile(0),mTree(0),mMetaTree(0),mEntries(0),mLastInput(-1),mReady(false),mWritten(false) {}
  bool Init(TFile* file,const std::vector<std::string>& species,const std::string& selectionProvenance,
      const std::string& mixingProvenance,const std::string& producerId,Int_t compressionLevel,Long64_t autoFlush) {
    if(mReady||!mError.empty())return Fail("reduced-tree writer already initialized/failed");
    if(!file||!file->IsOpen()||file->IsZombie()||!file->IsWritable())return Fail("reduced-tree output is not writable");
    if(!femto_candidate_store::ValidSpecies(species,mError))return false;
    if(selectionProvenance.empty()||mixingProvenance.empty()||producerId.empty())return Fail("reduced-tree provenance/producer ID is required");
    if(compressionLevel<0||compressionLevel>9||autoFlush==0)return Fail("invalid reduced-tree compression/autoFlush");
    if(file->GetKey(kTreeName)||file->GetList()->FindObject(kTreeName)||file->GetKey(kMetadataName)||file->GetList()->FindObject(kMetadataName))
      return Fail("refusing existing reduced-tree output");
    mSpecies=species;mFile=file;file->SetCompressionLevel(compressionLevel);file->cd();
    std::string hex=TUUID().AsString();hex.erase(std::remove(hex.begin(),hex.end(),'-'),hex.end());
    if(hex.size()!=32)return Fail("cannot create persistent segment identity");
    std::istringstream high(hex.substr(0,16)),low(hex.substr(16,16));high>>std::hex>>mFrame.segmentHigh;low>>std::hex>>mFrame.segmentLow;
    if(!high||!low)return Fail("cannot parse persistent segment identity");
    mTree=new TTree(kTreeName,"reconstructed Lambda/helium events; pair analysis is downstream");mTree->SetDirectory(file);mTree->SetAutoFlush(autoFlush);
    if(!femto_candidate_store::Bind(mTree,mEvent,mCandidates,true,mError,true)||!BindExtra(mTree,mFrame,mExtra,true,mError))return false;
    mMetaTree=new TTree(kMetadataName,"completed producer segments independent of transport ROOT files");mMetaTree->SetDirectory(file);
    mMeta.segmentHigh=mFrame.segmentHigh;mMeta.segmentLow=mFrame.segmentLow;*mMeta.producerId=producerId;
    *mMeta.layout=Layout(species);*mMeta.selection=selectionProvenance;*mMeta.mixing=mixingProvenance;
    if(!BindMetadata(mMetaTree,mMeta,true,mError))return false;mReady=true;return true;
  }
  bool Fill(const Event& event,const FemtoCandidateStore& physical,const FemtoNuclearPidPairQaStore& qa,const Payload& payload) {
    if(!mReady||mWritten||!mError.empty())return Fail("reduced-tree writer is not ready");
    if(!femto_candidate_store::ValidEvent(event,mError,true))return false;
    if(event.inputIndex<=mLastInput||!mIdentities.Insert(event.runId,event.eventId))return Fail("duplicate event or non-increasing source input index");
    mEvent=event;mCandidates.Clear();
    for(FemtoCandidateStore::const_iterator it=physical.begin();it!=physical.end();++it)
      if(SpeciesIndex(it->first)<0&&it->first!="deuteron"&&it->first!="triton")return Fail("unsupported physical species "+it->first);
    for(FemtoNuclearPidPairQaStore::const_iterator it=qa.begin();it!=qa.end();++it)
      if(SpeciesIndex(it->first)<0&&it->first!="deuteron"&&it->first!="triton")return Fail("unsupported QA species "+it->first);
    for(size_t i=0;i<mSpecies.size();++i) {
      FemtoCandidateStore::const_iterator p=physical.find(mSpecies[i]);FemtoNuclearPidPairQaStore::const_iterator q=qa.find(mSpecies[i]);
      mCandidates.physicalPresent->push_back(p!=physical.end()?1:0);mCandidates.qaPresent->push_back(q!=qa.end()?1:0);
      if(p!=physical.end())for(size_t j=0;j<p->second.size();++j) {
        if(p->second[j].speciesKey!=mSpecies[i])return Fail("physical species/key mismatch");
        if(!femto_candidate_store::ValidCandidate(p->second[j],event,mError))return false;
        femto_candidate_store::Append(mCandidates,p->second[j],Int_t(i),femto_candidate_store::kPhysical,false);
      }
      if(q!=qa.end())for(size_t j=0;j<q->second.size();++j) {
        if(q->second[j].candidate.speciesKey!=mSpecies[i])return Fail("QA species/key mismatch");
        if(!femto_candidate_store::ValidCandidate(q->second[j].candidate,event,mError))return false;
        femto_candidate_store::Append(mCandidates,q->second[j].candidate,Int_t(i),femto_candidate_store::kPidQa,q->second[j].selected);
      }
    }
    mFrame.ordinal=mEntries;
    if(!EncodePayload(event,physical,payload,mFrame,mExtra,mError))return false;
    if(mTree->Fill()<=0)return Fail("reduced-tree TTree Fill failed");++mEntries;mLastInput=event.inputIndex;return true;
  }
  bool Write() {
    if(!mReady||mWritten||!mError.empty()||!mFile||!mFile->IsOpen())return Fail("reduced-tree writer cannot Write");
    mMeta.entries=mEntries;mMeta.completed=kTRUE;mFile->cd();
    if(mMetaTree->Fill()<=0||mTree->Write("",TObject::kOverwrite)<=0||mMetaTree->Write("",TObject::kOverwrite)<=0)return Fail("reduced-tree metadata/tree Write failed");
    mWritten=true;return true;
  }
  Long64_t Entries()const{return mEntries;}
  const std::string& LastError()const{return mError;}
 private:
  TFile* mFile;TTree* mTree;TTree* mMetaTree;Event mEvent;Frame mFrame;ExtraBuffers mExtra;
  femto_candidate_store::Buffers mCandidates;Metadata mMeta;std::vector<std::string> mSpecies;
  Long64_t mEntries;Int_t mLastInput;bool mReady,mWritten;std::string mError;EventIdentities mIdentities;
  bool Fail(const std::string& error){if(mError.empty())mError=error;return false;}
  Int_t SpeciesIndex(const std::string& key)const{for(size_t i=0;i<mSpecies.size();++i)if(mSpecies[i]==key)return Int_t(i);return -1;}
  Writer(const Writer&);Writer& operator=(const Writer&);
};

class Reader {
 public:
  Reader() : mFile(0),mTree(0),mFileIndex(-1),mEntry(0),mEntries(0),mTotal(0),mReady(false),
      mFinished(false),mNewSegment(false),mNewFile(false),mRequireMixing(true),mHaveSegment(false),mLastInput(-1) {}
  ~Reader(){Close();}
  bool Init(const char* input,const std::vector<std::string>& species,const std::string& selectionProvenance,
      const std::string& mixingProvenance,bool requireOriginalMixing) {
    if(mReady||!mFiles.empty()||!mError.empty())return Fail("reduced-tree reader already initialized/failed");
    if(!input||!*input||!femto_candidate_store::ValidSpecies(species,mError)||selectionProvenance.empty()||mixingProvenance.empty())
      return Fail("reduced-tree input/species/provenance is required");
    mSpecies=species;mSelection=selectionProvenance;mMixing=mixingProvenance;mRequireMixing=requireOriginalMixing;
    const std::string path=Absolute(input,gSystem->WorkingDirectory());
    if(path.size()>=5&&path.compare(path.size()-5,5,".root")==0)mFiles.push_back(path);
    else {
      std::ifstream in(path.c_str());if(!in)return Fail("cannot read reduced-tree input list "+path);
      const size_t slash=path.find_last_of('/');const std::string base=slash==std::string::npos?gSystem->WorkingDirectory():path.substr(0,slash);
      std::string line;while(std::getline(in,line)) {
        const size_t comment=line.find('#');if(comment!=std::string::npos)line.erase(comment);
        std::istringstream fields(line);std::string item,extra;if(!(fields>>item))continue;
        if(fields>>extra)return Fail("reduced-tree list requires one path per line");mFiles.push_back(Absolute(item,base));
      }
      if(in.bad())return Fail("reduced-tree list read failed");
    }
    if(mFiles.empty())return Fail("reduced-tree input list is empty");
    std::set<std::string> paths;std::set<SegmentKey> segments;
    for(size_t i=0;i<mFiles.size();++i) {
      if(!paths.insert(mFiles[i]).second)return Fail("duplicate reduced-tree input path");
      TFile* file=TFile::Open(mFiles[i].c_str(),"READ");
      if(!file||file->IsZombie()){if(file)delete file;return Fail("cannot open reduced-tree input "+mFiles[i]);}
      TTree* tree=0;std::map<SegmentKey,Long64_t> records;bool good=CheckFile(file,tree,records);
      if(good)for(std::map<SegmentKey,Long64_t>::const_iterator j=records.begin();j!=records.end();++j)
        if(!segments.insert(j->first).second){good=false;Fail("duplicate producer segment across reduced-tree inputs");break;}
      const Long64_t n=good?tree->GetEntries():-1;
      if(good&&(n<0||n>std::numeric_limits<Long64_t>::max()-mTotal)){good=false;Fail("invalid reduced-tree entry count");}
      if(good){mFileEntries.push_back(n);mTotal+=n;}
      if(tree)tree->ResetBranchAddresses();file->Close();delete file;if(!good)return false;
    }
    mReady=true;return true;
  }
  bool Next(Event& event,FemtoCandidateStore& physical,FemtoNuclearPidPairQaStore& qa,Payload& payload) {
    physical.clear();qa.clear();payload=Payload();mNewSegment=false;mNewFile=false;
    if(!mError.empty()||mFinished)return false;if(!mReady)return Fail("reduced-tree reader is not initialized");
    while(!mTree||mEntry>=mEntries) {
      if(mTree&&!FinishFile())return false;
      if(mFileIndex+1>=Int_t(mFiles.size())){Close();mFinished=true;return false;}
      if(!Open(mFileIndex+1))return false;if(mEntries==0)continue;
    }
    const bool first=mEntry==0;
    if(mTree->GetEntry(mEntry)<=0)return Fail("reduced-tree entry read failed");
    if(!femto_candidate_store::ValidEvent(mEvent,mError,true)||!ValidateCandidates())return false;
    const SegmentKey key(mFrame.segmentHigh,mFrame.segmentLow);
    std::map<SegmentKey,Long64_t>::const_iterator meta=mFileSegments.find(key);
    if(meta==mFileSegments.end()||mFrame.ordinal!=mReadCounts[key]||mFrame.ordinal>=meta->second)return Fail("invalid/duplicate reduced-tree segment ordinal");
    if(!mHaveSegment||key!=mCurrentSegment) {
      if(mHaveSegment&&mReadCounts[mCurrentSegment]!=mFileSegments[mCurrentSegment])return Fail("interleaved/incomplete reduced-tree segment");
      if(!mStartedSegments.insert(key).second)return Fail("repeated/noncontiguous reduced-tree segment");
      mCurrentSegment=key;mHaveSegment=true;mLastInput=-1;mNewSegment=true;
    }
    if(mEvent.inputIndex<=mLastInput||!mIdentities.Insert(mEvent.runId,mEvent.eventId))return Fail("duplicate physical event or reordered source input index");
    for(size_t i=0;i<mSpecies.size();++i) {
      if((*mCandidates.physicalPresent)[i])physical[mSpecies[i]];
      if((*mCandidates.qaPresent)[i])qa[mSpecies[i]];
    }
    for(size_t i=0;i<mCandidates.species->size();++i) {
      const Int_t si=(*mCandidates.species)[i],collection=(*mCandidates.collection)[i],selected=(*mCandidates.qaSelected)[i];
      if(si<0||si>=Int_t(mSpecies.size())||(collection!=femto_candidate_store::kPhysical&&collection!=femto_candidate_store::kPidQa)||
         (selected!=0&&selected!=1)||(collection==femto_candidate_store::kPhysical&&selected!=0))return Fail("invalid reduced-tree candidate collection/species/flag");
      if((collection==femto_candidate_store::kPhysical&&!(*mCandidates.physicalPresent)[si])||
         (collection==femto_candidate_store::kPidQa&&!(*mCandidates.qaPresent)[si]))return Fail("candidate belongs to absent reduced-tree collection");
      FemtoCandidate c;if(!femto_candidate_store::Decode(mCandidates,i,mSpecies[si],c,mError)||!femto_candidate_store::ValidCandidate(c,mEvent,mError))return false;
      if(collection==femto_candidate_store::kPhysical)physical[mSpecies[si]].push_back(c);
      else{FemtoNuclearPidPairQaCandidate item;item.candidate=c;item.selected=selected!=0;qa[mSpecies[si]].push_back(item);}
    }
    if(!DecodePayload(mEvent,physical,mFrame,mExtra,payload,mError))return false;
    event=mEvent;++mEntry;++mReadCounts[key];mLastInput=mEvent.inputIndex;mNewFile=first;return true;
  }
  bool Finished()const{return mFinished&&mError.empty();}
  bool NewSegment()const{return mNewSegment;}
  bool NewFile()const{return mNewFile;}
  Long64_t EntriesTotal()const{return mTotal;}
  const std::string& LastError()const{return mError;}
 private:
  TFile* mFile;TTree* mTree;Int_t mFileIndex;Long64_t mEntry,mEntries,mTotal;
  Event mEvent;Frame mFrame;ExtraBuffers mExtra;femto_candidate_store::Buffers mCandidates;
  std::vector<std::string> mFiles,mSpecies;std::vector<Long64_t> mFileEntries;
  std::string mSelection,mMixing,mError;bool mReady,mFinished,mNewSegment,mNewFile,mRequireMixing,mHaveSegment;
  SegmentKey mCurrentSegment;Int_t mLastInput;EventIdentities mIdentities;
  std::map<SegmentKey,Long64_t> mFileSegments,mReadCounts;std::set<SegmentKey> mStartedSegments;
  bool Fail(const std::string& error){if(mError.empty())mError=error;return false;}
  static std::string Absolute(const std::string& path,const std::string& base){return (!path.empty()&&path[0]=='/')||path.find("://")!=std::string::npos?path:base+"/"+path;}
  void Close(){if(mTree)mTree->ResetBranchAddresses();if(mFile){mFile->Close();delete mFile;}mFile=0;mTree=0;mEntry=mEntries=0;}
  bool CheckFile(TFile* file,TTree*& tree,std::map<SegmentKey,Long64_t>& records) {
    TNamed* status=dynamic_cast<TNamed*>(file->Get(femto_candidate_store::kProducerStatusName));
    if(!status||std::string(status->GetTitle())!="completed")return Fail("reduced-tree producer is incomplete/missing");
    tree=dynamic_cast<TTree*>(file->Get(kTreeName));TTree* meta=dynamic_cast<TTree*>(file->Get(kMetadataName));
    if(!tree||!meta||meta->GetEntries()<=0)return Fail("missing reduced-tree event/segment metadata");
    Metadata row;if(!BindMetadata(meta,row,false,mError)){meta->ResetBranchAddresses();return false;}
    Long64_t total=0;
    for(Long64_t i=0;i<meta->GetEntries();++i) {
      if(meta->GetEntry(i)<=0){meta->ResetBranchAddresses();return Fail("cannot read reduced-tree segment metadata");}
      const SegmentKey key(row.segmentHigh,row.segmentLow);
      if(!row.producerId||!row.layout||!row.selection||!row.mixing||(!row.segmentHigh&&!row.segmentLow)||row.entries<0||!row.completed||row.producerId->empty()||
          *row.layout!=Layout(mSpecies)||*row.selection!=mSelection||row.mixing->empty()||
          (mRequireMixing&&*row.mixing!=mMixing)||!records.insert(std::make_pair(key,row.entries)).second||
          row.entries>std::numeric_limits<Long64_t>::max()-total) {
        meta->ResetBranchAddresses();return Fail("reduced-tree segment schema/provenance/identity mismatch");
      }
      total+=row.entries;
    }
    meta->ResetBranchAddresses();if(total!=tree->GetEntries())return Fail("reduced-tree entry count differs from completed segment metadata");
    mEvent=Event();mFrame=Frame();
    return femto_candidate_store::Bind(tree,mEvent,mCandidates,false,mError,true)&&BindExtra(tree,mFrame,mExtra,false,mError);
  }
  bool Open(Int_t index) {
    Close();mFileIndex=index;mFile=TFile::Open(mFiles[index].c_str(),"READ");
    if(!mFile||mFile->IsZombie())return Fail("cannot reopen reduced-tree input");
    mFileSegments.clear();mReadCounts.clear();mHaveSegment=false;
    if(!CheckFile(mFile,mTree,mFileSegments))return false;
    mEntries=mTree->GetEntries();if(mEntries!=mFileEntries[index])return Fail("reduced-tree input changed after preflight");return true;
  }
  bool FinishFile() {
    for(std::map<SegmentKey,Long64_t>::const_iterator i=mFileSegments.begin();i!=mFileSegments.end();++i)
      if(mReadCounts[i->first]!=i->second)return Fail("reduced-tree producer segment was not completely read");
    return true;
  }
  bool ValidateCandidates() {
    if(!mCandidates.species||!mCandidates.collection||!mCandidates.qaSelected||!mCandidates.quality||
       !mCandidates.floats||!mCandidates.ints||!mCandidates.physicalPresent||!mCandidates.qaPresent)
      return Fail("null reduced-tree canonical vector");
    const size_t n=mCandidates.species->size();
    if(n>std::numeric_limits<size_t>::max()/femto_candidate_store::kFloatFields||n>std::numeric_limits<size_t>::max()/femto_candidate_store::kIntFields||
      mCandidates.collection->size()!=n||mCandidates.qaSelected->size()!=n||mCandidates.quality->size()!=n||
      mCandidates.floats->size()!=n*femto_candidate_store::kFloatFields||mCandidates.ints->size()!=n*femto_candidate_store::kIntFields||
      mCandidates.physicalPresent->size()!=mSpecies.size()||mCandidates.qaPresent->size()!=mSpecies.size())return Fail("inconsistent reduced-tree canonical vector lengths");
    for(size_t i=0;i<mSpecies.size();++i)if(((*mCandidates.physicalPresent)[i]!=0&&(*mCandidates.physicalPresent)[i]!=1)||
      ((*mCandidates.qaPresent)[i]!=0&&(*mCandidates.qaPresent)[i]!=1))return Fail("invalid reduced-tree collection presence flag");return true;
  }
  Reader(const Reader&);Reader& operator=(const Reader&);
};

} // namespace femto_lambda_tree
#endif
