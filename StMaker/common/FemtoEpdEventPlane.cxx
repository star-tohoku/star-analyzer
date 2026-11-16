#include "FemtoEpdEventPlane.h"
#include "FemtoCandidateEventStoreConfig.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoEpdHit.h"
#include "StEpdUtil/StEpdGeom.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TVector3.h"
#include <algorithm>
#include <cstdio>
#include <fstream>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>

namespace femto_epd {
namespace {
void Require(bool condition,const std::string& message) { if(!condition)throw std::runtime_error("EPD: "+message); }
bool Finite(double x) { return std::isfinite(x); }
bool HexHash(const std::string& s) {
  if(s.size()!=64)return false;
  for(size_t i=0;i<s.size();++i)if(!((s[i]>='0'&&s[i]<='9')||(s[i]>='a'&&s[i]<='f')))return false;
  return true;
}
std::string Quote(const std::string& p) {
  std::string s="'";for(size_t i=0;i<p.size();++i)s+=p[i]=='\''?"'\\''":std::string(1,p[i]);return s+"'";
}
template<class T> std::vector<T> Vector(const YAML::Node& n,const std::string& key) {
  Require(n[key].IsSequence(),key+" must be a sequence");std::vector<T> v;
  for(size_t i=0;i<n[key].size();++i)v.push_back(n[key][i].as<T>());return v;
}
template<class T> std::vector<std::vector<T> > Matrix(const YAML::Node& n,const std::string& key) {
  Require(n[key].IsSequence(),key+" must be a sequence");std::vector<std::vector<T> > v;
  for(size_t i=0;i<n[key].size();++i){Require(n[key][i].IsSequence(),key+" rows must be sequences");std::vector<T> row;for(size_t j=0;j<n[key][i].size();++j)row.push_back(n[key][i][j].as<T>());v.push_back(row);}return v;
}
double Polynomial(const std::vector<double>& p,double x) {
  double y=0.;for(size_t i=0;i<p.size();++i)y+=p[i]*std::pow(x,double(i));return y;
}
TProfile2D* Profile2(TFile* f,const std::string& name,size_t runs) {
  TProfile2D* p=dynamic_cast<TProfile2D*>(f->Get(name.c_str()));
  Require(p&&p->GetNbinsX()>=9&&p->GetNbinsY()>=int(runs),"missing/wrong calibration axes: "+name);return p;
}
TProfile3D* Profile3(TFile* f,const std::string& name,size_t runs,int terms) {
  TProfile3D* p=dynamic_cast<TProfile3D*>(f->Get(name.c_str()));
  Require(p&&p->GetNbinsX()>=9&&p->GetNbinsY()>=terms&&p->GetNbinsZ()>=int(runs),"missing/wrong calibration axes: "+name);return p;
}
double Value(TProfile2D* p,int cent,int run) {
  const int b=p->GetBin(cent+1,run+1);const double v=p->GetBinContent(b);
  Require(p->GetBinEntries(b)>0&&Finite(v),std::string("empty/nonfinite calibration bin: ")+p->GetName());return v;
}
double Value(TProfile3D* p,int cent,int term,int run) {
  const int b=p->GetBin(cent+1,term+1,run+1);const double v=p->GetBinContent(b);
  Require(p->GetBinEntries(b)>0&&Finite(v),std::string("empty/nonfinite calibration bin: ")+p->GetName());return v;
}
bool Nonzero(double x,double y) { return Finite(x)&&Finite(y)&&x!=0.&&y!=0.; }
struct RandomScope {
  TRandom* old;
  explicit RandomScope(TRandom* r):old(gRandom){gRandom=r;}
  ~RandomScope(){gRandom=old;}
};
UInt_t EventSeed(UInt_t seed,Int_t run,Int_t event) {
  ULong64_t z=(ULong64_t(UInt_t(run))<<32)^UInt_t(event)^seed;
  z+=0x9e3779b97f4a7c15ULL;z=(z^(z>>30))*0xbf58476d1ce4e5b9ULL;
  z=(z^(z>>27))*0x94d049bb133111ebULL;z^=z>>31;
  const UInt_t value=UInt_t(z)^UInt_t(z>>32);return value?value:seed;
}
}

std::string Sha256File(const std::string& path) {
  Require(!path.empty(),"empty file path for SHA256");
  const std::string command="/usr/bin/sha256sum -- "+Quote(path);
  FILE* pipe=popen(command.c_str(),"r");Require(pipe!=0,"cannot run SHA256");
  char line[4096];const bool ok=fgets(line,sizeof(line),pipe)!=0;const int code=pclose(pipe);
  Require(ok&&code==0,"cannot checksum "+path);
  const std::string hash=std::string(line).substr(0,64);Require(HexHash(hash),"invalid SHA256 output for "+path);return hash;
}

void Settings::Load(const std::string& mainconf) {
  const YAML::Node main=YAML::LoadFile(mainconf);
  if(!main["eventPlane"].IsDefined())return;
  Require(main["eventPlane"].IsScalar(),"mainconf.eventPlane must be an explicit YAML path");
  path=femto_candidate_store::ConfigReference(mainconf,main["eventPlane"].as<std::string>());
  const YAML::Node n=YAML::LoadFile(path);femto_candidate_store::UniqueKeys(n,"EPD config");
  const char* keys[]={"writeDetailedLedger","source","harmonic","bins","invalidPolicy","recenterFile","recenterSha256","shiftFile","shiftSha256","runListFile","runListSha256","randomSeed","centralityThresholds","minNMip","maxTileWeight","shiftTerms","epdRowGroups","triggerIds","minVz","maxVz","vertexCenterX","vertexCenterY","maxVertexRadius","pileupMinCoefficients","pileupMaxCoefficients","tpcMinPt","tpcMaxPt","tpcEtaRanges","tpcMinNHitsFit","tpcMinFitRatio","tpcMaxDca"};
  std::set<std::string> required;
  for(size_t i=0;i<sizeof(keys)/sizeof(keys[0]);++i){required.insert(keys[i]);Require(n[keys[i]].IsDefined(),std::string("missing setting ")+keys[i]);}
  for(YAML::const_iterator i=n.begin();i!=n.end();++i)Require(required.count(i->first.as<std::string>()),"unknown setting "+i->first.as<std::string>());
  Require(n["source"].as<std::string>()=="xialeiEpdPsi1"&&n["invalidPolicy"].as<std::string>()=="reject","unsupported source/invalid policy");
  writeDetailedLedger=n["writeDetailedLedger"].as<bool>();
  harmonic=n["harmonic"].as<int>();bins=n["bins"].as<int>();shiftTerms=n["shiftTerms"].as<int>();
  randomSeed=n["randomSeed"].as<UInt_t>();
  Require(harmonic==1&&bins==8&&shiftTerms>0&&randomSeed!=0,"invalid harmonic/bins/shiftTerms/randomSeed");
  recenterFile=femto_candidate_store::ConfigReference(mainconf,n["recenterFile"].as<std::string>());
  shiftFile=femto_candidate_store::ConfigReference(mainconf,n["shiftFile"].as<std::string>());
  runListFile=femto_candidate_store::ConfigReference(mainconf,n["runListFile"].as<std::string>());
  recenterSha256=n["recenterSha256"].as<std::string>();shiftSha256=n["shiftSha256"].as<std::string>();runListSha256=n["runListSha256"].as<std::string>();
  Require(HexHash(recenterSha256)&&HexHash(shiftSha256)&&HexHash(runListSha256),"explicit lowercase SHA256 fingerprints are required");
#define EPD_DOUBLE(name) name=n[#name].as<double>();Require(Finite(name),"nonfinite " #name)
  EPD_DOUBLE(minNMip);EPD_DOUBLE(maxTileWeight);EPD_DOUBLE(minVz);EPD_DOUBLE(maxVz);
  EPD_DOUBLE(vertexCenterX);EPD_DOUBLE(vertexCenterY);EPD_DOUBLE(maxVertexRadius);
  EPD_DOUBLE(tpcMinPt);EPD_DOUBLE(tpcMaxPt);EPD_DOUBLE(tpcMinFitRatio);EPD_DOUBLE(tpcMaxDca);
#undef EPD_DOUBLE
  tpcMinNHitsFit=n["tpcMinNHitsFit"].as<int>();
  Require(minNMip>=0&&maxTileWeight>=minNMip&&minVz<maxVz&&maxVertexRadius>0&&tpcMinPt>=0&&tpcMinPt<tpcMaxPt&&tpcMinFitRatio>=0&&tpcMinFitRatio<=1&&tpcMaxDca>=0&&tpcMinNHitsFit>=0,"invalid cut ranges");
  centralityThresholds=Vector<int>(n,"centralityThresholds");triggerIds=Vector<int>(n,"triggerIds");
  Require(centralityThresholds.size()==9&&!triggerIds.empty(),"nine calibration centrality thresholds and triggers are required");
  for(size_t i=0;i<centralityThresholds.size();++i)Require(centralityThresholds[i]>=0&&(i==0||centralityThresholds[i]>centralityThresholds[i-1]),"unordered calibration centrality thresholds");
  epdRowGroups=Matrix<int>(n,"epdRowGroups");tpcEtaRanges=Matrix<double>(n,"tpcEtaRanges");
  Require(epdRowGroups.size()==4&&tpcEtaRanges.size()==2,"four EPD row groups and two TPC eta groups are required");
  std::set<int> rows;
  for(size_t i=0;i<epdRowGroups.size();++i){Require(epdRowGroups[i].size()==2&&epdRowGroups[i][0]>0&&epdRowGroups[i][0]<=epdRowGroups[i][1]&&epdRowGroups[i][1]<=16,"invalid EPD row group");for(int r=epdRowGroups[i][0];r<=epdRowGroups[i][1];++r)Require(rows.insert(r).second,"overlapping EPD row groups");}
  for(size_t i=0;i<tpcEtaRanges.size();++i)Require(tpcEtaRanges[i].size()==2&&Finite(tpcEtaRanges[i][0])&&Finite(tpcEtaRanges[i][1])&&tpcEtaRanges[i][0]<tpcEtaRanges[i][1],"invalid TPC eta range");
  pileupMinCoefficients=Vector<double>(n,"pileupMinCoefficients");pileupMaxCoefficients=Vector<double>(n,"pileupMaxCoefficients");
  Require(pileupMinCoefficients.size()==5&&pileupMaxCoefficients.size()==5,"pileup polynomials require five coefficients");
  for(size_t i=0;i<5;++i)Require(Finite(pileupMinCoefficients[i])&&Finite(pileupMaxCoefficients[i]),"nonfinite pileup coefficient");
  const std::string method="EPD East AB corrected first-order plane; normalized tile-weight Q1; recenter+shift; original all-plane validity; independent event-keyed SplitMix64/TRandom3 geometry sampling; exact Float_t downstream binning\n";
  definition=method+YAML::Dump(n);
  YAML::Node pairing;
  for(YAML::const_iterator i=n.begin();i!=n.end();++i)if(i->first.as<std::string>()!="writeDetailedLedger")pairing[i->first.as<std::string>()]=i->second;
  pairingDefinition=method+YAML::Dump(pairing);
  enabled=true;
}

Provider::Provider():mRecenter(0),mShift(0),mShiftSin(0),mShiftCos(0) {
  for(int h=0;h<2;++h){for(int g=0;g<2;++g)for(int c=0;c<2;++c)mTpc[h][g][c]=0;for(int g=0;g<7;++g)for(int c=0;c<2;++c)mEpd[h][g][c]=0;}
}
Provider::~Provider(){if(mShift){mShift->Close();delete mShift;}if(mRecenter){mRecenter->Close();delete mRecenter;}}
void Provider::Init(const Settings& settings,bool loadCalibration) {
  mSettings=settings;Require(settings.enabled,"provider initialized without configuration");
  if(!loadCalibration)return;
  Require(Sha256File(settings.recenterFile)==settings.recenterSha256,"recenter SHA256 mismatch");
  Require(Sha256File(settings.shiftFile)==settings.shiftSha256,"shift SHA256 mismatch");
  Require(Sha256File(settings.runListFile)==settings.runListSha256,"run-list SHA256 mismatch");
  std::ifstream in(settings.runListFile.c_str());Require(bool(in),"unreadable run list");std::string line;std::set<int> seen;
  while(std::getline(in,line)){
    const size_t at=line.find('#');if(at!=std::string::npos)line.erase(at);
    if(line.find_first_not_of(" \t\r\n")==std::string::npos)continue;
    std::istringstream s(line);int run;std::string extra;
    Require(bool(s>>run)&&!(s>>extra)&&run>0&&seen.insert(run).second,"invalid/duplicate run-list row");
    mRuns.push_back(run);
  }
  Require(!in.bad()&&!mRuns.empty(),"empty or failed run-list read");
  TDirectory::TContext directory(gDirectory);
  mRecenter=TFile::Open(settings.recenterFile.c_str(),"READ");mShift=TFile::Open(settings.shiftFile.c_str(),"READ");
  Require(mRecenter&&!mRecenter->IsZombie()&&!mRecenter->TestBit(TFile::kRecovered)&&mShift&&!mShift->IsZombie()&&!mShift->TestBit(TFile::kRecovered),"invalid calibration ROOT file");
  for(int h=0;h<2;++h)for(int g=0;g<2;++g)for(int c=0;c<2;++c)mTpc[h][g][c]=Profile2(mRecenter,Form("TPCQ%c_%c_recen%s",c?'y':'x',g?'B':'A',h?"_psi2":""),mRuns.size());
  for(int h=0;h<2;++h)for(int g=0;g<7;++g)for(int c=0;c<2;++c)mEpd[h][g][c]=Profile2(mRecenter,Form("EPDQ%c_ring%d_recen%s",c?'y':'x',g,h?"_psi2":""),mRuns.size());
  mShiftSin=Profile3(mShift,"EPDshiftpar_sin_ring4",mRuns.size(),settings.shiftTerms);
  mShiftCos=Profile3(mShift,"EPDshiftpar_cos_ring4",mRuns.size(),settings.shiftTerms);
}

Result Provider::Calculate(StPicoDst* pico) const {
  Require(mRecenter&&mShift,"calculation requires loaded calibration");
  Result r;if(!pico||!pico->event())return r;const StPicoEvent* event=pico->event();
  const Settings& s=mSettings;bool trigger=false;for(size_t i=0;i<s.triggerIds.size();++i)trigger=trigger||event->isTrigger(s.triggerIds[i]);
  if(!trigger){r.status=1;return r;}
  const TVector3 pv=event->primaryVertex();const Float_t vx=pv.X(),vy=pv.Y(),vz=pv.Z();
  if(!Finite(vx)||!Finite(vy)||!Finite(vz)||vz<s.minVz||vz>s.maxVz||(vx-s.vertexCenterX)*(vx-s.vertexCenterX)+(vy-s.vertexCenterY)*(vy-s.vertexCenterY)>s.maxVertexRadius*s.maxVertexRadius){r.status=2;return r;}
  const std::vector<int>::const_iterator run=std::find(mRuns.begin(),mRuns.end(),event->runId());Require(run!=mRuns.end(),"run not covered by calibration run list");r.runIndex=int(run-mRuns.begin());
  // The calibration multiplicity is the unqualified primary-track count,
  // independently of the centrality used for the Lambda pair pool.
  int primary=0;for(unsigned i=0;i<pico->numberOfTracks();++i){const StPicoTrack* t=pico->track(i);if(t&&t->isPrimary())++primary;}
  const double tof=event->nBTOFMatch();if(primary<Polynomial(s.pileupMinCoefficients,tof)||primary>Polynomial(s.pileupMaxCoefficients,tof)){r.status=3;return r;}
  for(int i=8;i>=0;--i)if(primary>=s.centralityThresholds[i]){r.calibrationCent=i;break;}
  if(r.calibrationCent<0){r.status=4;return r;}
  // Preserve the source event validity, including both TPC harmonics and
  // all EPD groups, although the reported angle uses EPD-AB Q1 only.
  double tpc[2][2][2]={{{0}}},tw[2]={0.,0.};
  for(unsigned i=0;i<pico->numberOfTracks();++i){
    const StPicoTrack* t=pico->track(i);if(!t||!t->isPrimary())continue;
    const Float_t pt=t->pMom().Perp(),phi=t->pMom().Phi(),eta=t->pMom().PseudoRapidity(),dca=t->gDCA(pv).Mag();
    const int fit=t->nHitsFit(),maximum=t->nHitsMax();const Float_t quality=Float_t(fit)/Float_t(maximum);
    if(!Finite(pt)||!Finite(phi)||!Finite(eta)||!Finite(dca)||std::fabs(dca)>s.tpcMaxDca||fit<s.tpcMinNHitsFit||quality<s.tpcMinFitRatio||!(pt>s.tpcMinPt&&pt<s.tpcMaxPt))continue;
    int group=-1;for(int g=0;g<2;++g)if(eta>s.tpcEtaRanges[g][0]&&eta<s.tpcEtaRanges[g][1]){group=g;break;}if(group<0)continue;
    for(int h=0;h<2;++h){tpc[h][group][0]+=pt*std::cos(double(h+1)*phi);tpc[h][group][1]+=pt*std::sin(double(h+1)*phi);}++tw[group];
  }
  for(int h=0;h<2;++h)for(int g=0;g<2;++g){if(tw[g]==0||!Nonzero(tpc[h][g][0],tpc[h][g][1])){r.status=5;return r;}for(int c=0;c<2;++c)tpc[h][g][c]=tpc[h][g][c]/tw[g]-Value(mTpc[h][g][c],r.calibrationCent,r.runIndex);if(!Nonzero(tpc[h][g][0],tpc[h][g][1])){r.status=5;return r;}}
  double q[2][7][2]={{{0}}},weight[7]={0.};
  TRandom3 random(EventSeed(s.randomSeed,event->runId(),event->eventId()));RandomScope randomScope(&random);
  // StEpdGeom captures gRandom in its constructor. Its lifetime must stay
  // inside this scope so neither sampling nor destruction touches the
  // caller's random generator.
  StEpdGeom geometry;
  for(unsigned i=0;i<pico->numberOfEpdHits();++i){
    const StPicoEpdHit* hit=pico->epdHit(i);if(!hit||hit->id()>0)continue;
    const int row=hit->row();int group=-1;for(int g=0;g<4;++g)if(row>=s.epdRowGroups[g][0]&&row<=s.epdRowGroups[g][1]){group=g;break;}if(group<0)continue;
    // The source samples geometry before rejecting sub-threshold MIP hits.
    const TVector3 point=geometry.RandomPointOnTile(hit->id())-pv;const double phi=point.Phi(),mip=hit->nMIP();
    if(!Finite(phi)||!Finite(mip)){r.status=6;return r;}if(mip<s.minNMip)continue;const double w=std::min(mip,s.maxTileWeight);
    const int groups[]={group,group<2?4:5,6};
    for(int j=0;j<3;++j){const int g=groups[j];weight[g]+=w;for(int h=0;h<2;++h){q[h][g][0]+=w*std::cos((h+1)*phi);q[h][g][1]+=w*std::sin((h+1)*phi);}}
    if(group<2)++r.hits;
  }
  for(int h=0;h<2;++h)for(int g=0;g<7;++g){
    if(weight[g]==0||!Nonzero(q[h][g][0],q[h][g][1])){r.status=6;return r;}
    for(int c=0;c<2;++c)q[h][g][c]/=weight[g];
    if(h==0&&g==4){r.rawQx=q[h][g][0];r.rawQy=q[h][g][1];r.weight=weight[g];}
    for(int c=0;c<2;++c)q[h][g][c]-=Value(mEpd[h][g][c],r.calibrationCent,r.runIndex);
    if(!Nonzero(q[h][g][0],q[h][g][1])){r.status=7;return r;}
  }
  r.recenteredQx=q[0][4][0];r.recenteredQy=q[0][4][1];const double recentered=std::atan2(r.recenteredQy,r.recenteredQx);r.psi1=recentered;
  for(int i=0;i<s.shiftTerms;++i)r.psi1+=(2.0/(i+1))*(-Value(mShiftSin,r.calibrationCent,i,r.runIndex)*std::cos((i+1)*recentered)+Value(mShiftCos,r.calibrationCent,i,r.runIndex)*std::sin((i+1)*recentered));
  // The source wraps with a single correction in each direction.
  const double pi=TMath::Pi();if(r.psi1 < -pi)r.psi1+=2*pi;if(r.psi1>pi)r.psi1-=2*pi;
  r.bin=XialeiBin(r.psi1);r.valid=r.bin>=0;r.status=r.valid?0:8;return r;
}
} // namespace femto_epd
