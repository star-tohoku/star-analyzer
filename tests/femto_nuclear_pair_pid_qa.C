int femto_nuclear_pair_pid_qa(const char*,const char*,const char*);
#ifndef __CINT__
#include "StMaker/common/FemtoLambdaLegacy.h"
#include "StMaker/common/FemtoNuclearPairPidSparse.h"
#include "FemtoNuclearPidPairQaCandidate.h"
#include "FemtoCandidate.h"
#include "HistManager.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"
#include "yaml-cpp/yaml.h"
#include "TClonesArray.h"
#include "TH1.h"
#include "TH2.h"
#include "TAxis.h"
#include "TFile.h"
#include "TNamed.h"
#include "TSystem.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <sys/resource.h>

namespace NuclearPairQaTest {
unsigned long long checks=0;
int Cells(const TH1& h){int n=h.GetNbinsX()+2;if(h.GetDimension()>1)n*=h.GetNbinsY()+2;if(h.GetDimension()>2)n*=h.GetNbinsZ()+2;return n;}
void Need(bool ok,const std::string& message){++checks;if(!ok)throw std::runtime_error(message);}
std::string Num(int v){std::ostringstream s;s<<v;return s.str();}
std::string Base(const std::string& path){const size_t p=path.rfind("/config/");Need(p!=std::string::npos,"explicit main under config");return path.substr(0,p)+"/config/";}
std::string Ref(const std::string& base,const YAML::Node& node,const char* key){const std::string s=node[key].as<std::string>();return s[0]=='/'?s:base+s;}
std::vector<std::string> Keys(const std::string& path){const YAML::Node d=YAML::LoadFile(path);std::vector<std::string> out;for(YAML::const_iterator i=d["histograms"].begin();i!=d["histograms"].end();++i)out.push_back(i->first.as<std::string>());return out;}
void SameHist(const TH1& a,const TH1& b,const std::string& label){
  Need(std::string(a.ClassName())==b.ClassName()&&a.GetDimension()==b.GetDimension()&&Cells(a)==Cells(b),label+" shape");
  Need(a.GetEntries()==b.GetEntries(),label+" entries");
  for(int cell=0;cell<Cells(a);++cell){Need(a.GetBinContent(cell)==b.GetBinContent(cell),label+" exact content incl flow");Need(a.GetBinError(cell)==b.GetBinError(cell),label+" exact errors incl flow");}
}
void OldEqual(FemtoLambdaLegacy& on,FemtoLambdaLegacy& off,const std::string& base,const YAML::Node& main){
  const char* refs[]={"femtoHist","nuclearHist"};
  for(int which=0;which<2;++which){const std::vector<std::string> keys=Keys(Ref(base,main,refs[which]));HistManager* a=which?on.NuclearHistograms():on.RootHistograms();HistManager* b=which?off.NuclearHistograms():off.RootHistograms();
    for(unsigned k=0;k<keys.size();++k){TH1* x=a->Get(keys[k].c_str());TH1* y=b->Get(keys[k].c_str());Need(x&&y,"old histogram present");SameHist(*x,*y,keys[k]);}}
}
double Kstar(const FemtoCandidate& a,const FemtoCandidate& b){TLorentzVector x=a.P4();x.Boost(-(a.P4()+b.P4()).BoostVector());return x.P();}
void Projection(TH2& qa,const TH2& canonical,int cent){
  TH1D* p=qa.ProjectionX("nuclear_pair_qa_x",0,qa.GetNbinsY()+1,"e");Need(p!=0,"QA ProjectionX");p->SetDirectory(0);
  double entries=0;for(int x=0;x<=qa.GetNbinsX()+1;++x){const int ref=canonical.GetBin(x,cent+1);entries+=canonical.GetBinContent(ref);Need(p->GetBinContent(x)==canonical.GetBinContent(ref),"QA to existing CF closure incl flow");Need(std::fabs(std::pow(p->GetBinError(x),2)-std::pow(canonical.GetBinError(ref),2))<1.e-9,"QA to CF variance closure");}Need(qa.GetEntries()==entries,"QA unit-count entries closure including flow");delete p;
}
double QaEntries(HistManager& hist,const std::vector<std::string>& keys){double count=0;for(unsigned i=0;i<keys.size();++i)count+=hist.Get(keys[i].c_str())->GetEntries();return count;}
void Half(const TH2& all,const TH2& low,const TH2& high){
  Need(all.GetEntries()==low.GetEntries()+high.GetEntries(),"QA signal entries partition");for(int b=0;b<Cells(all);++b){Need(all.GetBinContent(b)==low.GetBinContent(b)+high.GetBinContent(b),"QA signal content partition");Need(std::fabs(std::pow(all.GetBinError(b),2)-std::pow(low.GetBinError(b),2)-std::pow(high.GetBinError(b),2))<1.e-9,"QA signal variance partition");}
}
void SparseSlice(THnSparseD& h,const TH2& qa,int axis,int cent,int repr,int mode,int stage,int region){
  for(int d=0;d<3;++d)h.GetAxis(d)->SetRange(0,h.GetAxis(d)->GetNbins()+1);
  h.GetAxis(3)->SetRange(cent+1,cent+1);h.GetAxis(4)->SetRange(repr+1,repr+1);h.GetAxis(5)->SetRange(mode+1,mode+1);h.GetAxis(7)->SetRange(stage+1,stage+1);
  const int lo[]={1,3,0,1,2},hi[]={2,3,0,1,2};h.GetAxis(6)->SetRange(lo[region]+1,hi[region]+1);
  TH2D* p=h.Projection(axis,0,"EA");Need(p!=0,"joint slice");p->SetDirectory(0);
  Need(p->GetNbinsX()==qa.GetNbinsX()&&p->GetNbinsY()==qa.GetNbinsY(),"joint/dense axes shape");
  for(int b=0;b<Cells(qa);++b){Need(p->GetBinContent(b)==qa.GetBinContent(b),"joint/dense flow content closure");Need(std::fabs(std::pow(p->GetBinError(b),2)-std::pow(qa.GetBinError(b),2))<1.e-9,"joint/dense variance closure");}delete p;
}
struct Pico {
  TClonesArray events,tracks,tofs;TClonesArray* arrays[StPicoArrays::NAllPicoArrays];StPicoDst dst;
  Pico():events("StPicoEvent"),tracks("StPicoTrack"),tofs("StPicoBTofPidTraits"){for(int i=0;i<StPicoArrays::NAllPicoArrays;++i)arrays[i]=0;arrays[StPicoArrays::Event]=&events;arrays[StPicoArrays::Track]=&tracks;arrays[StPicoArrays::BTofPidTraits]=&tofs;new(events[0])StPicoEvent();StPicoDst::set(arrays);}
  ~Pico(){StPicoDst::unset();}
  StPicoTrack* Add(FemtoLambdaLegacy& helper,int type,double r,double pull,bool tof,double m2){
    const int index=tracks.GetEntriesFast();StPicoTrack* t=new(tracks[index])StPicoTrack();t->setId(index+100);t->setNHitsFit(25);t->setNHitsDedx(20);t->setPrimaryMomentum(r,0,0);t->setGlobalMomentum(r,0,0);
    const double p=t->pMom().Mag(),a=helper.NuclearDedxPull(type,p,1.),b=helper.NuclearDedxPull(type,p,2.);Need(std::isfinite(a)&&b>a,"synthetic inverse PID response");t->setDedx((1.+(pull-a)/(b-a))*1.e-6);
    if(tof){const int ti=tofs.GetEntriesFast();StPicoBTofPidTraits* b=new(tofs[ti])StPicoBTofPidTraits();b->setTrackIndex(index);b->setBTofMatchFlag(1);b->setBeta(float(1./std::sqrt(1.+m2/(p*p))));t->setBTofPidTraitsIndex(ti);}return t;
  }
};
bool Approved(int type,double r,bool tof,double m2){
  bool domain=type==2?(r>=.75&&r<2.5):(r>.61&&r<5.);
  if(!std::isfinite(r)||r<=0)return false;
  if(r>=2){const double lo=type==2?1.78285:3.268773,hi=type==2?2.06485:3.693147;domain=domain&&tof&&std::isfinite(m2)&&m2>lo&&m2<hi;}return domain;
}
void Collect(FemtoLambdaLegacy& on,FemtoLambdaLegacy& off){
  Pico fixture;
  fixture.Add(on,3,1.02,0,false,0);fixture.Add(on,3,1.02,4,false,0);fixture.Add(on,2,1.02,0,false,0);
  fixture.Add(on,3,.65,0,false,0);fixture.Add(on,3,.9,0,false,0);fixture.Add(on,3,2.2,0,false,0);
  fixture.Add(on,3,2.2,0,true,3.48);fixture.Add(on,3,2.2,0,true,1.92);fixture.Add(on,3,3.3,0,true,3.48);
  fixture.Add(on,3,1.02,0,false,0)->setNHitsDedx(14);fixture.Add(on,3,1.02,0,false,0)->setGlobalMomentum(.05,0,1.);
  FemtoCandidateStore a,b;FemtoNuclearPidPairQaStore qa,disabled;
  on.CollectNuclei(&fixture.dst,42,8,a,&qa);off.CollectNuclei(&fixture.dst,42,8,b,&disabled);Need(disabled.empty(),"off does not create QA candidates");
  const char* species[]={"deuteron","triton","he3","he4"};
  for(int s=0;s<4;++s){Need(a[species[s]].size()==b[species[s]].size(),"collection ON/OFF count invariant");for(unsigned i=0;i<a[species[s]].size();++i){const FemtoCandidate& x=a[species[s]][i];const FemtoCandidate& y=b[species[s]][i];Need(x.trk.trackIndex==y.trk.trackIndex&&x.px==y.px&&x.py==y.py&&x.pz==y.pz&&x.energy==y.energy,"collection ON/OFF candidate invariant");}}
  unsigned droppedOwn=0,droppedWinner=0;
  for(int s=2;s<4;++s){const std::vector<FemtoNuclearPidPairQaCandidate>& rows=qa[species[s]];std::map<int,bool> got;for(unsigned i=0;i<rows.size();++i){const FemtoCandidate& c=rows[i].candidate;Need(got.insert(std::make_pair(c.trk.trackIndex,rows[i].selected)).second,"one QA candidate per target/track");Need(c.eventIndex==42,"source event provenance");Need(std::fabs(c.P4().P()/2.-fixture.dst.track(c.trk.trackIndex)->pMom().Mag())<1.e-6,"QA raw R from one Z correction");bool selected=false;for(unsigned j=0;j<a[species[s]].size();++j)selected=selected||a[species[s]][j].trk.trackIndex==c.trk.trackIndex;Need(selected==rows[i].selected,"original-double selected flag agrees with actual selected store");const double pull=s==2?c.trk.nSigmaHe3:c.trk.nSigmaHe4;if(!selected&&std::fabs(pull)>2)++droppedOwn;if(!selected)++droppedWinner;}
    for(unsigned i=0;i<fixture.dst.numberOfTracks();++i){StPicoTrack* t=fixture.dst.track(i);bool tof=false;double m2=0;if(t->isTofTrack()&&t->bTofPidTraitsIndex()>=0){const double beta=fixture.dst.btofPidTraits(t->bTofPidTraitsIndex())->btofBeta();tof=std::isfinite(beta)&&beta>0;if(tof)m2=t->pMom().Mag2()*(1./(beta*beta)-1.);}
      const bool expected=t->nHitsDedx()>=15&&t->gPt()>=.1&&Approved(s,t->pMom().Mag(),tof,m2);Need((got.count(i)>0)==expected,"prePID retains quality/domain/TOF and removes own-pull/winner gates");}
  }
  Need(droppedOwn>0&&droppedWinner>0,"prePID has rejected-own/winner diagnostic population");
  const double edges[]={.61,.7,.75,.85,.95,1.75,1.8,2.,2.5,3.2,5.};
  for(int s=2;s<4;++s)for(unsigned i=0;i<sizeof(edges)/sizeof(edges[0]);++i)for(int side=-1;side<=1;++side)for(int tof=0;tof<2;++tof){const double r=side<0?std::nextafter(edges[i],-std::numeric_limits<double>::infinity()):side>0?std::nextafter(edges[i],std::numeric_limits<double>::infinity()):edges[i];const double m2=s==2?1.92:3.48;Need(on.PassNuclearPidQaDomain(s,r,tof,m2)==Approved(s,r,tof,m2),"prePID exact/nextafter boundary oracle");}
}
}

int femto_nuclear_pair_pid_qa(const char* onMain,const char* offMain,const char* outputDir){
  using namespace NuclearPairQaTest;
  try{
    checks=0;Need(onMain&&*onMain&&offMain&&*offMain&&outputDir&&*outputDir,"explicit ON/OFF/output paths");
    const YAML::Node main=YAML::LoadFile(onMain),old=YAML::LoadFile(offMain);Need(bool(main["nuclearPidPairQaHist"])&&!old["nuclearPidPairQaHist"],"ON/OFF explicit feature toggle");
    for(YAML::const_iterator it=main.begin();it!=main.end();++it){const std::string key=it->first.as<std::string>();if(key!="nuclearPidPairQaHist")Need(old[key]&&old[key].as<std::string>()==it->second.as<std::string>(),"only QA reference differs");}
    Need(gSystem->AccessPathName(outputDir),"new output directory required");Need(gSystem->mkdir(outputDir,true)==0,"output mkdir");
    FemtoLambdaLegacy on,off;Need(on.Init(onMain,std::cerr)&&off.Init(offMain,std::cerr),"ON/OFF helpers Init");
    Need(on.NuclearPidPairQaHistograms()&&!off.NuclearPidPairQaHistograms(),"separate optional QA managers");Need(on.NuclearPidPairSparse()&&!off.NuclearPidPairSparse(),"separate optional joint helpers");
    const std::string base=Base(onMain);const YAML::Node maker=YAML::LoadFile(Ref(base,main,"maker"));const double mean=maker["lambdaSignalMean"].as<double>(),width=maker["lambdaSignalSigma"].as<double>()*maker["lambdaSignalNSigma"].as<double>();
    FemtoCandidate l,n;l.speciesKey="lambda";l.source=kFemtoCandResonance;l.eventIndex=10;l.reso.dau1EventIndex=10;l.reso.dau2EventIndex=10;l.reso.dau1Index=2;l.reso.dau2Index=3;TLorentzVector lp;lp.SetXYZM(.32,.14,.8,maker["lambdaPairMass"].as<double>());l.SetP4(lp);
    const float masses[]={float(mean-.5*width),float(mean+.5*width),float(mean-2*width),float(mean+2*width),float(mean+5*width)};
    unsigned pairCalls=0;
    for(int mode=0;mode<2;++mode)for(int cent=0;cent<9;++cent)for(int type=2;type<4;++type)for(int v=0;v<3;++v){
      n.speciesKey=type==2?"he3":"he4";n.source=kFemtoCandTrack;n.eventIndex=10+mode;n.trk.trackIndex=mode?2:9;n.charge=2;n.SetP4(on.NuclearP4(TVector3(v==2?7.:.45,v==2?0.:.17,.28),type));
      n.trk.nSigmaHe3=type==2?(v==0?-15.f:v==1?-1.25f:15.f):7.375f;n.trk.nSigmaHe4=type==3?(v==0?-16.f:v==1?1.375f:16.f):-8.625f;
      for(unsigned mass=0;mass<sizeof(masses)/sizeof(masses[0]);++mass){l.reso.invMass=masses[mass];const double ks=v==0?-.05:v==1?.217:1.07,qlab=(l.P4().Vect()-n.P4().Vect()).Mag();Need(on.FillPair(l,n,ks,qlab,cent,mode!=0)&&off.FillPair(l,n,ks,qlab,cent,mode!=0),"same accepted selected pair");const bool filled=on.FillPrePidNuclearPairQa(l,n,ks,cent,mode!=0);Need(on.MassRegion(l.reso.invMass)==0||filled,"prePID independent stage fill");++pairCalls;}
    }
    HistManager* qa=on.NuclearPidPairQaHistograms();THnSparseD* joint=on.NuclearPidPairSparse()->Histogram();Need(joint!=0,"joint exists");
    const char* labels[]={"4He","3He","3HeAs4He"},*canonical[]={"he4","he3","he3_as_he4"},*regions[]={"","_SBPos","_SBNeg","_signalLow","_signalHigh"},*rootRegions[]={"_signal","_rightSB","_leftSB","_signalLow","_signalHigh"};
    unsigned maps=0;
    for(int stage=0;stage<2;++stage)for(int mode=0;mode<2;++mode)for(int repr=0;repr<3;++repr)for(int cent=0;cent<9;++cent)for(int variable=0;variable<2;++variable){TH2* split[5]={0,0,0,0,0};
      for(int region=0;region<5;++region){const std::string key=std::string(variable?"hKstarPOverQ":"hKstarNSigma")+(stage?"PrePid":"")+"_"+(mode?"Mixed_":"")+labels[repr]+regions[region]+"_CentBin"+Num(cent);TH2* h=dynamic_cast<TH2*>(qa->Get(key.c_str()));Need(h&&std::string(h->ClassName())=="TH2D","all QA keys TH2D");Need(h->GetNbinsX()==200&&h->GetXaxis()->GetXmin()==0&&h->GetXaxis()->GetXmax()==1&&h->GetNbinsY()==120&&h->GetYaxis()->GetXmin()==(variable?0:-12)&&h->GetYaxis()->GetXmax()==(variable?6:12),"QA axes independent contract");split[region]=h;
        TH2* cf=dynamic_cast<TH2*>(on.RootHistograms()->Get((std::string("hKstar")+(mode?"ME":"SE")+"VsCent_lambda_"+canonical[repr]+rootRegions[region]).c_str()));Need(cf!=0,"old CF reference");Projection(*h,*cf,cent);SparseSlice(*joint,*h,variable?1:2,cent,repr,mode,stage,region);++maps;}
      Half(*split[0],*split[3],*split[4]);
    }
    // Disjoint fixture pulls distinguish the selected original He3 identity from
    // the artificial He4 mass representation. Neither momentum nor PID is remapped.
    const char* physical[]={"3He","3HeAs4He"};for(int i=0;i<2;++i){TH2* h=dynamic_cast<TH2*>(qa->Get((std::string("hKstarNSigma_")+physical[i]+"_CentBin0").c_str()));Need(h->Integral(0,h->GetNbinsX()+1,h->GetYaxis()->FindBin(-1.25),h->GetYaxis()->FindBin(-1.25))==2,"mapped He3 retains original He3 pull");Need(h->Integral(0,h->GetNbinsX()+1,h->GetYaxis()->FindBin(-8.625),h->GetYaxis()->FindBin(-8.625))==0,"mapped He3 never silently uses He4 pull");}
    const double rawOracle=std::sqrt(.45*.45+.17*.17+.28*.28);
    for(int stage=0;stage<2;++stage)for(int repr=0;repr<3;++repr){const std::string tail=std::string(stage?"PrePid_":"_")+labels[repr]+"_CentBin0";TH2* r=dynamic_cast<TH2*>(qa->Get(("hKstarPOverQ"+tail).c_str()));const int expectedR=r->GetYaxis()->FindBin(rawOracle),doubleZ=r->GetYaxis()->FindBin(2.*rawOracle);Need(expectedR!=doubleZ,"R fixture distinguishes double Z correction");Need(r->Integral(0,r->GetNbinsX()+1,expectedR,expectedR)==4,"independent raw-R oracle in dense maps");Need(r->Integral(0,r->GetNbinsX()+1,doubleZ,doubleZ)==0,"raw R is not physical He momentum");}
    TH2* he4Ns=dynamic_cast<TH2*>(qa->Get("hKstarNSigma_4He_CentBin0"));Need(he4Ns->Integral(0,he4Ns->GetNbinsX()+1,he4Ns->GetYaxis()->FindBin(1.375),he4Ns->GetYaxis()->FindBin(1.375))==2,"He4 uses own He4 pull");Need(he4Ns->Integral(0,he4Ns->GetNbinsX()+1,he4Ns->GetYaxis()->FindBin(7.375),he4Ns->GetYaxis()->FindBin(7.375))==0,"He4 never uses competing He3 pull");
    const double before=joint->GetEntries();l.reso.invMass=float(mean);n.speciesKey="he3";n.eventIndex=10;n.trk.trackIndex=2;
    Need(!on.FillPair(l,n,.2,.3,8,false)&&!on.FillPrePidNuclearPairQa(l,n,.2,8,false),"shared-track rejects both stages");
    n.trk.trackIndex=9;Need(!on.FillPair(l,n,.2,.3,8,true)&&!on.FillPrePidNuclearPairQa(l,n,.2,8,true),"same-event ME rejects both stages");
    Need(!on.FillPair(l,n,std::numeric_limits<double>::quiet_NaN(),.3,8,false)&&!on.FillPrePidNuclearPairQa(l,n,std::numeric_limits<double>::quiet_NaN(),8,false),"nonfinite pair rejects both stages");Need(joint->GetEntries()==before,"rejected pairs never partially fill joint");
    const std::vector<std::string> qaKeys=Keys(Ref(base,main,"nuclearPidPairQaHist"));const double denseBefore=QaEntries(*qa,qaKeys),jointBefore=joint->GetEntries();
    n.trk.nSigmaHe3=std::numeric_limits<float>::quiet_NaN();Need(on.FillPair(l,n,.2,.3,8,false)&&off.FillPair(l,n,.2,.3,8,false),"invalid QA pull cannot reject or change regular CF pair");on.FillPrePidNuclearPairQa(l,n,.2,8,false);Need(QaEntries(*qa,qaKeys)==denseBefore&&joint->GetEntries()==jointBefore,"invalid QA tuple atomically skips both dense projections and joint");
    Collect(on,off);OldEqual(on,off,base,main);
    for(int d=0;d<8;++d)joint->GetAxis(d)->SetRange(0,0);
    {TFile file((std::string(outputDir)+"/qa_fixture.root").c_str(),"CREATE");Need(!file.IsZombie()&&on.Write(&file),"QA file write");}
    {TFile file((std::string(outputDir)+"/qa_fixture.root").c_str(),"READ");Need(file.Get("hNuclearPidPairJoint")&&file.Get("NuclearPidPairJointMetadata"),"joint + metadata readback");for(int mode=0;mode<2;++mode){const std::string key=std::string("hKstarNSigma_")+(mode?"Mixed_":"")+"3HeAs4He_CentBin8";Need(file.Get((std::string(mode?"mix/":"true/")+key).c_str())&&!file.Get((std::string(mode?"true/":"mix/")+key).c_str()),"QA SE/ME directory readback");}}
    struct rusage usage;getrusage(RUSAGE_SELF,&usage);std::ofstream summary((std::string(outputDir)+"/summary.json").c_str());summary<<"{\"passed\":true,\"checks\":"<<checks<<",\"qaMaps\":"<<maps<<",\"pairCalls\":"<<pairCalls<<",\"peakRssKiB\":"<<usage.ru_maxrss<<"}\n";summary.close();Need(bool(summary),"summary write");
    std::cout<<"NUCLEAR_PAIR_PID_QA_PASS checks="<<checks<<" maps="<<maps<<" pairCalls="<<pairCalls<<" peakRssKiB="<<usage.ru_maxrss<<std::endl;return 0;
  }catch(const std::exception& e){std::cerr<<"NUCLEAR_PAIR_PID_QA_FAIL "<<e.what()<<std::endl;return 1;}
}
#endif
