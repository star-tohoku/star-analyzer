int femto_nuclear_pair_pid_pool(const char*,const char*,const char*);
#ifndef __CINT__
// Only the test changes access. Include dependencies first so the access macro
// affects StFemtoMaker's test boundary, not standard/ROOT dependency classes.
#include "StMaker.h"
#include "StPhiKKReconstruction.h"
#include "FemtoCandidate.h"
#include "FemtoNuclearPidPairQaCandidate.h"
#include "FemtoPhiMixSampler.h"
#include "cuts/FemtoConfig.h"
#include "cuts/MixingConfig.h"
#include "StarClassLibrary/StPhysicalHelixD.hh"
#include "TString.h"
#include <deque>
#include <map>
#include <string>
#include <vector>
#include <sstream>
#define private public
#include "StMaker/StFemtoMaker/StFemtoMaker.h"
#undef private
#include "StMaker/common/FemtoLambdaLegacy.h"
#include "StMaker/common/FemtoNuclearPairPidSparse.h"
#include "ConfigManager.h"
#include "HistManager.h"
#include "StChain.h"
#include "yaml-cpp/yaml.h"
#include "TH1.h"
#include "TAxis.h"
#include "TSystem.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <sys/resource.h>

namespace PairPidPoolTest {
unsigned long long checks=0;
int Cells(const TH1& h){int n=h.GetNbinsX()+2;if(h.GetDimension()>1)n*=h.GetNbinsY()+2;if(h.GetDimension()>2)n*=h.GetNbinsZ()+2;return n;}
void Need(bool ok,const std::string& why){++checks;if(!ok)throw std::runtime_error(why);}
typedef std::vector<int> Cell;
typedef std::map<Cell,double> Counts;
struct Event {int id,cent;FemtoCandidateStore physical;FemtoNuclearPidPairQaStore qa;};
double Ks(const FemtoCandidate& a,const FemtoCandidate& b){TLorentzVector x=a.P4();x.Boost(-(a.P4()+b.P4()).BoostVector());return x.P();}
FemtoCandidate Lambda(int event,double mass){FemtoCandidate c;c.speciesKey="lambda";c.source=kFemtoCandResonance;c.eventIndex=event;c.reso.dau1EventIndex=event;c.reso.dau2EventIndex=event;c.reso.dau1Index=2;c.reso.dau2Index=3;c.reso.invMass=mass;TLorentzVector p;p.SetXYZM(.3+.01*event,.17,.8,1.115683);c.SetP4(p);return c;}
FemtoCandidate Nucleus(FemtoLambdaLegacy& helper,int event,int type,int track,double ns){FemtoCandidate c;c.speciesKey=type==2?"he3":"he4";c.source=kFemtoCandTrack;c.eventIndex=event;c.trk.trackIndex=track;c.charge=2;c.SetP4(helper.NuclearP4(TVector3(.8+.02*event,.1,.3),type));c.trk.nSigmaHe3=type==2?ns:8.2;c.trk.nSigmaHe4=type==3?ns:-7.4;return c;}
void Qa(Event& e,const FemtoCandidate& n,bool selected){FemtoNuclearPidPairQaCandidate q;q.candidate=n;q.selected=selected;e.qa[n.speciesKey].push_back(q);}
bool Nonempty(const FemtoCandidateStore& s){for(FemtoCandidateStore::const_iterator i=s.begin();i!=s.end();++i)if(!i->second.empty())return true;return false;}
void OraclePair(Counts& expected,THnSparseD& axes,const FemtoCandidate& l,const FemtoCandidate& n,int cent,bool mixed,double mean,double width){
  if(FemtoCandidatesShareTrack(l,n)||(mixed&&l.eventIndex==n.eventIndex)||(!mixed&&l.eventIndex!=n.eventIndex))return;
  const double delta=double(l.reso.invMass)-mean;int mass=-1;if(std::fabs(delta)<=width)mass=delta<0?1:2;else if(delta< -width&&delta>=-4.*width)mass=0;else if(delta>width&&delta<=4.*width)mass=3;if(mass<0)return;
  const int repr=n.speciesKey=="he4"?0:1;
  for(int mapped=0;mapped<(repr==1?2:1);++mapped){FemtoCandidate representation=n;if(mapped){TLorentzVector p;p.SetVectM(n.P4().Vect(),3.72742);representation.SetP4(p);}
    const double point[]={Ks(l,representation),n.P4().P()/2.,repr==0?n.trk.nSigmaHe4:n.trk.nSigmaHe3,double(cent),double(mapped?2:repr),double(mixed?1:0),double(mass),1.};Cell cell(8);for(int d=0;d<8;++d)cell[d]=axes.GetAxis(d)->FindBin(point[d]);expected[cell]+=1.;}
}
void OracleQa(Counts& expected,THnSparseD& axes,const Event& current,const std::deque<Event>& buffer,int maxMix,double mean,double width){
  FemtoCandidateStore::const_iterator lc=current.physical.find("lambda");
  if(lc!=current.physical.end())for(unsigned a=0;a<lc->second.size();++a)for(FemtoNuclearPidPairQaStore::const_iterator target=current.qa.begin();target!=current.qa.end();++target)for(unsigned b=0;b<target->second.size();++b)OraclePair(expected,axes,lc->second[a],target->second[b].candidate,current.cent,false,mean,width);
  for(unsigned depth=0;depth<buffer.size()&&depth<unsigned(maxMix);++depth){const Event& old=buffer[buffer.size()-1-depth];
    if(lc!=current.physical.end())for(unsigned a=0;a<lc->second.size();++a)for(FemtoNuclearPidPairQaStore::const_iterator target=old.qa.begin();target!=old.qa.end();++target)for(unsigned b=0;b<target->second.size();++b)OraclePair(expected,axes,lc->second[a],target->second[b].candidate,current.cent,true,mean,width);
    FemtoCandidateStore::const_iterator lb=old.physical.find("lambda");if(lb!=old.physical.end())for(unsigned a=0;a<lb->second.size();++a)for(FemtoNuclearPidPairQaStore::const_iterator target=current.qa.begin();target!=current.qa.end();++target)for(unsigned b=0;b<target->second.size();++b)OraclePair(expected,axes,lb->second[a],target->second[b].candidate,current.cent,true,mean,width);
  }
}
void OldHistEqual(StFemtoMaker& on,StFemtoMaker& off,const YAML::Node& main,const std::string& base){
  const char* key[]={"femtoHist","nuclearHist"};for(int which=0;which<2;++which){std::string ref=main[key[which]].as<std::string>();if(ref[0]!='/')ref=base+ref;const YAML::Node hist=YAML::LoadFile(ref);for(YAML::const_iterator i=hist["histograms"].begin();i!=hist["histograms"].end();++i){const std::string name=i->first.as<std::string>();TH1* a=(which?on.mLambdaLegacy->NuclearHistograms():on.mLambdaLegacy->RootHistograms())->Get(name.c_str());TH1* b=(which?off.mLambdaLegacy->NuclearHistograms():off.mLambdaLegacy->RootHistograms())->Get(name.c_str());Need(a&&b&&Cells(*a)==Cells(*b)&&a->GetEntries()==b->GetEntries(),"regular histogram structure/entries invariant");for(int cell=0;cell<Cells(*a);++cell){Need(a->GetBinContent(cell)==b->GetBinContent(cell),"QA loops leave regular contents unchanged");Need(a->GetBinError(cell)==b->GetBinError(cell),"QA loops leave regular errors unchanged");}}}
}
}

int femto_nuclear_pair_pid_pool(const char* onMain,const char* offMain,const char* outputDir){
  using namespace PairPidPoolTest;
  try{
    checks=0;Need(onMain&&offMain&&outputDir,"explicit paths");Need(gSystem->AccessPathName(outputDir)&&gSystem->mkdir(outputDir,true)==0,"new output dir");
    Need(ConfigManager::GetInstance().LoadConfig(onMain),"actual production configuration loaded");
    const YAML::Node main=YAML::LoadFile(onMain);const std::string path=onMain,base=path.substr(0,path.rfind("/config/"))+"/config/";std::string makerRef=main["maker"].as<std::string>();if(makerRef[0]!='/')makerRef=base+makerRef;const YAML::Node maker=YAML::LoadFile(makerRef);const double mean=maker["lambdaSignalMean"].as<double>(),width=maker["lambdaSignalSigma"].as<double>()*maker["lambdaSignalNSigma"].as<double>();std::string mixRef=main["mixing"].as<std::string>();if(mixRef[0]!='/')mixRef=base+mixRef;const YAML::Node mixing=YAML::LoadFile(mixRef);const double vz=.5*(mixing["minVz"].as<double>()+mixing["maxVz"].as<double>());
    StChain* chain=new StChain("pairPidPoolRegression");StFemtoMaker* on=new StFemtoMaker("poolOn",0,"unused-on.root");StFemtoMaker* off=new StFemtoMaker("poolOff",0,"unused-off.root");StFemtoMaker* pair[]={on,off};
    for(int i=0;i<2;++i){pair[i]->mLambdaLegacy=new FemtoLambdaLegacy;Need(pair[i]->mLambdaLegacy->Init(i?offMain:onMain,std::cerr),"helper init");pair[i]->m_histManager=pair[i]->mLambdaLegacy->RootHistograms();pair[i]->mLambdaMixMinVz=mixing["minVz"].as<double>();pair[i]->mLambdaMixMaxVz=mixing["maxVz"].as<double>();pair[i]->mLambdaMaxMixEvents=2;pair[i]->m_cent9=4;}
    THnSparseD* joint=on->mLambdaLegacy->NuclearPidPairSparse()->Histogram();Need(joint!=0,"joint exists");Need(std::fabs(on->mLambdaLegacy->NuclearP4(TVector3(),3).M()-3.72742)<1.e-12,"independent mapped-mass oracle equals configured He4 mass");
    std::vector<Event> events(8);for(int i=0;i<8;++i){events[i].id=i;events[i].cent=i==7?5:4;}
    Qa(events[0],Nucleus(*on->mLambdaLegacy,0,3,9,6.1),false);
    events[1].physical["lambda"].push_back(Lambda(1,mean-.5*width));Qa(events[1],Nucleus(*on->mLambdaLegacy,1,3,9,-1.1),false);
    FemtoCandidate n2=Nucleus(*on->mLambdaLegacy,2,3,7,.5);events[2].physical["he4"].push_back(n2);Qa(events[2],n2,true);Qa(events[2],Nucleus(*on->mLambdaLegacy,2,2,8,-2.5),false);
    events[3].physical["lambda"]; // Nonempty map, but no physical candidates.
    Qa(events[3],Nucleus(*on->mLambdaLegacy,3,3,9,4.1),false);
    events[4].physical["lambda"].push_back(Lambda(4,mean+.5*width));FemtoCandidate n4=Nucleus(*on->mLambdaLegacy,4,2,9,-.7);events[4].physical["he3"].push_back(n4);Qa(events[4],n4,true);Qa(events[4],Nucleus(*on->mLambdaLegacy,4,3,8,5.2),false);
    events[5].physical["lambda"].push_back(Lambda(5,mean-2*width));Qa(events[5],Nucleus(*on->mLambdaLegacy,5,3,2,3.3),false);
    FemtoCandidate n6=Nucleus(*on->mLambdaLegacy,6,3,3,-.3);events[6].physical["he4"].push_back(n6);Qa(events[6],n6,true);
    events[7].physical["lambda"].push_back(Lambda(7,mean+2*width));Qa(events[7],Nucleus(*on->mLambdaLegacy,7,3,9,2.7),false);
    std::map<int,std::deque<Event> > oraclePool;Counts expected;unsigned refused=0;
    const FemtoConfig& config=ConfigManager::GetInstance().GetFemtoConfig();
    for(unsigned e=0;e<events.size();++e){const Event& current=events[e];const int bin=on->GetMixingBin(vz,current.cent,0.);OracleQa(expected,*joint,current,oraclePool[bin],2,mean,width);
      for(int state=0;state<2;++state){StFemtoMaker& m=*pair[state];m.m_eventCandidates=current.physical;m.mLambdaNuclearPidQaCandidates=state?FemtoNuclearPidPairQaStore():current.qa;m.m_cent9=current.cent;
        for(unsigned channel=0;channel<config.channels.size();++channel){const FemtoConfig::ChannelDef& ch=config.channels[channel];if(ch.name!="lambda_he4"&&ch.name!="lambda_he3")continue;m.FillSameEventPairs(ch);m.FillMixedEventPairs(ch,channel,vz,current.cent,0.);}
        m.FillNuclearPidQaPairs(vz,current.cent,0.);m.StoreEventForMixing(vz,current.cent,0.);
      }
      if(Nonempty(current.physical)){oraclePool[bin].push_back(current);while(oraclePool[bin].size()>unsigned(ConfigManager::GetInstance().GetMixingConfig().bufferSize))oraclePool[bin].pop_front();}else ++refused;
      Need(on->m_mixingPool.size()==off->m_mixingPool.size(),"QA-only event never creates a physical pool key");
      for(std::map<int,std::deque<StFemtoMaker::FemtoMixingEvent> >::const_iterator it=on->m_mixingPool.begin();it!=on->m_mixingPool.end();++it){Need(off->m_mixingPool[it->first].size()==it->second.size()&&it->second.size()==oraclePool[it->first].size(),"physical admission and eviction unchanged");}
    }
    Counts actual;for(Long64_t bin=0;bin<joint->GetNbins();++bin){Int_t coords[8];const double value=joint->GetBinContent(bin,coords);if(coords[7]!=2||!value)continue;Cell key(coords,coords+8);actual[key]+=value;Need(std::fabs(std::pow(joint->GetBinError(bin),2)-value)<1.e-9,"prePID unit-count sparse variance");}
    Need(actual==expected,"prePID SE/ME exact tuple oracle: source pulls/R, cent, mass, both directions, last-two physical events");Need(refused==2,"two QA-only events exercised");OldHistEqual(*on,*off,main,base);
    const double regularSE=on->mLambdaLegacy->RootHistograms()->Get("hKstarSE_lambda_he3_signal")->GetEntries(),regularME=on->mLambdaLegacy->RootHistograms()->Get("hKstarME_lambda_he4_signal")->GetEntries();Need(regularSE>0&&regularME>0,"explicitly nonempty regular SE/ME controls rule out vacuous invariance");
    struct rusage usage;getrusage(RUSAGE_SELF,&usage);std::ofstream summary((std::string(outputDir)+"/summary.json").c_str());summary<<"{\"passed\":true,\"checks\":"<<checks<<",\"events\":"<<events.size()<<",\"qaOnlyRefused\":"<<refused<<",\"prePidNonzeroCells\":"<<actual.size()<<",\"regularSEHe3Signal\":"<<regularSE<<",\"regularMEHe4Signal\":"<<regularME<<",\"peakRssKiB\":"<<usage.ru_maxrss<<"}\n";summary.close();Need(bool(summary),"summary saved");delete chain;
    std::cout<<"NUCLEAR_PAIR_PID_POOL_PASS checks="<<checks<<" refusedQaOnly="<<refused<<" nonzeroCells="<<actual.size()<<std::endl;return 0;
  }catch(const std::exception& e){std::cerr<<"NUCLEAR_PAIR_PID_POOL_FAIL "<<e.what()<<std::endl;return 1;}
}
#endif
