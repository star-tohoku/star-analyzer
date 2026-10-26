int test_he4_acceptance_samples(const char*,const char*);
#ifndef __CINT__
#include "test_femto_nuclear_pid_domain.C"
#include "TH2.h"
#include <vector>
namespace {
const char* SampleVariables[]={"BetaVsEta","BetaVsYLab","PtVsYLab"};
const char* SampleStages[]={"lowR_noTof","lowR_tofMatch","lowR_tofM2","highR_tofM2"};
std::string SampleKey(int v,int s){return std::string("hNucleus_")+SampleVariables[v]+"_he4_"+SampleStages[s];}
struct SampleTrack{double r;int tofMode;double m2;int group;SampleTrack(double a,int b,double c,int d):r(a),tofMode(b),m2(c),group(d){}};
void SampleAdd(DomainPicoFixture& f,FemtoLambdaLegacy& helper,const SampleTrack& row){
 const int index=f.tracks.GetEntriesFast();StPicoTrack* t=new(f.tracks[index])StPicoTrack();t->setId(index+20);t->setNHitsFit(25);t->setNHitsDedx(20);t->setPrimaryMomentum(row.r,0.,0.);t->setGlobalMomentum(row.r,0.,0.);
 const double r=t->pMom().Mag(),a=helper.NuclearDedxPull(3,r,1.),b=helper.NuclearDedxPull(3,r,2.);t->setDedx((1.-a/(b-a))*1.e-6);
 if(row.tofMode){const int k=f.tofs.GetEntriesFast();StPicoBTofPidTraits* v=new(f.tofs[k])StPicoBTofPidTraits();v->setTrackIndex(index);v->setBTofMatchFlag(1);v->setBeta(row.tofMode==2?0.:1./std::sqrt(1.+row.m2/(r*r)));t->setBTofPidTraitsIndex(k);}
}
}
#endif
int test_he4_acceptance_samples(const char* source,const char* output){
#ifndef __CINT__
 try{DomainChecks=0;DomainRequire(source&&*source&&output&&*output,"explicit inputs");DomainRequire(gSystem->AccessPathName(output)&&gSystem->mkdir(output,kTRUE)==0,"new output");FemtoLambdaLegacy helper;DomainRequire(helper.Init(source,std::cerr),"actual mainconf Init");DomainPicoFixture f;
 std::vector<SampleTrack> rows;
 // group: -1 rejected, 0 low noTOF only, 1 low TOF outside m2, 2 low TOF+m2, 3 high TOF+m2.
 rows.push_back(SampleTrack(.60,0,0,-1));
 rows.push_back(SampleTrack(.90,0,0,0));rows.push_back(SampleTrack(.95,1,3.48096,2));rows.push_back(SampleTrack(1.0,1,7.76906,1));rows.push_back(SampleTrack(1.1,2,0,0));
 rows.push_back(SampleTrack(std::nextafter(2.f,0.f),0,0,0));
 rows.push_back(SampleTrack(2.0,0,0,-1));rows.push_back(SampleTrack(2.0,1,3.48096,3));
 rows.push_back(SampleTrack(std::nextafter(2.f,3.f),1,7.76906,-1));
 rows.push_back(SampleTrack(3.3,0,0,-1));rows.push_back(SampleTrack(3.3,1,3.48096,3));rows.push_back(SampleTrack(5.0,1,3.48096,-1));
 for(unsigned i=0;i<rows.size();++i)SampleAdd(f,helper,rows[i]);
 FemtoCandidateStore store;helper.CollectNuclei(&f.dst,0,8,store);DomainRequire(store["he4"].size()==7,"selection unchanged: five low and two high He4 candidates");
 for(int v=0;v<3;++v){
  TH2* inclusive=dynamic_cast<TH2*>(helper.RootHistograms()->Get((std::string("hNucleus_")+SampleVariables[v]+"_he4").c_str()));DomainRequire(inclusive&&inclusive->GetEntries()==7,"inclusive He4 acceptance remains complete");TH2* sample[4];TH2* expected[4];const double counts[]={5,2,1,2};
  for(int s=0;s<4;++s){sample[s]=dynamic_cast<TH2*>(helper.RootHistograms()->Get(SampleKey(v,s).c_str()));DomainRequire(sample[s]&&sample[s]->GetEntries()==counts[s],"expected independent sample entries");expected[s]=static_cast<TH2*>(sample[s]->Clone(TString::Format("expected_%d_%d",v,s)));expected[s]->Reset();expected[s]->SetDirectory(0);}
  for(unsigned i=0;i<rows.size();++i){if(rows[i].group<0)continue;const StPicoTrack* t=f.dst.track(i);TLorentzVector p4=helper.NuclearP4(t->pMom(),3);FemtoCandidate saved;saved.SetP4(p4);p4=saved.P4();const double x=v==0?p4.Eta():p4.Rapidity(),y=v==2?p4.Pt():p4.P()/p4.E();if(rows[i].group<3)expected[0]->Fill(x,y);if(rows[i].group==1||rows[i].group==2)expected[1]->Fill(x,y);if(rows[i].group==2)expected[2]->Fill(x,y);if(rows[i].group==3)expected[3]->Fill(x,y);}
  for(int bin=0;bin<(inclusive->GetNbinsX()+2)*(inclusive->GetNbinsY()+2);++bin){
   for(int s=0;s<4;++s)DomainRequire(sample[s]->GetBinContent(bin)==expected[s]->GetBinContent(bin),"each sample matches independent per-track expected bin, including flows");
   DomainRequire(sample[0]->GetBinContent(bin)+sample[3]->GetBinContent(bin)==inclusive->GetBinContent(bin),"binwise low noTOF + high TOF+m2 closes to inclusive");DomainRequire(sample[2]->GetBinContent(bin)<=sample[1]->GetBinContent(bin)&&sample[1]->GetBinContent(bin)<=sample[0]->GetBinContent(bin),"binwise low m2 subset of matched subset of no requirement");
  }
  for(int s=0;s<4;++s)delete expected[s];
 }
 const std::string path=std::string(output)+"/he4_acceptance_samples.root";{TFile file(path.c_str(),"CREATE");DomainRequire(!file.IsZombie()&&helper.Write(&file),"write ROOT");}{TFile file(path.c_str(),"READ");for(int v=0;v<3;++v)for(int s=0;s<4;++s)DomainRequire(dynamic_cast<TH2*>(file.Get(SampleKey(v,s).c_str()))!=0,"all twelve TH2 saved at top level");TNamed* meta=dynamic_cast<TNamed*>(file.Get("FemtoHe4AcceptanceSamples"));DomainRequire(meta&&std::string(meta->GetTitle()).find("thresholdPOverQ: 2")!=std::string::npos&&std::string(meta->GetTitle()).find("do not sum low-R samples")!=std::string::npos,"sample definitions and boundary recorded");}
 std::cout<<"PASS He4 acceptance samples: "<<DomainChecks<<" checks"<<std::endl;return 0;
 }catch(const std::exception& e){std::cerr<<"FAIL He4 acceptance samples after "<<DomainChecks<<" checks: "<<e.what()<<std::endl;return 1;}
#else
 return 1;
#endif
}
