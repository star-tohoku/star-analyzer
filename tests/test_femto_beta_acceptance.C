int test_femto_beta_acceptance(const char*,const char*,const char*,const char*);
#ifndef __CINT__
#include "test_femto_nuclear_pid_domain.C"
#include "FemtoLambdaProvider.h"
#include "TH2.h"
namespace {
const char* BetaNames[]={"hLambda_BetaVsEta_signal","hLambda_BetaVsYLab_signal","hNucleus_BetaVsEta_he4","hNucleus_BetaVsYLab_he4"};
void CheckBetaBin(FemtoLambdaLegacy& h,int species,const TLorentzVector& v,double entries) {
  for(int axis=0;axis<2;++axis){TH2* map=dynamic_cast<TH2*>(h.RootHistograms()->Get(BetaNames[2*species+axis]));DomainRequire(map&&map->GetEntries()==entries,"one entry per accepted candidate");const double x=axis==0?v.Eta():v.Rapidity();DomainRequire(map->GetBinContent(map->FindBin(x,v.P()/v.E()))>=1,"X=eta or lab Y, vertical beta=P/E from candidate physical P4");}
}
}
#endif
int test_femto_beta_acceptance(const char* source,const char* unused1,const char* unused2,const char* output){
#ifndef __CINT__
 try{DomainChecks=0;DomainRequire(source&&*source&&output&&*output,"explicit test arguments");DomainRequire(gSystem->AccessPathName(output)&&gSystem->mkdir(output,kTRUE)==0,"new output directory");
 FemtoLambdaLegacy helper;DomainRequire(helper.Init(source,std::cerr),"production Init with beta profile");
 for(int i=0;i<4;++i){TH2* h=dynamic_cast<TH2*>(helper.RootHistograms()->Get(BetaNames[i]));DomainRequire(h&&h->GetDimension()==2,"four beta maps are TH2");DomainRequire(h->GetYaxis()->GetXmin()==0&&h->GetYaxis()->GetXmax()==1&&h->GetNbinsY()==200,"beta axis from YAML");DomainRequire(h->GetXaxis()->GetXmin()==-5&&h->GetXaxis()->GetXmax()==5,"wide lab eta/Y axes");}
 DomainPicoFixture fixture;
 fixture.Add(helper,.60,false);fixture.Add(helper,.90,false);fixture.Add(helper,3.3,false);fixture.Add(helper,3.3,true);
 fixture.dst.btofPidTraits(0)->setBeta(1./std::sqrt(1.+3.48096/(3.3*3.3)));
 FemtoCandidateStore store;helper.CollectNuclei(&fixture.dst,0,8,store);
 DomainRequire(store["he4"].size()==2,"accepted low-noTOF and high-TOF He4 only");
 for(unsigned i=0;i<2;++i)CheckBetaBin(helper,1,store["he4"][i].P4(),2);
 // The high-R candidate has physical P=2R, so raw-rigidity beta would land elsewhere.
 const TLorentzVector nucleus=store["he4"][1].P4();const double wrongBeta=(nucleus.P()/2.)/std::sqrt(std::pow(nucleus.P()/2.,2)+nucleus.M2());TH2* hn=dynamic_cast<TH2*>(helper.RootHistograms()->Get(BetaNames[2]));DomainRequire(hn->GetBinContent(hn->FindBin(nucleus.Eta(),wrongBeta))==0,"He4 beta uses physical momentum, not rigidity");
 const int protonIndex=fixture.tracks.GetEntriesFast();StPicoTrack* p=new(fixture.tracks[protonIndex])StPicoTrack();p->setGlobalMomentum(.3,.4,1.2);
 const int pionIndex=fixture.tracks.GetEntriesFast();StPicoTrack* pi=new(fixture.tracks[pionIndex])StPicoTrack();pi->setGlobalMomentum(-.2,.15,-.35);
 FemtoLambdaCandidate c;c.protonIndex=protonIndex;c.pionIndex=pionIndex;c.candidate.source=kFemtoCandResonance;c.candidate.eventIndex=0;c.candidate.speciesKey="lambda";c.candidate.reso.dau1EventIndex=c.candidate.reso.dau2EventIndex=0;c.candidate.reso.dau1Index=protonIndex;c.candidate.reso.dau2Index=pionIndex;
 for(int i=0;i<2;++i){TLorentzVector v;v.SetXYZM(.4,.2,i==0?2.:-2.,1.115683);c.candidate.SetP4(v);c.mass=1.11596;c.candidate.reso.invMass=c.mass;DomainRequire(helper.FillLambdaAcceptance(c,&fixture.dst),"signal Lambda acceptance fills");CheckBetaBin(helper,0,c.candidate.P4(),i+1);}
 c.mass=1.095;c.candidate.reso.invMass=c.mass;DomainRequire(helper.FillLambdaAcceptance(c,&fixture.dst),"sideband acceptance call succeeds without filling");
 for(int k=0;k<4;++k)DomainRequire(helper.FillPair(c.candidate,store["he4"][0],.1,.2,8,false,&fixture.dst),"pair filling remains valid");
 for(int i=0;i<4;++i)DomainRequire(helper.RootHistograms()->Get(BetaNames[i])->GetEntries()==2,"sideband and multiple pairs do not multiply single-candidate beta acceptance");
 const std::string path=std::string(output)+"/beta_acceptance.root";{TFile f(path.c_str(),"CREATE");DomainRequire(!f.IsZombie()&&helper.Write(&f),"write beta acceptance ROOT");}
 {TFile f(path.c_str(),"READ");for(int i=0;i<4;++i){TH2* h=dynamic_cast<TH2*>(f.Get(BetaNames[i]));DomainRequire(h&&h->GetEntries()==2,"four TH2 maps saved at ROOT top level");}}
 std::cout<<"PASS beta acceptance: "<<DomainChecks<<" checks"<<std::endl;return 0;
 }catch(const std::exception& e){std::cerr<<"FAIL beta acceptance after "<<DomainChecks<<" checks: "<<e.what()<<std::endl;return 1;}
#else
 return 1;
#endif
}
