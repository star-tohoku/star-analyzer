// ACLiC test: ROOT5 + libStarAnaConfig/libStCommon. No KF library required.
#include "StMaker/common/FemtoLambdaLegacy.h"
#include "FemtoCandidate.h"
#include "HistManager.h"
#include "TH1.h"
#include "TH2.h"
#include "TMemFile.h"
#include "TSystem.h"
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {
void CheckLegacy(bool condition,const char* text) { if(!condition) throw std::runtime_error(text); }
double ToyKstar(const FemtoCandidate& a,const FemtoCandidate& b) {
  TLorentzVector q=a.P4(),pair=q+b.P4(); q.Boost(-pair.BoostVector()); return q.P();
}
}
int femto_lambda_legacy(const char* mainconf) {
 try {
  FemtoLambdaLegacy helper;
  CheckLegacy(helper.Init(mainconf,std::cerr),"helper Init");
  CheckLegacy(helper.MassRegion(1.11596)==1,"signal");
  CheckLegacy(helper.MassRegion(1.095)==2,"left sideband");
  CheckLegacy(helper.MassRegion(1.14)==3,"right sideband");
  CheckLegacy(helper.MassRegion(1.20)==0,"outside narrow windows");
  const double mean=1.11596,width=.00296*3.;
  const double edges[]={mean-width,mean+width,mean-4*width,mean+4*width};
  for(int i=0;i<4;++i) for(int j=-1;j<=1;++j) {
    const double mass=edges[i]+j*1.e-8,delta=mass-mean;
    const int expected=std::fabs(delta)<=width?1:(delta < -width && delta>=-4*width?2:(delta>width && delta<=4*width?3:0));
    CheckLegacy(helper.MassRegion(mass)==expected,"legacy boundary inequality");
  }
  double pulls[4]={.1,.1,999.,999.};
  CheckLegacy(helper.SelectNuclearType(1.,pulls,false,-999.)==0,"tie chooses d before t");
  CheckLegacy(helper.SelectNuclearType(2.5,pulls,false,-999.)==-1,"rigidity upper boundary excluded");
  CheckLegacy(helper.SelectNuclearType(2.4999,pulls,false,-999.)==0,"rigidity below boundary accepted");
  const TVector3 raw(.3,.4,.5);
  const double masses[]={1.87561,2.80892,2.80839,3.72742};
  for(int s=0;s<4;++s) {
    double p[4]={999.,999.,999.,999.}; p[s]=0.;
    CheckLegacy(helper.SelectNuclearType(1.,p,false,-999.)==s,"each nuclear hypothesis including rare He4");
    const TLorentzVector q=helper.NuclearP4(raw,s);
    CheckLegacy(std::fabs(q.P()-raw.Mag()*(s>=2?2.:1.))<1.e-12,"Z correction once");
    CheckLegacy(std::fabs(q.M()-masses[s])<1.e-12,"nuclear fixed mass");
  }
  FemtoCandidate l,n;
  l.eventIndex=10; l.source=kFemtoCandResonance; l.speciesKey="lambda";
  l.reso.dau1EventIndex=10;l.reso.dau1Index=7;l.reso.dau2EventIndex=10;l.reso.dau2Index=8;
  TLorentzVector lp;lp.SetVectM(TVector3(.2,.1,1.),1.115683);l.SetP4(lp);
  n.eventIndex=10;n.source=kFemtoCandTrack;n.speciesKey=helper.SpeciesKey();n.trk.trackIndex=9;
  n.SetP4(helper.NuclearP4(raw,3));
  const double k=ToyKstar(l,n),ql=(l.P4().Vect()-n.P4().Vect()).Mag();
  const double mass[]={1.11596,1.095,1.14,1.20};
  for(int i=0;i<4;++i) {
    l.reso.invMass=mass[i];
    CheckLegacy(std::fabs(ToyKstar(l,n)-k)<1.e-12,"raw mass must not change kstar");
    CheckLegacy(helper.FillPair(l,n,k,ql,2,false),"SE pair");
    n.eventIndex=11;n.trk.trackIndex=7; // Reused local daughter index in a DIFFERENT event.
    CheckLegacy(helper.FillPair(l,n,k,ql,2,true),"ME reused local index");
    n.eventIndex=10;n.trk.trackIndex=9;
  }
  n.trk.trackIndex=7;
  CheckLegacy(!helper.FillPair(l,n,k,ql,2,false),"same-event daughter overlap rejected");
  n.trk.trackIndex=9;
  CheckLegacy(!helper.FillPair(l,n,k,ql,2,true),"same event forbidden in ME");
  n.eventIndex=12;
  CheckLegacy(!helper.FillPair(l,n,k,ql,2,false),"cross event forbidden in SE");
  const std::string channel="lambda_"+helper.SpeciesKey();
  for(int m=0;m<2;++m) {
    const std::string mode=m?"ME":"SE";
    CheckLegacy(helper.RootHistograms()->Get(("hKstar"+mode+"_"+channel).c_str())->GetEntries()==4,"base exactly one fill per pair");
    const char* suffix[]={"_signal","_leftSB","_rightSB"};
    for(int r=0;r<3;++r) CheckLegacy(helper.RootHistograms()->Get(("hKstar"+mode+"_"+channel+suffix[r]).c_str())->GetEntries()==1,"mass-region fanout once");
  }
  TMemFile file("femto_lambda_legacy_test","RECREATE");
  CheckLegacy(helper.Write(&file),"directory write");
  const char* legacy=helper.SpeciesKey()=="deuteron"?"d":helper.SpeciesKey()=="triton"?"t":helper.SpeciesKey()=="he3"?"3He":"4He";
  const std::string sp=legacy;
  TH2* se=dynamic_cast<TH2*>(file.Get(("true/hKstarMass_"+sp+"_CentBin2").c_str()));
  TH2* me=dynamic_cast<TH2*>(file.Get(("mix/hKstarMass_Mixed_"+sp+"_CentBin2").c_str()));
  CheckLegacy(se && me && se->GetEntries()==4 && me->GetEntries()==4,"both full-mass directories once");
  CheckLegacy(!file.Get(("true/hKstarMass_Mixed_"+sp+"_CentBin2").c_str()),"no empty ME artifact in true/");
  CheckLegacy(se->GetNbinsX()==200 && se->GetNbinsY()==200 &&
              std::fabs(se->GetXaxis()->GetXmin())<1.e-12 &&
              std::fabs(se->GetXaxis()->GetXmax()-1.)<1.e-12 &&
              std::fabs(se->GetYaxis()->GetXmin()-1.05)<1.e-12 &&
              std::fabs(se->GetYaxis()->GetXmax()-1.25)<1.e-12,"legacy full-mass axes");
  std::cout<<"PASS femto_lambda_legacy "<<helper.SpeciesKey()<<std::endl;
  return 0;
 } catch(const std::exception& e) { std::cerr<<"FAIL femto_lambda_legacy "<<e.what()<<std::endl;return 1; }
}

