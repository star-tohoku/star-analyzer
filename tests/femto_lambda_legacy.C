// ACLiC test: ROOT5 + libStarAnaConfig/libStCommon. No KF library required.
// rootcint only needs the callable entry point. The compiled fixture deliberately
// uses C++11 calibration helpers that ROOT5's interpreter cannot parse.
int femto_lambda_legacy(const char* mainconf);
#ifndef __CINT__
#include "StMaker/common/FemtoLambdaLegacy.h"
#include "FemtoCandidate.h"
#include "FemtoLambdaProvider.h"
#include "StMaker/common/NuclearIdDeDxVsMom.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "TClonesArray.h"
#include "HistManager.h"
#include "TH1.h"
#include "TH2.h"
#include "TMemFile.h"
#include "TSystem.h"
#include <cmath>
#include <iostream>
#include <limits>
#include <algorithm>
#include <vector>
#include <stdexcept>

namespace {
void CheckLegacy(bool condition,const char* text) { if(!condition) throw std::runtime_error(text); }
double ToyKstar(const FemtoCandidate& a,const FemtoCandidate& b) {
  TLorentzVector q=a.P4(),pair=q+b.P4(); q.Boost(-pair.BoostVector()); return q.P();
}
int LegacyCells(const TH1& h) {
  return (h.GetNbinsX()+2)*(h.GetDimension()>1?h.GetNbinsY()+2:1);
}
void CheckPartition(const TH1* all,const TH1* low,const TH1* high) {
  CheckLegacy(all && low && high,"split histograms exist");
  CheckLegacy(LegacyCells(*all)==LegacyCells(*low) && LegacyCells(*all)==LegacyCells(*high),
              "split cell shape");
  CheckLegacy(all->GetEntries()==low->GetEntries()+high->GetEntries(),"split entry partition");
  for(int bin=0;bin<LegacyCells(*all);++bin) {
    CheckLegacy(all->GetBinContent(bin)==low->GetBinContent(bin)+high->GetBinContent(bin),
                "split content partition including under/overflow");
    const double va=std::pow(all->GetBinError(bin),2);
    const double vh=std::pow(low->GetBinError(bin),2)+std::pow(high->GetBinError(bin),2);
    CheckLegacy(std::fabs(va-vh)<1.e-10*std::max(1.,va),"split variance partition");
  }
}
void CheckSplitSignal(const char* mainconf,const FemtoCandidate& lambda,const FemtoCandidate& nucleus) {
  FemtoLambdaLegacy helper;
  CheckLegacy(helper.Init(mainconf,std::cerr),"split helper Init");
  const double mean=1.11596,width=.00296*3.;
  CheckLegacy(helper.SignalHalf(mean)==2,"exact mean belongs to high half only");
  CheckLegacy(helper.SignalHalf(std::numeric_limits<double>::quiet_NaN())==0 &&
              helper.SignalHalf(std::numeric_limits<double>::infinity())==0 &&
              helper.SignalHalf(-std::numeric_limits<double>::infinity())==0,"nonfinite mass has no half");
  const double edges[]={mean-width,mean,mean+width};
  std::vector<float> masses;
  for(int i=0;i<3;++i) {
    const float stored=static_cast<float>(edges[i]);
    masses.push_back(std::nextafter(stored,-std::numeric_limits<float>::infinity()));
    masses.push_back(stored);
    masses.push_back(std::nextafter(stored,std::numeric_limits<float>::infinity()));
    for(int j=-1;j<=1;++j) {
      const double mass=edges[i]+j*1.e-8;
      const int expected=helper.MassRegion(mass)==1?(mass<mean?1:2):0;
      CheckLegacy(helper.SignalHalf(mass)==expected,"double boundary follows original signal exactly");
    }
  }
  masses.push_back(static_cast<float>(mean-width/2));
  masses.push_back(static_cast<float>(mean+width/2));
  masses.push_back(1.095f);masses.push_back(1.14f);masses.push_back(1.20f);
  const std::string channel="lambda_"+helper.SpeciesKey();
  TH1* inclusive=helper.RootHistograms()->Get(("hKstarSE_"+channel+"_signal").c_str());
  const double lo=inclusive->GetXaxis()->GetXmin(),hi=inclusive->GetXaxis()->GetXmax();
  const double kstars[]={lo-.1,lo,lo+.25*(hi-lo),hi-1.e-8,hi,hi+.1};
  int expected[3]={0,0,0};
  for(size_t i=0;i<masses.size();++i) {
    FemtoCandidate l=lambda,n=nucleus;l.reso.invMass=masses[i];
    const int half=helper.MassRegion(l.reso.invMass)==1?(double(l.reso.invMass)<mean?1:2):0;
    CheckLegacy(helper.SignalHalf(l.reso.invMass)==half,"stored Float_t boundary classification");
    for(int cent=0;cent<9;++cent) for(unsigned k=0;k<sizeof(kstars)/sizeof(kstars[0]);++k) {
      if(half)++expected[half];
      n.eventIndex=l.eventIndex;n.trk.trackIndex=9;
      CheckLegacy(helper.FillPair(l,n,kstars[k],.2,cent,false),"split SE accepted");
      n.eventIndex=l.eventIndex+1;n.trk.trackIndex=l.reso.dau1Index;
      CheckLegacy(helper.FillPair(l,n,kstars[k],.2,cent,true),"split ME accepted with cross-event local index");
    }
  }
  CheckLegacy(expected[1]>0 && expected[2]>0,"both halves exercised");
  for(int mode=0;mode<2;++mode)for(int cent=0;cent<2;++cent) {
    const std::string stem=std::string("hKstar")+(mode?"ME":"SE")+(cent?"VsCent":"")+"_"+channel;
    TH1* signal=helper.RootHistograms()->Get((stem+"_signal").c_str());
    TH1* low=helper.RootHistograms()->Get((stem+"_signalLow").c_str());
    TH1* high=helper.RootHistograms()->Get((stem+"_signalHigh").c_str());
    CheckPartition(signal,low,high);
    CheckLegacy(low->GetEntries()==expected[1] && high->GetEntries()==expected[2],"split exact expected entries");
  }
  CheckLegacy(helper.RootHistograms()->Get("hLambda_PtVsYLab_signal")->GetEntries()==0,
              "pair filling never weights single-particle acceptance");
}
void CheckPairMt(const char* mainconf,const FemtoCandidate& lambda,const FemtoCandidate& nucleus) {
  FemtoLambdaLegacy helper;
  CheckLegacy(helper.Init(mainconf,std::cerr),"mT helper Init");
  const int type=helper.SpeciesKey()=="deuteron"?0:helper.SpeciesKey()=="triton"?1:helper.SpeciesKey()=="he3"?2:3;
  const char* shortName[]={"d","t","3He","4He"};
  const char* suffix[]={"","_SBPos","_SBNeg","_signalLow","_signalHigh"};
  const char* canonical[]={"_signal","_rightSB","_leftSB","_signalLow","_signalHigh"};
  const double lambdaMass=1.115683; // Explicit fixture matches production maker YAML.
  const double nuclearMass=helper.NuclearP4(TVector3(),type).M();
  const double averageMass=.5*(lambdaMass+nuclearMass);
  const float masses[]={1.112f,1.119f,1.095f,1.14f,1.20f};
  const double kstars[]={-.1,0.,.217,.999999,1.,1.1};
  // Non-collinear transverse vectors, then a large pair-pT exercising mT overflow.
  const double vectors[3][4]={{.7,-.3,-.2,.8},{3.,4.,-.9,.3},{20.,15.,6.,-3.}};
  for(int mode=0;mode<2;++mode)for(int cent=0;cent<9;++cent) {
    TH2* maps[5]={0,0,0,0,0};TH2* expected[5]={0,0,0,0,0};
    const std::string stem=std::string(mode?"Mixed_":"")+shortName[type];
    const std::string tail="_CentBin"+std::to_string(cent);
    for(int region=0;region<5;++region) {
      const std::string name="hKstarMt_"+stem+suffix[region]+tail;
      maps[region]=dynamic_cast<TH2*>(helper.NuclearHistograms()->Get(name.c_str()));
      CheckLegacy(maps[region] && std::string(maps[region]->ClassName())=="TH2D","all 90 mT keys are TH2D");
      CheckLegacy(maps[region]->GetNbinsX()==200 && maps[region]->GetNbinsY()==200 &&
                  maps[region]->GetXaxis()->GetXmin()==0. && maps[region]->GetXaxis()->GetXmax()==1. &&
                  maps[region]->GetYaxis()->GetXmin()==0. && maps[region]->GetYaxis()->GetXmax()==10.,"mT axes");
      expected[region]=static_cast<TH2*>(maps[region]->Clone(("expected_"+name).c_str()));
      expected[region]->SetDirectory(0);expected[region]->Reset();
    }
    for(int vector=0;vector<3;++vector)for(int longitudinal=0;longitudinal<2;++longitudinal) {
      FemtoCandidate l=lambda,n=nucleus;
      TLorentzVector lp,np;
      lp.SetXYZM(vectors[vector][0],vectors[vector][1],longitudinal?17.:2.1,lambdaMass);
      np.SetXYZM(vectors[vector][2],vectors[vector][3],longitudinal?-23.:-.9,nuclearMass);
      l.SetP4(lp);n.SetP4(np);
      n.eventIndex=l.eventIndex+(mode?1:0);n.trk.trackIndex=mode?l.reso.dau1Index:9;
      // Independent oracle: vector-sum pair pT / 2 and configured mean mass.
      // Read back Float_t components exactly as production sees them.
      const double px=double(l.px)+double(n.px),py=double(l.py)+double(n.py);
      const double mt=std::sqrt(averageMass*averageMass+.25*(px*px+py*py));
      if(vector==0) {
        const double wrongScalar=.5*(l.P4().Pt()+n.P4().Pt());
        const double wrongMt=std::sqrt(averageMass*averageMass+wrongScalar*wrongScalar);
        CheckLegacy(maps[0]->GetYaxis()->FindBin(mt)!=maps[0]->GetYaxis()->FindBin(wrongMt),
                    "fixture distinguishes vector sum from scalar pT sum");
        CheckLegacy(maps[0]->GetYaxis()->FindBin(mt)!=maps[0]->GetYaxis()->FindBin((l.P4()+n.P4()).Mt()),
                    "fixture distinguishes average-mass mT from pair total transverse mass");
      } else if(vector==2) {
        CheckLegacy(maps[0]->GetYaxis()->FindBin(mt)==maps[0]->GetNbinsY()+1,"fixture reaches mT overflow");
      }
      for(unsigned mass=0;mass<sizeof(masses)/sizeof(masses[0]);++mass) {
        l.reso.invMass=masses[mass];
        const int classification=helper.MassRegion(l.reso.invMass);
        for(unsigned k=0;k<sizeof(kstars)/sizeof(kstars[0]);++k) {
          CheckLegacy(helper.FillPair(l,n,kstars[k],.2,cent,mode!=0),"mT SE/ME accepted pair");
          if(classification==1) {
            expected[0]->Fill(kstars[k],mt);
            expected[l.reso.invMass<1.11596?3:4]->Fill(kstars[k],mt);
          } else if(classification==2) expected[2]->Fill(kstars[k],mt);
          else if(classification==3) expected[1]->Fill(kstars[k],mt);
        }
      }
    }
    for(int region=0;region<5;++region) {
      CheckLegacy(maps[region]->GetEntries()==expected[region]->GetEntries(),"mT exact entry count");
      for(int bin=0;bin<LegacyCells(*maps[region]);++bin) {
        CheckLegacy(maps[region]->GetBinContent(bin)==expected[region]->GetBinContent(bin),"mT exact expected bin including flows");
        CheckLegacy(maps[region]->GetBinError(bin)==expected[region]->GetBinError(bin),"mT exact expected Sumw2");
      }
      TH2* byCent=dynamic_cast<TH2*>(helper.RootHistograms()->Get(
        (std::string("hKstar")+(mode?"ME":"SE")+"VsCent_lambda_"+helper.SpeciesKey()+canonical[region]).c_str()));
      CheckLegacy(byCent!=0,"mT canonical centrality reference");
      TH1D* projection=maps[region]->ProjectionX("toy_mt_projection",0,maps[region]->GetNbinsY()+1,"e");
      projection->SetDirectory(0);
      for(int x=0;x<=maps[region]->GetNbinsX()+1;++x) {
        const int reference=byCent->GetBin(x,cent+1);
        CheckLegacy(projection->GetBinContent(x)==byCent->GetBinContent(reference),"mT ProjectionX reproduces old kstar cent slice");
        CheckLegacy(std::fabs(std::pow(projection->GetBinError(x),2)-std::pow(byCent->GetBinError(reference),2))<1.e-10,
                    "mT ProjectionX squared-error closure");
      }
      delete projection;delete expected[region];
    }
    CheckPartition(maps[0],maps[3],maps[4]);
  }
  TMemFile file("femto_lambda_mt_test","RECREATE");
  CheckLegacy(helper.Write(&file),"mT output write");
  for(int mode=0;mode<2;++mode)for(int cent=0;cent<9;++cent)for(int region=0;region<5;++region) {
    const std::string name=std::string("hKstarMt_")+(mode?"Mixed_":"")+shortName[type]+suffix[region]+
                           "_CentBin"+std::to_string(cent);
    CheckLegacy(file.Get((std::string(mode?"mix/":"true/")+name).c_str())!=0,"mT correct SE/ME directory");
    CheckLegacy(file.Get((std::string(mode?"true/":"mix/")+name).c_str())==0,"mT no opposite-directory artifact");
  }
  CheckLegacy(helper.RootHistograms()->Get("hLambda_PtVsYLab_signal")->GetEntries()==0,
              "mT pair maps never multiply candidate acceptance");
}

struct LambdaPicoFixture {
  TClonesArray event,tracks;
  TClonesArray* arrays[StPicoArrays::NAllPicoArrays];
  StPicoDst dst;
  LambdaPicoFixture():event("StPicoEvent"),tracks("StPicoTrack") {
    for(int i=0;i<StPicoArrays::NAllPicoArrays;++i)arrays[i]=0;
    arrays[StPicoArrays::Event]=&event;arrays[StPicoArrays::Track]=&tracks;
    new(event[0]) StPicoEvent();StPicoDst::set(arrays);
  }
  ~LambdaPicoFixture() { StPicoDst::unset(); }
};
void CheckAcceptance(const char* mainconf) {
  FemtoLambdaLegacy helper;
  CheckLegacy(helper.Init(mainconf,std::cerr),"acceptance helper Init");
  LambdaPicoFixture fixture;
  StPicoTrack* proton=new(fixture.tracks[0])StPicoTrack();
  StPicoTrack* pion=new(fixture.tracks[1])StPicoTrack();
  proton->setId(101);proton->setNHitsFit(25);proton->setGlobalMomentum(.3,.4,1.2);
  pion->setId(202);pion->setNHitsFit(-25);pion->setGlobalMomentum(-.2,.15,-.35);
  // Both daughters intentionally have zero pMom: primary momenta cannot supply this QA.
  CheckLegacy(proton->pMom().Mag2()==0 && pion->pMom().Mag2()==0,"secondary-track fixture");
  FemtoLambdaCandidate c;c.pdg=3122;c.protonIndex=0;c.pionIndex=1;c.protonId=101;c.pionId=202;
  TLorentzVector lambda;lambda.SetXYZM(.4,.2,2.,1.115683);c.candidate.SetP4(lambda);
  c.candidate.speciesKey="lambda";c.candidate.eventIndex=0;
  c.candidate.reso.dau1EventIndex=0;c.candidate.reso.dau2EventIndex=0;
  c.candidate.reso.dau1Index=0;c.candidate.reso.dau2Index=1;
  const float masses[]={1.11596f-.00296f,1.11596f,1.11596f+.00296f,1.095f,1.14f};
  for(unsigned i=0;i<sizeof(masses)/sizeof(masses[0]);++i) {
    c.mass=masses[i];c.candidate.reso.invMass=c.mass;
    helper.FillLambda(c,2,100.,TVector3());
    CheckLegacy(helper.FillLambdaAcceptance(c,&fixture.dst),"candidate acceptance fill");
  }
  CheckLegacy(helper.RootHistograms()->Get("hLambda_InvMass")->GetEntries()==5,
              "old full-mass Lambda QA remains unrestricted");
  const char* names[]={"hLambda_PtVsYLab_signal","hLambdaProton_PtVsYLab_signal","hLambdaPion_PtVsYLab_signal"};
  TLorentzVector p4[3];p4[0]=c.candidate.P4();
  p4[1].SetVectM(proton->gMom(),.9382720813);p4[2].SetVectM(pion->gMom(),.13957039);
  for(int i=0;i<3;++i) {
    TH2* h=dynamic_cast<TH2*>(helper.RootHistograms()->Get(names[i]));
    CheckLegacy(h && h->GetEntries()==3,"one entry per signal candidate, no sideband");
    CheckLegacy(h->GetBinContent(h->FindBin(p4[i].Rapidity(),p4[i].Pt()))==3,
                "acceptance X=lab rapidity Y=pT with correct mass/global momentum");
  }
  c.mass=1.11596f;
  CheckLegacy(!helper.FillLambdaAcceptance(c,0),"signal acceptance fails missing Pico source");
  c.protonIndex=100;
  CheckLegacy(!helper.FillLambdaAcceptance(c,&fixture.dst),"signal acceptance fails invalid daughter index");
  c.protonIndex=0;
  proton->setGlobalMomentum(std::numeric_limits<double>::quiet_NaN(),.4,1.2);
  CheckLegacy(!helper.FillLambdaAcceptance(c,&fixture.dst),"signal acceptance fails nonfinite daughter");
  proton->setGlobalMomentum(.3,.4,1.2);
  c.mass=1.095f;
  CheckLegacy(helper.FillLambdaAcceptance(c,0),"non-signal needs no acceptance source");
  for(int i=0;i<3;++i)
    CheckLegacy(helper.RootHistograms()->Get(names[i])->GetEntries()==3,"invalid source causes no partial acceptance fill");

  // Four selected hypotheses, no Lambda pairing: target-nucleus acceptance is track based.
  for(int s=0;s<4;++s) {
    StPicoTrack* n=new(fixture.tracks[s+2])StPicoTrack();
    n->setId(300+s);n->setNHitsFit(s%2?-25:25);n->setNHitsDedx(20);
    n->setPrimaryMomentum(.6,0.,.8);n->setGlobalMomentum(.6,0.,.8);
    // StPicoTrack setter takes GeV/cm; stored/getter calibration uses keV/cm.
    n->setDedx(1.e-6*NuclearIdDeDxVsMom::GetMean(static_cast<NuclearIdDeDxVsMom::SpeciesIndex>(s),n->pMom().Mag()));
  }
  FemtoCandidateStore store;
  helper.CollectNuclei(&fixture.dst,0,2,store);
  const char* species[]={"deuteron","triton","he3","he4"};
  for(int s=0;s<4;++s)CheckLegacy(store[species[s]].size()==1,"each final nuclear species selected once");
  TH2* nucleus=dynamic_cast<TH2*>(helper.RootHistograms()->Get(("hNucleus_PtVsYLab_"+helper.SpeciesKey()).c_str()));
  const FemtoCandidate& selected=store[helper.SpeciesKey()][0];
  CheckLegacy(nucleus && nucleus->GetEntries()==1,"target nucleus filled once without Lambda pairs");
  CheckLegacy(nucleus->GetBinContent(nucleus->FindBin(selected.P4().Rapidity(),selected.P4().Pt()))==1,
              "selected nuclear P4 includes Z correction exactly once");
  CheckLegacy(helper.RootHistograms()->Get(("hKstarSE_lambda_"+helper.SpeciesKey()).c_str())->GetEntries()==0,
              "acceptance does not create pairs");
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
  CheckSplitSignal(mainconf,l,n);
  CheckAcceptance(mainconf);
  CheckPairMt(mainconf,l,n);
  std::cout<<"PASS femto_lambda_legacy "<<helper.SpeciesKey()<<std::endl;
  return 0;
 } catch(const std::exception& e) { std::cerr<<"FAIL femto_lambda_legacy "<<e.what()<<std::endl;return 1; }
}
#endif // !__CINT__
