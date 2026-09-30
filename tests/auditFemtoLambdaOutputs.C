#include "TFile.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TKey.h"
#include "TNamed.h"
#include "TParameter.h"
#include "TAxis.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TSystem.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <sstream>
#include <string>
#include <vector>

// Read-only end-to-end audit. Original anaLambdaNuclearId is the legacy reference.
// Histograms and trees are never modified. Only the CSV report is written.
// Fully qualified vector element names are required by the ROOT5 CINT dictionary.
namespace FemtoLambdaAudit {
void Require(bool ok, const std::string& what) { if (!ok) throw std::runtime_error(what); }
std::string Number(int n) { std::ostringstream s; s << n; return s.str(); }
TH1* Hist(TFile& f, const std::string& key) {
  TH1* h = dynamic_cast<TH1*>(f.Get(key.c_str()));
  Require(h != 0, std::string(f.GetName()) + " missing " + key);
  return h;
}
int CellCount(const TH1& h) {
  return (h.GetNbinsX()+2)*(h.GetDimension()>1?h.GetNbinsY()+2:1)*(h.GetDimension()>2?h.GetNbinsZ()+2:1);
}
bool SameAxis(const TAxis& a, const TAxis& b) {
  if(a.GetNbins()!=b.GetNbins() || a.GetXmin()!=b.GetXmin() || a.GetXmax()!=b.GetXmax())return false;
  for(int i=1;i<=a.GetNbins()+1;++i)
    if(a.GetBinLowEdge(i)!=b.GetBinLowEdge(i))return false;
  return true;
}
bool SameShape(const TH1& a,const TH1& b) {
  if(std::string(a.ClassName())!=b.ClassName() || a.GetDimension()!=b.GetDimension() ||
     !SameAxis(*a.GetXaxis(),*b.GetXaxis()))return false;
  if(a.GetDimension()>1 && !SameAxis(*a.GetYaxis(),*b.GetYaxis()))return false;
  if(a.GetDimension()>2 && !SameAxis(*a.GetZaxis(),*b.GetZaxis()))return false;
  return true;
}
bool SameHist(const TH1& a, const TH1& b) {
  if(!SameShape(a,b) || a.GetEntries()!=b.GetEntries())return false;
  for(int i=0;i<CellCount(a);++i)
    if(a.GetBinContent(i)!=b.GetBinContent(i) || a.GetBinError(i)!=b.GetBinError(i))return false;
  return true;
}
void AxisContract(const TAxis& a,int bins,double lo,double hi,const std::string& what) {
  Require(a.GetNbins()==bins && std::fabs(a.GetXmin()-lo)<1.e-12 &&
      std::fabs(a.GetXmax()-hi)<1.e-12,what+" axis contract");
  for(int i=1;i<=bins+1;++i)
    Require(std::fabs(a.GetBinLowEdge(i)-(lo+(hi-lo)*(i-1)/bins))<1.e-12,what+" bin edges");
}
bool SameVariance(double a,double b) {
  // Only derived sums of squared errors need rounding tolerance. Contents,
  // original QA errors and direct aliases remain exact comparisons.
  return std::fabs(a-b)<=1.e-12*std::max(1.,std::max(std::fabs(a),std::fabs(b)));
}
double Variance(const TH1& h,int bin) { const double e=h.GetBinError(bin);return e*e; }
void Branch(TTree& t,const char* key,void* address,const char* type) {
  TLeaf* leaf=t.GetLeaf(key);
  Require(leaf && std::string(leaf->GetTypeName())==type,std::string(t.GetName())+" missing/wrong branch "+key);
  Require(t.SetBranchAddress(key,address)>=0,std::string(t.GetName())+" cannot bind "+key);
}
Long64_t Parameter(TFile& f,const char* name) {
  TParameter<Long64_t>* p=dynamic_cast<TParameter<Long64_t>*>(f.Get(name));
  Require(p!=0,std::string("missing count parameter ")+name);return p->GetVal();
}
struct Event {
  int input,run,event,status,cent,selected;
};
std::vector<FemtoLambdaAudit::Event> Events(TFile& f,bool requireSelected) {
  TTree* t=dynamic_cast<TTree*>(f.Get("eventLedger"));
  Require(t!=0,"event ledger missing");Event e;
  Branch(*t,"inputIndex",&e.input,"Int_t");Branch(*t,"runId",&e.run,"Int_t");
  Branch(*t,"eventId",&e.event,"Int_t");Branch(*t,"status",&e.status,"Int_t");
  Branch(*t,"cent9",&e.cent,"Int_t");
  const bool hasSelected=t->GetLeaf("selectedLambda")!=0;
  Require(!requireSelected || hasSelected,"event selectedLambda missing");
  if(hasSelected)Branch(*t,"selectedLambda",&e.selected,"Int_t");
  std::vector<FemtoLambdaAudit::Event> out;
  for(Long64_t i=0;i<t->GetEntries();++i){
    e.selected=-1;
    Require(t->GetEntry(i)>0,"event ledger read failed");
    Require(e.input==i+1,"event inputIndex not consecutive / one based");
    Require(e.status>=0 && e.status<=4,"event ledger contains incomplete/error status");
    Require(e.status!=4 || (e.cent>=0 && e.cent<9),"reconstructed event has invalid cent9");
    if(hasSelected)Require(e.selected>=0 && (e.status==4 || e.selected==0),"event selected count/status mismatch");
    out.push_back(e);
  }
  t->ResetBranchAddresses();return out;
}
void SameEvents(const std::vector<FemtoLambdaAudit::Event>& a,const std::vector<FemtoLambdaAudit::Event>& b,bool selected) {
  Require(a.size()==b.size(),"event ledger size mismatch");
  for(size_t i=0;i<a.size();++i)
    Require(a[i].input==b[i].input && a[i].run==b[i].run && a[i].event==b[i].event &&
        a[i].status==b[i].status && a[i].cent==b[i].cent && (!selected || a[i].selected==b[i].selected),
        "event ledger input/run/event/status/cent9 mismatch at input "+Number(i+1));
}
struct Candidate {
  int run,event,proton,pion,input;
  float mass,px,py,pz;
  bool operator<(const Candidate& o) const {
    if(run!=o.run)return run<o.run;if(event!=o.event)return event<o.event;
    if(proton!=o.proton)return proton<o.proton;if(pion!=o.pion)return pion<o.pion;
    if(mass!=o.mass)return mass<o.mass;if(px!=o.px)return px<o.px;
    if(py!=o.py)return py<o.py;return pz<o.pz;
  }
};
std::vector<FemtoLambdaAudit::Candidate> Candidates(TFile& f,bool standalone,const std::vector<FemtoLambdaAudit::Event>* events=0) {
  TTree* t=dynamic_cast<TTree*>(f.Get(standalone?"KfLambdaCandidates":"selectedLambda"));
  Require(t!=0,"candidate ledger missing");Candidate c;Bool_t selected=true;
  Branch(*t,"runId",&c.run,"Int_t");Branch(*t,"eventId",&c.event,"Int_t");
  Branch(*t,"protonIndex",&c.proton,"Int_t");Branch(*t,"pionIndex",&c.pion,"Int_t");
  Branch(*t,"mass",&c.mass,"Float_t");Branch(*t,"px",&c.px,"Float_t");
  Branch(*t,"py",&c.py,"Float_t");Branch(*t,"pz",&c.pz,"Float_t");
  if(standalone)Branch(*t,"selected",&selected,"Bool_t");
  else Branch(*t,"inputIndex",&c.input,"Int_t");
  std::vector<int> counts(events?events->size():0,0);
  std::vector<FemtoLambdaAudit::Candidate> out;
  for(Long64_t i=0;i<t->GetEntries();++i){
    Require(t->GetEntry(i)>0,"candidate read failed");if(!selected)continue;
    Require(TMath::Finite(c.mass)&&TMath::Finite(c.px)&&TMath::Finite(c.py)&&TMath::Finite(c.pz),
        "nonfinite selected candidate");
    if(events){
      Require(c.input>0 && size_t(c.input)<=events->size(),"candidate inputIndex outside ledger");
      const Event& e=(*events)[c.input-1];
      Require(e.run==c.run && e.event==c.event && e.status==4,"candidate/event ledger identity mismatch");
      ++counts[c.input-1];
    }
    out.push_back(c);
  }
  if(events)for(size_t i=0;i<events->size();++i)
    Require(counts[i]==(*events)[i].selected,"per-event selected candidate count mismatch");
  t->ResetBranchAddresses();std::sort(out.begin(),out.end());return out;
}
void SameCandidates(const std::vector<FemtoLambdaAudit::Candidate>& a,const std::vector<FemtoLambdaAudit::Candidate>& b) {
  Require(a.size()==b.size(),"selected KF candidate count mismatch");
  for(size_t i=0;i<a.size();++i)
    Require(!(a[i]<b[i]) && !(b[i]<a[i]),"selected KF run/event/index/mass/momentum mismatch");
}
double LabeledCount(TH1& h,const char* label) {
  int found=0;double value=0;
  for(int i=1;i<=h.GetNbinsX();++i)if(std::string(h.GetXaxis()->GetBinLabel(i))==label){
    ++found;value=h.GetBinContent(i);
  }
  Require(found==1,std::string("missing/duplicate event label ")+label);return value;
}
TH2* CentHist(TFile& f,const std::string& name) {
  TH2* h=dynamic_cast<TH2*>(Hist(f,name));Require(h!=0,name+" is not TH2");
  AxisContract(*h->GetXaxis(),200,0.,1.,name+" kstar");
  AxisContract(*h->GetYaxis(),9,-.5,8.5,name+" cent9");
  for(int x=0;x<=h->GetNbinsX()+1;++x)
    Require(h->GetBinContent(x,0)==0 && h->GetBinContent(x,10)==0 &&
        h->GetBinError(h->GetBin(x,0))==0 && h->GetBinError(h->GetBin(x,10))==0,
        name+" centrality under/overflow not empty");
  return h;
}
void InclusiveCentClosure(const TH1& inclusive,const TH2& byCent,const std::string& what) {
  Require(inclusive.GetDimension()==1 && SameAxis(*inclusive.GetXaxis(),*byCent.GetXaxis()),what+" axes");
  Require(inclusive.GetEntries()==byCent.GetEntries(),what+" entries");
  for(int x=0;x<=inclusive.GetNbinsX()+1;++x){
    double count=0,var=0;
    for(int y=1;y<=9;++y){count+=byCent.GetBinContent(x,y);var+=Variance(byCent,byCent.GetBin(x,y));}
    Require(count==inclusive.GetBinContent(x) && SameVariance(var,Variance(inclusive,x)),what+" contents/errors");
  }
}
void SliceClosure(const TH1& slice,const TH2& byCent,int cent,const std::string& what) {
  Require(slice.GetDimension()==1 && SameAxis(*slice.GetXaxis(),*byCent.GetXaxis()),what+" axes");
  double count=0;
  for(int x=0;x<=slice.GetNbinsX()+1;++x){
    const int bin=byCent.GetBin(x,cent+1);
    Require(slice.GetBinContent(x)==byCent.GetBinContent(bin) &&
        slice.GetBinError(x)==byCent.GetBinError(bin),what+" contents/errors");
    count+=slice.GetBinContent(x);
  }
  Require(count==slice.GetEntries(),what+" unit-weight entries");
}
void MassSliceClosure(const TH2& mass,const TH2& byCent,int cent,const std::string& what) {
  AxisContract(*mass.GetXaxis(),200,0.,1.,what+" kstar");
  AxisContract(*mass.GetYaxis(),200,1.05,1.25,what+" mass");
  Require(SameAxis(*mass.GetXaxis(),*byCent.GetXaxis()),what+" canonical kstar axes");
  double entries=0;
  for(int x=0;x<=mass.GetNbinsX()+1;++x){
    double count=0,var=0;
    for(int y=0;y<=mass.GetNbinsY()+1;++y){count+=mass.GetBinContent(x,y);var+=Variance(mass,mass.GetBin(x,y));}
    const int bin=byCent.GetBin(x,cent+1);
    Require(count==byCent.GetBinContent(bin) && SameVariance(var,Variance(byCent,bin)),what+" projection contents/errors");
    entries+=count;
  }
  Require(entries==mass.GetEntries(),what+" unit-weight entries");
}
void PairAliases(TFile& f,TFile& legacy,const std::string& shortName,const std::string& base,bool mixed) {
  const std::string dir=mixed?"mix/":"true/",prefix=mixed?"Mixed_":"",mode=mixed?"ME":"SE";
  TH1* inclusive=Hist(f,"hKstar"+mode+"_"+base);
  TH2* cent=CentHist(f,"hKstar"+mode+"VsCent_"+base);
  InclusiveCentClosure(*inclusive,*cent,mode+" all-mass centrality");
  double total=0;
  for(int c=0;c<9;++c){
    const std::string key=dir+"hKstarMass_"+prefix+shortName+"_CentBin"+Number(c);
    TH2* mass=dynamic_cast<TH2*>(Hist(f,key));Require(mass!=0,key+" is not TH2");
    Require(SameShape(*mass,*Hist(legacy,key)),key+" original class/axes mismatch");
    MassSliceClosure(*mass,*cent,c,key);total+=mass->GetEntries();
  }
  Require(total==inclusive->GetEntries(),mode+" all-mass alias entry mismatch");
  const char* sb[]={"","_SBNeg","_SBPos"};
  const char* channel[]={"_signal","_leftSB","_rightSB"};
  for(int r=0;r<3;++r){
    const std::string old=dir+"hKstar_"+prefix+shortName+sb[r];
    const std::string canon="hKstar"+mode+"_"+base+channel[r];
    Require(SameHist(*Hist(f,old),*Hist(f,canon)),canon+" inclusive alias class/axes/content/error mismatch");
    Require(SameShape(*Hist(f,old),*Hist(legacy,old)),old+" original class/axes mismatch");
    TH2* regionCent=CentHist(f,"hKstar"+mode+"VsCent_"+base+channel[r]);
    InclusiveCentClosure(*Hist(f,old),*regionCent,canon+" centrality");
    for(int c=0;c<9;++c){
      const std::string key=old+"_CentBin"+Number(c);
      Require(SameShape(*Hist(f,key),*Hist(legacy,key)),key+" original class/axes mismatch");
      SliceClosure(*Hist(f,key),*regionCent,c,key);
    }
    const std::string qlab=dir+"hQlab_"+prefix+shortName+sb[r];
    Require(SameShape(*Hist(f,qlab),*Hist(legacy,qlab)),qlab+" original class/axes mismatch");
    Require(Hist(f,qlab)->GetEntries()==Hist(f,old)->GetEntries(),qlab+" entry mismatch");
  }
}
void MixingCounters(TFile& f,double fullME,double& forward,double& reverse) {
  TH2* h=dynamic_cast<TH2*>(Hist(f,"hMixSamplerQA"));Require(h!=0,"mix QA not TH2");
  AxisContract(*h->GetXaxis(),15,-.5,14.5,"mix QA status");
  AxisContract(*h->GetYaxis(),4,-.5,3.5,"mix QA channel");
  // Enum indices from include/FemtoMixingSampler.h; do not include its C++11
  // sampler in this ROOT5-compatible audit macro.
  double q[15]={0};
  for(int x=0;x<=16;++x)for(int y=0;y<=5;++y){
    const double value=h->GetBinContent(x,y);
    Require(value>=0 && TMath::Finite(value) && value==std::floor(value),"invalid mix QA count");
    if(x==0 || x==16 || y==0 || y==5)Require(value==0,"mix QA flow not empty");
    else q[x-1]+=value; // Sum channel slots; do not assume base is channel zero.
  }
  Require(q[1]==q[11]+q[12],"ME eligible forward/reverse mismatch");
  Require(q[2]==q[9]+q[10] && q[2]==fullME,"ME filled forward/reverse/histogram mismatch");
  Require(q[0]==q[2]+q[5] && q[1]==q[0]+q[4] && q[3]==q[1]-q[2],"ME attempted/skipped mismatch");
  Require(q[5]==q[13]+q[14]+q[8],"ME skipped reason mismatch");
  Require(q[4]==0 && q[8]==0,"bufferAll cap loss / unresolved pair direction");
  Require(q[9]<=q[11] && q[10]<=q[12],"ME filled exceeds eligible direction");
  // A direction can legitimately be zero for rare nuclei in 10000 events.
  forward=q[9];reverse=q[10];
}
void AuditEnergy(const std::string& dir,const std::string& energy,std::ostream& report) {
  TFile legacy((dir+"/"+energy+"/legacy.root").c_str(),"READ");
  TFile standalone((dir+"/"+energy+"/standalone.root").c_str(),"READ");
  Require(!legacy.IsZombie()&&!standalone.IsZombie(),"reference ROOT input missing");
  TNamed* kfStatus=dynamic_cast<TNamed*>(standalone.Get("KFRunStatus"));
  Require(kfStatus && std::string(kfStatus->GetTitle())=="completed","standalone KF run not completed");
  const std::vector<FemtoLambdaAudit::Candidate> reference=Candidates(standalone,true);
  TH1* kfEvents=Hist(standalone,"hKfEventSelection");
  Require(LabeledCount(*kfEvents,"read events")==10000,"standalone input count not 10000");
  TTree* standaloneLedger=dynamic_cast<TTree*>(standalone.Get("eventLedger"));
  const bool hasStandaloneLedger=standaloneLedger!=0;
  std::vector<FemtoLambdaAudit::Event> standaloneEvents;
  if(hasStandaloneLedger)standaloneEvents=Events(standalone,false);
  else std::cout<<"NOTE "<<energy<<": standalone has no all-event ledger; full rejected-event identity closure is unavailable. Candidate tuples and read/reconstructed counts are checked.\n";
  const char* shortName[]={"d","t","3He","4He"};
  const char* species[]={"deuteron","triton","he3","he4"};
  std::vector<FemtoLambdaAudit::Event> deuteronEvents;
  for(int s=0;s<4;++s){
    TFile f((dir+"/"+energy+"/"+shortName[s]+".root").c_str(),"READ");
    Require(!f.IsZombie(),"new ROOT input missing");
    TNamed* status=dynamic_cast<TNamed*>(f.Get("FemtoLambdaRunStatus"));
    Require(status&&std::string(status->GetTitle())=="completed","run not completed");
    const std::vector<FemtoLambdaAudit::Event> events=Events(f,true);
    Require(events.size()==10000,"not exactly 10000 input events");
    if(s==0)deuteronEvents=events;else SameEvents(deuteronEvents,events,true);
    if(hasStandaloneLedger)SameEvents(standaloneEvents,events,false);
    Long64_t reconstructed=0,selected=0;
    for(size_t i=0;i<events.size();++i){if(events[i].status==4)++reconstructed;selected+=events[i].selected;}
    Require(Parameter(f,"inputEvents")==10000 && Parameter(f,"reconstructedEvents")==reconstructed &&
        Parameter(f,"selectedLambdaCandidates")==selected,"event ledger metadata count mismatch");
    Require(LabeledCount(*kfEvents,"reconstructed")==reconstructed,"standalone reconstructed event count mismatch");
    const std::vector<FemtoLambdaAudit::Candidate> candidates=Candidates(f,false,&events);
    SameCandidates(reference,candidates);
    Require(selected==Long64_t(candidates.size()),"candidate/event selected total mismatch");
    Require(SameHist(*Hist(f,"hLambda_InvMass"),*Hist(standalone,"hLambda_InvMass")),"Lambda mass class/axes/content/error mismatch");
    TDirectory* oldNuclear=legacy.GetDirectory("true");Require(oldNuclear!=0,"legacy true directory missing");
    TIter next(oldNuclear->GetListOfKeys());TKey* key=0;int checked=0;
    while((key=static_cast<TKey*>(next()))){
      const std::string n=key->GetName();
      if(n.find("hDedxP")!=0 && n.find("hM2P")!=0 && n.find("hPvs")!=0)continue;
      // Do not hide differences by normalizing or relaxing this comparison.
      // In particular, a different centrality RNG/event population is a mismatch,
      // not evidence of successful nuclear-ID migration.
      Require(SameHist(*Hist(f,"true/"+n),*Hist(legacy,"true/"+n)),
          energy+"/"+shortName[s]+" original nuclear QA class/axes/content/error mismatch: "+n);
      ++checked;
    }
    Require(checked>0,"no nuclear QA audited");
    const std::string base="lambda_"+std::string(species[s]);
    PairAliases(f,legacy,shortName[s],base,false);
    PairAliases(f,legacy,shortName[s],base,true);
    const double fullSE=Hist(f,"hKstarSE_"+base)->GetEntries(),fullME=Hist(f,"hKstarME_"+base)->GetEntries();
    double forward=0,reverse=0;MixingCounters(f,fullME,forward,reverse);
    report<<energy<<","<<shortName[s]<<","<<events.size()<<","
          <<Hist(legacy,"hLambda_InvMass")->GetEntries()<<","<<reference.size()<<","
          <<Hist(legacy,"true/hKstar_"+std::string(shortName[s]))->GetEntries()<<","
          <<Hist(f,"true/hKstar_"+std::string(shortName[s]))->GetEntries()<<","
          <<fullSE<<","<<fullME<<","<<checked<<","<<reconstructed<<","<<forward<<","<<reverse
          <<","<<events.size()<<","<<(hasStandaloneLedger?"checked":"unavailable")<<"\n";
    report.flush();
    std::cout<<"PASS "<<energy<<"/"<<shortName[s]<<": ledger, KF tuples, nuclear QA, pair aliases and ME counters\n";
  }
}
}
void auditFemtoLambdaOutputs(const char* directory="rootfile/femto_lambda_kf_validation_20260930",
                            const char* energy="") {
  try {
    const std::string dir(directory),one(energy?energy:"");
    FemtoLambdaAudit::Require(one.empty() || one=="13p5" || one=="3p85","energy must be empty, 13p5 or 3p85");
    const std::string reportPath=dir+(one.empty()?"/audit.csv":"/audit_"+one+".csv");
    std::ofstream report(reportPath.c_str());FemtoLambdaAudit::Require(bool(report),"cannot create audit CSV");
    report<<"energy,species,input_events,legacy_helix_lambda,kf_lambda,legacy_signal_SE,kf_signal_SE,kf_fullmass_SE,kf_fullmass_ME,nuclear_QA_equal,reconstructed_events,ME_forward,ME_reverse,ledger_events_checked,standalone_all_event_ledger\n";
    if(one.empty() || one=="13p5")FemtoLambdaAudit::AuditEnergy(dir,"13p5",report);
    if(one.empty() || one=="3p85")FemtoLambdaAudit::AuditEnergy(dir,"3p85",report);
    report.close();
    std::cout<<"PASS: four-species event identity/status/cent9 closure, exact standalone KF candidates, original nuclear QA, histogram classes/axes/errors, per-cent aliases and ME directional bookkeeping. See explicit standalone-ledger limitation above.\n";
    gSystem->Exit(0);
  } catch(const std::exception& e) {
    std::cerr<<"FAIL: "<<e.what()<<std::endl;gSystem->Exit(1);
  }
}

