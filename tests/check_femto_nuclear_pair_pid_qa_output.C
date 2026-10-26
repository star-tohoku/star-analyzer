void check_femto_nuclear_pair_pid_qa_output(const char*,const char*,const char*,int,int);
#ifndef __CINT__
#include <TFile.h>
#include <TDirectory.h>
#include <TKey.h>
#include <TH1.h>
#include <TH2.h>
#include <THnSparse.h>
#include <TTree.h>
#include <TLeaf.h>
#include <TSystem.h>
#include <TError.h>
#include <fstream>
#include <algorithm>
#include <stdexcept>
#include <string>
#include <cmath>
#include <iostream>

namespace NuclearPidPairOutputCheck {
long long checks=0,histograms=0,trees=0,newHistograms=0;
void Need(bool value,const std::string& label) {
  ++checks;if(!value)throw std::runtime_error(label);
}
int Cells(TH1* h) {
  int n=h->GetNbinsX()+2;
  if(h->GetDimension()>1)n*=h->GetNbinsY()+2;
  if(h->GetDimension()>2)n*=h->GetNbinsZ()+2;
  return n;
}
void SameHistogram(TH1* a,TH1* b,const std::string& path) {
  Need(b&&std::string(a->ClassName())==b->ClassName(),"Missing/changed old histogram: "+path);
  Need(a->GetDimension()==b->GetDimension()&&Cells(a)==Cells(b),"Old histogram shape: "+path);
  for(int d=0;d<a->GetDimension();++d) {
    TAxis* x=d==0?a->GetXaxis():d==1?a->GetYaxis():a->GetZaxis();
    TAxis* y=d==0?b->GetXaxis():d==1?b->GetYaxis():b->GetZaxis();
    Need(x->GetNbins()==y->GetNbins(),"Old axis count: "+path);
    for(int i=1;i<=x->GetNbins()+1;++i)Need(x->GetBinLowEdge(i)==y->GetBinLowEdge(i),"Old axis edge: "+path);
  }
  Need(a->GetEntries()==b->GetEntries(),"Old entries changed: "+path);
  for(int i=0;i<Cells(a);++i) {
    Need(a->GetBinContent(i)==b->GetBinContent(i),"Old content changed: "+path);
    Need(a->GetBinError(i)==b->GetBinError(i),"Old error changed: "+path);
  }
  ++histograms;
}
void SameTree(TTree* a,TTree* b,const std::string& path) {
  Need(b&&a->GetEntries()==b->GetEntries(),"Old tree entries changed: "+path);
  TObjArray* leaves=a->GetListOfLeaves();
  Need(leaves->GetEntries()==b->GetListOfLeaves()->GetEntries(),"Old tree schema changed: "+path);
  for(Long64_t entry=0;entry<a->GetEntries();++entry) {
    a->GetEntry(entry);b->GetEntry(entry);
    for(int j=0;j<leaves->GetEntries();++j) {
      TLeaf* x=static_cast<TLeaf*>(leaves->At(j));TLeaf* y=b->GetLeaf(x->GetName());
      Need(y&&std::string(x->GetTypeName())==y->GetTypeName()&&x->GetNdata()==y->GetNdata(),"Old leaf changed: "+path);
      for(int k=0;k<x->GetNdata();++k) {
        const double p=x->GetValue(k),q=y->GetValue(k);
        Need(p==q||(std::isnan(p)&&std::isnan(q)),"Old numeric leaf changed: "+path+"/"+x->GetName());
      }
    }
  }
  ++trees;
}
void Compare(TDirectory* source,TFile& on,const std::string& prefix) {
  TIter next(source->GetListOfKeys());TKey* key=0;
  while((key=static_cast<TKey*>(next()))) {
    const std::string path=prefix+key->GetName();TObject* a=source->Get(key->GetName());
    if(dynamic_cast<TDirectory*>(a))Compare(static_cast<TDirectory*>(a),on,path+"/");
    else if(dynamic_cast<TH1*>(a))SameHistogram(static_cast<TH1*>(a),dynamic_cast<TH1*>(on.Get(path.c_str())),path);
    else if(dynamic_cast<TTree*>(a))SameTree(static_cast<TTree*>(a),dynamic_cast<TTree*>(on.Get(path.c_str())),path);
  }
}
void InspectNew(TDirectory* source,TFile& off,const std::string& prefix) {
  TIter next(source->GetListOfKeys());TKey* key=0;
  while((key=static_cast<TKey*>(next()))) {
    const std::string path=prefix+key->GetName();TObject* a=source->Get(key->GetName());
    if(dynamic_cast<TDirectory*>(a))InspectNew(static_cast<TDirectory*>(a),off,path+"/");
    else if(dynamic_cast<TH1*>(a)&&!off.Get(path.c_str())) {
      TH2* h=dynamic_cast<TH2*>(a);Need(h&&std::string(h->ClassName())=="TH2D","Unexpected new histogram type: "+path);
      Need(h->GetNbinsX()==200&&h->GetXaxis()->GetXmin()==0&&h->GetXaxis()->GetXmax()==1,"New kstar axis: "+path);
      Need(h->GetNbinsY()==120,"New diagnostic Y bins: "+path);
      const bool ns=path.find("hKstarNSigma")!=std::string::npos;
      Need(ns||path.find("hKstarPOverQ")!=std::string::npos,"Unknown new TH2: "+path);
      Need(h->GetYaxis()->GetXmin()==(ns?-12.:0.)&&h->GetYaxis()->GetXmax()==(ns?12.:6.),"New Y range: "+path);
      double count=0.;
      for(int cell=0;cell<Cells(h);++cell) {
        const double n=h->GetBinContent(cell),error=h->GetBinError(cell);count+=n;
        Need(std::isfinite(n)&&n>=0&&std::fabs(error*error-n)<=1.e-10*std::max(1.,n),"Non-unit-count diagnostic errors: "+path);
      }
      Need(count==h->GetEntries(),"Raw diagnostic count/entries mismatch: "+path);++newHistograms;
    }
  }
}
}

void check_femto_nuclear_pair_pid_qa_output(const char* offPath,const char* onPath,const char* jsonPath,int expectedNew,int expectedEvents) {
  using namespace NuclearPidPairOutputCheck;
  try {
    checks=histograms=trees=newHistograms=0;
    Need(offPath&&onPath&&jsonPath&&expectedNew>0&&expectedEvents>0,"Explicit off/on/output/count/events required");
    TFile off(offPath,"READ"),on(onPath,"READ");
    Need(!off.IsZombie()&&!on.IsZombie(),"Cannot read off/on ROOT");
    Compare(&off,on,"");InspectNew(&on,off,"");Need(newHistograms==expectedNew,"Unexpected added histogram count");
    TTree* ledger=dynamic_cast<TTree*>(on.Get("eventLedger"));
    Need(ledger&&ledger->GetEntries()==expectedEvents,"Input event ledger differs from requested test events");
    Need(on.Get("FemtoLambdaNuclearPairPidQaDefinition")&&on.Get("NuclearPidPairJointMetadata"),"Missing new QA definitions");
    Need(!off.Get("hNuclearPidPairJoint"),"QA-off unexpectedly has joint sparse");
    THnSparseD* joint=dynamic_cast<THnSparseD*>(on.Get("hNuclearPidPairJoint"));
    Need(joint&&joint->GetNdimensions()==8,"Missing 8D THnSparseD joint");
    const int bins[]={200,120,120,9,3,2,4,2};
    for(int d=0;d<8;++d)Need(joint->GetAxis(d)->GetNbins()==bins[d],"Joint axis bin mismatch");
    Need(std::isfinite(joint->GetEntries())&&joint->GetEntries()>=0,"Invalid joint entries");
    std::ofstream output(jsonPath);Need(output.good(),"Cannot create explicit readback JSON");
    output<<"{\n  \"complete\": true,\n  \"oldHistogramsExactlyEqual\": "<<histograms
      <<",\n  \"oldTreesExactlyEqual\": "<<trees<<",\n  \"newTH2\": "<<newHistograms
      <<",\n  \"jointEntries\": "<<joint->GetEntries()<<",\n  \"jointAllocatedBins\": "<<joint->GetNbins()
      <<",\n  \"checks\": "<<checks<<"\n}\n";
    output.close();Need(!output.fail(),"Readback JSON write failed");
    std::cout<<"NUCLEAR_PID_PAIR_OUTPUT_PASS oldHist="<<histograms<<" oldTrees="<<trees<<" newTH2="<<newHistograms
      <<" jointEntries="<<joint->GetEntries()<<" checks="<<checks<<std::endl;
  } catch(const std::exception& e) {Error("check_femto_nuclear_pair_pid_qa_output","%s",e.what());gSystem->Exit(1);}
}
#endif
