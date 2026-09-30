// Validation-only control macro: compare original StNuclearIdMaker used by
// anaLambdaNuclearId.C against the new common builder on identical input events.
// Load StNuclearIdMaker in addition to the normal config/common/Pico libraries.
#include "StChain.h"
#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StMaker/StNuclearIdMaker/StNuclearIdMaker.h"
#include "StMaker/common/FemtoLambdaLegacy.h"
#include "StMaker/common/CentralityHelper.h"
#include "ConfigManager.h"
#include "cuts/CentralityCutConfig.h"
#include "TChain.h"
#include "TMemFile.h"
#include "TDirectory.h"
#include "TKey.h"
#include "TH1.h"
#include "TRandom.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {
struct NucleusClosureEntry { int id,type; TVector3 momentum; };
bool NucleusClosureLess(const NucleusClosureEntry& a,const NucleusClosureEntry& b) {
 if(a.id!=b.id)return a.id<b.id;
 if(a.type!=b.type)return a.type<b.type;
 if(a.momentum.X()!=b.momentum.X())return a.momentum.X()<b.momentum.X();
 if(a.momentum.Y()!=b.momentum.Y())return a.momentum.Y()<b.momentum.Y();
 return a.momentum.Z()<b.momentum.Z();
}
}
int compare_femto_lambda_nuclear(const char* input,const char* mainconf,int nEvents=1000) {
 StChain* chain=0;
 try {
  if(!input || !*input || !mainconf || !*mainconf || nEvents<=0) throw std::runtime_error("explicit arguments required");
  if(!ConfigManager::GetInstance().LoadConfig(mainconf)) throw std::runtime_error("LoadConfig");
  FemtoLambdaLegacy helper;
  if(!helper.Init(mainconf,std::cerr)) throw std::runtime_error("helper Init");
  CentralityHelper cent;
  const CentralityCutConfig& cuts=ConfigManager::GetInstance().GetCentralityCuts();
  if(!cent.Init(cuts)) throw std::runtime_error("test centrality Init");
  chain=new StChain("nuclearClosure");
  StPicoDstMaker* reader=new StPicoDstMaker(StPicoDstMaker::IoRead,input,"picoDst");
  reader->SetStatus("*",0);reader->SetStatus("Event",1);reader->SetStatus("Track",1);
  reader->SetStatus("BTofPidTraits",1);reader->SetStatus("BTofHit",1);
  // Name deliberately resolves nuclearHist in the new mainconf; implementation is
  // the unchanged ORIGINAL StNuclearIdMaker, not StLambdaNuclearMaker.
  StNuclearIdMaker* original=new StNuclearIdMaker("nuclear",reader,"");
  if(chain->Init()!=kStOK || !reader->chain() || reader->chain()->GetEntries()<nEvents)
   throw std::runtime_error("chain Init/input shortage");
  unsigned long long candidates=0;
  for(int e=0;e<nEvents;++e) {
   chain->Clear();
   // Replay original multiplicity-smearing RNG for test-only duplicate centrality.
   TRandom* replay=static_cast<TRandom*>(gRandom->Clone());
   if(chain->Make(e)!=kStOK) {delete replay;throw std::runtime_error("input Make failure");}
   StPicoDst* dst=reader->picoDst();StPicoEvent* event=dst?dst->event():0;
   if(!event)throw std::runtime_error("missing Pico event");
   int c9=-1,c16=-1;double ref=-1,weight=1.;CentralityRejectReason reason=kCentralityOk;
   int raw=cuts.mode=="fxtmult"?event->fxtMult():event->refMult();
   bool accepted=true;
   TRandom* productionRandom=gRandom;gRandom=replay;
   if(cent.IsEnabled()) accepted=cent.CheckBadRun(event->runId(),reason) &&
    cent.CheckPileup(raw,event->nBTOFMatch(),event->primaryVertex().Z(),reason) &&
    cent.ComputeBins(event,raw,event->primaryVertex().Z(),c9,c16,ref,weight,reason) &&
    cent.AcceptCentBin(c9,ref,reason);
   gRandom=productionRandom;delete replay;
   FemtoCandidateStore store;
   if(accepted)helper.CollectNuclei(dst,e,c9,store);
   std::vector<NucleusClosureEntry> a,b;
   const std::vector<Int_t>& ids=original->GetNuclearIdList();
   const std::vector<Int_t>& types=original->GetNuclearTypeList();
   const std::vector<TVector3>& moms=original->GetNuclearMomList();
   for(size_t i=0;i<ids.size();++i) {NucleusClosureEntry v={ids[i],types[i],moms[i]};a.push_back(v);}
   const char* keys[]={"deuteron","triton","he3","he4"};
   for(int s=0;s<4;++s)for(size_t i=0;i<store[keys[s]].size();++i) {
    const FemtoCandidate& n=store[keys[s]][i];
    NucleusClosureEntry v={dst->track(n.trk.trackIndex)->id(),s,n.P4().Vect()};b.push_back(v);
   }
   std::sort(a.begin(),a.end(),NucleusClosureLess);std::sort(b.begin(),b.end(),NucleusClosureLess);
   if(a.size()!=b.size())throw std::runtime_error("selected nuclear multiplicity mismatch");
   for(size_t i=0;i<a.size();++i) {
    if(a[i].id!=b[i].id || a[i].type!=b[i].type ||
      (a[i].momentum-b[i].momentum).Mag()>1.e-6)
       throw std::runtime_error("selected ID/species/Z momentum mismatch");
   }
   candidates+=a.size();
  }
  TMemFile out("femto_nuclear_closure","RECREATE");
  TDirectory* legacy=out.mkdir("original");legacy->cd();original->WriteHistograms();
  TDirectory* modern=out.mkdir("new");if(!helper.Write(modern))throw std::runtime_error("helper Write");
  int histograms=0;
  TIter next(legacy->GetListOfKeys());TKey* key=0;
  while((key=static_cast<TKey*>(next()))) {
   const std::string name=key->GetName();
   if(name.find("hDedxP")!=0 && name.find("hPvs")!=0 && name!="hM2P")continue;
   TH1* old=dynamic_cast<TH1*>(legacy->Get(name.c_str()));
   TH1* now=dynamic_cast<TH1*>(modern->Get(("true/"+name).c_str()));
   if(!old || !now || old->GetNbinsX()!=now->GetNbinsX() || old->GetNbinsY()!=now->GetNbinsY() || old->GetEntries()!=now->GetEntries())
    throw std::runtime_error("nuclear QA histogram dimensions/entries mismatch: "+name);
   for(int bin=0;bin<(old->GetNbinsX()+2)*(old->GetNbinsY()+2);++bin)
    if(old->GetBinContent(bin)!=now->GetBinContent(bin) || old->GetBinError(bin)!=now->GetBinError(bin))
     throw std::runtime_error("nuclear QA bin mismatch: "+name);
   ++histograms;
  }
  if(histograms!=27) throw std::runtime_error("expected all 27 configured common nuclear QA histograms");
  chain->Finish();delete chain;chain=0;
  std::cout<<"PASS original anaLambdaNuclearId nuclear closure: events="<<nEvents
   <<" candidates="<<candidates<<" histograms="<<histograms<<std::endl;
  return 0;
 }catch(const std::exception& e){
  std::cerr<<"FAIL original nuclear closure: "<<e.what()<<std::endl;
  delete chain;return 1;
 }
}


