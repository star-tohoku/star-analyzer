// Read-only old/new output audit after adding Lambda acceptance and signal halves.
// Example (parameters are the maker YAML values):
// root4star -b -q 'tests/check_femto_lambda_acceptance.C("before.root","after.root","deuteron",1.11596,.00296,3.,1.115683)'
#include "auditFemtoLambdaOutputs.C"
#include "TLorentzVector.h"
#include <set>

namespace FemtoLambdaAcceptanceAudit {
void HistPaths(TDirectory& directory,const std::string& prefix,std::set<std::string>& paths) {
  TIter next(directory.GetListOfKeys());TKey* key=0;
  while((key=static_cast<TKey*>(next()))) {
    TObject* object=directory.Get(key->GetName());
    const std::string name=prefix+key->GetName();
    TDirectory* nested=dynamic_cast<TDirectory*>(object);
    if(nested)HistPaths(*nested,name+"/",paths);
    else if(dynamic_cast<TH1*>(object))paths.insert(name);
  }
}
void Partition(TFile& file,const std::string& stem) {
  using namespace FemtoLambdaAudit;
  TH1* signal=Hist(file,stem+"_signal");
  TH1* low=Hist(file,stem+"_signalLow");
  TH1* high=Hist(file,stem+"_signalHigh");
  Require(SameShape(*signal,*low) && SameShape(*signal,*high),stem+" split class/axes");
  Require(signal->GetEntries()==low->GetEntries()+high->GetEntries(),stem+" split entries");
  for(int bin=0;bin<CellCount(*signal);++bin)
    Require(signal->GetBinContent(bin)==low->GetBinContent(bin)+high->GetBinContent(bin) &&
            SameVariance(Variance(*signal,bin),Variance(*low,bin)+Variance(*high,bin)),
            stem+" split bin/variance including flow");
}
void UnitWeightMap(TH1& h,const std::string& what) {
  using namespace FemtoLambdaAudit;
  Require(h.GetDimension()==2,what+" not TH2");
  double sum=0;
  for(int bin=0;bin<CellCount(h);++bin) {
    const double content=h.GetBinContent(bin);
    Require(TMath::Finite(content) && content>=0 && content==std::floor(content),
            what+" invalid count");
    Require(SameVariance(content,Variance(h,bin)),what+" non-unit-weight variance");
    sum+=content;
  }
  Require(sum==h.GetEntries(),what+" entries including flow");
}
}

void check_femto_lambda_acceptance(const char* oldFile,const char* newFile,const char* species,
                                  double mean,double sigma,double nSigma,double pairMass) {
  try {
    using namespace FemtoLambdaAudit;
    using namespace FemtoLambdaAcceptanceAudit;
    Require(oldFile && *oldFile && newFile && *newFile && species && *species,
            "explicit old/new ROOT paths and species required");
    const std::string nuclear(species);
    Require(nuclear=="deuteron" || nuclear=="triton" || nuclear=="he3" || nuclear=="he4",
            "unknown target nuclear species");
    Require(TMath::Finite(mean) && TMath::Finite(sigma) && TMath::Finite(nSigma) &&
            TMath::Finite(pairMass) && sigma>0 && nSigma>0 && pairMass>0,
            "explicit finite maker window/mass values required");
    TFile old(oldFile,"READ"),now(newFile,"READ");
    Require(!old.IsZombie() && !now.IsZombie(),"old/new ROOT unavailable");
    TNamed* oldStatus=dynamic_cast<TNamed*>(old.Get("FemtoLambdaRunStatus"));
    TNamed* nowStatus=dynamic_cast<TNamed*>(now.Get("FemtoLambdaRunStatus"));
    Require(oldStatus && nowStatus && std::string(oldStatus->GetTitle())=="completed" &&
            std::string(nowStatus->GetTitle())=="completed","incomplete old/new run");
    std::set<std::string> oldHist,newHist;
    HistPaths(old,"",oldHist);HistPaths(now,"",newHist);
    Require(oldHist.size()==180 && newHist.size()==192,"expected 180 old / 192 new histogram keys");
    for(std::set<std::string>::const_iterator i=oldHist.begin();i!=oldHist.end();++i)
      Require(SameHist(*Hist(old,*i),*Hist(now,*i)),"existing histogram changed: "+*i);
    std::set<std::string> expectedNew;
    const char* maps[]={"hLambda_PtVsYLab_signal","hLambdaProton_PtVsYLab_signal",
                       "hLambdaPion_PtVsYLab_signal"};
    for(int i=0;i<3;++i)expectedNew.insert(maps[i]);
    expectedNew.insert("hNucleus_PtVsYLab_"+nuclear);
    const std::string channel="lambda_"+nuclear;
    for(int mode=0;mode<2;++mode)for(int dimension=0;dimension<2;++dimension) {
      const std::string stem=std::string("hKstar")+(mode?"ME":"SE")+(dimension?"VsCent":"")+"_"+channel;
      expectedNew.insert(stem+"_signalLow");expectedNew.insert(stem+"_signalHigh");
      Partition(now,stem);
    }
    for(std::set<std::string>::const_iterator i=newHist.begin();i!=newHist.end();++i)
      if(!oldHist.count(*i))Require(expectedNew.erase(*i)==1,"unexpected new histogram: "+*i);
    Require(expectedNew.empty(),"missing new histogram");
    for(int mode=0;mode<2;++mode) {
      const std::string stem=std::string("hKstar")+(mode?"ME":"SE");
      const char* halves[]={"_signalLow","_signalHigh"};
      for(int half=0;half<2;++half) {
        TH1* one=Hist(now,stem+"_"+channel+halves[half]);
        TH2* two=CentHist(now,stem+"VsCent_"+channel+halves[half]);
        InclusiveCentClosure(*one,*two,stem+halves[half]);
      }
    }

    const std::vector<FemtoLambdaAudit::Event> oldEvents=Events(old,true),newEvents=Events(now,true);
    SameEvents(oldEvents,newEvents,true);
    const char* parameters[]={"inputEvents","reconstructedEvents","selectedLambdaCandidates"};
    for(int i=0;i<3;++i)Require(Parameter(old,parameters[i])==Parameter(now,parameters[i]),
                               std::string("metadata count changed: ")+parameters[i]);
    const std::vector<FemtoLambdaAudit::Candidate> before=Candidates(old,false,&oldEvents);
    const std::vector<FemtoLambdaAudit::Candidate> after=Candidates(now,false,&newEvents);
    SameCandidates(before,after);
    TH2* lambdaMap=dynamic_cast<TH2*>(Hist(now,maps[0]));
    Require(lambdaMap!=0,"Lambda acceptance not TH2");
    TH2* expected=static_cast<TH2*>(lambdaMap->Clone("audit_expected_lambda_acceptance"));
    expected->SetDirectory(0);expected->Reset();
    const double width=sigma*nSigma;
    Long64_t signals=0;
    for(size_t i=0;i<after.size();++i) {
      const Candidate& c=after[i];
      // Exactly the existing unbinned MassRegion predicate, NOT a binned mass integral.
      if(std::fabs(double(c.mass)-mean)>width)continue;
      ++signals;
      TLorentzVector p4;p4.SetXYZM(c.px,c.py,c.pz,pairMass);
      // FemtoCandidate stores the fixed-mass energy as Float_t before P4() is read.
      p4.SetE(static_cast<Float_t>(p4.E()));
      expected->Fill(p4.Rapidity(),p4.Pt());
    }
    Require(SameHist(*lambdaMap,*expected),"Lambda lab rapidity/pT differs from candidate ledger");
    delete expected;
    for(int i=0;i<3;++i) {
      TH1* map=Hist(now,maps[i]);
      UnitWeightMap(*map,maps[i]);
      Require(SameShape(*lambdaMap,*map) && map->GetEntries()==signals,
              std::string(maps[i])+" signal-candidate count/axes mismatch");
    }
    TH1* nuclei=Hist(now,"hNucleus_PtVsYLab_"+nuclear);
    UnitWeightMap(*nuclei,"selected nucleus acceptance");
    Require(SameShape(*lambdaMap,*nuclei),"nucleus/Lambda lab axes differ");
    std::cout<<"PASS acceptance/split audit "<<species
             <<": old_histograms_exact="<<oldHist.size()
             <<", new_histograms="<<newHist.size()
             <<", events_exact="<<newEvents.size()
             <<", selected_candidates_exact="<<after.size()
             <<", signal_lambda_and_each_daughter="<<signals
             <<", selected_nucleus_entries="<<nuclei->GetEntries()<<std::endl;
    std::cout<<"NOTE: nucleus multiplicity and daughter momentum-source semantics are verified by "
                "the in-memory Pico fixture in femto_lambda_legacy.C; old ROOT has no nuclear/daughter ledger."
             <<std::endl;
    gSystem->Exit(0);
  } catch(const std::exception& e) {
    std::cerr<<"FAIL acceptance/split audit: "<<e.what()<<std::endl;gSystem->Exit(1);
  }
}
