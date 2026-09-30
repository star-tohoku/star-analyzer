// Read-only audit for the additive centrality/mass-region k* versus pair-mT maps.
// root4star -b -q 'tests/check_femto_lambda_mt.C+("old.root","new.root","deuteron")'
#include "auditFemtoLambdaOutputs.C"
#include <set>

namespace FemtoLambdaMtAudit {
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
TH2* MtHist(TFile& file,const std::string& name) {
  using namespace FemtoLambdaAudit;
  TH2* h=dynamic_cast<TH2*>(Hist(file,name));
  Require(h && std::string(h->ClassName())=="TH2D",name+" not TH2D");
  AxisContract(*h->GetXaxis(),200,0.,1.,name+" kstar");
  AxisContract(*h->GetYaxis(),200,0.,10.,name+" pair mT");
  return h;
}
void ProjectionClosure(TH2& mt,const TH2& canonical,int cent,const std::string& what) {
  using namespace FemtoLambdaAudit;
  Require(SameAxis(*mt.GetXaxis(),*canonical.GetXaxis()),what+" kstar axes differ");
  // Explicitly include both mT flow bins: the TH2 bin range is not a pair cut.
  TH1D* projection=mt.ProjectionX(("audit_mt_projection_"+Number(cent)).c_str(),0,mt.GetNbinsY()+1,"e");
  Require(projection!=0,what+" ProjectionX failed");
  projection->SetDirectory(0);
  double entries=0;
  for(int x=0;x<=mt.GetNbinsX()+1;++x) {
    const int bin=canonical.GetBin(x,cent+1);
    double count=0,variance=0;
    for(int y=0;y<=mt.GetNbinsY()+1;++y) {
      const int cell=mt.GetBin(x,y);
      const double value=mt.GetBinContent(cell);
      Require(TMath::Finite(value) && value>=0 && value==std::floor(value),what+" invalid unit count");
      Require(SameVariance(value,Variance(mt,cell)),what+" non-unit-weight error");
      count+=value;variance+=Variance(mt,cell);
    }
    Require(count==canonical.GetBinContent(bin) &&
            SameVariance(variance,Variance(canonical,bin)),what+" sum-to-canonical mismatch");
    Require(projection->GetBinContent(x)==canonical.GetBinContent(bin) &&
            SameVariance(Variance(*projection,x),Variance(canonical,bin)),what+" ProjectionX mismatch");
    entries+=count;
  }
  Require(entries==mt.GetEntries(),what+" entries including flow mismatch");
  delete projection;
}
void Partition(const TH2& signal,const TH2& low,const TH2& high,const std::string& what) {
  using namespace FemtoLambdaAudit;
  Require(SameShape(signal,low) && SameShape(signal,high),what+" split shape mismatch");
  Require(signal.GetEntries()==low.GetEntries()+high.GetEntries(),what+" split entries mismatch");
  for(int bin=0;bin<CellCount(signal);++bin)
    Require(signal.GetBinContent(bin)==low.GetBinContent(bin)+high.GetBinContent(bin) &&
            SameVariance(Variance(signal,bin),Variance(low,bin)+Variance(high,bin)),
            what+" split content/variance including flow mismatch");
}
}

void check_femto_lambda_mt(const char* oldFile,const char* newFile,const char* species) {
  try {
    using namespace FemtoLambdaAudit;
    using namespace FemtoLambdaMtAudit;
    Require(oldFile && *oldFile && newFile && *newFile && species && *species,
            "explicit old/new ROOT paths and target species required");
    const std::string nuclear(species);
    Require(nuclear=="deuteron" || nuclear=="triton" || nuclear=="he3" || nuclear=="he4",
            "unknown target nuclear species");
    const std::string shortName=nuclear=="deuteron"?"d":nuclear=="triton"?"t":nuclear=="he3"?"3He":"4He";
    TFile old(oldFile,"READ"),now(newFile,"READ");
    Require(!old.IsZombie() && !now.IsZombie(),"old/new ROOT unavailable");
    TNamed* oldStatus=dynamic_cast<TNamed*>(old.Get("FemtoLambdaRunStatus"));
    TNamed* nowStatus=dynamic_cast<TNamed*>(now.Get("FemtoLambdaRunStatus"));
    Require(oldStatus && nowStatus && std::string(oldStatus->GetTitle())=="completed" &&
            std::string(nowStatus->GetTitle())=="completed","incomplete old/new run");
    TNamed* definition=dynamic_cast<TNamed*>(now.Get("FemtoLambdaPairMtDefinition"));
    Require(definition && std::string(definition->GetTitle()).find("average_mass")==0,
            "missing/wrong pair-mT definition metadata");
    std::set<std::string> oldHist,newHist;
    HistPaths(old,"",oldHist);HistPaths(now,"",newHist);
    Require(oldHist.size()==192 && newHist.size()==282,"expected 192 old / 282 new histogram keys");
    for(std::set<std::string>::const_iterator i=oldHist.begin();i!=oldHist.end();++i)
      Require(SameHist(*Hist(old,*i),*Hist(now,*i)),"existing histogram changed: "+*i);
    std::set<std::string> expectedNew;
    const char* suffix[]={"","_SBPos","_SBNeg","_signalLow","_signalHigh"};
    const char* canonical[]={"_signal","_rightSB","_leftSB","_signalLow","_signalHigh"};
    for(int mode=0;mode<2;++mode)for(int cent=0;cent<9;++cent) {
      const std::string prefix=mode?"mix/hKstarMt_Mixed_":"true/hKstarMt_";
      TH2* region[5]={0,0,0,0,0};
      for(int r=0;r<5;++r) {
        const std::string name=prefix+shortName+suffix[r]+"_CentBin"+Number(cent);
        expectedNew.insert(name);region[r]=MtHist(now,name);
        TH2* reference=CentHist(now,std::string("hKstar")+(mode?"ME":"SE")+
                               "VsCent_lambda_"+nuclear+canonical[r]);
        ProjectionClosure(*region[r],*reference,cent,name);
        if(r<3) {
          const std::string oldKey=std::string(mode?"mix/hKstar_Mixed_":"true/hKstar_")+
                                   shortName+suffix[r]+"_CentBin"+Number(cent);
          SliceClosure(*Hist(now,oldKey),*reference,cent,oldKey+" canonical alias");
        }
      }
      Partition(*region[0],*region[3],*region[4],prefix+shortName+"_CentBin"+Number(cent));
    }
    for(std::set<std::string>::const_iterator i=newHist.begin();i!=newHist.end();++i)
      if(!oldHist.count(*i))Require(expectedNew.erase(*i)==1,"unexpected new histogram: "+*i);
    Require(expectedNew.empty(),"missing new mT histogram");
    const std::vector<FemtoLambdaAudit::Event> oldEvents=Events(old,true),newEvents=Events(now,true);
    SameEvents(oldEvents,newEvents,true);
    const char* parameters[]={"inputEvents","reconstructedEvents","selectedLambdaCandidates"};
    for(int i=0;i<3;++i)
      Require(Parameter(old,parameters[i])==Parameter(now,parameters[i]),
              std::string("metadata count changed: ")+parameters[i]);
    const std::vector<FemtoLambdaAudit::Candidate> before=Candidates(old,false,&oldEvents);
    const std::vector<FemtoLambdaAudit::Candidate> after=Candidates(now,false,&newEvents);
    SameCandidates(before,after);
    std::cout<<"PASS mT audit "<<species<<": old_histograms_exact="<<oldHist.size()
             <<", new_histograms="<<newHist.size()<<", new_mT_maps=90"
             <<", events_exact="<<newEvents.size()<<", selected_candidates_exact="<<after.size()
             <<", ProjectionX_and_variance_closure=90, split_2D_closure=18"<<std::endl;
    std::cout<<"NOTE: pair-mT numerical definition is tested with non-collinear P4 fixtures in "
                "femto_lambda_legacy.C; event output has no nuclear pair ledger."<<std::endl;
    gSystem->Exit(0);
  } catch(const std::exception& e) {
    std::cerr<<"FAIL mT audit: "<<e.what()<<std::endl;gSystem->Exit(1);
  }
}
