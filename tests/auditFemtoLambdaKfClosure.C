#include "auditFemtoLambdaOutputs.C"

// Intermediate read-only closure; deliberately does NOT claim legacy/pair audit.
// The full auditFemtoLambdaOutputs remains mandatory after legacy.root finishes.
void auditFemtoLambdaKfClosure(const char* directory="rootfile/femto_lambda_kf_validation_20260930",
                              const char* energy="13p5") {
  try {
    const std::string dir(directory),e(energy?energy:"");
    FemtoLambdaAudit::Require(e=="13p5" || e=="3p85","energy must be 13p5 or 3p85");
    TFile standalone((dir+"/"+e+"/standalone.root").c_str(),"READ");
    FemtoLambdaAudit::Require(!standalone.IsZombie(),"standalone ROOT missing");
    TNamed* status=dynamic_cast<TNamed*>(standalone.Get("KFRunStatus"));
    FemtoLambdaAudit::Require(status && std::string(status->GetTitle())=="completed","standalone incomplete");
    const std::vector<FemtoLambdaAudit::Candidate> reference=FemtoLambdaAudit::Candidates(standalone,true);
    TH1* kfEvents=FemtoLambdaAudit::Hist(standalone,"hKfEventSelection");
    FemtoLambdaAudit::Require(FemtoLambdaAudit::LabeledCount(*kfEvents,"read events")==10000,
        "standalone input count not 10000");
    const double kfReconstructed=FemtoLambdaAudit::LabeledCount(*kfEvents,"reconstructed");
    TTree* standaloneLedger=dynamic_cast<TTree*>(standalone.Get("eventLedger"));
    std::vector<FemtoLambdaAudit::Event> standaloneEvents;
    if(standaloneLedger)standaloneEvents=FemtoLambdaAudit::Events(standalone,false);
    else std::cout<<"NOTE: standalone has no all-event ledger; its rejected-event identities cannot be compared.\n";
    const std::string path=dir+"/audit_kf_"+e+".csv";
    std::ofstream report(path.c_str());
    FemtoLambdaAudit::Require(bool(report),"cannot create KF audit CSV");
    report<<"energy,species,input_events,reconstructed_events,selected_kf_candidates,ledger_events_checked,standalone_all_event_ledger,legacy_and_pair_audit\n";
    const char* names[]={"d","t","3He","4He"};
    std::vector<FemtoLambdaAudit::Event> referenceEvents;
    for(int s=0;s<4;++s) {
      TFile f((dir+"/"+e+"/"+names[s]+".root").c_str(),"READ");
      FemtoLambdaAudit::Require(!f.IsZombie(),"Femto ROOT missing");
      TNamed* runStatus=dynamic_cast<TNamed*>(f.Get("FemtoLambdaRunStatus"));
      FemtoLambdaAudit::Require(runStatus && std::string(runStatus->GetTitle())=="completed","Femto run incomplete");
      const std::vector<FemtoLambdaAudit::Event> events=FemtoLambdaAudit::Events(f,true);
      FemtoLambdaAudit::Require(events.size()==10000,"Femto input count not 10000");
      if(s==0)referenceEvents=events;else FemtoLambdaAudit::SameEvents(referenceEvents,events,true);
      if(standaloneLedger)FemtoLambdaAudit::SameEvents(standaloneEvents,events,false);
      Long64_t reconstructed=0,selected=0;
      for(size_t i=0;i<events.size();++i) {
        if(events[i].status==4)++reconstructed;
        selected+=events[i].selected;
      }
      FemtoLambdaAudit::Require(FemtoLambdaAudit::Parameter(f,"inputEvents")==10000 &&
          FemtoLambdaAudit::Parameter(f,"reconstructedEvents")==reconstructed &&
          FemtoLambdaAudit::Parameter(f,"selectedLambdaCandidates")==selected,"ledger/metadata count mismatch");
      FemtoLambdaAudit::Require(kfReconstructed==reconstructed,"standalone reconstructed count mismatch");
      const std::vector<FemtoLambdaAudit::Candidate> candidates=FemtoLambdaAudit::Candidates(f,false,&events);
      FemtoLambdaAudit::SameCandidates(reference,candidates);
      FemtoLambdaAudit::Require(selected==Long64_t(candidates.size()),"selected total mismatch");
      FemtoLambdaAudit::Require(FemtoLambdaAudit::SameHist(*FemtoLambdaAudit::Hist(f,"hLambda_InvMass"),
          *FemtoLambdaAudit::Hist(standalone,"hLambda_InvMass")),"KF mass histogram class/axes/content/error mismatch");
      report<<e<<","<<names[s]<<","<<events.size()<<","<<reconstructed<<","<<selected<<","
            <<events.size()<<","<<(standaloneLedger?"checked":"unavailable")<<",not_run\n";
      report.flush();
      std::cout<<"PASS "<<e<<"/"<<names[s]<<": "<<events.size()<<" inputs, "<<reconstructed
               <<" reconstructed, "<<selected<<" selected KF tuples; four-species event ledger closure\n";
    }
    report.close();
    std::cout<<"PASS: intermediate KF-only closure. Original anaLambdaNuclearId nuclear QA / pair aliases / ME audit NOT run here.\n";
    gSystem->Exit(0);
  } catch(const std::exception& error) {
    std::cerr<<"FAIL: "<<error.what()<<std::endl;
    gSystem->Exit(1);
  }
}

