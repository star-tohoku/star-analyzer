int test_femto_he4_pid_adoption(const char*,const char*,const char*,const char*);
#ifndef __CINT__
#include "test_femto_nuclear_pid_domain.C"
#include "StMaker/common/FemtoNuclearDedxCalibration.h"
#include "StMaker/common/NuclearIdDeDxVsMom.h"
namespace {
void AdoptSave(const std::string& path,const YAML::Node& node){std::ofstream o(path.c_str());o<<node;o.close();DomainRequire(o.good(),"write explicit fixture YAML");}
std::string AdoptFixture(const char* source,const std::string& output,YAML::Node nuclear){
  std::string base=gSystem->WorkingDirectory();base+="/config/";
  YAML::Node main=YAML::LoadFile(source);
  for(YAML::const_iterator i=main.begin();i!=main.end();++i){std::string key=i->first.as<std::string>(),value=i->second.as<std::string>();if(value[0]!='/')main[key]=base+value;}
  const std::string config=output+"/config";gSystem->mkdir(config.c_str(),kTRUE);
  AdoptSave(config+"/cuts.yaml",nuclear);main["nuclearid"]=config+"/cuts.yaml";AdoptSave(config+"/main.yaml",main);return config+"/main.yaml";
}
void AdoptTrack(DomainPicoFixture& f,FemtoLambdaLegacy& helper,int s,double r,bool tof,double mass2){
  int index=f.tracks.GetEntriesFast();StPicoTrack* t=new(f.tracks[index])StPicoTrack();
  t->setId(index+20);t->setNHitsFit(25);t->setNHitsDedx(20);t->setPrimaryMomentum(r,0.,0.);t->setGlobalMomentum(r,0.,0.);
  const double p=t->pMom().Mag(),a=helper.NuclearDedxPull(s,p,1.),b=helper.NuclearDedxPull(s,p,2.);
  t->setDedx((1.-a/(b-a))*1.e-6);
  if(tof){int k=f.tofs.GetEntriesFast();StPicoBTofPidTraits* v=new(f.tofs[k])StPicoBTofPidTraits();v->setTrackIndex(index);v->setBTofMatchFlag(1);v->setBeta(1./std::sqrt(1.+mass2/(p*p)));t->setBTofPidTraitsIndex(k);}
}
}
#endif
int test_femto_he4_pid_adoption(const char* mainconf,const char* calibration,const char* expected,const char* outputDirectory){
#ifndef __CINT__
 try{
  DomainChecks=0;const std::string out=std::string(gSystem->WorkingDirectory())+"/"+outputDirectory;
  DomainRequire(gSystem->AccessPathName(out.c_str())&&gSystem->mkdir(out.c_str(),kTRUE)==0,"new output directory");
  FemtoNuclearDedxCalibration c;DomainRequire(c.Load(calibration,std::cerr),"accepted polynomial loads");
  FemtoLambdaLegacy helper;DomainRequire(helper.Init(mainconf,std::cerr),"actual production mainconf initializes");
  std::ifstream reference(expected);DomainRequire(bool(reference),"independent coefficient reference readable");std::string line;std::getline(reference,line);
  int s;double r,mean,sigma;unsigned rows=0;
  while(reference>>s>>r>>mean>>sigma){++rows;
    DomainRequire(std::fabs(c.Mean(s,r)-mean)<1.e-11&&std::fabs(c.Sigma(s,r)-sigma)<1.e-12,"accepted-study mean/sigma agreement");
    DomainRequire(std::fabs(helper.NuclearDedxPull(s,r,mean+.3*sigma)-.3)<1.e-10,"production pulls use accepted parameters");
    if(s<3){DomainRequire(!c.UsesCalibration(s,r),"d/t/He3 not overridden");DomainRequire(c.Mean(s,r)==NuclearIdDeDxVsMom::GetMean(static_cast<NuclearIdDeDxVsMom::SpeciesIndex>(s),r)&&c.Sigma(s,r)==NuclearIdDeDxVsMom::GetSigma(static_cast<NuclearIdDeDxVsMom::SpeciesIndex>(s),r),"bit-exact Git mean/sigma for d/t/He3");}
  }
  DomainRequire(reference.eof()&&rows==4004,"all independent coefficient grid rows consumed");
  const double mass[]={3.48096,7.76906,1.92385,3.48096};
  const double lo[]={3.268773,7.2564275,1.78285,3.268773},hi[]={3.693147,8.2816925,2.06485,3.693147};
  const double inf=std::numeric_limits<double>::infinity();
  for(int species=0;species<4;++species){double pulls[]={999,999,999,999};pulls[species]=.1;
    DomainRequire(helper.SelectNuclearType(1.6,pulls,false,0)==species,"all species low R no TOF");
    DomainRequire(helper.PassNuclearPidQaDomain(species,1.6,false,0),"pair QA low R no TOF");
    const double probes[]={2.,std::nextafter(2.,inf),2.2};
    for(unsigned i=0;i<3;++i){double p=probes[i];DomainRequire(helper.SelectNuclearType(p,pulls,false,mass[species])==-1,"all species high R rejects missing TOF");DomainRequire(helper.SelectNuclearType(p,pulls,true,mass[species])==species,"all species high R accepts correct TOF");DomainRequire(!helper.PassNuclearPidQaDomain(species,p,false,mass[species])&&helper.PassNuclearPidQaDomain(species,p,true,mass[species]),"pair QA shares all-species TOF gate");
      const double edges[]={lo[species],hi[species]};for(int j=0;j<2;++j)DomainRequire(helper.SelectNuclearType(p,pulls,true,edges[j])==-1,"strict mass2 boundary");
      DomainRequire(helper.SelectNuclearType(p,pulls,true,std::nextafter(lo[species],hi[species]))==species,"mass2 immediately inside accepted");
    }
    // Test the TOF gate itself: production d has a separate historical veto just below 2.
    DomainRequire(helper.PassNuclearHighRigidityTof(species,std::nextafter(2.,0.),false,0),"immediately below 2 TOF not required for any species");
  }
  double h4[]={999,999,999,.1};
  DomainRequire(helper.SelectNuclearType(.61,h4,false,0)==-1&&!helper.PassNuclearPidQaDomain(3,.61,false,0),"He4 minimum strictly exclusive in selection and QA");
  DomainRequire(helper.SelectNuclearType(std::nextafter(.61,inf),h4,false,0)==3,"immediately above He4 minimum accepted");
  const double oldVeto[]={.9,1.775,3.3,4.5};for(unsigned i=0;i<4;++i)DomainRequire(helper.SelectNuclearType(oldVeto[i],h4,true,mass[3])==3,"accepted study He4 domain has no old veto");
  DomainRequire(helper.SelectNuclearType(5.,h4,true,mass[3])==-1,"maximum remains strict");
  const double competing[]={.2,999,999,.1};DomainRequire(helper.SelectNuclearType(.60,competing,false,0)==-1,"failed He4 winner never reassigned");
  YAML::Node doc=YAML::LoadFile(calibration);
  for(int variant=0;variant<7;++variant){YAML::Node bad=YAML::Clone(doc);
    if(variant==0)bad["species"]["he4"].remove("mean");
    if(variant==1)bad["species"]["he4"]["sigma"]=YAML::Load("[1, 2]");
    if(variant==2)bad["species"]["he4"]["mean"][0]=std::numeric_limits<double>::quiet_NaN();
    if(variant==3)bad["species"]["he4"]["sigma"]=YAML::Load("[-1, 0, 0]");
    if(variant==4)bad["species"]["he3"]=bad["species"]["he4"];
    if(variant==5)bad["source"]["year"]=2019;
    // Positive endpoints but negative at an interior stationary point.
    if(variant==6)bad["species"]["he4"]["sigma"]=YAML::Load("[-0.01, 0, 1]");
    std::string path=out+"/bad"+TString::Format("%d",variant).Data()+".yaml";AdoptSave(path,bad);FemtoNuclearDedxCalibration reject;DomainRequire(!reject.Load(path,std::cerr),"invalid polynomial fails loudly");
  }
  {FemtoNuclearDedxCalibration reject;DomainRequire(!reject.Load(out+"/missing.yaml",std::cerr),"explicit missing calibration fails");}
  YAML::Node main=YAML::LoadFile(mainconf);std::string cuts=main["nuclearid"].as<std::string>();YAML::Node nuclear=YAML::LoadFile("config/"+cuts);
  nuclear["dedxCalibration"]["calibrationFile"]=out+"/missing.yaml";std::string fixture=AdoptFixture(mainconf,out,nuclear);{FemtoLambdaLegacy reject;std::ostringstream error;DomainRequire(!reject.Init(fixture.c_str(),error)&&error.str().find("Cannot read calibration")!=std::string::npos,"missing explicit calibration makes provider Init fail for intended reason");}
  nuclear["dedxCalibration"]["enabled"]=false;nuclear.remove("speciesMinPOverQ");fixture=AdoptFixture(mainconf,out,nuclear);{FemtoLambdaLegacy reject;std::ostringstream error;DomainRequire(!reject.Init(fixture.c_str(),error)&&error.str().find("requires speciesMinPOverQ")!=std::string::npos,"exclusive policy without minima makes Init fail for intended reason");}
  {DomainPicoFixture f;
   AdoptTrack(f,helper,3,.60,false,0);AdoptTrack(f,helper,3,.62,false,0);AdoptTrack(f,helper,3,.90,false,0);AdoptTrack(f,helper,3,1.775,false,0);AdoptTrack(f,helper,3,3.3,false,mass[3]);AdoptTrack(f,helper,3,3.3,true,mass[3]);
   AdoptTrack(f,helper,2,1.6,false,0);AdoptTrack(f,helper,2,2.2,false,mass[2]);AdoptTrack(f,helper,2,2.2,true,mass[2]);AdoptTrack(f,helper,2,2.2,true,mass[3]);
   FemtoCandidateStore store;FemtoNuclearPidPairQaStore qa;helper.CollectNuclei(&f.dst,0,8,store,&qa);
   DomainRequire(store["he4"].size()==4,"He4 collection applies new minimum, no veto and TOF");DomainRequire(store["he3"].size()==2,"physical He3 collection low noTOF/high species TOF only");
   for(size_t k=0;k<store["he3"].size();++k){const FemtoCandidate& n=store["he3"][k];DomainRequire(helper.IsPairNucleus(2),"He3 remains enabled pair nucleus");FemtoCandidate l;l.eventIndex=0;l.source=kFemtoCandResonance;l.speciesKey="lambda";l.px=.3;l.energy=1.2;l.reso.invMass=1.1157;l.reso.dau1EventIndex=l.reso.dau2EventIndex=0;l.reso.dau1Index=100;l.reso.dau2Index=101;DomainRequire(helper.FillPair(l,n,.1,.2,8,false,&f.dst),"accepted He3 fills physical and mass-hypothesis pair spectra");}
   DomainRequire(helper.RootHistograms()->Get("hKstarSE_lambda_he3")->GetEntries()==2,"physical He3 kstar has only accepted candidates");DomainRequire(helper.RootHistograms()->Get("hKstarSE_lambda_he3_as_he4")->GetEntries()==2,"He3-as-He4 background kstar has same accepted candidates");
   for(size_t k=0;k<store["he4"].size();++k){const FemtoCandidate& n=store["he4"][k];const StPicoTrack* t=f.dst.track(n.trk.trackIndex);DomainRequire(std::fabs(n.trk.nSigmaHe4-helper.NuclearDedxPull(3,t->pMom().Mag(),t->dEdx()))<1.e-5,"stored He4 pull uses new coefficients");}
  }
  {TFile file((out+"/adopted_profile.root").c_str(),"CREATE");DomainRequire(!file.IsZombie()&&helper.Write(&file),"write adopted metadata");}
  {TFile file((out+"/adopted_profile.root").c_str(),"READ");TNamed* d=dynamic_cast<TNamed*>(file.Get("FemtoNuclearDedxCalibration"));TNamed* gates=dynamic_cast<TNamed*>(file.Get("FemtoNuclearPidDomain"));DomainRequire(d&&std::string(d->GetTitle()).find("explicit_yaml_log10_polynomial")!=std::string::npos,"ROOT metadata records polynomial");DomainRequire(gates&&std::string(gates->GetTitle()).find("exclusiveMinimum=1")!=std::string::npos&&std::string(gates->GetTitle()).find("speciesPOverQVeto.he4: []")!=std::string::npos,"ROOT metadata records strict minimum and empty He4 veto");}
  std::cout<<"PASS He4 production PID adoption: "<<DomainChecks<<" checks"<<std::endl;return 0;
 }catch(const std::exception& e){std::cerr<<"FAIL He4 PID adoption after "<<DomainChecks<<" checks: "<<e.what()<<std::endl;return 1;}
#else
 return 1;
#endif
}
