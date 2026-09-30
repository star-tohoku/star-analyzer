// Compile with ACLiC after loading core libraries only. No concrete KF provider
// is included or linked: this verifies the optional-dependency boundary.
#include "StMaker/StFemtoMaker/StFemtoMaker.h"
#include "FemtoLambdaProvider.h"
#include "ConfigManager.h"
#include "YamlParser.h"
#include "StChain.h"
#include "TSystem.h"
#include "TFile.h"
#include "TNamed.h"
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <unistd.h>

namespace femto_lambda_maker_test {
void Require(bool ok, const std::string& why) {
  if (!ok) throw std::runtime_error(why);
}
std::string Read(const std::string& path) {
  std::ifstream input(path.c_str());
  Require(bool(input), "cannot read fixture source " + path);
  std::ostringstream text; text << input.rdbuf(); return text.str();
}
std::string Scalar(const std::string& text, const std::string& key, const std::string& value, bool remove=false) {
  std::istringstream input(text); std::ostringstream output; std::string line;
  bool found=false;
  while (std::getline(input,line)) {
    if (line.compare(0,key.size()+1,key+":")==0) {
      Require(!found, "duplicate fixture key "+key); found=true;
      if (!remove) output << key << ": " << value << '\n';
    } else output << line << '\n';
  }
  Require(found,"fixture key absent "+key);
  return output.str();
}
struct Fixtures {
  std::string dir, maker, mixing;
  std::map<std::string,std::string> references;
  std::vector<std::string> files;
  explicit Fixtures(const char* source) {
    char pattern[]="/tmp/star_femto_maker_checks_XXXXXX";
    char* created=mkdtemp(pattern);
    Require(created!=0,"mkdtemp failed");
    dir=created;
    Require(gSystem->mkdir((dir+"/config/mainconf").c_str(),true)==0,"fixture directory creation failed");
    std::string main=source ? source : "";
    Require(!main.empty(),"mainconf argument is required");
    if(main[0]!='/') main=std::string(gSystem->WorkingDirectory())+"/"+main;
    const size_t at=main.rfind("/config/");
    Require(at!=std::string::npos,"source mainconf must be below config");
    Require(YamlParser::ParseFile(main.c_str(),references),"source mainconf cannot be parsed");
    for(std::map<std::string,std::string>::iterator i=references.begin();i!=references.end();++i)
      if(!i->second.empty() && i->second[0]!='/') i->second=main.substr(0,at)+"/config/"+i->second;
    maker=Read(references["maker"]); mixing=Read(references["mixing"]);
    // ConfigManager currently accepts config-relative references only. Preserve
    // an isolated copied subtree instead of teaching the test a different resolver.
    for(std::map<std::string,std::string>::iterator i=references.begin();i!=references.end();++i) {
      const std::string name="base_"+i->first+".yaml";
      Write("config/"+name,Read(i->second));
      i->second=name;
    }
  }
  ~Fixtures() {
    for(size_t i=0;i<files.size();++i) unlink(files[i].c_str());
    rmdir((dir+"/config/mainconf").c_str()); rmdir((dir+"/config").c_str()); rmdir(dir.c_str());
  }
  std::string Write(const std::string& name,const std::string& contents) {
    const std::string path=dir+"/"+name;
    std::ofstream out(path.c_str()); out << contents; out.close();
    Require(bool(out),"cannot write fixture "+path); files.push_back(path); return path;
  }
  std::string Main(const std::string& name,const std::string& makerText="",const std::string& mixText="",
                   const std::string& histText="",const std::string& removeReference="") {
    std::map<std::string,std::string> refs=references;
    if(!makerText.empty()) { Write("config/"+name+"_maker.yaml",makerText); refs["maker"]=name+"_maker.yaml"; }
    if(!mixText.empty()) { Write("config/"+name+"_mix.yaml",mixText); refs["mixing"]=name+"_mix.yaml"; }
    if(!histText.empty()) { Write("config/"+name+"_hist.yaml",histText); refs["femtoHist"]=name+"_hist.yaml"; }
    if(!removeReference.empty()) refs.erase(removeReference);
    std::ostringstream main;
    for(std::map<std::string,std::string>::const_iterator i=refs.begin();i!=refs.end();++i)
      main << i->first << ": " << i->second << '\n';
    return Write("config/mainconf/"+name+".yaml",main.str());
  }
  std::string Output(const std::string& name) {
    const std::string path=dir+"/"+name+".root"; files.push_back(path); return path;
  }
};
class MockProvider : public FemtoLambdaProvider {
public:
  static int destroyed;
  explicit MockProvider(bool result=true):initResult(result) {}
  virtual ~MockProvider(){++destroyed;}
  virtual bool Init(const char*,std::ostream&){return initResult;}
  virtual bool AcceptEvent(const StPicoEvent&,int)const{return true;}
  virtual double ComputeVr(double x,double y)const{return std::sqrt(x*x+y*y);}
  virtual bool Process(StPicoDst*,int){return true;}
  virtual void Clear(){}
  virtual const std::vector<FemtoLambdaCandidate>& Candidates()const{return candidates;}
  virtual const FemtoLambdaEventStats& Stats()const{return stats;}
  virtual const std::string& LastError()const{return error;}
private:
  bool initResult;
  std::vector<FemtoLambdaCandidate> candidates;
  FemtoLambdaEventStats stats;
  std::string error;
};
int MockProvider::destroyed=0;

void Case(Fixtures& f,const std::string& name,const std::string& main,
          bool expected,bool provider=true,bool providerInit=true,bool existingOutput=false) {
  const std::string output=f.Output(name);
  if(existingOutput) { std::ofstream preserved(output.c_str()); preserved << "do-not-overwrite"; }
  if (!ConfigManager::GetInstance().LoadConfig(main.c_str())) {
    Require(!expected,"valid fixture configuration failed to load: "+name);
    std::cout << "PASS configuration rejection " << name << std::endl;
    return;
  }
  const int before=MockProvider::destroyed;
  StChain* chain=new StChain;
  StFemtoMaker* maker=new StFemtoMaker("femto",0,output.c_str());
  maker->SetMainConfigPath(main.c_str());
  if(provider) maker->SetLambdaProvider(new MockProvider(providerInit));
  const int result=maker->Init();
  const bool accepted=result==kStOK;
  delete chain;
  Require(accepted==expected, name+(expected?": valid Init rejected":": invalid Init accepted"));
  Require(MockProvider::destroyed-before==(provider?1:0),name+": provider ownership leak/double-delete");
  if(existingOutput) Require(Read(output)=="do-not-overwrite","existing output was changed");
  std::cout << "PASS maker Init " << name << " status=" << result << std::endl;
}

void NoKfLibraries() {
  const std::string libs=gSystem->GetLibraries();
  Require(libs.find("libStKfParticleCommon.so")==std::string::npos &&
          libs.find("libKFParticle.so")==std::string::npos,
          "test process loaded optional KF libraries");
}
}

int femto_lambda_maker_checks(const char* mainconf) {
  using namespace femto_lambda_maker_test;
  try {
    NoKfLibraries();
    Fixtures f(mainconf);
    Case(f,"valid_mock",f.Main("valid_mock"),true);
    Case(f,"provider_absent",f.Main("provider_absent"),false,false);
    Case(f,"provider_init_failure",f.Main("provider_init_failure"),false,true,false);
    Case(f,"existing_output",f.Main("existing_output"),false,true,true,true);
    Case(f,"unsupported_builder",f.Main("unsupported_builder",
         Scalar(f.maker,"species_lambda_builderType","track")),false);
    Case(f,"missing_builder",f.Main("missing_builder",
         Scalar(f.maker,"species_lambda_particleKey","",true)),false);
    Case(f,"missing_channel_field",f.Main("missing_channel_field",
         Scalar(f.maker,"channel_1_signalMax","",true)),false);
    Case(f,"duplicate_channel",f.Main("duplicate_channel",
         Scalar(f.maker,"channel_1_name","lambda_deuteron")),false);
    Case(f,"unsupported_closepair",f.Main("unsupported_closepair",
         Scalar(f.maker,"closePairEnabled","true")),false);
    Case(f,"missing_kf_reference",f.Main("missing_kf_reference","","","","kf"),false);
    const char* mixKeys[]={"nVzBins","nCentralityBins","nEventPlaneBins","bufferSize",
                          "mixingMode","mixBothDirections","minVz","maxVz","maxMixEvents","vzOutOfRangePolicy"};
    for(size_t i=0;i<sizeof(mixKeys)/sizeof(mixKeys[0]);++i) {
      const std::string name=std::string("missing_mix_")+mixKeys[i];
      Case(f,name,f.Main(name,"",Scalar(f.mixing,mixKeys[i],"",true)),false);
    }
    Case(f,"malformed_mix_bins",f.Main("malformed_mix_bins","",Scalar(f.mixing,"nVzBins","not_a_number")),false);
    Case(f,"malformed_mix_bool",f.Main("malformed_mix_bool","",Scalar(f.mixing,"mixBothDirections","not_a_boolean")),false);
    Case(f,"wrong_cent_bins",f.Main("wrong_cent_bins","",Scalar(f.mixing,"nCentralityBins","8")),false);
    Case(f,"unsupported_sampler",f.Main("unsupported_sampler","",Scalar(f.mixing,"mixingMode","randomSample")),false);
    Case(f,"incomplete_hist",f.Main("incomplete_hist","","",
         "histograms:\n  hDeliberatelyIncomplete:\n    title: test\n    nBins: 1\n    min: 0\n    max: 1\n"),false);
    NoKfLibraries();
    std::cout << "PASS core-only StFemtoMaker load + Lambda Init failure/ownership cases" << std::endl;
    return 0;
  } catch(const std::exception& error) {
    std::cerr << "FAIL maker checks: " << error.what() << std::endl; return 1;
  }
}
