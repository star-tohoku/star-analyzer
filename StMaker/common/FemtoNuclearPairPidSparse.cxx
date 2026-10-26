#include "FemtoNuclearPairPidSparse.h"
#include "TAxis.h"
#include "TDirectory.h"
#include "TNamed.h"
#include "yaml-cpp/yaml.h"
#include <cmath>
#include <fstream>
#include <iterator>
#include <ostream>
#include <sstream>
#include <stdexcept>

namespace {
template<class T> T Required(const YAML::Node& node,const char* key) {
  if(!node[key]) throw std::runtime_error(std::string("Missing jointSparse key: ")+key);
  return node[key].as<T>();
}
void Known(const YAML::Node& node,const char* const* keys,unsigned size,const char* label) {
  if(!node || !node.IsMap()) throw std::runtime_error(std::string(label)+" must be a map");
  for(YAML::const_iterator i=node.begin();i!=node.end();++i) {
    const std::string key=i->first.as<std::string>();
    bool found=false;
    for(unsigned j=0;j<size;++j) found=found || key==keys[j];
    if(!found) throw std::runtime_error(std::string("Unknown ")+label+" key: "+key);
  }
}
}

FemtoNuclearPairPidSparse::FemtoNuclearPairPidSparse()
  :mHistogram(0),mAccepted(0),mRejectedNonfinite(0),mRejectedCategory(0) {}
FemtoNuclearPairPidSparse::~FemtoNuclearPairPidSparse() { delete mHistogram; }

bool FemtoNuclearPairPidSparse::Init(const char* yamlPath,std::ostream& diagnostics) {
  delete mHistogram;mHistogram=0;mMetadata.clear();
  mAccepted=0;mRejectedNonfinite=0;mRejectedCategory=0;
  try {
    if(!yamlPath || !*yamlPath) throw std::runtime_error("Empty jointSparse config path");
    std::ifstream input(yamlPath);
    if(!input) throw std::runtime_error("Cannot read jointSparse configuration");
    const std::string contents((std::istreambuf_iterator<char>(input)),std::istreambuf_iterator<char>());
    if(input.bad() || contents.empty()) throw std::runtime_error("Empty/unreadable jointSparse configuration");
    const YAML::Node node=YAML::Load(contents)["jointSparse"];
    const char* topKeys[]={"schemaVersion","name","title","axes"};
    Known(node,topKeys,4,"jointSparse");
    if(Required<int>(node,"schemaVersion")!=1)
      throw std::runtime_error("Unsupported jointSparse schemaVersion");
    const std::string name=Required<std::string>(node,"name");
    const std::string title=Required<std::string>(node,"title");
    if(name!="hNuclearPidPairJoint" || title.empty())
      throw std::runtime_error("Invalid jointSparse name/title");
    const char* names[]={"Kstar","POverQ","NSigma","Cent9","Representation",
                         "MixingMode","MassClass","PidStage"};
    const YAML::Node axes=node["axes"];
    Known(axes,names,8,"jointSparse axes");
    Int_t bins[8];Double_t minimum[8],maximum[8];std::string titles[8];
    const char* axisKeys[]={"bins","min","max","title"};
    for(unsigned i=0;i<8;++i) {
      const YAML::Node axis=axes[names[i]];
      Known(axis,axisKeys,4,names[i]);
      bins[i]=Required<int>(axis,"bins");
      minimum[i]=Required<double>(axis,"min");maximum[i]=Required<double>(axis,"max");
      titles[i]=Required<std::string>(axis,"title");
      if(bins[i]<=0 || !std::isfinite(minimum[i]) || !std::isfinite(maximum[i]) ||
         !(maximum[i]>minimum[i]) || titles[i].empty())
        throw std::runtime_error(std::string("Invalid jointSparse axis: ")+names[i]);
    }
    // Category meanings are part of the schema, not adjustable physics cuts.
    const int categories[]={9,3,2,4,2};
    for(unsigned i=3;i<8;++i)
      if(bins[i]!=categories[i-3] || minimum[i]!=-0.5 || maximum[i]!=categories[i-3]-0.5)
        throw std::runtime_error(std::string("Invalid jointSparse category mapping: ")+names[i]);
    mHistogram=new THnSparseD(name.c_str(),title.c_str(),8,bins,minimum,maximum);
    mHistogram->Sumw2();
    for(unsigned i=0;i<8;++i) {
      mHistogram->GetAxis(i)->SetName(names[i]);
      mHistogram->GetAxis(i)->SetTitle(titles[i].c_str());
    }
    const char* repr[]={"he4","he3","he3as4"};
    const char* mode[]={"SE","ME"};
    const char* mass[]={"SBNeg","signalLow","signalHigh","SBPos"};
    const char* stage[]={"selected","prePID"};
    for(int i=0;i<3;++i) mHistogram->GetAxis(4)->SetBinLabel(i+1,repr[i]);
    for(int i=0;i<2;++i) {
      mHistogram->GetAxis(5)->SetBinLabel(i+1,mode[i]);
      mHistogram->GetAxis(7)->SetBinLabel(i+1,stage[i]);
    }
    for(int i=0;i<4;++i) mHistogram->GetAxis(6)->SetBinLabel(i+1,mass[i]);
    std::ostringstream metadata;
    metadata << "schemaVersion: 1\nresolvedConfiguration: " << yamlPath
      << "\naxisOrder: [Kstar,POverQ,NSigma,Cent9,Representation,MixingMode,MassClass,PidStage]\n"
      << "representation: {0: he4, 1: he3, 2: he3as4}\n"
      << "mixingMode: {0: SE, 1: ME}\n"
      << "exclusiveMassClass: {0: SBNeg, 1: signalLow, 2: signalHigh, 3: SBPos}\n"
      << "pidStage: {0: selected, 1: prePID}\n"
      << "signalProjection: sum(signalLow,signalHigh); do not count signal again\n"
      << "sourcePull: he4 uses nSigmaHe4; he3 and he3as4 use nSigmaHe3\n"
      << "rigidity: source physical momentum magnitude/2; unchanged for he3as4\n"
      << "weights: unit pair counts; Sumw2 enabled\n"
      << "flowPolicy: finite Kstar/POverQ/NSigma flow bins retained; projections must include flows\n"
      << "prePID: removes own TPC nSigma and best-species competition only; other species-domain and TOF requirements remain\n"
      << "mixing: caller supplies original physical pool admission and both-direction provenance; QA does not select events\n"
      << "resolution: explicit YAML binning, not original native low-R study resolution\n"
      << "----- exact consumed QA YAML -----\n" << contents;
    mMetadata=metadata.str();
    diagnostics << "[FemtoNuclearPairPidSparse] " << name << " loaded from " << yamlPath
                << "; 8 dimensions; unit counts and retained finite flow bins" << std::endl;
    return true;
  } catch(const std::exception& error) {
    delete mHistogram;mHistogram=0;mMetadata.clear();
    diagnostics << "[FemtoNuclearPairPidSparse] Init failed: " << error.what() << std::endl;
    return false;
  }
}

void FemtoNuclearPairPidSparse::Fill(double kstar,double rigidity,double pull,int cent9,
    int representation,bool mixed,int massClass,bool prePid) {
  if(!mHistogram) return;
  if(!std::isfinite(kstar) || !std::isfinite(rigidity) || !std::isfinite(pull)) {
    ++mRejectedNonfinite;return;
  }
  if(cent9<0 || cent9>8 || representation<0 || representation>2 || massClass<0 || massClass>3) {
    ++mRejectedCategory;return;
  }
  Double_t values[]={kstar,rigidity,pull,static_cast<double>(cent9),
    static_cast<double>(representation),mixed?1.:0.,static_cast<double>(massClass),prePid?1.:0.};
  mHistogram->Fill(values,1.);++mAccepted;
}

bool FemtoNuclearPairPidSparse::Write(TDirectory* directory) const {
  if(!mHistogram || !directory) return false;
  TDirectory* previous=gDirectory;
  if(!directory->cd()) return false;
  bool ok=mHistogram->Write(mHistogram->GetName(),TObject::kOverwrite)>0;
  std::ostringstream contents;
  contents << mMetadata << "\n----- QA fill counters -----\naccepted: " << mAccepted
    << "\nrejectedNonfinite: " << mRejectedNonfinite
    << "\nrejectedCategory: " << mRejectedCategory << "\n";
  TNamed metadata("NuclearPidPairJointMetadata",contents.str().c_str());
  if(metadata.Write(metadata.GetName(),TObject::kOverwrite)<=0) ok=false;
  if(previous) previous->cd();
  return ok;
}
