#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>

Int_t epdEarly81Headers(const char* inputs,const char* output,const char* job,
                       Long64_t events,const char* mainconf) {
  if(!inputs || !output || !job || !mainconf || events!=-1 ||
     !gSystem->AccessPathName(output)) return 2;
  std::ifstream list(inputs);
  std::ofstream out(output);
  if(!list || !out) return 2;
  std::string path;
  int files=0;
  while(std::getline(list,path)) {
    std::istringstream fields(path);
    std::string file,extra;
    if(!(fields>>file) || (fields>>extra) || file[0]!='/') return 2;
    std::unique_ptr<TFile> root(TFile::Open(file.c_str(),"READ"));
    Long64_t entries=-1;
    const char* state="unreadable";
    if(root.get() && !root->IsZombie() && !root->TestBit(TFile::kRecovered)) {
      TTree* pico=dynamic_cast<TTree*>(root->Get("PicoDst"));
      if(pico && pico->GetBranch("Event") && pico->GetBranch("Track") && pico->GetBranch("EpdHit")) {
        entries=pico->GetEntries();
        state=entries>0?"has_events":"empty_tree";
      } else state="missing_tree_or_branches";
    }
    out<<file.substr(file.find_last_of('/')+1)<<'\t'<<entries<<'\t'<<state
       <<'\t'<<(root.get()?root->GetBytesRead():0)<<std::endl;
    if(!out.good()) return 2;
    ++files;
  }
  if(!list.eof() || !files) return 2;
  out.close();
  if(!out.good()) return 2;
  std::cout<<"EARLY_PICO_HEADERS_PASS files="<<files<<" job="<<job<<std::endl;
  return 0;
}
