int femto_nuclear_pair_pid_sparse(const char*,const char*);
#ifndef __CINT__
#include "StMaker/common/FemtoNuclearPairPidSparse.h"
#include "yaml-cpp/yaml.h"
#include "TH2D.h"
#include "TAxis.h"
#include "THnSparse.h"
#include "TFile.h"
#include "TNamed.h"
#include "TSystem.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace PairPidSparseTest {
unsigned long long checks=0;
void Need(bool ok,const std::string& reason){++checks;if(!ok)throw std::runtime_error(reason);}
struct Row {double k,r,ns;int cent,repr,mode,mass,stage;};
std::string Resolve(const std::string& base,const std::string& ref){return !ref.empty()&&ref[0]=='/'?ref:base+ref;}
void Axes(THnSparseD& h){
  const int bins[]={200,120,120,9,3,2,4,2};
  const double lo[]={0,0,-12,-.5,-.5,-.5,-.5,-.5},hi[]={1,6,12,8.5,2.5,1.5,3.5,1.5};
  Need(h.GetNdimensions()==8,"joint dimension contract");
  for(int d=0;d<8;++d)Need(h.GetAxis(d)->GetNbins()==bins[d]&&h.GetAxis(d)->GetXmin()==lo[d]&&h.GetAxis(d)->GetXmax()==hi[d],"joint axis contract");
}
void Select(THnSparseD& h,int cent,int repr,int mode,int stage,int massLo,int massHi){
  for(int d=0;d<3;++d)h.GetAxis(d)->SetRange(0,h.GetAxis(d)->GetNbins()+1);
  const int value[]={cent,repr,mode,0,stage};
  for(int d=3;d<8;++d)if(d!=6)h.GetAxis(d)->SetRange(value[d-3]+1,value[d-3]+1);
  h.GetAxis(6)->SetRange(massLo+1,massHi+1);
}
void Same(const TH2& a,const TH2& b,const std::string& name){
  Need(a.GetNbinsX()==b.GetNbinsX()&&a.GetNbinsY()==b.GetNbinsY(),name+" shape");
  for(int x=0;x<=a.GetNbinsX()+1;++x)for(int y=0;y<=a.GetNbinsY()+1;++y){
    const int bin=a.GetBin(x,y);
    Need(a.GetBinContent(bin)==b.GetBinContent(bin),name+" flow-inclusive content");
    Need(std::fabs(std::pow(a.GetBinError(bin),2)-std::pow(b.GetBinError(bin),2))<1.e-9,name+" variance");
  }
}
}

int femto_nuclear_pair_pid_sparse(const char* mainconf,const char* outputDir){
  using namespace PairPidSparseTest;
  try{
    checks=0;Need(mainconf&&*mainconf&&outputDir&&*outputDir,"explicit paths required");
    std::string mainPath=mainconf;if(mainPath[0]!='/')mainPath=std::string(gSystem->WorkingDirectory())+"/"+mainPath;
    const size_t at=mainPath.rfind("/config/");Need(at!=std::string::npos,"mainconf below config");
    const YAML::Node main=YAML::LoadFile(mainPath);Need(bool(main["nuclearPidPairQaHist"]),"explicit QA histogram reference");
    const std::string yaml=Resolve(mainPath.substr(0,at)+"/config/",main["nuclearPidPairQaHist"].as<std::string>());
    Need(gSystem->AccessPathName(outputDir),"refuse existing output directory");Need(gSystem->mkdir(outputDir,true)==0,"make test output");
    {FemtoNuclearPairPidSparse invalid;std::ostringstream diag;Need(!invalid.Init("",diag),"empty configuration rejected");}
    FemtoNuclearPairPidSparse helper;Need(helper.Init(yaml.c_str(),std::cerr),"sparse Init");
    THnSparseD* h=helper.Histogram();Need(h!=0,"sparse allocated");Axes(*h);
    const double physical[4][3]={{-.03,-.10,-13.},{.217,.93,-1.31},{.999,5.99,11.9},{1.07,6.10,13.}};
    std::vector<Row> rows;
    for(int c=0;c<9;++c)for(int r=0;r<3;++r)for(int m=0;m<2;++m)for(int q=0;q<4;++q)for(int s=0;s<2;++s)for(int v=0;v<4;++v){
      Row row={physical[v][0],physical[v][1],physical[v][2],c,r,m,q,s};rows.push_back(row);
      helper.Fill(row.k,row.r,row.ns,c,r,m!=0,q,s!=0);
    }
    Need(h->GetEntries()==rows.size(),"each exclusive joint tuple filled once");
    const double before=h->GetEntries(),nan=std::numeric_limits<double>::quiet_NaN();
    helper.Fill(nan,1.,0.,0,0,false,0,false);helper.Fill(.2,nan,0.,0,0,false,0,false);helper.Fill(.2,1.,nan,0,0,false,0,false);
    helper.Fill(.2,1.,0.,-1,0,false,0,false);helper.Fill(.2,1.,0.,9,0,false,0,false);
    helper.Fill(.2,1.,0.,0,-1,false,0,false);helper.Fill(.2,1.,0.,0,3,false,0,false);
    helper.Fill(.2,1.,0.,0,0,false,-1,false);helper.Fill(.2,1.,0.,0,0,false,4,false);
    Need(h->GetEntries()==before,"invalid category/nonfinite never filled");
    unsigned projections=0;
    for(int c=0;c<9;++c)for(int r=0;r<3;++r)for(int m=0;m<2;++m)for(int s=0;s<2;++s)for(int region=0;region<5;++region){
      const int low[]={0,1,2,3,1},high[]={0,1,2,3,2};
      Select(*h,c,r,m,s,low[region],high[region]);
      for(int axis=1;axis<=2;++axis){
        // A preserves the original target axes even when their flow bins are
        // explicitly selected. Without A ROOT5 expands the projected axes.
        TH2D* actual=h->Projection(axis,0,"EA");Need(actual!=0,"joint 2D projection");actual->SetDirectory(0);
        TH2D expected("expected_pair_pid_sparse","",200,0.,1.,120,axis==1?0.:-12.,axis==1?6.:12.);expected.SetDirectory(0);expected.Sumw2();
        for(unsigned i=0;i<rows.size();++i){const Row& row=rows[i];if(row.cent==c&&row.repr==r&&row.mode==m&&row.stage==s&&row.mass>=low[region]&&row.mass<=high[region])expected.Fill(row.k,axis==1?row.r:row.ns);}
        Same(*actual,expected,"independent tuple oracle");
        double entries=0.;for(int x=0;x<=actual->GetNbinsX()+1;++x)for(int y=0;y<=actual->GetNbinsY()+1;++y)entries+=actual->GetBinContent(x,y);
        Need(entries==expected.GetEntries(),"flow-inclusive projected counts");delete actual;++projections;
      }
    }
    for(int d=0;d<8;++d)h->GetAxis(d)->SetRange(0,0);
    const std::string output=std::string(outputDir)+"/sparse_fixture.root";
    {TFile file(output.c_str(),"CREATE");Need(!file.IsZombie(),"output CREATE");Need(helper.Write(&file),"sparse Write success");file.Write();}
    {TFile file(output.c_str(),"READ");THnSparseD* stored=dynamic_cast<THnSparseD*>(file.Get("hNuclearPidPairJoint"));Need(stored!=0,"joint ROOT readback");Axes(*stored);Need(stored->GetEntries()==before,"joint entries persisted");Need(dynamic_cast<TNamed*>(file.Get("NuclearPidPairJointMetadata"))!=0,"joint metadata persisted");}
    std::ofstream summary((std::string(outputDir)+"/summary.json").c_str());summary<<"{\"passed\":true,\"checks\":"<<checks<<",\"projections\":"<<projections<<",\"jointRows\":"<<rows.size()<<"}\n";summary.close();Need(bool(summary),"write summary");
    std::cout<<"PAIR_PID_SPARSE_PASS checks="<<checks<<" projections="<<projections<<" rows="<<rows.size()<<std::endl;return 0;
  }catch(const std::exception& e){std::cerr<<"PAIR_PID_SPARSE_FAIL "<<e.what()<<std::endl;return 1;}
}
#endif
