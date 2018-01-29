#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TSystem.h"
//#include <sys/types.h>
//#include <dirent.h>
//#include <libgen.h>
#include <string>
#include <iostream>
#include <algorithm>
using namespace std;

void PlotTime(const char* dirname, const char* name, TH1F* timeplot) {
  TTree* rdtree = new TTree(name,name);
  string fullname = std::string(dirname) + std::string(name);
  rdtree->ReadFile(fullname.c_str(),"evt/I:time/F:mod/C",',');
  rdtree->Project(timeplot->GetName(),"1000.0*time");
}


void PlotAllTimes(const char* dirname, double maxtime) {
  // find all the csv files in this directory
  static const string csv(".csv");
  vector<TH1F*> plots;
//  auto dir`p = opendir(dirname);
  void* dirp = gSystem->OpenDirectory(dirname);
  std::vector<string> fnames;
  if(dirp != 0){
    const char* dent(0);
    while ( (dent = gSystem->GetDirEntry(dirp)) != 0){
      string fname(dent);
      // compare the suffix
      if(fname.size() >= csv.size() && 0 == fname.compare (fname.length() - csv.length(), csv.length(), csv)){
	cout << "found entry " << fname << endl;
	fnames.push_back(fname);
      }
    }
    std::sort(fnames.begin(),fnames.end());
    for(auto const& fname : fnames ){
      string name = fname.substr(1,fname.find_last_of(".")-1);
      string title = fname.substr(1,fname.find_last_of(".")-1) + std::string(" Execution Time;Time (msec)");
      TH1F* plot = new TH1F(name.c_str(),title.c_str(),500,0.0,maxtime);
      PlotTime(dirname,fname.c_str(),plot);
      plots.push_back(plot);
    }
    TCanvas* acan = new TCanvas("alltimes","times",1200,1200);
    int nxcel = (int)ceil(sqrt(plots.size()));
    int nycel = (int)ceil(plots.size()/nxcel);
    acan->Divide(nxcel,nycel);
    for(unsigned iplot=0;iplot<plots.size();++iplot){
      acan->cd(iplot+1);
      plots[iplot]->Draw();
//      TCanvas* can = new TCanvas(plots[iplot]->GetName(),plots[iplot]->GetName(),400,400);
//      can->Divide(1,1);
//      can->cd(1);
//      plots[iplot]->Draw();
//      char rname[100];
//      snprintf(rname,100,"%s%s",plots[iplot]->GetName(),".root");
//      can->SaveAs(rname);
    }
  }
}
