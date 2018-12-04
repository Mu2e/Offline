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
  float ttime(0.0);
  float noRSDtime(0.0);
  static const string csv(".csv");
  vector<TH1F*> plots;
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
    bool first(true);
    float norm(0.0);

    std::sort(fnames.begin(),fnames.end());
    for(auto const& fname : fnames ){
      string name = fname.substr(1,fname.find_last_of(".")-1);
      string title = fname.substr(1,fname.find_last_of(".")-1) + std::string(" Execution Time;Time (msec)");
      TH1F* plot = new TH1F(name.c_str(),title.c_str(),500,0.0,maxtime);
      PlotTime(dirname,fname.c_str(),plot);
      plots.push_back(plot);
      if(first){
	first = false;
	norm = plot->GetEntries();
      }
      float mtime = plot->GetMean()*plot->GetEntries()/norm;
      cout << "module " << name << " time = " << mtime << endl;
      ttime += mtime;
      if(name != "RSD")noRSDtime += plot->GetMean()*plot->GetEntries()/norm;
    }
    TCanvas* atcan = new TCanvas("atcan","times",800,800);
    int nxcel = (int)ceil(sqrt(plots.size()));
    int nycel = (int)ceil(plots.size()/float(nxcel));
    atcan->Divide(nxcel,nycel);
    for(unsigned iplot=0;iplot<plots.size();++iplot){
      atcan->cd(iplot+1);
      plots[iplot]->Draw();
    }
    cout << "Total time = " << ttime << endl;
    cout << "No RSD time = " << noRSDtime << endl;
  }
}
