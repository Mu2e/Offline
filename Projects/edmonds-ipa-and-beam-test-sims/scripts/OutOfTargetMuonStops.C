#include "TFile.h"
#include "TTree.h"
#include "TH2.h"

void OutOfTargetMuonStops(std::string filename) {
  
  TFile* file = new TFile(filename.c_str(), "READ");
  TTree* stops = (TTree*) file->Get("outOfTargetDumper/stops");

  double bin_width = 10;
  double min_x = -3000;
  double max_x = -5000;
  int n_bins_x = (max_x - min_x) / bin_width;
  double min_y = -500;
  double max_y = 500;
  int n_bins_y = (max_y - min_y) / bin_width;
  double min_z = 6300;
  double max_z = 7500;
  int n_bins_z = (max_z - min_z) / bin_width;
  double bin_width_r = 1;
  double min_r = 335;
  double max_r = 380;
  int n_bins_r = (max_r - min_r) / bin_width_r;
  TH2F* hStopsRZ = new TH2F("hStopsRZ", "", n_bins_z,min_z,max_z, n_bins_r,min_r,max_r);

  stops->Draw("sqrt(y*y + (x+3904)*(x+3904)):z>>hStopsRZ", "y<500 && y>-500 && x<-3000 && x>-5000 && z>6300 && z<7500", "COLZ");
}
