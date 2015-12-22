#include "TFile.h"
#include "TTree.h"
#include "TH2.h"

#include <sstream>

void OutOfTargetMuonStops(std::string filename) {
  
  TFile* file = new TFile(filename.c_str(), "READ");
  TTree* stops = (TTree*) file->Get("outOfTargetDumper/stops");

  double bin_width = 10;
  double central_line = -3904;
  double max_r = 355;
  double min_x = central_line-max_r;
  double max_x = central_line+max_r;

  int n_bins_x = (max_x - min_x) / bin_width;
  double min_y = -max_r;
  double max_y = +max_r;
  int n_bins_y = (max_y - min_y) / bin_width;
  double min_z = 6300;
  double max_z = 7500;
  int n_bins_z = (max_z - min_z) / bin_width;
  double bin_width_r = 1;
  double min_r = 330;
  int n_bins_r = (max_r - min_r) / bin_width_r;
  TH2F* hStopsRZ = new TH2F("hStopsRZ", "", n_bins_z,min_z,max_z, n_bins_r,min_r,max_r);
  TH2F* hStopsYX = new TH2F("hStopsYX", "", n_bins_x,min_x,max_x, n_bins_y,min_y,max_y);

  std::stringstream geomcut;
  geomcut << "y<" << max_y << " && y>" << min_y << " && x<" << max_x << " && x>" << min_x << " && z>" << min_z << " && z<" << max_z;
  stops->Draw("sqrt(y*y + (x+3904)*(x+3904)):z>>hStopsRZ", geomcut.str().c_str(), "COLZ");
  stops->Draw("y:x>>hStopsYX", geomcut.str().c_str(), "COLZ");
}
