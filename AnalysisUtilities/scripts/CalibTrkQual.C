#include <iostream>
#include <fstream>
#include <math.h>

#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TCut.h"

using namespace std;

double EfficiencyToTrkQual(TTree* inpt_tree, const char* leaf_name, double effic, TCut signal_cut, int subdivs=1000000);
std::string CalibTrkQual_Eff(TTree* tree, std::string train_name, std::string leaf_name, TCut signal_cut);

void CalibTrkQual(TTree* tree, std::string train_name, std::string leaf_name, TCut signal_cut) {

  std::cout << "Calibrating " << train_name << "..." << std::endl;

  std::string xmlfilename = "AnalysisConditions/weights/" + train_name + ".weights.xml";
  std::ifstream xmlfile(xmlfilename, std::ofstream::in);

  std::string outfilename = "AnalysisConditions/weights/" + train_name + ".weights.xml_calibrating";
  std::ofstream outfile(outfilename, std::ofstream::out);

  const int n_chars = 256;
  char xmlline[n_chars];
  if (!xmlfile.good()) {
    std::cout << "Problem with xml file " << xmlfilename << std::endl;
    return;
  }
  while (xmlfile.good()) {
    xmlfile.getline(xmlline, n_chars);

    std::string xml_string(xmlline);

    if (xml_string.find("<Calibration>") != std::string::npos) {
      std::cout << "Calibration block already exists! Skipping to next weights file..." << std::endl;
      return;
    }

    if (xml_string.find("</MethodSetup>") != std::string::npos) {
      // Add the calibration block just before the end
      outfile << "  <Calibration>" << std::endl;

      const auto& calibtable_eff = CalibTrkQual_Eff(tree, train_name, leaf_name, signal_cut);
      outfile << calibtable_eff;

      // TODO: if we ever have a definition of background rejection for TrkQual, we can calibrate it here
      //      const auto& calibtable_bkg_rej = CalibTrkQual_BkgRej(tree, train_name);
      //      outfile << calibtable_bkg_rej;

      outfile << "  </Calibration>" << std::endl;
    }    
    outfile << xml_string << std::endl;
  }
  std::rename(outfilename.c_str(), xmlfilename.c_str());
  outfile.close();
  std::cout << "Done!" << std::endl;
}

std::string CalibTrkQual_Eff(TTree* tree, std::string train_name, std::string leaf_name, TCut signal_cut) {

  double min_eff = 0.0;
  double max_eff = 1.0;
  double eff_step = 0.01;
  std::stringstream calibtable;
  calibtable.str("");

  int counter = 0;
  for (double i_eff = min_eff; i_eff <= max_eff; i_eff += eff_step) {
    double trkqual_cut = EfficiencyToTrkQual(tree, leaf_name.c_str(), i_eff, signal_cut);
    calibtable << "    <Calib Index=\"" << counter << "\" CalibVal=\"" << 1-i_eff << "\" Val=\"" << trkqual_cut << "\"/>" << std::endl;
    ++counter;
  }
  calibtable << "    <Calib Index=\"" << counter << "\" CalibVal=\"0\" Val=\"0.0\"/>" << std::endl;
  return calibtable.str();
}

double EfficiencyToTrkQual(TTree* inpt_tree, const char* leaf_name, double effic, TCut signal_cut, int subdivs) {
  // Given a target efficiency and a TrkAna tree with TrkQual training, function returns a TrkQual cut value
  // that approximately achieves the target efficiency
  //
  // Returned cut value's significant figures grow roughly as log10(subdivs)

  string leaf_str(leaf_name);
  string leaf_hist_str = "dequal."+leaf_str+">>hist("+to_string(subdivs)+",0,1)";
  TCut std_cuts = signal_cut;

  inpt_tree->Draw(leaf_hist_str.c_str(), std_cuts, "goff");
  TH1F *hist = (TH1F*) gDirectory->Get("hist");
  TH1F *cumhist = (TH1F*)hist->GetCumulative();
  double target = (hist->GetEntries())*(1-effic);
  float N = floor(log2(subdivs));
 
  int i = 0;
  int ind = subdivs/2;
  int step = subdivs/4;
  double cur_val = 0.;
  while (i < N) {
    cur_val = cumhist->GetBinContent(ind);

    if (cur_val < target) {
      ind += step;
    } else {
      ind -= step;
    }
    step = step/2;
    i += 1;
  }

  double approx_trkqual = ((double) ind) / subdivs;
  // if (approx_trkqual < 10.0/subdivs) {
  //   approx_trkqual = 0;
  // }
  //  string add_on = " && dequal.";
  //  //  cout << "Relative Efficiency " << inpt_tree->GetEntries((std_cuts+add_on+train_str+" > "+to_string(approx_trkqual)).c_str())/hist->GetEntries() << "\n" << endl;
  return approx_trkqual;  
}
