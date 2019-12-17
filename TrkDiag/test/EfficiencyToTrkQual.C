#include <iostream>
#include <fstream>
#include <math.h>

#include "TFile.h"
#include "TH1.h"
#include "TTree.h"

using namespace std;

double EfficiencyToTrkQual(TTree* inpt_tree, const char* train_name, double effic, int subdivs=10000);

void CalibrateTrkQualEff(TTree* tree, std::string train_name) {


  double min_eff = 0.0;
  double max_eff = 1.0;
  double eff_step = 0.1;

  std::cout << "Creating calib table for " << train_name << "..." << std::endl;
  std::string outfilename = "TrkDiag/test/" + train_name + "Calib.txt";

  std::ofstream outfile(outfilename, std::ofstream::out);
  
  outfile << "TABLE " << train_name << "Calib" << std::endl;
  outfile << "# Training file = " << tree->GetCurrentFile()->GetName() << std::endl;
  outfile << "# Training tree = " << tree->GetDirectory()->GetName() << "/" << tree->GetName() << std::endl;
  outfile << "# Training name = " << train_name << std::endl;
  
  int counter = 0;
  for (double i_eff = min_eff; i_eff <= max_eff; i_eff += eff_step) {
    double trkqual_cut = EfficiencyToTrkQual(tree, train_name.c_str(), i_eff);
    outfile << counter << ", " << i_eff << ", " << trkqual_cut << std::endl;
    ++counter;
  }
  outfile.close();

  std::cout << "Done!" << std::endl;
}

double EfficiencyToTrkQual(TTree* inpt_tree, const char* train_name, double effic, int subdivs) {
  // Given a target efficiency and a TrkAna tree with TrkQual training, function returns a TrkQual cut value
  // that approximately achieves the target efficiently
  //
  // Returned cut value's significant figures grow roughly as log10(subdivs)

  string train_str(train_name);
  string leaf_hist_str = "dequal."+train_str+">>hist("+to_string(subdivs)+",0,1)";
  const char* std_cuts = "de.t0>700 && de.t0<1695 && deent.td > 0.577350 && deent.td < 1.000 && deent.d0>-80 && "
                         "de.d0<105 && detrkqual.MaxRadius>450 && detrkqual.MaxRadius<680";


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
  string add_on = " && dequal.";
  //  cout << "Relative Efficiency " << inpt_tree->GetEntries((std_cuts+add_on+train_str+" > "+to_string(approx_trkqual)).c_str())/hist->GetEntries() << "\n" << endl;
  return approx_trkqual;  
}
