#include "AnalysisUtilities/scripts/CalibTrkQual.C+"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
void RunCalibTrkQual(){

  TFile* file = new TFile("../Offline/out/trkana-CeplusEndpoint-mix.trkqual-dev.root", "READ");
  TTree* tree = (TTree*) file->Get("TrkAnaPos/trkana");

  const int n_trkquals = 1;
  std::string train_names[n_trkquals] = {"TrkQualPos"};//{"TrkQual", "TrkQualNeg"};
  std::string leaf_names[n_trkquals] = {"TrkQual"};//{"TrkQual", "TrkQualNeg"};
  for (int i_trkqual = 0; i_trkqual < n_trkquals; ++i_trkqual) {
    std::string train_name = train_names[i_trkqual];
    std::string leaf_name = leaf_names[i_trkqual];

    TCut time_cut("de.t0>700 && de.t0<1695");
    TCut tandip_cut("deent.td > 0.577350 && deent.td < 1.000");
    TCut d0_cut("deent.d0>-80 && de.d0<105");
    TCut max_radius_cut("detrkqual.MaxRadius>450 && detrkqual.MaxRadius<680");
    TCut signal_cut = time_cut + tandip_cut + d0_cut + max_radius_cut;

    CalibTrkQual(tree, train_name, leaf_name, signal_cut);
  }
  std::cout << "All calibrations finished!" << std::endl;
}

