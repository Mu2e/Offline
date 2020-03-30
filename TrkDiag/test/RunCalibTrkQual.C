//#include "TrkDiag/test/TrainTrkQual.C+"
#include "TrkDiag/test/CalibTrkQual.C+"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
void RunCalibTrkQual(){

  TFile* file = new TFile("../out/trkana-CeEndpoint-mix.TrkQualNew.root", "READ");
  TTree* tree = (TTree*) file->Get("TrkAnaNeg/trkana");

  const int n_trkquals = 2;
  std::string train_names[n_trkquals] = {"TrkQual", "TrkQualNew"};
  for (int i_trkqual = 0; i_trkqual < n_trkquals; ++i_trkqual) {
    std::string train_name = train_names[i_trkqual];

    TCut time_cut("de.t0>700 && de.t0<1695");
    TCut tandip_cut("deent.td > 0.577350 && deent.td < 1.000");
    TCut d0_cut("deent.d0>-80 && de.d0<105");
    TCut max_radius_cut("detrkqual.MaxRadius>450 && detrkqual.MaxRadius<680");
    TCut signal_cut = time_cut + tandip_cut + d0_cut + max_radius_cut;

    CalibTrkQual(tree, train_name, signal_cut);
  }
}

