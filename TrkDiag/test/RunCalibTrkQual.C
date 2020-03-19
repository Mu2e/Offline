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
    CalibTrkQual(tree, train_name);
  }
}

