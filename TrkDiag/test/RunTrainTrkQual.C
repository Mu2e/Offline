#include "TrkDiag/test/TrainTrkQual.C+"
#include "TrkDiag/test/EfficiencyToTrkQual.C+"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
void RunTrainTrkQual(){

  TFile* file = new TFile("out/trkana-CeEndpoint-mix.TrkQualNew.root", "READ");
  TTree* tree = (TTree*) file->Get("TrkAnaNeg/trkana");
  std::string train_name = "TrkQual";

  //  TrainTrkQual(tree, train_name);

  CalibrateTrkQualEff(tree, train_name);
}

