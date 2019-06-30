#include "TrkDiag/test/TrainPID.C+"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
{
  TChain* cetree = new TChain();
  cetree->Add("/data/TAR_FeM.root/TrkAnaNeg/trkana");
  TChain* mutree = new TChain();
  mutree->Add("/data/TAR_CRY.root/TrkAnaNeg/trkana");
  mutree->Add("/data/TAR_CRY-p.root/TrkAnaPos/trkana");
  mutree->Add("/data/TAR_DSC.root/TrkAnaNeg/trkana");
  mutree->Add("/data/TAR_DSC-p.root/TrkAnaPos/trkana");
  TrainPID(cetree,mutree);
//  cetree->Print();
//  mutree->Print();
}

