#include "TrkDiag/test/TrainPID.C+"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
{
  TChain* cetree = new TChain();
  cetree->Add("/data/TARCeE.root/TrkAnaNeg/trkana");
  TChain* mutree = new TChain();
  mutree->Add("/data/TAR_CRY.root/TrkAnaNeg/trkana");
  mutree->Add("/data/TAR_CRY-p.root/TrkAnaPos/trkana");
  mutree->Add("/data/TARDS-cosmic.root/TrkAnaNeg/trkana");
  mutree->Add("/data/TARDS-cosmic-p.root/TrkAnaPos/trkana");
  TrainPID(cetree,mutree);
//  cetree->Print();
//  mutree->Print();
}

