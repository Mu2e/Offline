#include "TrkDiag/test/TrainPID.C+"
{
  TFile cefile("/data/TARCeEMix.root");
  TTree* cetree = (TTree*)cefile.Get("TrkAnaNeg/trkana");
  TFile mufile("/data/TARCRY.root");
  TTree* mutree = (TTree*)mufile.Get("TrkAnaNeg/trkana");
  TrainPID(cetree,mutree);
}

