#include "LoadChain.C" 
#include "TrkDiag/test/TrkRecoDiag.C+"
#include "TCanvas.h"
#include <string>
void DrawTRD(const char* dir){
  string files = string(dir) + string("/files.txt");
  TChain* t = LoadChain(files.c_str(),"TrkRecoDiag/trdiag");
  TrkRecoDiag trd(t,100000);
  trd.Loop();
  TCanvas* trdcan = new TCanvas("trdcan","TRD Canvas",1000,500);
  trdcan->Divide(2,1);
  trdcan->cd(1);
  trd._eff->Draw();
  trdcan->cd(2);
  gPad->SetLogy();
  trd._rej->Draw();
}
