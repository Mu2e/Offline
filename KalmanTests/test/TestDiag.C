// simple script to test diagnostics
#include "TH1.h"
#include <vector>
#include "TCanvas.h"
#include <iostream>
#include <TROOT.h>
#include <TStyle.h>

void TestDiag(std::vector<TH1*>& plots,unsigned nrows=2,unsigned ncols=2) {
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat("emruo");
  unsigned ican(0);
  unsigned npads = nrows*ncols;
  TCanvas* can(0);
  unsigned ipad = npads+1;
  for(unsigned iplot=0;iplot<plots.size();++iplot){
    if(ipad>npads){
      char canname[80];
      snprintf(canname,80,"TestCan%i",ican);
      can = new TCanvas(canname,"Test",800,800);
      can->Divide(ncols,nrows);
      ++ican;
      ipad=1;
    }
    can->cd(ipad++);
    plots[iplot]->Draw();
  }
}
