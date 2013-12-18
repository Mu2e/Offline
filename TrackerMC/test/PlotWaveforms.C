#include "TDirectory.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TArc.h"
#include "TH2F.h"
#include "TString.h"
#include "TList.h"
#include "TRegexp.h"
#include <iostream>
#include <math.h>

void PlotWaveforms(TDirectory* tdir,unsigned nmax=20, unsigned nps=4,const char* name=0){
  gStyle->SetOptStat(0);
  int ican(0);
  int iplot(0);
  char cname[20];
  TCanvas* can(0);
  std::vector<TCanvas*> cans;
  TString match = "SWF*";
  std::cout << "Match = " << match << std::endl;
  TRegexp re(match, kTRUE);
  tdir->ReadAll();
  TList* dlist = tdir->GetList();
  std::cout << "Directory has " << dlist->GetSize() << std::endl;
  dlist->First()->Print();
  TIter nextobj(dlist);
  TObject* obj(0);
  unsigned ndiv = max(1,static_cast<int>(ceil(sqrt(nps))));
  std::cout << " ndiv = " << ndiv << std::endl;
  int mps = ndiv*ndiv;
  while(iplot < nmax && (obj = (TObject*) nextobj())){
    TString s = obj->GetName();
    std::cout << "Object name " << s << std::endl;
    if (s.Index(re) == kNPOS) continue;
    div_t divide = div(iplot,mps);
     std::cout << "divide " << iplot << " by " << nps << " gives  quot " << divide.quot << " rem " << divide.rem << std::endl;
    if(divide.rem == 0){
      ++ican;
      snprintf(cname,20,"can_%i",ican);
      std::cout << "Creating canvas " << cname << std::endl;
      can = new TCanvas(cname,cname,800,800);
      std::cout << "canvas " << can->GetName() << std::endl;
      cans.push_back(can);
      can->Clear();
      can->Divide(ndiv,ndiv);
    }
    can->cd(divide.rem+1);
    std::cout << "Object pointer " << obj << std::endl;
    std::cout << "Drawing object " << obj->GetName() << std::endl;
    obj->Draw();
    ++iplot;
  }
  if(name != 0){
    char fname[100];
    std::cout << "Found Canvases " << cans.size() << std::endl;
    for(size_t jcan=0;jcan<cans.size();jcan++){
      snprintf(fname,100,"%s_%s.png",name,cans[jcan]->GetTitle());
      std::cout << "Saving to file " << fname << std::endl;
      cans[jcan]->SaveAs(fname);
    }
  }
}
