#include "TDirectory.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TArc.h"
#include "TH2F.h"
#include <iostream>
#include <string>
#include <vector>

// define a global vector of names. This can be overwritten on invocation

class PlotHelices {
  public:
  PlotHelices(TDirectory* tdir) : _drawfz(true), _csize(250), _tdir(tdir)
  {
  // standard names
    _names.push_back("trk_sh");
    _names.push_back("notselected_sh");
    _names.push_back("selected_sh");
    _names.push_back("mctsh");
    _names.push_back("tc_sh");
    _names.push_back("bkg_used_sh");
    _names.push_back("pri_notused_sh");
    _names.push_back("pri_used_sh");
  }
  PlotHelices(TDirectory* tdir, std::vector<std::string>const& pnames) : _drawfz(true), _csize(250) ,_tdir(tdir),  _names(pnames) {}
  
  void plot(int nmax=20, int nps=3,const char* canname="hcan");

  bool _drawfz;
  int _csize;
  TDirectory* _tdir; // directory
  std::vector<TCanvas*> cans;
  std::vector<std::string> _names;  // names of all the plots
};


void PlotHelices::plot(int nmax, int nps,const char* canname){
  gStyle->SetOptStat(0);
  int ican(-1);
  int ny = 1;
  if(_drawfz)ny=2;

  TCanvas* cans[100];
  bool drawlegend(true);
  TLegend* leg = new TLegend(0.8,0.6,1.0,1.0);
  for(int iplot=1;iplot<=nmax;++iplot){
    std::vector<TH2F*> xyplots, fzplots;
    for(auto name : _names) {
      char xyname[100], fzname[100];

      snprintf(xyname,100,"%sxy%i",name.c_str(),iplot);
      snprintf(fzname,100,"%sfz%i",name.c_str(),iplot);
      TH2F* xyplot = (TH2F*)_tdir->Get(xyname);
      if(xyplot != 0){
	xyplot->SetStats(0);
	xyplots.push_back(xyplot);
      }
      TH2F* fzplot = (TH2F*)_tdir->Get(fzname);
      if(fzplot != 0){
	fzplot->SetStats(0);
	fzplots.push_back(fzplot);
      }
      if(iplot ==1){
	if(xyplot != 0) leg->AddEntry(xyplot,xyplot->GetTitle(),"l");
      }
    }
    if(xyplots.size() > 0){
      div_t divide = div(iplot-1,nps);
      if(divide.rem == 0){
	if(canname != 0 && ican>=0){
	  char fname[100];
	  snprintf(fname,100,"%s.pdf",cans[ican]->GetTitle());
	  cans[ican]->SaveAs(fname);
	}
	++ican;
	char cname[50];
	snprintf(cname,20,"%s_%i",canname,ican);
	cans[ican] = new TCanvas(cname,cname,_csize*nps,_csize*ny);
	cans[ican]->Clear();
	cans[ican]->Divide(nps,ny);
      }
      cans[ican]->cd(divide.rem+1);
      for(size_t ixy = 0;ixy < xyplots.size(); ++ixy){
	if(ixy == 0)
	  xyplots[ixy]->Draw();
	else
	  xyplots[ixy]->Draw("same");
      }
      leg->Draw();

      if(_drawfz && fzplots.size() > 0){
	cans[ican]->cd(divide.rem+nps+1);
	for(size_t ifz = 0;ifz < fzplots.size(); ++ifz){
	  if(ifz == 1)
	    fzplots[ifz]->Draw();
	  else
	    fzplots[ifz]->Draw("same");
	}
      }
    }
  }
  // save last canvas
  if(ican > 0){
    char fname[100];
    snprintf(fname,100,"%s.pdf",cans[ican]->GetTitle());
    cans[ican]->Update();
    cans[ican]->SaveAs(fname);
  }
}

