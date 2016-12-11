#include "TDirectory.h"
#include "TLegend.h"
#include "TSpectrum2.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TArc.h"
#include "TH2F.h"
#include <iostream>
void PlotTimeSpectra2D(TDirectory* tdir,double sigma=2.0,double minn=5,unsigned nmax=20, unsigned nps=2,const char* name=0){
  gStyle->SetOptStat(0);
  int ican(-1);
  unsigned iplot(0);
  TCanvas* cans[100];
  char canname[100];
  TSpectrum2 tp2(100);
  TH1F* dummy = new TH1F("dummy","dummy",10,0,1);
  dummy->SetMarkerStyle(23);
  dummy->SetMarkerColor(kRed);
  dummy->SetMarkerSize(1);
  for(size_t ievt=0;ievt<10000;++ievt){
    bool first(true);
    char rname[100];
    char sname[100];
    char cname[100];
      snprintf(rname,100,"rawptspectrum%lu",ievt);
      snprintf(sname,100,"looseptspectrum%lu",ievt);
      snprintf(cname,100,"convptspectrum%lu",ievt);
      TH2F* rh = (TH2F*)tdir->Get(rname);
      TH2F* th = (TH2F*)tdir->Get(sname);
      TH2F* ch = (TH2F*)tdir->Get(cname);
      if(rh != 0 && th != 0 && ch != 0){
	div_t divide = div(iplot,nps*nps);
	//      std::cout << "divide " << iplot << " by " << nps << " gives  quot " << divide.quot << " rem " << divide.rem << std::endl;
	if(divide.rem == 0){
	  ++ican;
	  snprintf(canname,20,"can_%i",ican);
	  cans[ican] = new TCanvas(canname,canname,800,600);
	  cans[ican]->Clear();
	  cans[ican]->Divide(nps,2);
	}
	cans[ican]->cd(divide.rem+1);
	rh->SetStats(0);
	th->SetStats(0);
	ch->SetStats(0);
	th->SetLineColor(kGreen);
	th->SetMarkerColor(kGreen);
	th->SetMarkerStyle(2);
	// make an absolute threshold
	double thresh(0.99);
	double maxn = th->GetMaximum();
	if(maxn > minn) thresh = minn/maxn;
	std::cout << "threshold = " << thresh << std::endl;
	tp2.Search(th,sigma,"nobackgroundnomarkovnodraw",thresh);
// 	th->SetTitle("Time Spectrum;nsec;#phi");
//        th->Draw("box");
	th->Draw();
//       rh->GetXaxis()->SetTitle("nsec");
  //      rh->GetYaxis()->SetTitle("#phi");
        rh->Draw("same");
	ch->SetMarkerColor(kRed);
        ch->SetMarkerStyle(5);
        ch->Draw("same");
        if(first){
          first = false;
          TLegend* leg = new TLegend(0.0,0.7,0.3,0.9);
          leg->AddEntry(rh,"All hits","P");
          leg->AddEntry(th,"Selected hits","P");
          leg->AddEntry(ch,"Conversion hits","P");
          leg->AddEntry(dummy,"TSpectrum Peak","P");
          leg->Draw();
        }
	++iplot;
	if(iplot > nmax)break;
      }
      if(iplot > nmax)break;
  }
  if(name != 0){
    char fname[100];
    for(size_t jcan=0;jcan<=ican;jcan++){
      snprintf(fname,100,"%s_%s.png",name,cans[jcan]->GetTitle());
      cans[jcan]->SaveAs(fname);
    }
  }
}
