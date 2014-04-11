#include "TDirectory.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TArc.h"
#include "TH2F.h"
#include <iostream>

void PlotTimeSpectra(TDirectory* tdir,unsigned nmax=20, int nps=3, const char* name=0){
  gStyle->SetOptStat(0);
  TCanvas* cans[100];
  TH1F* dummy = new TH1F("dummy","dummy",10,0,1);
  dummy->SetMarkerStyle(23);
  dummy->SetMarkerColor(kRed);
//  dummy->SetMarkerSize(1);
  int ican(-1);
  int iplot(0);
  bool first(true);
  for(size_t ievt=0;ievt<10000;++ievt){
    char canname[100];
    char rname[100];
    char tname[100];
//    char lname[100];
    char cname[100];
    char pname[100];
    snprintf(rname,100,"rawtspectrum%lu",ievt);
    snprintf(tname,100,"tightnodeltatspectrum%lu",ievt);
  //  snprintf(lname,100,"loosetspectrum%lu",ievt);
    snprintf(cname,100,"convtspectrum%lu",ievt);
    snprintf(pname,100,"protontspectrum%i",ievt);
    TH1F* rh = (TH1F*)tdir->Get(rname);
    TH1F* th = (TH1F*)tdir->Get(tname);
//    TH1F* lh = (TH1F*)tdir->Get(lname);
    TH1F* ch = (TH1F*)tdir->Get(cname);
    TH1F* ph = (TH1F*)tdir->Get(pname);
    if(rh != 0 && th != 0 && ch != 0 && ph != 0){
      div_t divide = div(iplot,nps*nps);
      //      std::cout << "divide " << iplot << " by " << nps << " gives  quot " << divide.quot << " rem " << divide.rem << std::endl;
      if(divide.rem == 0){
	++ican;
	snprintf(canname,20,"can_%i",ican);
	cans[ican] = new TCanvas(canname,canname,800,600);
	cans[ican]->Clear();
	cans[ican]->Divide(nps,2);
      }
      unsigned ipave = divide.rem+1;
      cans[ican]->cd(ipave);
      rh->SetStats(0);
      th->SetStats(0);
//      lh->SetStats(0);
      ch->SetStats(0);
      ph->SetStats(0);
//
      char title[100];
      snprintf(title,100,"Time Spectrum event %lu",ievt);
      rh->SetTitle(title);
      rh->GetXaxis()->SetTitle("nsec");
      rh->GetYaxis()->SetTitle("# hits");
      rh->Draw();
      //        th->SetLineWidth(2);
//      lh->SetFillColor(kYellow);
//      lh->Draw("same");
      th->SetFillColor(kGreen);
      th->Draw("same");
      ch->SetLineStyle(1);
      ch->SetLineWidth(2);
      ch->Draw("same");
      ph->Draw("same");
      TLegend* leg(0);
      if(ipave==1){
	leg = new TLegend(0.1,0.7,0.5,0.9);
	leg->AddEntry(rh,"All hits","l");
	leg->AddEntry(th,"Selected hits","F");
//	leg->AddEntry(lh,"Loose Selected hits","F");
	leg->AddEntry(ch,"Conversion hits","LF");
	leg->AddEntry(ph,"Proton hits","LF");
	leg->AddEntry(dummy,"Reconstructed Time Peak","P");
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
