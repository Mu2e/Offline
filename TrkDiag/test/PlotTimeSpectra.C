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
  TH1F* tchits = new TH1F("tchits","tchits",10,0,1);
  tchits->SetMarkerStyle(23);
  tchits->SetMarkerColor(kRed);
  TH1F* tccalo = new TH1F("tccalo","tccalo",10,0,1);
  tccalo->SetMarkerStyle(34);
  tccalo->SetMarkerColor(kRed);
  int ican(-1);
  unsigned iplot(0);
  bool first(true);
  for(size_t ievt=0;ievt<10000;++ievt){
    char canname[100];
    char rname[100];
    char tname[100];
    char lname[100];
    char cname[100];
    char pname[100];
    char caname[100];
    snprintf(rname,100,"rawtspectrum%lu",ievt);
    snprintf(tname,100,"clusttspectrum%lu",ievt);
    snprintf(lname,100,"seltspectrum%lu",ievt);
    snprintf(cname,100,"allconvtspectrum%lu",ievt);
    snprintf(pname,100,"clustconvtspectrum%i",ievt);
    snprintf(caname,100,"calotspectrum%i",ievt);
    TH1F* rh = (TH1F*)tdir->Get(rname);
    TH1F* th = (TH1F*)tdir->Get(tname);
    TH1F* lh = (TH1F*)tdir->Get(lname);
    TH1F* ch = (TH1F*)tdir->Get(cname);
    TH1F* ph = (TH1F*)tdir->Get(pname);
    TH1F* ca = (TH1F*)tdir->Get(caname);
    if(rh != 0 && th != 0 &&lh != 0 && ch != 0 && ph != 0){
      div_t divide = div(iplot,nps*nps);
      //      std::cout << "divide " << iplot << " by " << nps << " gives  quot " << divide.quot << " rem " << divide.rem << std::endl;
      if(divide.rem == 0){
	++ican;
	snprintf(canname,20,"can_%i",ican);
	cans[ican] = new TCanvas(canname,canname,1200,1000);
	cans[ican]->Clear();
	cans[ican]->Divide(nps,nps);
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
      rh->GetYaxis()->SetTitle("# hits/20 nsec");
 //     rh->SetFillColor(kGreen);
      rh->Draw();
      th->SetTitle(title);
      th->GetXaxis()->SetTitle("nsec");
      th->GetYaxis()->SetTitle("# hits/20 nsec");
      th->SetFillColor(kGreen);
      th->Draw("same");
      //        th->SetLineWidth(2);
//      lh->SetFillColor(kYellow);
//      lh->Draw("same");
//      rh->Draw("same");
      ch->SetLineStyle(1);
      ch->SetLineWidth(2);
      ch->SetFillColor(0);
      ch->SetLineColor(kRed);
      ch->Draw("same");
      ph->Draw("same");
      if(ca != 0){
	ca->SetFillStyle(3001);
	ca->SetFillColor(kMagenta);
	ca->Draw("same");
      }
      TLegend* leg(0);
      if(ipave==1){
	leg = new TLegend(0.1,0.7,0.4,0.9);
	leg->AddEntry(rh,"All hits","l");
	leg->AddEntry(lh,"Selected e^{-} hits","F");
	leg->AddEntry(th,"Clustered hits","F");
	leg->AddEntry(ch,"All CE hits","LF");
	leg->AddEntry(ph,"Clustered CE hits","LF");
	if(ca)leg->AddEntry(ca,"CaloCluster","LF");
	leg->AddEntry(tchits,"Hits only Cluster Time","P");
	leg->AddEntry(tccalo,"Hits+Calo Cluster Time","P");
	leg->Draw();
      }
      ++iplot;
      if(iplot > nmax)break;
    }
    if(iplot > nmax)break;
  }
  if(name != 0){
    char fname[100];
    for(int jcan=0;jcan<=ican;jcan++){
      snprintf(fname,100,"%s_%s.png",name,cans[jcan]->GetTitle());
      cans[jcan]->SaveAs(fname);
    }
  }
}
