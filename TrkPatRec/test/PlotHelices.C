#include "TDirectory.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TArc.h"
#include "TH2F.h"
#include <iostream>

void PlotHelices(TDirectory* tdir,unsigned nmax=20, int nps=3,const char* name=0){
  gStyle->SetOptStat(0);
  int ican(-1);
  unsigned iplot(0);
  char gxyname[100];
  char cxyname[100];
  char cbxyname[100];
  char bxyname[100];
  char sxyname[100];
  char bsxyname[100];
  char mctxyname[100];
  char gfzname[100];
  char cfzname[100];
  char cbfzname[100];
  char mctfzname[100];
  char bfzname[100];
  char sfzname[100];
  char bsfzname[100];
  TCanvas* cans[100];
  for(size_t ievt=0;ievt<10000;++ievt){
    for(size_t itrk=0;itrk<100;++itrk){
      int jplot = 100*ievt + itrk;
      snprintf(gxyname,100,"gshxy%i",jplot);
      snprintf(cxyname,100,"cshxy%i",jplot);
      snprintf(cbxyname,100,"cbshxy%i",jplot);
      snprintf(mctxyname,100,"mctshxy%i",jplot);
      snprintf(bxyname,100,"bshxy%i",jplot);
      snprintf(sxyname,100,"sshxy%i",jplot);
      snprintf(bsxyname,100,"bsshxy%i",jplot);
      snprintf(gfzname,100,"gshphiz%i",jplot);
      snprintf(cfzname,100,"cshphiz%i",jplot);
      snprintf(cbfzname,100,"cbshphiz%i",jplot);
      snprintf(mctfzname,100,"mctshphiz%i",jplot);
      snprintf(bfzname,100,"bshphiz%i",jplot);
      snprintf(sfzname,100,"sshphiz%i",jplot);
      snprintf(bsfzname,100,"bsshphiz%i",jplot);
      TH2F* gshxy = (TH2F*)tdir->Get(gxyname);
      TH2F* cshxy = (TH2F*)tdir->Get(cxyname);
      TH2F* cbshxy = (TH2F*)tdir->Get(cbxyname);
      TH2F* mctshxy = (TH2F*)tdir->Get(mctxyname);
      TH2F* bshxy = (TH2F*)tdir->Get(bxyname);
      TH2F* sshxy = (TH2F*)tdir->Get(sxyname);
      TH2F* bsshxy = (TH2F*)tdir->Get(bsxyname);
      TH2F* gshfz = (TH2F*)tdir->Get(gfzname);
      TH2F* cshfz = (TH2F*)tdir->Get(cfzname);
      TH2F* cbshfz = (TH2F*)tdir->Get(cbfzname);
      TH2F* mctshfz = (TH2F*)tdir->Get(mctfzname);
      TH2F* bshfz = (TH2F*)tdir->Get(bfzname);
      TH2F* sshfz = (TH2F*)tdir->Get(sfzname);
      TH2F* bsshfz = (TH2F*)tdir->Get(bsfzname);
      if(gshxy != 0 && gshfz != 0) {
	div_t divide = div(iplot,nps);
	//	std::cout << "divide " << iplot << " by " << nps << " gives  quot " << divide.quot << " rem " << divide.rem << std::endl;
	if(divide.rem == 0){
	  ++ican;
	  char cname[50];
	  snprintf(cname,20,"can_%i",ican);
	  cans[ican] = new TCanvas(cname,cname,400*nps,800);
	  cans[ican]->Clear();
	  cans[ican]->Divide(nps,2);
	}
	gshxy->SetStats(0);
	cans[ican]->cd(divide.rem+1);
	gshxy->Draw();
	if(cshxy != 0)cshxy->Draw("same");
	if(cbshxy != 0)cbshxy->Draw("same");
	if(mctshxy != 0)mctshxy->Draw("same");
	if(bshxy != 0)bshxy->Draw("same");
	if(sshxy != 0)sshxy->Draw("same");
	if(bsshxy != 0)bsshxy->Draw("same");
	cans[ican]->cd(divide.rem+nps+1);
	gshfz->SetStats(0);
	gshfz->Draw();        
	if(cshfz != 0)cshfz->Draw("same");
	if(cbshfz != 0)cbshfz->Draw("same");
	if(mctshfz != 0)mctshfz->Draw("same");
	if(bshfz != 0)bshfz->Draw("same");
	if(sshfz != 0)sshfz->Draw("same");
	if(bsshfz != 0)bsshfz->Draw("same");
	++iplot;
	if(iplot > nmax)break;
      }
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
