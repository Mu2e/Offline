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
  char ce_stereo_used_xyname[100];
  char ce_stereo_notused_xyname[100];
  char ce_notstereo_used_xyname[100];
  char ce_notstereo_notused_xyname[100];
  char bkg_stereo_used_xyname[100];
  char bkg_stereo_notused_xyname[100];
  char bkg_notstereo_used_xyname[100];
  char bkg_notstereo_notused_xyname[100];
  char mctxyname[100];
  char ce_stereo_used_phizname[100];
  char ce_stereo_notused_phizname[100];
  char ce_notstereo_used_phizname[100];
  char ce_notstereo_notused_phizname[100];
  char bkg_stereo_used_phizname[100];
  char bkg_stereo_notused_phizname[100];
  char bkg_notstereo_used_phizname[100];
  char bkg_notstereo_notused_phizname[100];
  char mctfzname[100];

  TCanvas* cans[100];
  for(size_t ievt=0;ievt<100;++ievt){
    for(size_t itrk=0;itrk<2;++itrk){
      int jplot = 100*ievt + itrk;
      snprintf(ce_stereo_used_xyname,100,"ce_stereo_used_shxy%i",jplot);
      snprintf(ce_stereo_notused_xyname,100,"ce_stereo_notused_shxy%i",jplot);
      snprintf(ce_notstereo_used_xyname,100,"ce_notstereo_used_shxy%i",jplot);
      snprintf(ce_notstereo_notused_xyname,100,"ce_notstereo_notused_shxy%i",jplot);
      snprintf(bkg_stereo_used_xyname,100,"bkg_stereo_used_shxy%i",jplot);
      snprintf(bkg_stereo_notused_xyname,100,"bkg_stereo_notused_shxy%i",jplot);
      snprintf(bkg_notstereo_used_xyname,100,"bkg_notstereo_used_shxy%i",jplot);
      snprintf(bkg_notstereo_notused_xyname,100,"bkg_notstereo_notused_shxy%i",jplot);
      snprintf(mctxyname,100,"mctshxy%i",jplot);
      snprintf(ce_stereo_used_phizname,100,"ce_stereo_used_shphiz%i",jplot);
      snprintf(ce_stereo_notused_phizname,100,"ce_stereo_notused_shphiz%i",jplot);
      snprintf(ce_notstereo_used_phizname,100,"ce_notstereo_used_shphiz%i",jplot);
      snprintf(ce_notstereo_notused_phizname,100,"ce_notstereo_notused_shphiz%i",jplot);
      snprintf(bkg_stereo_used_phizname,100,"bkg_stereo_used_shphiz%i",jplot);
      snprintf(bkg_stereo_notused_phizname,100,"bkg_stereo_notused_shphiz%i",jplot);
      snprintf(bkg_notstereo_used_phizname,100,"bkg_notstereo_used_shphiz%i",jplot);
      snprintf(bkg_notstereo_notused_phizname,100,"bkg_notstereo_notused_shphiz%i",jplot);
      snprintf(mctfzname,100,"mctshphiz%i",jplot);
      TH2F* ce_stereo_used_shxy = (TH2F*)tdir->Get(ce_stereo_used_xyname);
      TH2F* ce_stereo_notused_shxy = (TH2F*)tdir->Get(ce_stereo_notused_xyname);
      TH2F* ce_notstereo_used_shxy = (TH2F*)tdir->Get(ce_notstereo_used_xyname);
      TH2F* ce_notstereo_notused_shxy = (TH2F*)tdir->Get(ce_notstereo_notused_xyname);
      TH2F* bkg_stereo_used_shxy = (TH2F*)tdir->Get(bkg_stereo_used_xyname);
      TH2F* bkg_stereo_notused_shxy = (TH2F*)tdir->Get(bkg_stereo_notused_xyname);
      TH2F* bkg_notstereo_used_shxy = (TH2F*)tdir->Get(bkg_notstereo_used_xyname);
      TH2F* bkg_notstereo_notused_shxy = (TH2F*)tdir->Get(bkg_notstereo_notused_xyname);
      TH2F* mctshxy = (TH2F*)tdir->Get(mctxyname);
      TH2F* ce_stereo_used_shphiz = (TH2F*)tdir->Get(ce_stereo_used_phizname);
      TH2F* ce_stereo_notused_shphiz = (TH2F*)tdir->Get(ce_stereo_notused_phizname);
      TH2F* ce_notstereo_used_shphiz = (TH2F*)tdir->Get(ce_notstereo_used_phizname);
      TH2F* ce_notstereo_notused_shphiz = (TH2F*)tdir->Get(ce_notstereo_notused_phizname);
      TH2F* bkg_stereo_used_shphiz = (TH2F*)tdir->Get(bkg_stereo_used_phizname);
      TH2F* bkg_stereo_notused_shphiz = (TH2F*)tdir->Get(bkg_stereo_notused_phizname);
      TH2F* bkg_notstereo_used_shphiz = (TH2F*)tdir->Get(bkg_notstereo_used_phizname);
      TH2F* bkg_notstereo_notused_shphiz = (TH2F*)tdir->Get(bkg_notstereo_notused_phizname);
      TH2F* mctshfz = (TH2F*)tdir->Get(mctfzname);

      if(ce_stereo_used_shxy != 0 && ce_stereo_used_shphiz != 0) {
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
	ce_stereo_used_shxy->SetStats(0);
	cans[ican]->cd(divide.rem+1);
	ce_stereo_used_shxy->Draw();
	if(ce_stereo_notused_shxy != 0)ce_stereo_notused_shxy->Draw("same");
	if(ce_notstereo_used_shxy != 0)ce_notstereo_used_shxy->Draw("same");
	if(ce_notstereo_notused_shxy != 0)ce_notstereo_notused_shxy->Draw("same");
	if(bkg_stereo_used_shxy != 0)bkg_stereo_used_shxy->Draw("same");
	if(bkg_stereo_notused_shxy != 0)bkg_stereo_notused_shxy->Draw("same");
	if(bkg_notstereo_used_shxy != 0)bkg_notstereo_used_shxy->Draw("same");
	if(bkg_notstereo_notused_shxy != 0)bkg_notstereo_notused_shxy->Draw("same");
	cans[ican]->cd(divide.rem+nps+1);
	ce_stereo_used_shphiz->Draw();
	ce_stereo_used_shphiz->SetStats(0);
	if(ce_stereo_notused_shphiz != 0)ce_stereo_notused_shphiz->Draw("same");
	if(ce_notstereo_used_shphiz != 0)ce_notstereo_used_shphiz->Draw("same");
	if(ce_notstereo_notused_shphiz != 0)ce_notstereo_notused_shphiz->Draw("same");
	if(bkg_stereo_used_shphiz != 0)bkg_stereo_used_shphiz->Draw("same");
	if(bkg_stereo_notused_shphiz != 0)bkg_stereo_notused_shphiz->Draw("same");
	if(bkg_notstereo_used_shphiz != 0)bkg_notstereo_used_shphiz->Draw("same");
	if(bkg_notstereo_notused_shphiz != 0)bkg_notstereo_notused_shphiz->Draw("same");
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
