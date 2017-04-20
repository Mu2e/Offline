#include "TDirectory.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TArc.h"
#include "TH2F.h"
#include <iostream>

void PlotHelices(TDirectory* tdir,int nmax=20, int nps=3,const char* name=0){
  gStyle->SetOptStat(0);
  int ican(-1);
  char ce_used_xyname[100];
  char ce_notused_xyname[100];
  char bkg_used_xyname[100];
  char selected_xyname[100];
  char notselected_xyname[100];
  char tc_xyname[100];
  char trk_xyname[100];
  char mctxyname[100];
//
  char ce_used_fzname[100];
  char ce_notused_fzname[100];
  char bkg_used_fzname[100];
  char selected_fzname[100];
  char notselected_fzname[100];
  char tc_fzname[100];
  char trk_fzname[100];
  char mctfzname[100];
  TCanvas* cans[100];
  for(int iplot=1;iplot<nmax;++iplot){
    snprintf(ce_used_xyname,100,"ce_used_shxy%i",iplot);
    snprintf(ce_notused_xyname,100,"ce_notused_shxy%i",iplot);
    snprintf(bkg_used_xyname,100,"bkg_used_shxy%i",iplot);
    snprintf(selected_xyname,100,"selected_shxy%i",iplot);
    snprintf(notselected_xyname,100,"notselected_shxy%i",iplot);
    snprintf(tc_xyname,100,"tc_shxy%i",iplot);
    snprintf(trk_xyname,100,"trk_shxy%i",iplot);
    snprintf(mctxyname,100,"mctshxy%i",iplot);
    TH2F* ce_used_shxy = (TH2F*)tdir->Get(ce_used_xyname);
    if(ce_used_shxy == 0)break;
    TH2F* ce_notused_shxy = (TH2F*)tdir->Get(ce_notused_xyname);
    TH2F* bkg_used_shxy = (TH2F*)tdir->Get(bkg_used_xyname);
    TH2F* selected_shxy = (TH2F*)tdir->Get(selected_xyname);
    TH2F* notselected_shxy = (TH2F*)tdir->Get(notselected_xyname);
    TH2F* tc_shxy = (TH2F*)tdir->Get(tc_xyname);
    TH2F* trk_shxy = (TH2F*)tdir->Get(trk_xyname);
    TH2F* mctshxy = (TH2F*)tdir->Get(mctxyname);
//
    snprintf(ce_used_fzname,100,"ce_used_shfz%i",iplot);
    snprintf(ce_notused_fzname,100,"ce_notused_shfz%i",iplot);
    snprintf(bkg_used_fzname,100,"bkg_used_shfz%i",iplot);
    snprintf(selected_fzname,100,"selected_shfz%i",iplot);
    snprintf(notselected_fzname,100,"notselected_shfz%i",iplot);
    snprintf(tc_fzname,100,"tc_shfz%i",iplot);
    snprintf(trk_fzname,100,"trk_shfz%i",iplot);
    snprintf(mctfzname,100,"mctshfz%i",iplot);
    TH2F* ce_used_shfz = (TH2F*)tdir->Get(ce_used_fzname);
    TH2F* ce_notused_shfz = (TH2F*)tdir->Get(ce_notused_fzname);
    TH2F* bkg_used_shfz = (TH2F*)tdir->Get(bkg_used_fzname);
    TH2F* selected_shfz = (TH2F*)tdir->Get(selected_fzname);
    TH2F* notselected_shfz = (TH2F*)tdir->Get(notselected_fzname);
    TH2F* tc_shfz = (TH2F*)tdir->Get(tc_fzname);
    TH2F* trk_shfz = (TH2F*)tdir->Get(trk_fzname);
    TH2F* mctshfz = (TH2F*)tdir->Get(mctfzname);
    div_t divide = div(iplot-1,nps);
    if(divide.rem == 0){
      if(name != 0 && ican>=0){
	char fname[100];
	snprintf(fname,100,"%s_%s.png",name,cans[ican]->GetTitle());
	cans[ican]->SaveAs(fname);
      }
      ++ican;
      char cname[50];
      snprintf(cname,20,"hcan_%i",ican);
      cans[ican] = new TCanvas(cname,cname,400*nps,800);
      cans[ican]->Clear();
      cans[ican]->Divide(nps,2);
    }
    ce_used_shxy->SetStats(0);
    cans[ican]->cd(divide.rem+1);
    if(trk_shxy != 0)trk_shxy->Draw();
    if(notselected_shxy != 0)notselected_shxy->Draw("same");
    if(selected_shxy != 0)selected_shxy->Draw("same");
    if(mctshxy != 0)mctshxy->Draw("same");
    if(tc_shxy != 0)tc_shxy->Draw("same");
    if(ce_notused_shxy != 0)ce_notused_shxy->Draw("same");
    if(bkg_used_shxy != 0)bkg_used_shxy->Draw("same");
    if(ce_used_shxy != 0)ce_used_shxy->Draw("same");
    cans[ican]->cd(divide.rem+nps+1);
    if(trk_shfz != 0)trk_shfz->Draw();
    if(notselected_shfz != 0)notselected_shfz->Draw("same");
    if(selected_shfz != 0)selected_shfz->Draw("same");
    if(mctshfz != 0)mctshfz->Draw("same");
    if(tc_shfz != 0)tc_shfz->Draw("same");
    if(ce_notused_shfz != 0)ce_notused_shfz->Draw("same");
    if(bkg_used_shfz != 0)bkg_used_shfz->Draw("same");
    if(ce_used_shfz != 0)ce_used_shfz->Draw("same");
   } 
}
