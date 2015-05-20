#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"

void PanelU(TTree* pad) {
  TH1F* nhits = new TH1F("nhits","N hits/panel",5,-0.5,4.5);
  TH1F* nactive = new TH1F("nactive","N hits/panel",5,-0.5,4.5);
  nhits->SetLineColor(kBlack);
  nactive->SetLineColor(kRed);
  nhits->SetStats(0);
  nactive->SetStats(0);
  TH1F* nres = new TH1F("nres","log_{10} N States/panel",100,0,4);
  nres->SetStats(0);

  TH1F* hpat = new TH1F("hpat","Hit Pattern of best result",3,-0.5,2.5);
  unsigned ibin(1);
  hpat->GetXaxis()->SetBinLabel(ibin++,"Null");
  hpat->GetXaxis()->SetBinLabel(ibin++,"Same-sign");
  hpat->GetXaxis()->SetBinLabel(ibin++,"Opposite-sign");
  hpat->GetXaxis()->SetLabelSize(0.08);
  hpat->SetStats(0);
  
  pad->Project("nhits","nhits");
  pad->Project("nactive","nactive");
  pad->Project("nres","log10(nres)");

  pad->Project("hpat","results[0]._hpat","nhits>=2","",100000);

  TCanvas* ncan = new TCanvas("ncan","ncan",1200,800);
  ncan->Divide(2,2);
  ncan->cd(1);
//  gPad->SetLogy();
  nhits->Draw();
  nactive->Draw("same");
  TLegend* nleg = new TLegend(0.6,0.6,0.9,0.9);
  nleg->AddEntry(nhits,"All Hits","L");
  nleg->AddEntry(nactive,"Active Hits","L");
  nleg->Draw();
  ncan->cd(2);
  gPad->SetLogy();
  nres->Draw();
  ncan->cd(3);
  hpat->Draw();
}

void PanelChisq(TTree* pad) {
  TH1F* chisq1 = new TH1F("chisq1","Best Panel #chi^{2};u #chi^{2}",100,0,20.0);  
  TH1F* chisq2 = new TH1F("chisq2","Best Panel #chi^{2};u #chi^{2}",100,0,20.0);  
  TH1F* chisq3 = new TH1F("chisq3","Best Panel #chi^{2};u #chi^{2}",100,0,20.0);  


  TH1F* dchisq1 = new TH1F("dchisq1","Panel #chi^{2} Next Best - Best;#Delta u #chi^{2}",100,0,20.0);  
  TH1F* dchisq2 = new TH1F("dchisq2","Panel #chi^{2} Next Best - Best;#Delta u #chi^{2}",100,0,20.0);  
  TH1F* dchisq3 = new TH1F("dchisq3","Panel #chi^{2} Next Best - Best;#Delta u #chi^{2}",100,0,20.0);  

  chisq1->SetLineColor(kBlack);
  chisq2->SetLineColor(kGreen);
  chisq3->SetLineColor(kRed);
  chisq1->SetStats(0);
  chisq2->SetStats(0);
  chisq3->SetStats(0);

  dchisq1->SetLineColor(kBlack);
  dchisq2->SetLineColor(kGreen);
  dchisq3->SetLineColor(kRed);
  dchisq1->SetStats(0);
  dchisq2->SetStats(0);
  dchisq3->SetStats(0);

  pad->Project("chisq1","results[0]._chisq","nhits==1","",1000000);
  pad->Project("chisq2","results[0]._chisq","nhits==2","",1000000);
  pad->Project("chisq3","results[0]._chisq","nhits==3","",1000000);

  pad->Project("dchisq1","results[1]._chisq-results[0]._chisq","nhits==1","",1000000);
  pad->Project("dchisq2","results[1]._chisq-results[0]._chisq","nhits==2","",1000000);
  pad->Project("dchisq3","results[1]._chisq-results[0]._chisq","nhits==3","",1000000);

  TCanvas* chi2can = new TCanvas("chi2can","chisquared",1200,800);
  chi2can->Divide(2,2);
  TLegend* cleg = new TLegend(0.5,0.6,0.9,0.9);
  char tstring[50];
  snprintf(tstring,50,"1-hit Panels, <#chi^{2}>=%3.1f",chisq1->GetMean());
  cleg->AddEntry(chisq1,tstring,"L");
  snprintf(tstring,50,"2-hit Panels, <#chi^{2}>=%3.1f",chisq2->GetMean());
  cleg->AddEntry(chisq2,tstring,"L");
  snprintf(tstring,50,"3-hit Panels, <#chi^{2}>=%3.1f",chisq3->GetMean());
  cleg->AddEntry(chisq3,tstring,"L");
  chi2can->cd(1);
  chisq1->Draw();
  chisq2->Draw("same");
  chisq3->Draw("same");
  cleg->Draw();
  chi2can->cd(2);
  dchisq1->Draw();
  dchisq2->Draw("same");
  dchisq3->Draw("same");



}

void UProjections(TTree* pud) {
  TH1F* uhgpos = new TH1F("uhgpos","U Projection, hit position;U position WRT track (mm)",100,-6.0,6.0);
  TH1F* uhbpos = new TH1F("uhbpos","U Projection, hit position;U position WRT track (mm)",100,-6.0,6.0);
  TH1F* uhnpos = new TH1F("uhnpos","U Projection, hit position;U position WRT track (mm)",100,-6.0,6.0);

  TH1F* uherr = new TH1F("uherr","U Projection, hit error;Estimated error (mm)",100,0.0,0.25);
  TH1F* uterr = new TH1F("uterr","U Projection, track error;Projected fit error (mm)",100,0.0,2.0);

  TH1F* uhgpull = new TH1F("uhgpull","U Projection, hit pull",100,-10.0,10.0);
  TH1F* uhnpull = new TH1F("uhnpull","U Projection, hit pull",100,-10.0,10.0);

  uhgpos->SetStats(0);
  uhbpos->SetStats(0);
  uhnpos->SetStats(0);

  uhgpos->SetLineColor(kBlack);
  uhbpos->SetLineColor(kRed);
  uhnpos->SetLineColor(kBlue);

  uherr->SetStats(0);
  uterr->SetStats(0);

  uhgpull->SetStats(0);
  uhnpull->SetStats(0);

  uhgpull->SetLineColor(kBlack);
  uhnpull->SetLineColor(kBlue);

  pud->Project("uhgpos","uinfo._wcpos+uinfo._dr*uinfo._ambig-tupos","uinfo._active&&abs(uinfo._ambig) > 0");
  pud->Project("uhbpos","uinfo._wcpos-uinfo._dr*uinfo._ambig-tupos","uinfo._active&&abs(uinfo._ambig) > 0");
  pud->Project("uhnpos","uinfo._wcpos-tupos","uinfo._active&&uinfo._ambig == 0");

  pud->Project("uherr","uinfo._uerr");
  pud->Project("uterr","tuerr");

  pud->Project("uhgpull","(uinfo._wcpos+uinfo._dr*uinfo._ambig-tupos)/sqrt(_uerr*_uerr+tuerr*tuerr)","uinfo._active&&abs(uinfo._ambig) > 0");
  pud->Project("uhnpull","(uinfo._wcpos-tupos)/sqrt(_uerr*_uerr+tuerr*tuerr)","uinfo._active&&uinfo._ambig == 0");

  TCanvas* uproj = new TCanvas("uproj","uproj",1200,800);
  uproj->Divide(2,2);
  uproj->cd(1);
  uhgpos->Draw();
  uhbpos->Draw("same");
  uhnpos->Draw("same");
  TLegend* uleg = new TLegend(0.6,0.6,0.9,0.9);
  uleg->AddEntry(uhgpos,"Nominal non-null ambiguity","L");
  uleg->AddEntry(uhbpos,"Opposite non-null ambiguity","L");
  uleg->AddEntry(uhnpos,"Null ambiguity","L");
  uleg->Draw();
  uproj->cd(2);
  uherr->Draw();
  uproj->cd(3);
  uterr->Draw();
  uproj->cd(4);
  uhgpull->Draw();
  uhnpull->Draw("same");

}
