#include "TTree.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TCut.h"
#include "TPaveText.h"
#include "TLine.h"
#include <iostream>

void mu2e(TTree* dio, TTree* con, double diogenrange=5.0, double ndio=100000, double nconv=100000) {
// diogenrange is the momentum range over which the DIO events were generated
  double nstopped(7.56e17);
  double capfrac(0.609); 
  double decayfrac = 1.0 - capfrac;
  double ndecay = nstopped*decayfrac;
  double ncap = nstopped*capfrac;
  double conprob(1e-15);
  double momlow(103.3);
  double momhigh(104.7);
  double trueconvmom(105);

  unsigned nbins(100);
  double mmin(101);
  double mmax(106);
  double mevperbin = (mmax-mmin)/nbins;

  TH1F* timeshift = new TH1F("timeshift","T_{0} - muon conversion time;nsec",100,0,100);
  timeshift->SetStats(0);  
  con->Project("timeshift","t0-mct0","fitstatus>0");

  TProfile* truedio = new TProfile("truedio","MC true DIO spectrum;MeV/c;d#Gamma/(#Gamma dE) (MeV^{-1})",nbins,mmin,mmax,0,1);
  dio->Project("truedio","diowt:mcmom");
  truedio->SetMaximum(1e-12);
  truedio->SetMinimum(1e-24);
  TH1F* tdio = new TH1F("tdio","MC true DIO spectrum;MeV/c;d#Gamma/(#Gamma dE) (MeV^{-1})",nbins,mmin,mmax);
  tdio->Sumw2();
  dio->Project("tdio","mcmom","diowt");
  tdio->Scale(diogenrange/(ndio*mevperbin));
  tdio->SetMaximum(1e-12);
  tdio->SetMinimum(1e-24);
  truedio->SetLineColor(kRed);
  truedio->SetMarkerColor(kRed);
  tdio->SetLineColor(kBlue);
  tdio->SetMarkerColor(kBlue);
  truedio->SetStats(0);

  TCanvas* valid = new TCanvas("valid","validation",1200,800);
  valid->Divide(1,2);
  valid->cd(1);
  gPad->SetLogy();
  truedio->Draw();
  tdio->Draw("same");
  TLegend* vleg = new TLegend(0.7,0.8,0.9,0.9);
  vleg->AddEntry(truedio,"CTM formula","L");
  vleg->AddEntry(tdio,"Weighted sampling","L");
  vleg->Draw();

  valid->cd(2);
  timeshift->Draw();


  TH1F* diospec[4];
  TH1F* conspec[4];
// basic cuts
  TCut reco("fitstatus>0");
  double tdlow(0.57735027);
  double tdhigh(1.0);
  double t0min(710);
  char ctext[80];
  snprintf(ctext,80,"td>%f&&td<%f",tdlow,tdhigh);
  TCut pitch(ctext);
  snprintf(ctext,80,"t0>%f",t0min);
  TCut livegate(ctext);
  snprintf(ctext,80,"fitmom>%f&&fitmom<%f",momlow,momhigh);
  TCut momwin(ctext);
// cuts for different tightness of selection
  TCut ncuts[4], t0cuts[4], momcuts[4], fitcuts[4];
  ncuts[0] = "nactive>=20";
  ncuts[1] = "nactive>=20";
  ncuts[2] = "nactive>=25";
  ncuts[3] = "nactive>=30";
  t0cuts[0] = "";
  t0cuts[1] = "t0err<1.5";
  t0cuts[2] = "t0err<1.0";
  t0cuts[3] = "t0err<0.9";
  momcuts[0] = "";
  momcuts[1] = "fitmomerr<0.2";
  momcuts[2] = "fitmomerr<0.18";
  momcuts[3] = "fitmomerr<0.15";
  fitcuts[0] = "";
  fitcuts[1] = "fitcon>1e-4";
  fitcuts[2] = "fitcon>1e-3";
  fitcuts[3] = "fitcon>1e-2";


  for(unsigned ires=0;ires<4;ires++){
    char dioname[50];
    snprintf(dioname,50,"diospec%i",ires);
    char conname[50];
    snprintf(conname,50,"conspec%i",ires);
 
    diospec[ires] = new TH1F(dioname,"Reconstructed momentum;MeV/c;Events per experiment",nbins,mmin,mmax);
    diospec[ires]->SetStats(0);
    diospec[ires]->SetLineColor(kBlue);
    diospec[ires]->Sumw2();

    conspec[ires] = new TH1F(conname,"Reconstructed momentum;MeV/c;Events per experiment",nbins,mmin,mmax);
    conspec[ires]->SetStats(0);
    conspec[ires]->SetLineColor(kRed);
    conspec[ires]->Sumw2();

    TCut quality = ncuts[ires] && t0cuts[ires] && momcuts[ires] && fitcuts[ires];
    TCut final = (reco+pitch+livegate+quality);

    dio->Project(dioname,"fitmom","diowt"*final);
    diospec[ires]->Scale(ndecay*diogenrange/ndio);

    con->Project(conname,"fitmom",final);
    conspec[ires]->Scale(ncap*conprob/nconv);

  }

  TLegend* leg = new TLegend(0.7,0.8,0.9,0.9);
  leg->AddEntry(diospec[0],"DIO","L");
  leg->AddEntry(conspec[0],"Conversion","L");

  TPaveText* info = new TPaveText(0.4,0.8,0.7,0.9,"NDC");
  char text[80];
  snprintf(text,80,"%g stopped muons",nstopped);
  TString snstop(text);
  info->AddText(snstop);
  snprintf(text,80,"%g Conversion Rate",conprob);
  TString sconprob(text);
  info->AddText(sconprob);
  info->SetBorderSize(0);
  valid->SaveAs("mu2e_valid.png");

// log scale
  TCanvas* logcan = new TCanvas("mu2elog","mu2e log scale",1200,800);
  logcan->Clear();
  logcan->Divide(2,2);

  for(unsigned ires=0;ires<4;ires++){
    logcan->cd(ires+1);
    gPad->SetLogy();
    diospec[ires]->SetMinimum(1e-4);
    diospec[ires]->Draw();
    conspec[ires]->Draw("same");

    TPaveText* logtext = new TPaveText(0.1,0.2,0.4,0.5,"NDC");  
 
    char line[40];
    snprintf(line,80,"%4.3f<tan(#lambda)<%4.3f",tdlow,tdhigh);
    logtext->AddText(line);
    snprintf(line,80,"t0>%5.1f nsec",t0min);
    logtext->AddText(line);
    sprintf(line,"%s",ncuts[ires].GetTitle());
    logtext->AddText(line);
    sprintf(line,"%s",t0cuts[ires].GetTitle());
    logtext->AddText(line);
    sprintf(line,"%s",momcuts[ires].GetTitle());
    logtext->AddText(line);
    sprintf(line,"%s",fitcuts[ires].GetTitle());
    logtext->AddText(line);
    logtext->Draw();
    leg->Draw();
    info->Draw();
 
  }
  logcan->cd(0);
  logcan->SaveAs("mu2e_log.png"); 

// linear scale
  TCanvas* lincan = new TCanvas("mu2elin","mu2e linear scale",1200,800);
  lincan->Clear();
  lincan->Divide(2,2);
   for(unsigned ires=0;ires<4;ires++){
    lincan->cd(ires+1);
    TH1* diocopy = diospec[ires]->DrawCopy();
    diocopy->SetMinimum(-0.2);
    diocopy->SetMaximum(4);
    conspec[ires]->Draw("same");

    int istart = diospec[ires]->FindFixBin(momlow+0.5*mevperbin);
    int istop = diospec[ires]->FindFixBin(momhigh-0.5*mevperbin);
    cout << "Integration low edge " << diospec[ires]->GetBinLowEdge(istart) << " for cut at " << momlow << endl;
    cout << "Integration high edge " << diospec[ires]->GetBinLowEdge(istop)+mevperbin << " for cut at " << momhigh << endl;
    double dint_err, cint_err;
    double dint = diospec[ires]->IntegralAndError(istart,istop,dint_err);
    double cint = conspec[ires]->IntegralAndError(istart,istop,cint_err);

    TPaveText* inttext = new TPaveText(0.5,0.6,0.9,0.8,"NDC");
    char itext[50];
    snprintf(itext,50,"%4.2f MeV/c < P < %4.2f MeV/c",momlow,momhigh);
    inttext->AddText(itext);
    snprintf(itext,50,"DIO integral = %5.3f #pm %4.3f",dint,dint_err);
    inttext->AddText(itext);
    snprintf(itext,50,"Conv. integral = %5.3f #pm %4.3f",cint,cint_err);
    inttext->AddText(itext);
    inttext->Draw();


    TPaveText* lintext = new TPaveText(0.1,0.5,0.4,0.8,"NDC");  
    char line[40];
    snprintf(line,80,"%4.3f<tan(#lambda)<%4.3f",tdlow,tdhigh);
    lintext->AddText(line);
    snprintf(line,80,"t0>%5.1f nsec",t0min);
    lintext->AddText(line);
    sprintf(line,"%s",ncuts[ires].GetTitle());
    lintext->AddText(line);
    sprintf(line,"%s",t0cuts[ires].GetTitle());
    lintext->AddText(line);
    sprintf(line,"%s",momcuts[ires].GetTitle());
    lintext->AddText(line);
    sprintf(line,"%s",fitcuts[ires].GetTitle());
    lintext->AddText(line);
    lintext->Draw();

    TLine* momlowl = new TLine(momlow,0.0,momlow,1.5*conspec[ires]->GetBinContent(conspec[ires]->GetMaximumBin()));
    momlowl->SetLineColor(kBlack);
    momlowl->SetLineStyle(2);
    momlowl->SetLineWidth(2);
    momlowl->Draw();

    TLine* momhighl = new TLine(momhigh,0.0,momhigh,1.5*conspec[ires]->GetBinContent(conspec[ires]->GetMaximumBin()));
    momhighl->SetLineColor(kBlack);
    momhighl->SetLineStyle(2);
    momhighl->SetLineWidth(2);
    momhighl->Draw();
    leg->Draw();
    info->Draw();
 
  }
  lincan->cd(0);
  lincan->SaveAs("mu2e_lin.png");

  TCanvas* dioc = new TCanvas("dio","dio",1200,800);
  dioc->Divide(1,2);

  Double_t dmhi = trueconvmom;
  Double_t dmlow = trueconvmom - diogenrange;
  TH1F* diogen = new TH1F("diogen","True DIO momentum;MeV",100,dmlow,dmhi);
  TH1F* diowt = new TH1F("diowt","True DIO momentum;MeV",100,dmlow,dmhi);
//  diowt->Sumw2();
  dio->Project("diogen","mcmom");
  dio->Project("diowt","mcmom","diowt");
  double scale =ndecay*diogenrange/ndio;
  cout << "DIO scale factor = " << scale << endl;
  diowt->Scale(scale);
  diowt->SetLineColor(kBlue);
  diogen->SetLineColor(kRed);
  diowt->SetStats(0);
  diogen->SetStats(0);
  dioc->cd(1);
  gPad->SetLogy();
  diowt->Draw();
  diogen->Draw("same");
  TLegend* dioleg = new TLegend(.2,.4,.6,.6);
  dioleg->AddEntry(diogen,"Generated","l");
  dioleg->AddEntry(diowt,"Weighted","l");
  dioleg->Draw();


  dioc->cd(2);
  gPad->SetLogy();
  Int_t colors[4] = {kRed,kBlue,kGreen,kBlack};
  TH1F* diogenwin[4];
  TH1F* diogood[4];
  char* dopt[4] = {"","same","same","same"};
  char* cutset[4] = {"Cutset A","Cutset B","Cutset C","Cutset D"};
  TLegend* dgenwinleg = new TLegend(.5,.6,.7,.9);
  for(unsigned ires=0;ires<4;ires++){
    char dioname[50];
    snprintf(dioname,50,"diogenwin%i",ires);
    diogenwin[ires] = new TH1F(dioname,"True DIO momentum of fake signals;MeV",100,dmlow,dmhi);
    diogenwin[ires]->SetStats(0);
//   TH1F* diogood[ires] = new TH1F("diogood","True DIO momentum",100,dmlow,dmhi);
//    dio->Project("diogoodwt","mcmom",goodfit);


    TCut quality = ncuts[ires] && t0cuts[ires] && momcuts[ires] && fitcuts[ires];
    TCut final = (reco+pitch+livegate+quality);
    dio->Project(dioname,"mcmom",final+momwin);
    diogenwin[ires]->SetLineColor(colors[ires]);
    dgenwinleg->AddEntry(diogenwin[ires],cutset[ires],"l");
    diogenwin[ires]->Draw(dopt[ires]);
  }
  dgenwinleg->Draw();
  dioc->SaveAs("diocan.png");
}
