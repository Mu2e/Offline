#include "TTree.h"
#include "TH1F.h"
#include "TF1.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TCut.h"
#include "TPaveText.h"
#include "TLine.h"
#include <iostream>


// the following approximation is from Czarnecki etal, 'Muon decay in orbit:spectrum of high-energy electrons',
// for E>85 MeV
Double_t DIOCZ(Double_t *x, Double_t *par) {
  double ee = x[0];
  double norm = par[0];
  double mal(25133);
//    double mmu(105.654);
  double emu(105.194);
//    double emue(104.973);
//    double me(0.511);
  double a5(8.6434e-17);
  double a6(1.16874e-17);
  double a7(-1.87828e-19);
  double a8(9.16327e-20);
  double delta = emu - ee - ee*ee/(2*mal);
  return norm*(a5*pow(delta,5) + a6*pow(delta,6) + a7*pow(delta,7) + a8*pow(delta,8));
}

void mu2e(TTree* dio, TTree* con, double diogenrange, double ndio, double ncon,bool weightdio=true) {
// diogenrange is the momentum range over which the DIO events were generated
  double nstopped(7.56e17);
  double capfrac(0.609); 
  double decayfrac = 1.0 - capfrac;
  double ndecay = nstopped*decayfrac;
  double ncap = nstopped*capfrac;
  double conprob(1e-15);
  double momlow(103.3);
  double momhigh(104.7);
  double trueconvmom(104.973);

  unsigned nbins(100);
  double mmin(101);
  double mmax(106);
  double mevperbin = (mmax-mmin)/nbins;

  double conscale = ncap*conprob/ncon;
  cout << "Conversion scale factor =" << conscale << endl;
// dio spectrum
  TF1* diocz_f = new TF1("diocz_f",DIOCZ,85.0,105,1);
  diocz_f->SetLineColor(kGreen);
  diocz_f->SetParameter(0,1.0);
// integrate the DIO spectrum over the range specified.  This is relative to the free decay rate
  double dioint = diocz_f->Integral(trueconvmom-diogenrange,trueconvmom);
  double dioscale(1.0);
  if(weightdio){ 
    dioscale =ndecay*diogenrange/ndio;
  } else {
    dioscale = dioint*ndecay/ndio;
  }
  cout << "DIO scale factor = " << dioscale << endl;

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
  t0cuts[0] = "t0err<2";
  t0cuts[1] = "t0err<1.5";
  t0cuts[2] = "t0err<1.0";
  t0cuts[3] = "t0err<0.9";
  momcuts[0] = "fitmomerr<0.3";
  momcuts[1] = "fitmomerr<0.2";
  momcuts[2] = "fitmomerr<0.18";
  momcuts[3] = "fitmomerr<0.15";
  fitcuts[0] = "fitcon>1e-6";
  fitcuts[1] = "fitcon>1e-4";
  fitcuts[2] = "fitcon>1e-3";
  fitcuts[3] = "fitcon>1e-2";

  unsigned mu2ecut=2;


  for(unsigned icut=0;icut<4;icut++){
    char dioname[50];
    snprintf(dioname,50,"diospec%i",icut);
    char conname[50];
    snprintf(conname,50,"conspec%i",icut);
 
    diospec[icut] = new TH1F(dioname,"Reconstructed momentum;MeV/c;Events per experiment",nbins,mmin,mmax);
    diospec[icut]->SetStats(0);
    diospec[icut]->SetLineColor(kBlue);
    diospec[icut]->Sumw2();

    conspec[icut] = new TH1F(conname,"Reconstructed momentum;MeV/c;Events per experiment",nbins,mmin,mmax);
    conspec[icut]->SetStats(0);
    conspec[icut]->SetLineColor(kRed);
    conspec[icut]->Sumw2();

    TCut quality = ncuts[icut] && t0cuts[icut] && momcuts[icut] && fitcuts[icut];
    TCut final = (reco+pitch+livegate+quality);

    dio->Project(dioname,"fitmom","diowt"*final);
    diospec[icut]->Scale(dioscale);

    con->Project(conname,"fitmom",final);
    conspec[icut]->Scale(conscale);

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

// plot results
  TCanvas* mu2ecan = new TCanvas("mu2e","mu2e result",900,600);
  mu2ecan->Clear();
  mu2ecan->Divide(1,1);
  TCanvas* allcan = new TCanvas("mu2eall","mu2e results",1200,800);
  allcan->Clear();
  allcan->Divide(2,2);
  for(unsigned icut=0;icut<4;icut++){
    allcan->cd(icut+1);
    TH1* diocopy = diospec[icut]->DrawCopy();
    diocopy->SetMinimum(-0.2);
    diocopy->SetMaximum(5);
    conspec[icut]->Draw("same");

    int istart = diospec[icut]->FindFixBin(momlow+0.5*mevperbin);
    int istop = diospec[icut]->FindFixBin(momhigh-0.5*mevperbin);
//    cout << "Integration low edge " << diospec[icut]->GetBinLowEdge(istart) << " for cut at " << momlow << endl;
//    cout << "Integration high edge " << diospec[icut]->GetBinLowEdge(istop)+mevperbin << " for cut at " << momhigh << endl;
    double dint_err, cint_err;
    double dint = diospec[icut]->IntegralAndError(istart,istop,dint_err);
    double cint = conspec[icut]->IntegralAndError(istart,istop,cint_err);

    TPaveText* inttext = new TPaveText(0.5,0.65,0.9,0.8,"NDC");
    char itext[50];
    snprintf(itext,50,"%4.2f MeV/c < P < %4.2f MeV/c",momlow,momhigh);
    inttext->AddText(itext);
    snprintf(itext,50,"DIO integral = %5.3f #pm %4.3f",dint,dint_err);
    inttext->AddText(itext);
    snprintf(itext,50,"Conv. integral = %5.3f #pm %4.3f",cint,cint_err);
    inttext->AddText(itext);
    inttext->Draw();

    TPaveText* cuttext = new TPaveText(0.1,0.5,0.4,0.8,"NDC");  
    char line[40];
    snprintf(line,80,"%4.3f<tan(#lambda)<%4.3f",tdlow,tdhigh);
    cuttext->AddText(line);
    snprintf(line,80,"t0>%5.1f nsec",t0min);
    cuttext->AddText(line);
    sprintf(line,"%s",ncuts[icut].GetTitle());
    cuttext->AddText(line);
    sprintf(line,"%s",t0cuts[icut].GetTitle());
    cuttext->AddText(line);
    sprintf(line,"%s",momcuts[icut].GetTitle());
    cuttext->AddText(line);
    sprintf(line,"%s",fitcuts[icut].GetTitle());
    cuttext->AddText(line);
    cuttext->Draw();

    TLine* momlowl = new TLine(momlow,0.0,momlow,1.5*conspec[icut]->GetBinContent(conspec[icut]->GetMaximumBin()));
    momlowl->SetLineColor(kBlack);
    momlowl->SetLineStyle(2);
    momlowl->SetLineWidth(2);
    momlowl->Draw();

    TLine* momhighl = new TLine(momhigh,0.0,momhigh,1.5*conspec[icut]->GetBinContent(conspec[icut]->GetMaximumBin()));
    momhighl->SetLineColor(kBlack);
    momhighl->SetLineStyle(2);
    momhighl->SetLineWidth(2);
    momhighl->Draw();
    leg->Draw();
    info->Draw();


    if(icut == mu2ecut){
      mu2ecan->cd(0);
      diocopy->Draw();
      conspec[icut]->Draw("same");
      inttext->Draw();
      cuttext->Draw();
      momlowl->Draw();
      momhighl->Draw();
      leg->Draw();
      info->Draw();
    }

 


  }
  allcan->cd(0);
  allcan->SaveAs("mu2e_all.png");
  mu2ecan->SaveAs("mu2e.png");
  mu2ecan->SaveAs("mu2e.eps");

  TCanvas* dioc = new TCanvas("dio","dio",1200,800);
  dioc->Divide(2,2);

  Double_t dmhi = trueconvmom;
  Double_t dmlow = trueconvmom - diogenrange;
  TH1F* diogen = new TH1F("diogen","True DIO momentum;MeV",nbins,dmlow,dmhi);
  TH1F* diowt = new TH1F("diowt","True DIO momentum;MeV",nbins,dmlow,dmhi);
//  diowt->Sumw2();
  dio->Project("diogen","mcmom");
  dio->Project("diowt","mcmom","diowt");
  diowt->Scale(dioscale);
  diowt->SetLineColor(kBlue);
  diogen->SetLineColor(kRed);
  diowt->SetStats(0);
  diogen->SetStats(0);
  dioc->cd(1);
  gPad->SetLogy();
// dead-reconing on spectrum, accounting for bins
  double diofscale = ndecay*(dmhi-dmlow)/nbins;
  diocz_f->SetParameter(0,diofscale);
  diowt->Draw();
  diocz_f->Draw("same");
  diogen->Draw("same");
  TLegend* dioleg = new TLegend(.2,.4,.6,.6);
  dioleg->AddEntry(diogen,"Generated","l");
  dioleg->AddEntry(diowt,"Weighted","l");
  dioleg->AddEntry(diocz_f,"Czarnecki etal","l");
  dioleg->Draw();

  dioc->cd(3);
//  gPad->SetLogy();
  Int_t colors[4] = {kRed,kBlue,kGreen,kBlack};
  TH1F* diogenwin[4] = {0,0,0,0};
  const char* dopt[4] = {"","same","same","same"};
  const char* cutset[4] = {"Cutset A","Cutset B","Cutset C","Cutset D"};
  TLegend* dgenwinleg = new TLegend(.5,.6,.7,.9);
  diogenwin[0] = new TH1F("diogenwin_0","True momentum of DIO in signal box;MeV",100,dmlow,dmhi);
  for(unsigned icut=0;icut<4;icut++){
    char dioname[50];
    snprintf(dioname,50,"diogenwin%i",icut);
    diogenwin[icut] = new TH1F(dioname,"True momentum of DIO in signal box;MeV",100,dmlow,dmhi);
    diogenwin[icut]->SetStats(0);
//   TH1F* diogood[icut] = new TH1F("diogood","True DIO momentum",100,dmlow,dmhi);
//    dio->Project("diogoodwt","mcmom",goodfit);

    TCut quality = ncuts[icut] && t0cuts[icut] && momcuts[icut] && fitcuts[icut];
    TCut final = (reco+pitch+livegate+quality);
    dio->Project(dioname,"mcmom",final+momwin);
    diogenwin[icut]->SetFillColor(colors[icut]);
    dgenwinleg->AddEntry(diogenwin[icut],cutset[icut],"f");
    diogenwin[icut]->Draw(dopt[icut]);
  }
  dgenwinleg->Draw();

  dioc->cd(2);
  for(unsigned icut=0;icut<4;icut++){
    diospec[icut]->SetFillColor(colors[icut]);
    if(icut==0)
      diospec[icut]->Draw("Hist");
    else
      diospec[icut]->Draw("Histsame");
  }
  dgenwinleg->Draw();
  TLine* momlowl = new TLine(momlow,0.0,momlow,1.5*diospec[0]->GetBinContent(diospec[0]->GetMaximumBin()));
  momlowl->SetLineColor(kBlack);
  momlowl->SetLineStyle(2);
  momlowl->SetLineWidth(2);
  momlowl->Draw();

  TLine* momhighl = new TLine(momhigh,0.0,momhigh,1.5*diospec[0]->GetBinContent(diospec[0]->GetMaximumBin()));
  momhighl->SetLineColor(kBlack);
  momhighl->SetLineStyle(2);
  momhighl->SetLineWidth(2);
  momhighl->Draw();

  dioc->SaveAs("diocan.png");

}
