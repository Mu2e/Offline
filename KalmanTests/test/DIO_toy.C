// this macro tests the effect of resolution tails on the DIO spectrum in the signal region
// Dave Brown, 29 Jan 2013
//
#include "TH1F.h"
#include "TF1.h"
#include "TTree.h"
#include "TLegend.h"
#include "TPad.h"
#include "TPaveText.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TStyle.h"
#include "TCut.h"
#include "TMath.h"
#include "TProfile.h"
#include "TDirectory.h"
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

Double_t crystalball (Double_t *x, Double_t *par) {
  // par[0] : norm
  // par[1] : x0
  // par[2] : sigma
  // par[3] : n
  // par[4] : alpha
  // par[5] : fraction of exponential tail
  // par[6] : tail exponential lambda

  double dx = x[0]-par[1];
  if ( dx/fabs(par[2]) > -1.*par[4]) {
    double g = par[0]*TMath::Gaus(x[0], par[1], par[2]);
//    double g2 = par[5]*par[0]*TMath::Gaus(x[0], par[1], par[6]);
//    return g1+g2;
    double e = par[0]*par[5]*dx*exp(-(dx)/par[6])/(par[6]*par[6]);
    return g+e;
  }
  else {
    double A = pow(par[3]/fabs(par[4]), par[3])*exp(-0.5*par[4]*par[4]);
    double B = par[3]/fabs(par[4]) - fabs(par[4]);
    return par[0]*A*pow(B-dx/fabs(par[2]), -1.*par[3]);
  }
}

void
Dio_toy(unsigned ntrials,double corefactor=1.0,double ratefactor=1.0,double lambdafactor=1.0,unsigned nres=10) {
// fixed external numbers
  double momlow(95.0);
  double momhi(105.0);
  double intlow(103.4);
  double inthi(104.8);
  double plotlow(100.0);
  double plothi(105.0);
  double receff(0.11);
  double nstopped(6.7e17);
  double capfrac(0.609); 
  double decayfrac = 1.0 - capfrac;
  double ndecay = nstopped*decayfrac;
  size_t nbins(151);

// dio spectrum
  TF1* diocz_f = new TF1("diocz_f",DIOCZ,momlow,momhi,1);
  diocz_f->SetLineColor(kGreen);
  diocz_f->SetParameter(0,1.0);
// resolution function
  TF1* cball = new TF1("cball",crystalball,-10.0,5,7);
  cball->SetParName(0,"Norm");
  cball->SetParName(1,"x0");
  cball->SetParName(2,"sigma");
  cball->SetParName(3,"n");
  cball->SetParName(4,"alpha");
  cball->SetParName(5,"tailfrac");
  cball->SetParName(6,"taillambda");
// set parameters according to cutset 'C'
  cball->SetParameters(1.0,-0.8262,0.2416,3.551,0.3836,0.003,0.3259);
  cball->FixParameter(1,-0.8262);
  cball->FixParameter(2,0.2416*corefactor);
  cball->FixParameter(3,3.551);
  cball->FixParameter(4,0.3836);
  cball->FixParameter(5,0.003*ratefactor);
  cball->FixParameter(6,0.3259*lambdafactor);
  TH1D* rawdio = new TH1D("rawdio","DIO e^{-} momentum at decay;P_{decay} (MeV/c);N #mu2e",nbins,plotlow,plothi);
  TH1D* momres = new TH1D("momres","Reco e^{-} momentum resolution function;P_{reco}-P_{decay} (MeV/c)",nbins,-10,6);
  TH1D* recodio = new TH1D("recodio","Reco resolution smeared DIO e^{-} momentum;P_{reco} (MeV/c);N #mu2e",nbins,plotlow,plothi);
  recodio->Sumw2();
  rawdio->SetStats(0);
  recodio->SetStats(0);
  momres->SetStats(0);

// get normalization
  double dioint = diocz_f->Integral(momlow,momhi);
  double rawscale = ndecay*dioint/ntrials;
  cout << "DIO raw scale factor = " << rawscale << endl;


  for(unsigned itrial=0;itrial<ntrials;++itrial){
    double diomom = diocz_f->GetRandom();
    rawdio->Fill(diomom);
    for(unsigned jtrial=0;jtrial<nres;++jtrial){
      double mres = cball->GetRandom();
      momres->Fill(mres);
      recodio->Fill(diomom+mres); 
    }
  }
  rawdio->Scale(rawscale);
  recodio->Scale(rawscale*receff/nres);

  double mevperbin = (plothi-plotlow)/nbins;
  int istart = recodio->FindFixBin(intlow+0.5*mevperbin);
  int istop = recodio->FindFixBin(inthi-0.5*mevperbin);
  double dint_err;
  double dint = recodio->IntegralAndError(istart,istop,dint_err);

  TCanvas* can = new TCanvas("diocan","DIO Toy MC",1200,800);
  can->Divide(2,2);
  can->cd(1);
  gPad->SetLogy();
  rawdio->SetMinimum(0.5);
  rawdio->Draw();
  can->cd(2);
  gPad->SetLogy();
  momres->Fit(cball,"LIR");

  TPaveText* fit = new TPaveText(0.2,0.2,0.6,0.4,"NDC");
  fit->SetBorderSize(2);
  fit->SetTextSize(0.06);
  char ftext[80];
  snprintf(ftext,80,"Core #sigma = %3.2f MeV",cball->GetParameter(2));
  TString csig(ftext);
  fit->AddText(csig);
  snprintf(ftext,80,"Tail fraction = %4.4f",cball->GetParameter(5));
  TString tfrac(ftext);
  fit->AddText(tfrac);
  snprintf(ftext,80,"Tail #lambda = %4.4f MeV",cball->GetParameter(6));
  TString tlam(ftext);
  fit->AddText(tlam);

  fit->Draw();

  can->cd(3);
  gPad->SetLogy();
  recodio->SetMinimum(1.0e-4);
  recodio->Draw();

  TLine* intlowl = new TLine(intlow,0.0,intlow,1.5*recodio->GetBinContent(recodio->GetMaximumBin()));
  intlowl->SetLineColor(kBlack);
  intlowl->SetLineStyle(2);
  intlowl->SetLineWidth(2);
  intlowl->Draw();

  TLine* inthil = new TLine(inthi,0.0,inthi,1.5*recodio->GetBinContent(recodio->GetMaximumBin()));
  inthil->SetLineColor(kBlack);
  inthil->SetLineStyle(2);
  inthil->SetLineWidth(2);
  inthil->Draw();

  can->cd(4);

  TPaveText* info = new TPaveText(0.2,0.2,0.8,0.8,"NDC");
  info->SetBorderSize(2);
  info->SetTextSize(0.06);

  char text[80];
  snprintf(text,80,"%g stopped muons",nstopped);
  TString snstop(text);
  info->AddText(snstop);
  info->AddLine(0.0,1.0/6.0,1.0,1.0/6.0);
  snprintf(text,80,"%g reco + selection efficiency",receff);
  TString rece(text);
  info->AddText(rece);
  snprintf(text,80,"Core #sigma scale factor = %3.1f",corefactor);
  TString coref(text);
  info->AddText(coref);
  snprintf(text,80,"Exponential tail rate factor = %3.1f",ratefactor);
  TString ratef(text);
  info->AddText(ratef);
  snprintf(text,80,"Exponential #lambda factor = %3.1f",lambdafactor);
  TString lambdaf(text);
  info->AddText(lambdaf);
  snprintf(text,80,"DIO window integral = %4.2f #pm %3.2f",dint,dint_err);
  TString di(text);
  info->AddLine(0.0,5.0/6.0,1.0,5.0/6.0);
  info->AddText(di);
  info->Draw();

}
