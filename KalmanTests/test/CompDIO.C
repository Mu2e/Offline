#include "KalmanTests/test/DIOCZ.h"
#include "KalmanTests/test/DIOSC.h"
#include "TH1F.h"
#include "TF1.h"
#include "TLegend.h"
#include "TCanvas.h"
double emue(104.973);

Double_t  DIORatio(Double_t *x, Double_t *par) {
  double diocz = DIOCZ(x,par);
  double diosc = DIOSC(x,par);
  return diosc/diocz;
}


void CompDIO(double minrate=1e-20)
{
  double lowend(102.0);
  TF1* diocz = new TF1("_diocz",DIOCZ,lowend,emue,1);
  TF1* diosc = new TF1("_diosc",DIOSC,lowend,emue,1);
  diocz->SetParameter(0,1.0);
  diosc->SetParameter(0,1.0);
  diocz->SetLineColor(kRed);
  diosc->SetLineColor(kBlue);
  diocz->SetLineStyle(1);
  diosc->SetLineStyle(2);

  TH1F* diospec = new TH1F("diospec","DIO Spectrum Predictions;Electron Energy (MeV);1/#Gamma_{0} d#Gamma/dE", 100,lowend,emue);

  double x[1]; x[0] = lowend;
  double par[1]; par[0] = 1.0;
  double max = 1.2*DIOCZ(x,par);
  diospec->SetMaximum(max);
  diospec->SetMinimum(minrate);
  diospec->SetStats(0);
  TCanvas* diocan = new TCanvas("diocan","diocan",800,1200);
  diocan->Divide(1,2);
  diocan->cd(1);
  gPad->SetLogy();
  diospec->Draw();
  diocz->Draw("same");
  diosc->Draw("same");
  TLegend* leg = new TLegend(0.6,0.7,0.9,0.9);
  leg->AddEntry(diocz,"DIO 2011","L");
  leg->AddEntry(diosc,"DIO 2015","L");
  leg->Draw();

  TF1* ratio = new TF1("ratio",DIORatio,lowend,emue,1);
  ratio->SetParameter(0,1.0);
  ratio->SetLineColor(kBlack);
  TH1F* hratio = new TH1F("hratio","Ratio of 2015 prediction to 2011 prediction; Electron Energy (MeV)",100,lowend,emue);
  hratio->SetStats(0);
  hratio->SetMinimum(0.5);
  hratio->SetMaximum(1.0);
  diocan->cd(2);
  hratio->Draw();
  ratio->Draw("same");
}

