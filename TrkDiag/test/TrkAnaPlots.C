#include "TH1F.h"
#include "TF1.h"
#include "TTree.h"
#include "TLegend.h"
#include "TPad.h"
#include "TPaveText.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TStyle.h"
#include "TLine.h"
#include "TArrow.h"
#include "TCut.h"
#include "TBox.h"
#include "TMath.h"
#include "TProfile.h"
#include "TDirectory.h"
#include "Math/Math.h"
#include "THStack.h"


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


// The following is from Alexx Perloff, JetMetaAnalysis
double fnc_dscb(double*xx,double*pp) {
  double x   = xx[0];
// gaussian core
  double N   = pp[0];//norm
  double mu  = pp[1];//mean
  double sig = pp[2];//variance
  // transition parameters
  double a1  = pp[3];
  double p1  = pp[4];
  double a2  = pp[5];
  double p2  = pp[6];

  double u   = (x-mu)/sig;
  double A1  = TMath::Power(p1/TMath::Abs(a1),p1)*TMath::Exp(-a1*a1/2);
  double A2  = TMath::Power(p2/TMath::Abs(a2),p2)*TMath::Exp(-a2*a2/2);
  double B1  = p1/TMath::Abs(a1) - TMath::Abs(a1);
  double B2  = p2/TMath::Abs(a2) - TMath::Abs(a2);

  double result(N);
  if      (u<-a1) result *= A1*TMath::Power(B1-u,-p1);
  else if (u<a2)  result *= TMath::Exp(-u*u/2);
  else            result *= A2*TMath::Power(B2+u,-p2);
  return result;
}

void Draw(TTree* ta) {

  TH2F* evspep = new TH2F("evspep","Associated Cluster Energy vs Track Momentum;P (MeV/c);E (MeV)",50,0,200,50,0,200);
  TH2F* evspmp = new TH2F("evspmp","Associated Cluster Energy vs Track Momentum;P (MeV/c);E (MeV)",50,0,200,50,0,200);
  TH2F* evspem = new TH2F("evspem","Associated Cluster Energy vs Track Momentum;P (MeV/c);E (MeV)",50,0,200,50,0,200);
  TH2F* evspmm = new TH2F("evspmm","Associated Cluster Energy vs Track Momentum;P (MeV/c);E (MeV)",50,0,200,50,0,200);
  evspep->SetMaximum(50);
  evspmp->SetMaximum(50);
  evspep->SetStats(0);
  evspmp->SetStats(0);
  evspep->SetLineColor(kOrange);
  evspmp->SetLineColor(kCyan);
  evspem->SetMaximum(50);
  evspmm->SetMaximum(50);
  evspem->SetStats(0);
  evspmm->SetStats(0);
  evspem->SetLineColor(kRed);
  evspmm->SetLineColor(kBlue);

  TCanvas* can = new TCanvas("can","can",800,800);

  ta->Project("evspem","demc.eclust:dem.mom","dem.status>0&&tcnt.ndemc>0&&demmc.pdg==11");
  ta->Project("evspmm","demc.eclust:dem.mom","dem.status>0&&tcnt.ndemc>0&&demmc.pdg==13");
  ta->Project("evspep","demc.eclust:dem.mom","dem.status>0&&tcnt.ndemc>0&&demmc.pdg==-11");
  ta->Project("evspmp","demc.eclust:dem.mom","dem.status>0&&tcnt.ndemc>0&&demmc.pdg==-13");

  evspem->Draw();
  evspem->Draw("box");
  evspmm->Draw("boxsame");
  TLegend* leg = new TLegend(0.1,0.6,0.5,0.9);
  leg->AddEntry(evspem,"True e^{-}","L");
  leg->AddEntry(evspmm,"True #mu^{-}","L");
  leg->AddEntry(evspep,"True e^{+}","L");
  leg->AddEntry(evspmp,"True #mu^{+}","L");
  leg->Draw();
}

void MomResp(TTree* ta, double tqcut, double nmu) {
// cuts
  TCut reco("dem.status>0");
  char ctext[80];
  snprintf(ctext,80,"dem.trkqual>%f",tqcut);
  TCut goodfit(ctext);
  double tdlow(0.57735027);
  double tdhigh(1.0);
  double t0min(700.0);
  double t0max(1695.0);
  snprintf(ctext,80,"dem.td>%5.5f&&dem.td<%5.5f",tdlow,tdhigh);
  TCut rpitch = TCut(ctext);
  snprintf(ctext,80,"dem.t0>%f&&dem.t0<%f",t0min,t0max);
  TCut livegate = TCut(ctext);
  TCut cosmic = TCut("dem.d0<105 && dem.d0>-80 && (dem.d0+2/dem.om)>450 && (dem.d0+2/dem.om)<680");
  TCut rmomloose("dem.mom>100.0");
  TCut physics = rpitch+cosmic+livegate+rmomloose;

  TF1* dscb = new TF1("dscb",fnc_dscb,-10.0,5,7);
  dscb->SetParName(0,"Norm");
  dscb->SetParName(1,"x0");
  dscb->SetParName(2,"sigma");
  dscb->SetParName(3,"ANeg");
  dscb->SetParName(4,"PNeg");
  dscb->SetParName(5,"APos");
  dscb->SetParName(6,"PPos");

  TCanvas* rcan = new TCanvas("rcan","Momentum Resolution",600,600);
  rcan->Clear();
  gStyle->SetOptFit(111111);
  gStyle->SetOptStat("oumr");
  gPad->SetLogy();
  TH1F* momresp = new TH1F("momresp","momentum response at tracker;MeV/c",251,-10.0,4.0);
  momresp->Sumw2();
  TCut final = reco+goodfit;
  ta->Project("momresp","dem.mom-demmcgen.mom","evtinfo.evtwt"*final);
  //momresp->Scale(1.0/nmu);
  //    ta->Project(mname,"fit.mom-mcent.mom",final);
  double integral = momresp->GetEntries()*momresp->GetBinWidth(1);
  cout << "Integral = " << integral << " mean = " << momresp->GetMean() << " rms = " << momresp->GetRMS() << endl;
  dscb->SetParameters(0.05*integral,-0.6,0.3,0.7,3.0,3.0,3.0);
  dscb->SetNpx(1000);
  dscb->SetParLimits(3,0.0,50.0);
  dscb->SetParLimits(4,1.0,50.0);
  dscb->SetParLimits(5,0.0,50.0);
  dscb->SetParLimits(6,1.0,50.0);

  momresp->SetMinimum(0.5);
  momresp->Fit("dscb","LRQ");
  momresp->Fit("dscb","LRM");
}
void MomRes(TTree* ta, double tqcut) {
// cuts
  TCut reco("dem.status>0");
  char ctext[80];
  snprintf(ctext,80,"dem.trkqual>%f",tqcut);
  TCut goodfit(ctext);
  double tdlow(0.57735027);
  double tdhigh(1.0);
  double t0min(700.0);
  double t0max(1695.0);
  snprintf(ctext,80,"dem.td>%5.5f&&dem.td<%5.5f",tdlow,tdhigh);
  TCut rpitch = TCut(ctext);
  snprintf(ctext,80,"dem.t0>%f&&dem.t0<%f",t0min,t0max);
  TCut livegate = TCut(ctext);
  TCut cosmic = TCut("dem.d0<105 && dem.d0>-80 && (dem.d0+2/dem.om)>450 && (dem.d0+2/dem.om)<680");
  TCut rmomloose("dem.mom>100.0");
  TCut physics = rpitch+cosmic+livegate+rmomloose;

  TF1* dscb = new TF1("dscb",fnc_dscb,-2.0,2.5,7);
  dscb->SetParName(0,"Norm");
  dscb->SetParName(1,"x0");
  dscb->SetParName(2,"sigma");
  dscb->SetParName(3,"ANeg");
  dscb->SetParName(4,"PNeg");
  dscb->SetParName(5,"APos");
  dscb->SetParName(6,"PPos");

  TCanvas* rcan = new TCanvas("rcan","Momentum Resolution",600,600);
  rcan->Clear();
  gStyle->SetOptFit(111111);
  gStyle->SetOptStat("oumr");
  gPad->SetLogy();
  TH1F* momres = new TH1F("momres","momentum resolution at start of tracker;MeV/c",251,-4,4);
  momres->Sumw2();
  TCut final = reco+goodfit+physics;
  ta->Project("momres","dem.mom-demmcent.mom","evtinfo.evtwt"*final);
  //    ta->Project(mname,"fit.mom-mcent.mom",final);
  double integral = momres->GetEntries()*momres->GetBinWidth(1);
  cout << "Integral = " << integral << " mean = " << momres->GetMean() << " rms = " << momres->GetRMS() << endl;
  dscb->SetParameters(3*integral,momres->GetMean()+0.07,0.3*momres->GetRMS(),0.9,3.5,1.5,6.0);

  momres->SetMinimum(0.5);
  momres->Fit("dscb","LRQ");
  momres->Fit("dscb","LRM");

  TLine* zero = new TLine(0.0,0.0,0.0,momres->GetBinContent(momres->GetMaximumBin()));
  zero->SetLineStyle(2);
  zero->Draw();

  TPaveText* rtext = new TPaveText(0.1,0.5,0.4,0.9,"NDC");
  rtext->AddText("Reco Cuts");
  char line[40];
  snprintf(line,80,"%4.3f<tan(#lambda)<%4.3f",tdlow,tdhigh);
  rtext->AddText(line);
  snprintf(line,80,"t0>%5.1f nsec",t0min);
  rtext->AddText(line);
  sprintf(line,"%s",goodfit.GetTitle());
  rtext->AddText(line);
  rtext->Draw();

}

void Acc(TTree* ta, int ngen) {
  unsigned nbins(8);
  double bmax = nbins-0.5;

  TH1F* acc = new TH1F("acc","CE Acceptance #times Efficiency;;Cummulative a#times#epsilon",nbins,-0.5,bmax);
  TH1F* racc = new TH1F("racc","CE Acceptance #times Efficiency;;Relative a#times#epsilon",nbins,-0.5,bmax);
//  acc->Sumw2();
//  racc->Sumw2();
  unsigned ibin(1);
  acc->GetXaxis()->SetBinLabel(ibin++,"All CE");
  acc->GetXaxis()->SetBinLabel(ibin++,"MC Selection");
  acc->GetXaxis()->SetBinLabel(ibin++,"KF Track fit");
  acc->GetXaxis()->SetBinLabel(ibin++,"Fit Quality");
  acc->GetXaxis()->SetBinLabel(ibin++,"Livegate");
  acc->GetXaxis()->SetBinLabel(ibin++,"Reco pitch");
  acc->GetXaxis()->SetBinLabel(ibin++,"Cosmic Rejection");
  acc->GetXaxis()->SetBinLabel(ibin++,"Momentum window");


  ibin = 1;
  racc->GetXaxis()->SetBinLabel(ibin++,"All CE");
  racc->GetXaxis()->SetBinLabel(ibin++,"MC Selection");
  racc->GetXaxis()->SetBinLabel(ibin++,"KF Track fit");
  racc->GetXaxis()->SetBinLabel(ibin++,"Fit Quality");
  racc->GetXaxis()->SetBinLabel(ibin++,"Livegate");
  racc->GetXaxis()->SetBinLabel(ibin++,"Reco pitch");
  racc->GetXaxis()->SetBinLabel(ibin++,"Cosmic Rejection");
  racc->GetXaxis()->SetBinLabel(ibin++,"Momentum window");

  ibin = 0;
  const char* binnames[11] ={"0.0","1.0","2.0","3.0","4.0","5.0","6.0","7.0","8.0","9.0","10.0"};

  TCut mcsel = "demmc.ndigigood>15&&demmcent.mom>90.0";
  // &&demmcent.td>0.55&&demmcent.td<1.05";
  //&&fmod(demmcent.t0,1695.0)>500.0";
  TCut reco = "dem.status>0";
  TCut goodfit = "dem.trkqual>0.4";
  TCut livegate = "dem.t0>500.0&&dem.t0<1695";
  TCut rpitch = "dem.td>0.57735027&&dem.td<1.0";
  TCut cosmic = "dem.d0<105 && dem.d0>-80 && (dem.d0+2/dem.om)>450 && (dem.d0+2/dem.om)<680";
  TCut rmom = "dem.mom>100.0";

  ta->Project("acc",binnames[ibin++],"evtinfo.evtwt");
  ta->Project("+acc",binnames[ibin++],"evtinfo.evtwt"*mcsel);
  ta->Project("+acc",binnames[ibin++],"evtinfo.evtwt"*(mcsel+reco));
  ta->Project("+acc",binnames[ibin++],"evtinfo.evtwt"*(mcsel+reco+goodfit));
  ta->Project("+acc",binnames[ibin++],"evtinfo.evtwt"*(mcsel+reco+goodfit+livegate));
  ta->Project("+acc",binnames[ibin++],"evtinfo.evtwt"*(mcsel+reco+goodfit+livegate+rpitch));
  ta->Project("+acc",binnames[ibin++],"evtinfo.evtwt"*(mcsel+reco+goodfit+livegate+rpitch+cosmic));
  ta->Project("+acc",binnames[ibin++],"evtinfo.evtwt"*(mcsel+reco+goodfit+livegate+rpitch+cosmic+rmom));

  double all = acc->GetBinContent(1);
  double norm = ngen;
  if(ngen < 0)
    norm = all;
  double prev = norm;
  for(ibin=1;ibin<=nbins;ibin++){
    racc->SetBinContent(ibin,acc->GetBinContent(ibin)/prev);
    prev = acc->GetBinContent(ibin);
  }
  cout << "Found " << norm << "Entries." << endl;
  racc->SetMaximum(1.1);
  acc->Scale(1.0/(float)norm);
  acc->SetMaximum(1.1);
  acc->SetStats(0);
  racc->SetStats(0);
  acc->GetXaxis()->SetLabelSize(0.06);
  racc->GetXaxis()->SetLabelSize(0.06);
  acc->SetMarkerSize(2.0);
  racc->SetMarkerSize(2.0);
  acc->GetYaxis()->SetTitleSize(0.05);
  racc->GetYaxis()->SetTitleSize(0.05);

  gStyle->SetPaintTextFormat("5.4f");
  TCanvas* acan = new TCanvas("acan","Acceptance",1200,800);
  acan->Clear();
  acan->Divide(1,2);
  acan->cd(1);
  TPad* tp = (TPad*)acan->cd(1);
  tp->SetBottomMargin(0.15);
  acc->Draw("histtext0");
  acan->cd(2);
  tp = (TPad*)acan->cd(2);
  tp->SetBottomMargin(0.15);
  racc->Draw("histtext0");
}

