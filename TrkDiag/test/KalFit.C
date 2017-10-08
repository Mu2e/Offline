
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

// basic parameters
class KalFit {
  public:
  TTree* _tdiag;
  TCanvas* rcan;
  TCanvas* acan;
  double tdlow;
  double tdhigh;
  double t0min,t0max;
  double momlow;
  double momhigh;
  int minnhits;
  size_t icut;
  unsigned minnactive[4];
  double maxt0err[4];
  double maxmomerr[4];
  double minfitcon[4];
  TCut ncuts[4], t0cuts[4], momcuts[4], fitcuts[4], goodfit[4];
  TCut trkqualcut[4] = {"trkqual>0.3","trkqual>0.4","trkqual>0.5","trkqual>0.6"};
  TCut reco,cosmic,rmom,rpitch,livegate;
  TCut tpitch, tt0,tmom,tnhits,mcsel;
  TCut rmomloose;
  bool donecuts;

  KalFit(TTree*,size_t ic=2);
  void Cuts();
  void Hit();
  void T2d();
  void Trk();
  void AccPlots();
  void Acc(int ngen=-1);
  void Res(unsigned mincut=0,unsigned maxcut=3);
  void Res2(int ires);
  void Ambig(int acut=0);
  void Resid();
  void Con();
  void Error();
  void Drift();
  void NHits();
  void Nactive();
  void TrkQualRes();
  void Mom();
  void Rad();
  void MomTails(int iwt=-1);
  void MomTailsHits();
  void StrawMat();
};

KalFit::KalFit(TTree* trkdiag,size_t ic) : _tdiag(trkdiag), tdlow(0.57735027),
  tdhigh(1.0),
  t0min(700),t0max(1695),
  momlow(103.75),
  momhigh(105.0),
  minnhits(15),
  icut(ic),
  rmomloose("fit.mom>100.0"),
  donecuts(false)
{
  double t_minnactive[] = {15,22,25,30};
  double t_maxt0err[] = {10,2.0,1.5,1.3};
  double t_maxmomerr[]   = {1.0,0.28,0.25,0.22};
  double t_minfitcon[]  = {1e-15,1e-3,2e-3,1e-2};
  for(size_t jcut=0;jcut<4;++jcut){
    minnactive[jcut] = t_minnactive[jcut];
    maxt0err[jcut] = t_maxt0err[jcut];
    maxmomerr[jcut] = t_maxmomerr[jcut];
    minfitcon[jcut] = t_minfitcon[jcut];

  }
  Cuts();
}


void KalFit::Cuts() {
  for(size_t ic=0;ic<4;++ic){
    char cutstring[100];
    snprintf(cutstring,100,"nactive>=%i",minnactive[ic]);
    ncuts[ic] = TCut(cutstring);
    snprintf(cutstring,100,"t0err<%3.2f",maxt0err[ic]);
    t0cuts[ic] = TCut(cutstring);
    snprintf(cutstring,100,"fit.momerr<%4.3f",maxmomerr[ic]);
    momcuts[ic] = TCut(cutstring);
    snprintf(cutstring,100,"fit.con>%5.4f",minfitcon[ic]);
    fitcuts[ic] = TCut(cutstring);
  }
  char ctext[80];
  snprintf(ctext,80,"td>%5.5f&&td<%5.5f",tdlow,tdhigh);
  rpitch = TCut(ctext);
  snprintf(ctext,80,"t0>%f&&t0<%f",t0min,t0max);
  livegate = TCut(ctext);
  snprintf(ctext,80,"mcent.td>%4.3f&&mcent.td<%4.3f",tdlow-0.05,tdhigh+0.2);
  tpitch = TCut(ctext);
  snprintf(ctext,80,"fmod(mcmid.t0,1695.0)>%f",500.);
  tt0 = TCut(ctext);
  tmom = TCut("mcent.mom>100.0");
  snprintf(ctext,80,"mc.ndigigood>=%i",minnhits);
  tnhits = TCut(ctext);
  mcsel = tmom+tpitch+tt0+tnhits;

  reco = TCut("fit.status>0");
  cosmic = TCut("d0<105 && d0>-80 && (d0+2/om)>450 && (d0+2/om)<680");
  for(size_t jcut=0;jcut<4;++jcut){
//    goodfit[jcut] = ncuts[jcut]+momcuts[jcut]+fitcuts[jcut]+t0cuts[jcut];
    goodfit[jcut] = trkqualcut[jcut];
  }
  snprintf(ctext,80,"fit.mom>%f&&fit.mom<%f",momlow,momhigh);
  rmom = TCut(ctext);
  donecuts = true;
}

Double_t splitgaus(Double_t *x, Double_t *par) {
  Double_t retval;
  Double_t core;
  Double_t tail;
  Float_t xval = x[0];
  if(xval > par[1]) {
    core = exp(-0.5*pow((xval-par[1])/par[2],2))/par[2];
    tail = par[4]*exp(-0.5*pow((xval-par[1])/par[5],2))/par[5];
  } else {
    core = exp(-0.5*pow((xval-par[1])/par[3],2))/par[3];
    tail = (1/par[2]-1/par[3]+par[4]/par[5])*exp(-0.5*pow((xval-par[1])/par[6],2));
  }
  retval = par[0]*0.398942*(core+tail);
  // add a tail Gaussian
  return retval;
}

Double_t doublegausexp(Double_t *x, Double_t *par) {
  Double_t retval(0.0);
  Double_t core;
  Double_t gtail;
  Double_t etail;
  Float_t xval = x[0]-par[3];
//  double nu = par[7]+1.0;
  core = 0.398942*(1.0-par[1]-par[2])*exp(-0.5*pow(xval/par[4],2))/par[4];
  gtail = 0.398942*par[1]*exp(-0.5*pow(xval/par[5],2))/par[5];
  if(xval <= 0.0)
    etail = -par[2]*xval*exp(xval/par[6])/pow(par[6],2);
//    etail = 0.5*par[2]*xval*xval*exp(xval/par[6])/pow(par[6],3);
//    etail = TMath::Gamma(nu)*par[2]*pow((double)fabs(xval),par[7])*exp(xval/par[6])/pow(par[6],nu);
  else
    etail = 0;
  retval = par[0]*(core+gtail+etail);
  return retval;
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

void KalFit::Hit() {
  TCut goodhit = reco+ncuts[icut]+momcuts[icut];
  TCut active = "_active";

  TCanvas* dcan = new TCanvas("driftcan","driftcan",1200,800);
  TH1F* dres = new TH1F("dres","Drift radius resolution;mm",100,-1,1);
  TH1F* dpull = new TH1F("dpull","Drift radius pull",100,-10,10);


  TCanvas* rcan = new TCanvas("residcan","residcan",1200,800);
  TH1F* aresid = new TH1F("aresid","Residual;mm",100,-1.5,1.5);
  TH1F* iresid = new TH1F("iresid","Residual;mm",100,-1.5,1.5);
  TH1F* arpull = new TH1F("arpull","residual pull",100,-10,10);
  TH1F* irpull = new TH1F("irpull","residual pull",100,-10,10);


//  THStack* resid = new THStack("resid","TrkStrawHit Residual;mm");
  aresid->SetFillColor(kRed);
  aresid->SetStats(1);
  iresid->SetFillColor(kBlue);
  iresid->SetStats(0);

//  THStack* rpull = new THStack("rpull","TrkStrawHit Pull");
  arpull->SetFillColor(kRed);
  arpull->SetStats(1);
  irpull->SetFillColor(kBlue);
  irpull->SetStats(0);

  TCut cdrift("tshmc._dist>0.5");
  _tdiag->Project("dres","_rdrift-tshmc._dist",goodhit+active+cdrift);
  _tdiag->Project("dpull","(_rdrift-tshmc._dist)/_rdrifterr",goodhit+active+cdrift);

  _tdiag->Project("iresid","_resid",goodhit+!active);
  _tdiag->Project("aresid","_resid",goodhit+active);

  _tdiag->Project("irpull","_resid/_residerr",goodhit+!active);
  _tdiag->Project("arpull","_resid/_residerr",goodhit+active);

  dcan->Clear();
  dcan->Divide(2,1);
  dcan->cd(1);
  dres->Fit("gaus","","",-0.3,0.18);
  dcan->cd(2);
  dpull->Fit("gaus","","",-3,1.5);


  TLegend* rleg = new TLegend(0.1,0.8,0.4,0.9);
  rleg->AddEntry(aresid,"Active Hits","f");
  rleg->AddEntry(iresid,"Inactive Hits","f");
  rcan->Divide(2,1);
  rcan->cd(1);
  aresid->Fit("gaus");
  iresid->Draw("same");
  rleg->Draw();
  rcan->cd(2);
  arpull->Fit("gaus");
  irpull->Draw("same");

  TCanvas* t2dcan = new TCanvas("t2dcan","t2dcan",1200,800);
  TH2F* drad = new TH2F("drad","Drift radius;true drift radius (mm);reco drift radius (mm)",
      55,-0.1,2.6,55,-0.1,2.6);
  TProfile* pdresid = new TProfile("pdresid","Drift Residual vs drift radius;reco drift radius (mm);reco - true radius (mm)",25,0.0,2.5,0,2.5,"s");
  TProfile* pdresidt = new TProfile("pdresidt","Drift Residual vs true drift radius;MC true drift radius (mm);reco - true radius (mm)",25,0.0,2.5,0.,2.5,"s");
  TH2F* dresid = new TH2F("dresid","Drift Residual vs drift radius;reco drift radius (mm);reco - true radius (mm)",25,0.0,2.5,25,-0.5,0.8);
  TH2F* dresidt = new TH2F("dresidt","Drift Residual vs true drift radius;MC true drift radius (mm);reco - true radius (mm)",25,0.0,2.5,25,-0.5,0.8);
  _tdiag->Project("drad","_rdrift:tshmc._dist",goodhit+active);
  _tdiag->Project("pdresid","_rdrift-tshmc._dist:_rdrift",goodhit+active);
  _tdiag->Project("pdresidt","_rdrift-tshmc._dist:tshmc._dist",goodhit+active);
  _tdiag->Project("dresid","_rdrift-tshmc._dist:_rdrift",goodhit+active);
  _tdiag->Project("dresidt","_rdrift-tshmc._dist:tshmc._dist",goodhit+active);
  TH1F* pdresid_1 = new TH1F("pdresid_1","Mean drift residual vs drift radius;reco drift radius (mm); mean reco - true radius (mm)",25,0.0,2.5);
  TH1F* pdresid_2 = new TH1F("pdresid_2","Drift residual RMS vs drift radius;reco drift radius (mm); RMS reco - true radius (mm)",25,0.0,2.5);
  TH1F* pdresidt_1 = new TH1F("pdresidt_1","Mean drift residual vs true radius;true radius (mm); mean reco - true radius (mm)",25,0.0,2.5);
  TH1F* pdresidt_2 = new TH1F("pdresidt_2","Drift residual RMS vs true radius;true radius (mm); RMS reco - true radius (mm)",25,0.0,2.5);
  dresidt->SetLineColor(kBlue);
  dresidt->SetFillColor(kBlue);
  dresid->SetLineColor(kCyan);
  dresid->SetFillColor(kCyan);
  dresid->FitSlicesY();
  TH1D *dresid_1 = (TH1D*)gDirectory->Get("dresid_1");
  TH1D *dresid_2 = (TH1D*)gDirectory->Get("dresid_2");
  dresid_1->SetLineColor(kCyan);
  dresid_1->SetTitle("Drift Residual Mean vs Drift Radius;reco drift radius (mm);Residual mean (mm)");
  dresid_1->SetStats(0);
  dresid_2->SetLineColor(kCyan);
  dresid_2->SetTitle("Drift Residual Sigma vs Drift Radius;reco drift radius (mm);Residual #sigma (mm)");
  dresid_2->SetStats(0);
  drad->SetStats(0);
  dresid->SetStats(0);
  dresidt->SetStats(0);

  for(unsigned ibin=0;ibin<25;++ibin){
    pdresid_1->Fill(pdresid->GetBinCenter(ibin+1),pdresid->GetBinContent(ibin+1));
    pdresid_2->Fill(pdresid->GetBinCenter(ibin+1),pdresid->GetBinError(ibin+1));
    pdresidt_1->Fill(pdresidt->GetBinCenter(ibin+1),pdresidt->GetBinContent(ibin+1));
    pdresidt_2->Fill(pdresidt->GetBinCenter(ibin+1),pdresidt->GetBinError(ibin+1));
  }

  //dres->SetStats(0);
  t2dcan->Divide(3,2);
  t2dcan->cd(1);
  drad->Draw("colorZ");
  t2dcan->cd(2);
  dresid->Draw("box");
  t2dcan->cd(3);
  dresidt->Draw("box");
  t2dcan->cd(4);
  dresid_1->Fit("pol2");
  t2dcan->cd(5);
  dresid_2->Fit("pol2");
  t2dcan->cd(6);
  dres->Draw();

  TCanvas* drcan = new TCanvas("drcan","drcan",1200,800);
  drcan->Divide(2,2);
  drcan->cd(1);
  pdresid_1->Draw();
  drcan->cd(2);
  pdresid_2->Draw();
  drcan->cd(3);
  pdresidt_1->Draw();
  drcan->cd(4);
  pdresidt_2->Draw();

//  TCanvas* tcan = new TCanvas("ht0can","hit_t0can",1200,800);
//  TH1F* t0res = new TH1F("t0res","hit t0 resolution;nsec",100,-10,10);
//  TH1F* t0pull = new TH1F("t0pull","hit t0 pull",200,-30,30);
//  TH2F* dt0 = new TH2F("dt0","Hit t0;true t0 (nsec);reco t0 (nsec)",
//      100,0,2000,100,0,2000);
//  _tdiag->Project("t0res","_ht-_mcht","_active");
//  _tdiag->Project("t0pull","(_ht-_mcht)/_t0err","_active");
//  tcan->Clear();
//  tcan->Divide(2,2);
//  tcan->cd(1);
//  _tdiag->Draw("_ht:_mcht>>dt0","_active");
//  tcan->cd(2);
//  t0res->Fit("gaus");
//  tcan->cd(3);
//  t0pull->Fit("gaus");

  TCanvas* tdcan = new TCanvas("tdcan","tdcan",1200,800);
  TH1F* tdres = new TH1F("tdres","#Deltat V resolution;mm",100,-300,300);
  TH1F* tdpull = new TH1F("tdpull","#Deltat V pull",100,-15,15);
  TH2F* dtd = new TH2F("dtd","Hit V position;true V (mm);#Deltat V (mm)",
      100,-600,600,100,-600,600);
//  TH2F* pocatd = new TH2F("pocatd","Hit POCA V;true V (mm);POCA V (mm)",
//      100,-600,600,100,-600,600);

  TCut central("abs(tshmc._len)<50");
  _tdiag->Project("tdres","_tddist-tshmc._len","_active"+central);
  _tdiag->Project("tdpull","(_tddist-tshmc._len)/_tdderr","_active"+central);
  tdcan->Clear();
  tdcan->Divide(2,2);
  tdcan->cd(1);
  _tdiag->Draw("_tddist:tshmc._len>>dtd","_active");
  tdcan->cd(2);
  dtd->FitSlicesY(0,20,80);
  TH1D* dtd_1 = (TH1D*)gDirectory->Get("dtd_1");
  dtd_1->SetTitle("Average reco #Delta-t V vs true V;true V (mm);reco V (mm)");
  dtd_1->Fit("pol1");
  tdcan->cd(3);
  tdres->Fit("gaus","","",-100.,100.);
  tdcan->cd(4);
  tdpull->Fit("gaus","","",-3,3);

}

void KalFit::T2d(){
  TH2F* rt2d = new TH2F("rt2d","Reco drift distance vs #Delta t;Hit t - MC t_{0} (nsec);Reco Drift Distance (mm)",100,0,50.0,100,-0.05,2.55);
  TH2F* tt2d = new TH2F("tt2d","True drift distance vs #Delta t;Hit t - MC t_{0} (nsec);True Drift Distance (mm)",25,0,50.0,50,-0.05,2.55);
  TProfile* rt2dp = new TProfile("rt2dp","Reco drift distance vs #Delta t;Hit t - MC t_{0} (nsec);Reco Drift Distance (mm)",100,0,50.0,-0.05,2.55);
  TProfile* tt2dp = new TProfile("tt2dp","True drift distance vs #Delta t;Hit t - MC t_{0} (nsec);True Drift Distance (mm)",100,0,50.0,-0.05,2.55);
  TProfile* tt2dp2 = new TProfile("tt2dp2","True drift distance vs #Delta t;Hit t - MC t_{0} (nsec);True Drift Distance (mm)",100,0,50.0,-0.05,2.55,"s");
  _tdiag->Project("tt2d","tshmc._dist:_ht-mcmid.t0%1695","fit.status>0&&_active&&fit.con>1e-2");
  _tdiag->Project("rt2d","_rdrift:_ht-mcmid.t0%1695","fit.status>0&&_active&&fit.con>1e-2");
  _tdiag->Project("tt2dp","tshmc._dist:_ht-mcmid.t0%1695","fit.status>0&&_active&&fit.con>1e-2");
  _tdiag->Project("rt2dp","_rdrift:_ht-mcmid.t0%1695","fit.status>0&&_active&&fit.con>1e-2");
  _tdiag->Project("tt2dp2","tshmc._dist:_ht-mcmid.t0%1695","fit.status>0&&_active&&fit.con>1e-2");
  _tdiag->Project("rt2dp2","_rdrift:_ht-mcmid.t0%1695","fit.status>0&&_active&&fit.con>1e-2");

  TCanvas* t2dcan = new TCanvas("t2dcan","t2dcan",1200,800);
  t2dcan->Divide(2,2);
  t2dcan->cd(1);
  rt2d->Draw();
  t2dcan->cd(2);
  tt2d->Draw();
  t2dcan->cd(3);
  rt2dp->Fit("pol1","","",0,43);
  t2dcan->cd(4);
  tt2dp->Fit("pol3","","",0,43);

  TH1F* dresid = new TH1F("dresid","Drift residual RMS vs #Delta t;Hit t - MC t_{0};Drift RMS",100,0,50.0);
  for(unsigned ibin=0;ibin<100;++ibin){
    dresid->Fill(tt2dp2->GetBinCenter(ibin+1),tt2dp2->GetBinError(ibin+1));
  }
  tt2d->FitSlicesY();
  //TH1D *tt2d_1 = (TH1D*)gDirectory->Get("tt2d_1");
  TH1D *tt2d_2 = (TH1D*)gDirectory->Get("tt2d_2");

  TCanvas* t2drcan = new TCanvas("t2drcan","t2drcan",800,800);
  t2drcan->Divide(2,1);
  t2drcan->cd(1);
  dresid->Draw();
  t2drcan->cd(2);
  tt2d_2->Draw();


}

void KalFit::Trk () {
  TCanvas* tcan = new TCanvas("tt0can","trk_t0can",1200,800);
  TH1F* t0res = new TH1F("t0res","Final t0 resolution;nsec",100,-10,10);
  TH1F* t0pull = new TH1F("t0pull","Track t0 pull",100,-10,10);
  TH2F* dt0 = new TH2F("dt0","Track t0;true t0 (nsec);Initial t0 (nsec)",
      100,500,4000,100,500,4000);
  _tdiag->Project("t0res","fit.t0-mcmid.t0%1695","fit.status>0");
  _tdiag->Project("t0pull","(fit.t0-mcmid.t0%1695)/fit.t0err","fit.status>0");
  tcan->Clear();
  tcan->Divide(2,2);
  tcan->cd(1);
  _tdiag->Draw("fit.t0:mcmid.t0%1695>>dt0","fit.status>0");
  tcan->cd(3);
  t0res->Fit("gaus");
  tcan->cd(4);
  t0pull->Fit("gaus");


  TCanvas* pcan = new TCanvas("pullcan","pullcan",1200,800);
  TH1F* d0pull = new TH1F("d0pull","d0 pull",100,-10,10);
  TH1F* p0pull = new TH1F("p0pull","#phi0 pull",100,-10,10);
  TH1F* ompull = new TH1F("ompull","#omega pull",100,-10,10);
  TH1F* z0pull = new TH1F("z0pull","z0 pull",100,-10,10);
  TH1F* tdpull = new TH1F("tdpull","tan(#lambda) pull",100,-10,10);
  _tdiag->Project("d0pull","(fit.d0-mcent.d0)/fit.d0err","fit.status>0");
  _tdiag->Project("p0pull","(fit.p0-mcent.p0)/fit.p0err","fit.status>0");
  _tdiag->Project("ompull","(fit.om-mcent.om)/fit.omerr","fit.status>0");
  _tdiag->Project("z0pull","(fit.z0-mcent.z0)/fit.z0err","fit.status>0");
  _tdiag->Project("tdpull","(fit.td-mcent.td)/fit.tderr","fit.status>0");
  pcan->Clear();
  pcan->Divide(3,2);
  pcan->cd(1);
  d0pull->Fit("gaus");
  pcan->cd(2);
  p0pull->Fit("gaus");
  pcan->cd(3);
  ompull->Fit("gaus");
  pcan->cd(4);
  z0pull->Fit("gaus");
  pcan->cd(5);
  tdpull->Fit("gaus");

  TCanvas* fcan = new TCanvas("fitcan","fitcan",1200,800);
//  TH2F* mom = new TH2F("mom","momentum at first hit;true mom (MeV/c);reco mom (MeV/c)",
//    100,80,110,100,80,110);
  TH1F* mres = new TH1F("mres","momentum resolution at first hit;MeV/c",100,-2,2);
  TH1F* mpull = new TH1F("mpull","momentum pull at first hit",100,-10,10);
  TH1F* chisq = new TH1F("chisq","Chisq/NDof",100,0,10);
  _tdiag->Project("mres","fit.mom-mcent.mom","fit.status>0");
  _tdiag->Project("mpull","(fit.mom-mcent.mom)/fit.momerr","fit.status>0");
  _tdiag->Project("chisq","fit.chisq/fit.ndof","fit.status>0");
  fcan->Clear();
  fcan->Divide(2,2);
  fcan->cd(1);
  _tdiag->Draw("fit.mom:mcent.mom>>mom","fit.status>0");
  fcan->cd(2);
  mres->Draw();
  fcan->cd(3);
  mpull->Fit("gaus");
  fcan->cd(4);
  chisq->Draw();
}

void KalFit::AccPlots() {

  TH1F* nmc = new TH1F("nmc","N Straw Hits from CE;N straws",81,-0.5,80.5);
  TH1F* mcmom = new TH1F("mcmom","CE true momentum at tracker;CE momentum (MeV/c)",57,49,106);

  TH1F* fitcon = new TH1F("fitcon","log_{10} fit consistency",101,-8,0);
  TH1F* momerr = new TH1F("momerr","Fit momentum error;momentum error (MeV/c)",100,0,0.5);
  TH1F* t0err = new TH1F("t0err","Fit t_{0} error; t_{0} error (ns)",100,0,2.0);
  TH1F* na = new TH1F("na","Fit N active hits;N hits",71,-0.5,70.5);

  TH1F* t0 = new TH1F("t0","Track Fit t_{0};t_{0} (ns)",100,300,1800);
  TH1F* td = new TH1F("td","Track Fit tan(#lambda);tan(#lambda)",100,0.5,1.5);
  TH1F* d0 = new TH1F("d0","Track fit d_{0};d_{0} (mm)",100,-150,150);
  TH1F* rmax = new TH1F("rmax","Track fit rmax;d_{0}+2/#omega (mm)",100,300,900);

  TH1F* fitmom = new TH1F("fitmom","Track fit momentum;fit momentum (MeV/c)",100,98,107);

  _tdiag->Project("nmc","mc.ngood");
  _tdiag->Project("mcmom","mcent.mom",tnhits);

  _tdiag->Project("fitcon","log10(fit.con)",reco+tnhits+tmom);
  _tdiag->Project("momerr","fit.momerr",reco+tnhits+tmom);
  _tdiag->Project("t0err","t0err",reco+tnhits+tmom);
  _tdiag->Project("na","nactive",reco+tnhits+tmom);

  _tdiag->Project("t0","t0",reco+tnhits+tmom+goodfit[icut]);
  _tdiag->Project("td","td",reco+tnhits+tmom+goodfit[icut]+livegate);

  _tdiag->Project("d0","d0",reco+tnhits+tmom+goodfit[icut]+livegate);
  _tdiag->Project("rmax","d0+2.0/om",reco+tnhits+tmom+goodfit[icut]+livegate);

  _tdiag->Project("fitmom","fit.mom",reco+tnhits+tmom+goodfit[icut]+livegate+cosmic);

  TCanvas* pcan = new TCanvas("pcan","Pre-acceptance",1200,800);
  pcan->Clear();
  pcan->Divide(1,2);
  pcan->cd(1);
  gPad->SetLogy();
  nmc->Draw();
  TLine* nmccut = new TLine(20,0.0,20,nmc->GetMaximum());
  nmccut->SetLineColor(kBlack);
  nmccut->SetLineStyle(2);
  nmccut->SetLineWidth(2);
  nmccut->Draw();

  pcan->cd(2);
  gPad->SetLogy();
  mcmom->Draw();
  TLine* tmomcut = new TLine(100,0.0,100,mcmom->GetMaximum());
  tmomcut->SetLineColor(kBlack);
  tmomcut->SetLineStyle(2);
  tmomcut->SetLineWidth(2);
  tmomcut->Draw();

  TCanvas* fcan = new TCanvas("fcan","Fit Acceptance",1200,800);
  fcan->Clear();
  fcan->Divide(2,2);
  fcan->cd(1);
  fitcon->Draw();
  TLine* fitconcut = new TLine(log10(minfitcon[icut]),0.0,log10(minfitcon[icut]),fitcon->GetMaximum());
  fitconcut->SetLineColor(kBlack);
  fitconcut->SetLineStyle(2);
  fitconcut->SetLineWidth(2);
  fitconcut->Draw();

  fcan->cd(2);
  na->Draw();
  TLine* nacut = new TLine(minnactive[icut],0.0,minnactive[icut],na->GetMaximum());
  nacut->SetLineColor(kBlack);
  nacut->SetLineStyle(2);
  nacut->SetLineWidth(2);
  nacut->Draw();

  fcan->cd(3);
  t0err->Draw();
  TLine* t0errcut = new TLine(maxt0err[icut],0.0,maxt0err[icut],t0err->GetMaximum());
  t0errcut->SetLineColor(kBlack);
  t0errcut->SetLineStyle(2);
  t0errcut->SetLineWidth(2);
  t0errcut->Draw();

  fcan->cd(4);
  momerr->Draw();
  TLine* momerrcut = new TLine(maxmomerr[icut],0.0,maxmomerr[icut],momerr->GetMaximum());
  momerrcut->SetLineColor(kBlack);
  momerrcut->SetLineStyle(2);
  momerrcut->SetLineWidth(2);
  momerrcut->Draw();

  TCanvas* tcan = new TCanvas("tcan","Track parameter acceptance",1200,800);
  tcan->Divide(2,2);
  tcan->cd(1);
  t0->Draw();
  TLine* t0cut = new TLine(t0min,0.0,t0min,t0->GetMaximum());
  t0cut->SetLineColor(kBlack);
  t0cut->SetLineStyle(2);
  t0cut->SetLineWidth(2);
  t0cut->Draw();

  tcan->cd(2);
  td->Draw();
  TLine* tdcut_l = new TLine(tdlow,0.0,tdlow,td->GetMaximum());
  tdcut_l->SetLineColor(kBlack);
  tdcut_l->SetLineStyle(2);
  tdcut_l->SetLineWidth(2);
  tdcut_l->Draw();
  TLine* tdcut_h = new TLine(tdhigh,0.0,tdhigh,td->GetMaximum());
  tdcut_h->SetLineColor(kBlack);
  tdcut_h->SetLineStyle(2);
  tdcut_h->SetLineWidth(2);
  tdcut_h->Draw();

  tcan->cd(3);
  d0->Draw();
  TLine* d0cut = new TLine(105,0.0,105,d0->GetMaximum());
  d0cut->SetLineColor(kBlack);
  d0cut->SetLineStyle(2);
  d0cut->SetLineWidth(2);
  d0cut->Draw();

  tcan->cd(4);
  rmax->Draw();
  TLine* rmaxcut = new TLine(660,0.0,660,rmax->GetMaximum());
  rmaxcut->SetLineColor(kBlack);
  rmaxcut->SetLineStyle(2);
  rmaxcut->SetLineWidth(2);
  rmaxcut->Draw();

  TCanvas* mcan = new TCanvas("mcan","momentum",800,600);
  mcan->Divide(1,1);
  mcan->cd(1);
  fitmom->Draw();
  TLine* fitmomcut_l = new TLine(momlow,0.0,momlow,fitmom->GetMaximum());
  fitmomcut_l->SetLineColor(kBlack);
  fitmomcut_l->SetLineStyle(2);
  fitmomcut_l->SetLineWidth(2);
  fitmomcut_l->Draw();
  TLine* fitmomcut_h = new TLine(momhigh,0.0,momhigh,fitmom->GetMaximum());
  fitmomcut_h->SetLineColor(kBlack);
  fitmomcut_h->SetLineStyle(2);
  fitmomcut_h->SetLineWidth(2);
  fitmomcut_h->Draw();


}

void KalFit::Acc(int ngen) {
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
  _tdiag->Project("acc",binnames[ibin++],"evtwt");
  _tdiag->Project("+acc",binnames[ibin++],"evtwt"*mcsel);
  _tdiag->Project("+acc",binnames[ibin++],"evtwt"*(mcsel+reco));
  _tdiag->Project("+acc",binnames[ibin++],"evtwt"*(mcsel+reco+goodfit[icut]));
  _tdiag->Project("+acc",binnames[ibin++],"evtwt"*(mcsel+reco+goodfit[icut]+livegate));
  _tdiag->Project("+acc",binnames[ibin++],"evtwt"*(mcsel+reco+goodfit[icut]+livegate+rpitch));
  _tdiag->Project("+acc",binnames[ibin++],"evtwt"*(mcsel+reco+goodfit[icut]+livegate+rpitch+cosmic));
  _tdiag->Project("+acc",binnames[ibin++],"evtwt"*(mcsel+reco+goodfit[icut]+livegate+rpitch+cosmic+rmom));

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
  acan = new TCanvas("acan","Acceptance",1200,800);
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

void KalFit::Res(unsigned mincut,unsigned maxcut) {
//  TF1* sgau = new TF1("sgau",splitgaus,-1.5,1.5,7);
//  sgau->SetParName(0,"Norm");
//  sgau->SetParName(1,"Mean");
//  sgau->SetParName(2,"SigH");
//  sgau->SetParName(3,"SigL");
//  sgau->SetParName(4,"TFH");
//  sgau->SetParName(5,"TSigH");
//  sgau->SetParName(6,"TSigL");

//  TF1* degau = new TF1("degau",doublegausexp,-1.5,1.5,8);
//  TF1* degau = new TF1("degau",doublegausexp,-1.5,1.5,7);
//  degau->SetParName(0,"Norm");
//  degau->SetParName(1,"GTailFrac");
//  degau->SetParName(2,"ETailFrac");
//  degau->SetParName(3,"Mean");
//  degau->SetParName(4,"CoreSig");
//  degau->SetParName(5,"GTailSig");
//  degau->SetParName(6,"ETailLambda");
//  degau->SetParName(7,"ETailPower");

  TF1* cball = new TF1("cball",crystalball,-2.0,2.5,7);
  cball->SetParName(0,"Norm");
  cball->SetParName(1,"x0");
  cball->SetParName(2,"sigma");
  cball->SetParName(3,"n");
  cball->SetParName(4,"alpha");
  cball->SetParName(5,"tailfrac");
  cball->SetParName(6,"taillambda");

  TF1* dscb = new TF1("dscb",fnc_dscb,-2.0,2.5,7);
  dscb->SetParName(0,"Norm");
  dscb->SetParName(1,"x0");
  dscb->SetParName(2,"sigma");
  dscb->SetParName(3,"ANeg");
  dscb->SetParName(4,"PNeg");
  dscb->SetParName(5,"APos");
  dscb->SetParName(6,"PPos");

  TH1F* momres[4];
  TF1*  fitmomres[4];
  TH1F* effnorm = new TH1F("effnorm","effnorm",100,0,150);
//  _tdiag->Project("effnorm","mcent.mom");
  _tdiag->Project("effnorm","mcent.mom","evtwt"*mcsel);

  rcan = new TCanvas("rcan","Momentum Resolution",1200,800);
  rcan->Clear();
  unsigned ncan = maxcut-mincut+1;
  if(ncan==1)
    rcan->Divide(1,1);
  else if(ncan==2)
    rcan->Divide(2,1);
  else
    rcan->Divide(2,2);
  gStyle->SetOptFit(111111);
  gStyle->SetOptStat("oumr");
  for(unsigned ires=mincut;ires<maxcut+1;ires++){
    rcan->cd(ires+1);
    gPad->SetLogy();
    char mname[50];
    char fitname[50];
    snprintf(mname,50,"momres%i",ires);
    snprintf(fitname,50,"fitmomres%i",ires);
    momres[ires] = new TH1F(mname,"momentum resolution at start of tracker;MeV/c",251,-4,4);
    momres[ires]->Sumw2();
//  momres[ires]->SetStats(0);
    TCut final = reco+goodfit[ires]+rpitch+cosmic+livegate+rmomloose;
    _tdiag->Project(mname,"fit.mom-mcent.mom","evtwt"*final);
//    _tdiag->Project(mname,"fit.mom-mcent.mom",final);
    double integral = momres[ires]->GetEntries()*momres[ires]->GetBinWidth(1);
    cout << "Integral = " << integral << " mean = " << momres[ires]->GetMean() << " rms = " << momres[ires]->GetRMS() << endl;
    cball->SetParameters(3*integral,momres[ires]->GetMean()+0.07,0.3*momres[ires]->GetRMS(),3.0,1.0,0.02,0.2);
    cball->SetParLimits(5,0.001,0.4);
    cball->SetParLimits(6,0.1,momres[ires]->GetRMS());

    dscb->SetParameters(3*integral,momres[ires]->GetMean()+0.07,0.3*momres[ires]->GetRMS(),0.9,3.5,1.5,6.0);

    momres[ires]->SetMinimum(0.5);
//    momres[ires]->Fit("cball","LRQ");
//    momres[ires]->Fit("cball","LRM");
    momres[ires]->Fit("dscb","LRQ");
    momres[ires]->Fit("dscb","LRM");
//    fitmomres[ires] = new TF1(*cball);
    fitmomres[ires] = new TF1(*dscb);
    fitmomres[ires]->SetName(fitname);
    gDirectory->Append(fitmomres[ires]);

    TLine* zero = new TLine(0.0,0.0,0.0,momres[ires]->GetBinContent(momres[ires]->GetMaximumBin()));
    zero->SetLineStyle(2);
    zero->Draw();

    double keff = momres[ires]->GetEntries()/effnorm->GetEntries();
//    TPaveText* ttext = new TPaveText(0.1,0.75,0.4,0.9,"NDC");
//    ttext->AddText("Truth Cuts");
//    ttext->AddText(tnhits.GetTitle());
//    ttext->AddText(tmom.GetTitle());
//    ttext->AddText(tpitch.GetTitle());
//    ttext->Draw();

    TPaveText* rtext = new TPaveText(0.1,0.5,0.4,0.9,"NDC");
    rtext->AddText("Reco Cuts");
    char line[40];
    snprintf(line,80,"%4.3f<tan(#lambda)<%4.3f",tdlow,tdhigh);
    rtext->AddText(line);
    snprintf(line,80,"t0>%5.1f nsec",t0min);
    rtext->AddText(line);
    sprintf(line,"%s",goodfit[ires].GetTitle());
    rtext->AddText(line);
//    sprintf(line,"%s",ncuts[ires].GetTitle());
//    rtext->AddText(line);
//    sprintf(line,"%s",t0cuts[ires].GetTitle());
//    rtext->AddText(line);
//    sprintf(line,"%s",momcuts[ires].GetTitle());
//    rtext->AddText(line);
//    sprintf(line,"%s",fitcuts[ires].GetTitle());
//    rtext->AddText(line);
//    sprintf(line,"%s",rmomloose.GetTitle());
//    rtext->AddText(line);
    sprintf(line,"Eff=%4.4f",keff);
    rtext->AddText(line);
    rtext->Draw();

  }
  rcan->cd(0);
}

void KalFit::Res2(int ires) {
  TF1* cball2 = new TF1("cball2",crystalball,-5.0,2.0,7);
  cball2->SetParName(0,"Norm");
  cball2->SetParName(1,"x0");
  cball2->SetParName(2,"sigma");
  cball2->SetParName(3,"n");
  cball2->SetParName(4,"alpha");
  cball2->SetParName(5,"tailfrac");
  cball2->SetParName(6,"taillambda");

  unsigned nbins(250);
  TH1F* momres = new TH1F("momres","Tracker Momentum Resolution;P_{RECO}-P_{TRUE} (MeV/c);Arbitrary Units",nbins,-4.0,2.0);
  TCut final = goodfit[ires]+mcsel;
  _tdiag->Project("momres","fit.mom-mcent.mom",final);
  double integral = momres->GetEntries()*momres->GetBinWidth(1);
  cout << "Integral = " << integral << " mean = " << momres->GetMean() << " rms = " << momres->GetRMS() << endl;
  cball2->SetParameters(3*integral,momres->GetMean()+0.07,0.3*momres->GetRMS(),3.0,1.0,0.02,0.2);
  cball2->SetParLimits(5,0.001,0.4);
  cball2->SetParLimits(6,0.1,momres->GetRMS());

  momres->Fit("cball2","RNQ");
  momres->Fit("cball2","LRNQ");

  TH1F* effnorm = new TH1F("effnorm","effnorm",100,0,150);
  _tdiag->Project("effnorm","mcent.mom",mcsel);

  TCanvas* rcan2 = new TCanvas("rcan2","Momentum Resolution",800,800);
  rcan2->Clear();
  rcan2->cd(1);
  gPad->SetLogy();
  momres->SetMinimum(0.1);
  momres->SetStats(0);
  momres->Draw();
  cball2->SetNpx(nbins);
  cball2->Draw("same");
  TPaveText * rtext2 = new TPaveText(0.15,0.65,0.55,0.9,"NDC");
  char line[100];
//  snprintf(line,100,"Core Width = %3.1f#pm %2.1f KeV/c",1000.0*cball2->GetParameter(2),1000.0*cball2->GetParError(2));
//  rtext2->AddText(line);
//  snprintf(line,100,"High Tail Slope = %3.1f#pm %2.1f KeV/c",1000.0*cball2->GetParameter(6),1000.0*cball2->GetParError(6));
//  rtext2->AddText(line);
  snprintf(line,100,"Core Width = %3.0f KeV/c",1000.0*cball2->GetParameter(2));
  rtext2->AddText(line);
  snprintf(line,100,"High Tail Slope = %3.0f KeV/c",1000.0*cball2->GetParameter(6));
  rtext2->AddText(line);

  double tailsig(3.0);
  double taillow = cball2->GetParameter(1) + tailsig*cball2->GetParameter(2);
  int ilow = momres->FindBin(taillow)+1;
  double tailint = momres->Integral(ilow,nbins);
  double fullint = momres->GetEntries();
  snprintf(line,100,"High Tail Fraction = %2.1f %%",100*tailint/fullint);
  rtext2->AddText(line);
  snprintf(line,100,"Reco#timesSelection Efficiency = %2.0f %%",100*fullint/effnorm->GetEntries());
  rtext2->AddText(line);
  rtext2->Draw();
}

void KalFit::Ambig(int acut) {
  gStyle->SetOptStat(1111);

  TCut ghit("tshmc._rel==0");
  TCut delta("tshmc._rel>0");
  TCut bkg("tshmc._rel<0");
  TCut gambig("tshmc._ambig==_ambig");
  TCut bambig("tshmc._ambig!=_ambig&&_ambig!=0");
  TCut nambig("tsh._ambig==0");
  TCut active("tsh._active>0");
// apply requested cuts

  TCut goodtrk = (reco+goodfit[acut]+mcsel);

//  TCut goodtrk ="fit.status>0";

  TH1F* rdg = new TH1F("rdg","Hit fraction vs drift radius;true radius (mm);hit fraction",100,0.0,2.7);
  TH1F* rdn = new TH1F("rdn","Hit fraction vs drift radius;true radius (mm);hit fraction",100,0.0,2.7);
  TH1F* rdb = new TH1F("rdb","Hit fraction vs drift radius;true radius (mm);hit fraction",100,0.0,2.7);
  TH1F* rda = new TH1F("rda","Drift radius;true radius (mm)",100,0.0,2.7);
  TH1F* rdd = new TH1F("rdd","Drift radius;true radius (mm)",100,0.0,2.7);
  TH1F* rdf = new TH1F("rdf","Drift radius;true radius (mm)",100,0.0,2.7);
  TH1F* rdi = new TH1F("rdi","Drift radius;true radius (mm)",100,0.0,2.7);
  rdg->SetLineColor(kGreen);
  rdn->SetLineColor(kBlue);
  rdb->SetLineColor(kRed);
  rda->SetLineColor(kBlack);
  rdi->SetLineColor(kCyan);
  rdd->SetLineColor(kOrange);
  rdf->SetLineColor(kYellow);
  rdg->SetStats(0);
  rdn->SetStats(0);
//  rdb->SetStats(0);
  rdi->SetStats(0);
  rda->SetStats(0);
  rdd->SetStats(0);
  rdf->SetStats(0);
  rdg->Sumw2();
  rdn->Sumw2();
  rdb->Sumw2();
  rda->Sumw2();
  rdd->Sumw2();
  rdf->Sumw2();

  _tdiag->Project("rdg","tshmc._dist",goodtrk+active+gambig+ghit);
  _tdiag->Project("rdn","tshmc._dist",goodtrk+active+nambig+ghit);
  _tdiag->Project("rdb","tshmc._dist",goodtrk+active+bambig+ghit);
  _tdiag->Project("rda","tshmc._dist",goodtrk+active);
  _tdiag->Project("rdi","tshmc._dist",goodtrk+ghit+(!active));
  _tdiag->Project("rdd","tshmc._dist",goodtrk+active+delta);
  _tdiag->Project("rdf","tshmc._dist",goodtrk+active+bkg);
  Double_t ntotal = rda->GetEntries();
  Double_t nright = rdg->GetEntries();
  Double_t nneutral = rdn->GetEntries();
  Double_t nwrong = rdb->GetEntries();
  Double_t ndelta = rdd->GetEntries();
  Double_t nbkg = rdf->GetEntries();
  Double_t ninact = rdi->GetEntries();
  TH1F* rdgr = new TH1F(*rdg);
  TH1F* rdnr = new TH1F(*rdn);
  TH1F* rdbr = new TH1F(*rdb);
  rdgr->Divide(rda);
  rdnr->Divide(rda);
  rdbr->Divide(rda);

//  TH1F* frdg = new TH1F("frdg","True Drift radius, failed fits;radius (mm);N hits",100,-0.05,2.55);
//  TH1F* frdb = new TH1F("frdb","True Drift radius, failed fits;radius (mm);N hits",100,-0.05,2.55);
//  frdg->SetLineColor(kBlue);
//  frdb->SetLineColor(kRed);
//  frdg->SetStats(0);
//  frdb->SetStats(0);

//  _tdiag->Project("frdg","tshmc._dist",mcsel+active+ambig+(!goodfit));
//  _tdiag->Project("frdb","tshmc._dist",mcsel+active+(!ambig)+(!goodfit));

  TH1F* momres0 = new TH1F("momres0","Momentum resolution at start of tracker;p_{reco}-p_{true}(MeV/c)",151,-4,4);
  TH1F* momres1 = new TH1F("momres1","Momentum resolution at start of tracker;p_{reco}-p_{true}(MeV/c)",151,-4,4);
  TH1F* momres2 = new TH1F("momres2","Momentum resolution at start of tracker;p_{reco}-p_{true}(MeV/c)",151,-4,4);
  momres0->SetLineColor(kBlack);
  momres1->SetMarkerColor(kCyan);
  momres1->SetMarkerStyle(4);
  momres2->SetMarkerColor(kOrange);
  momres2->SetMarkerStyle(5);
 //  momres->SetStats(0);
  _tdiag->Project("momres0","fit.mom-mcent.mom",goodtrk);
  _tdiag->Project("momres1","fit.mom-mcent.mom",goodtrk&&"fit.status==1");
  _tdiag->Project("momres2","fit.mom-mcent.mom",goodtrk&&"fit.status==2");

  TH1F* afg = new TH1F("afg","Average hit fraction vs momentum resolution;p_{reco}-p_{true}(MeV/c);hit fraction",41,-4,4);
  TH1F* afn = new TH1F("afn","Average hit fraction vs momentum resolution;p_{reco}-p_{true}(MeV/c):hit fraction",41,-4,4);
  TH1F* afb = new TH1F("afb","Average hit fraction vs momentum resolution;p_{reco}-p_{true}(MeV/c);hit fraction",41,-4,4);
  TH1F* afa = new TH1F("afa","Average hit fraction vs momentum resolution;p_{reco}-p_{true}(MeV/c);hit fraction",41,-4,4);
  afg->SetStats(0);
  afn->SetStats(0);
  afb->SetStats(0);
  afa->SetStats(0);
  afg->SetLineColor(kGreen);
  afn->SetLineColor(kBlue);
  afb->SetLineColor(kRed);
  afg->Sumw2();
  afn->Sumw2();
  afb->Sumw2();
  afa->Sumw2();
  _tdiag->Project("afg","fit.mom-mcent.mom",goodtrk+active+gambig);
  _tdiag->Project("afn","fit.mom-mcent.mom",goodtrk+active+nambig);
  _tdiag->Project("afb","fit.mom-mcent.mom",goodtrk+active+bambig);
  _tdiag->Project("afa","fit.mom-mcent.mom",goodtrk+active);
  afg->Divide(afa);
  afn->Divide(afa);
  afb->Divide(afa);
  afg->SetMinimum(0.0);
  afg->SetMaximum(1.1);

  TCanvas* ambigcan = new TCanvas("ambigcan","Hit Ambiguity",1200,800);
  ambigcan->Divide(2,2);


  ambigcan->cd(1);
  rda->Draw();
  rdi->Draw("same");
  rdd->Draw("same");
  rdf->Draw("same");
  TLegend* drleg = new TLegend(0.15,0.15,0.55,0.5);
  char dtitle[100];
  snprintf(dtitle,100,"%4.0f Active hits",ntotal);
  drleg->AddEntry(rda,dtitle,"l");
  snprintf(dtitle,100,"%4.4f Delta-ray hits",ndelta/ntotal);
  drleg->AddEntry(rdd,dtitle,"l");
  snprintf(dtitle,100,"%4.4f Inactive good hits",ninact/ntotal);
  drleg->AddEntry(rdi,dtitle,"l");
  snprintf(dtitle,100,"%4.4f Background hits",nbkg/ntotal);
  drleg->AddEntry(rdf,dtitle,"l");
  drleg->Draw();

  ambigcan->cd(2);

//  TF1* ex = new TF1("ex","[0]*exp(-x/[1])/[1]+[2]*x+[3]",4);
//  ex->SetParameters(0.0,0.1,-0.01,0.04);
//  ex->SetParName(0,"ExpNorm");
//  ex->SetParName(1,"Lambda");
//  ex->SetParName(2,"Slope");
//  ex->SetParName(3,"Intercept");
//  rdbr->SetMaximum(1.1);
//  rdbr->SetMinimum(0.0);
//  rdbr->Fit(ex,"L");
  rdgr->Draw();
  rdnr->Draw("same");
  rdbr->Draw("same");

  TLegend* leg = new TLegend(0.4,0.35,0.9,0.6);
  char title[80];
  snprintf(title,80,"Correct ambiguity %4.3f",nright/ntotal);
  leg->AddEntry(rdgr,title,"l");
  snprintf(title,80,"Null ambiguity %4.3f",nneutral/ntotal);
  leg->AddEntry(rdnr,title,"l");
  snprintf(title,80,"Incorrect ambiguity %4.3f",nwrong/ntotal);
  leg->AddEntry(rdbr,title,"l");
  leg->Draw();

  ambigcan->cd(3);
//  TF1* degau = new TF1("degau",doublegausexp,-1.5,1.5,7);
//  degau->SetParName(0,"Norm");
//  degau->SetParName(1,"GTailFrac");
//  degau->SetParName(2,"ETailFrac");
//  degau->SetParName(3,"Mean");
//  degau->SetParName(4,"CoreSig");
//  degau->SetParName(5,"GTailSig");
//  degau->SetParName(6,"ETailLambda");
  double integral = momres0->GetEntries()*momres0->GetBinWidth(1);
//  degau->SetParameters(integral,0.1,0.2,0.0,0.3*momres0->GetRMS(),2*momres0->GetRMS(),0.5*momres0->GetRMS(),0.25);
//  degau->SetParLimits(1,0.02,0.3);
//  degau->SetParLimits(2,0.02,0.3);
//  degau->SetParLimits(4,0.05,momres0->GetRMS());
//  degau->SetParLimits(5,0.12,2*momres0->GetRMS());
//  degau->SetParLimits(6,0.1,momres0->GetRMS());
  TF1* cball = new TF1("cball",crystalball,-2.0,1.5,7);
  cball->SetParName(0,"Norm");
  cball->SetParName(1,"x0");
  cball->SetParName(2,"sigma");
  cball->SetParName(3,"n");
  cball->SetParName(4,"alpha");
  cball->SetParName(5,"tailfrac");
  cball->SetParName(6,"taillambda");
  cball->SetParameters(integral,0.0,0.1,1.0,1.0,0.05,0.5);
  cball->SetParLimits(5,0.001,0.4);
  cball->SetParLimits(6,0.1,momres0->GetRMS());

  gPad->SetLogy();
  momres0->Fit("cball","LIR");
  momres1->Draw("psame");
  momres2->Draw("psame");
  TLegend* mleg = new TLegend(0.13,0.6,0.43,0.85);
  snprintf(title,80,"%4.0f All Fits",momres0->GetEntries());
  mleg->AddEntry(momres0,title,"l");
  snprintf(title,80,"%4.0f Converged Fits",momres1->GetEntries());
  mleg->AddEntry(momres1,title,"p");
  snprintf(title,80,"%4.0f Unconverged Fits",momres2->GetEntries());
  mleg->AddEntry(momres2,title,"p");
  mleg->Draw();


  ambigcan->cd(4);
  afg->Draw();
  afn->Draw("same");
  afb->Draw("same");
  TLegend* fleg = new TLegend(0.16,0.35,0.625,0.6);
  fleg->AddEntry(afg,"Correct Ambiguity","l");
  fleg->AddEntry(afn,"No Ambiguity","l");
  fleg->AddEntry(afb,"Incorrect Ambiguity","l");
  fleg->Draw();

  ambigcan->cd(0);

}

void KalFit::Resid() {

  TCut delta("tshmc._proc==17");
  TCut primary("tshmc._proc==56");
  TCut gambig("tshmc._ambig==_ambig&&_ambig!=0");
  TCut bambig("tshmc._ambig!=_ambig&&_ambig!=0");
  TCut nambig("_ambig==0");
  TCut active("fit.status==1 && _active>0");

  TH1F* rdg = new TH1F("rdg","True Drift radius;radius (mm);N hits",100,-0.05,2.55);
  TH1F* rdb = new TH1F("rdb","True Drift radius;radius (mm);N hits",100,-0.05,2.55);
  TH1F* rdn = new TH1F("rdn","True Drift radius;radius (mm);N hits",100,-0.05,2.55);
  rdg->SetLineColor(kBlue);
  rdb->SetLineColor(kRed);
  rdn->SetLineColor(kGreen);
  rdg->SetStats(0);
  rdb->SetStats(0);
  rdn->SetStats(0);

  TH1F* rpullg = new TH1F("rpullg","Correct Ambiguity Residual Pull;Pull;N hits",100,-6,6);
  TH1F* rpullb = new TH1F("rpullb","Incorrect Ambiguity Residual Pull;Pull;N hits",100,-6,6);
  TH1F* rpulln = new TH1F("rpulln","No Assigned Ambiguity Residual Pull;Pull;N hits",100,-6,6);
  TH1F* rpulld = new TH1F("rpulld","#delta-ray Residual Pull;Pull;N hits",100,-6,6);
  rpullg->SetLineColor(kBlue);
  rpullb->SetLineColor(kRed);
  rpulln->SetLineColor(kGreen);
  rpulld->SetLineColor(kCyan);

  _tdiag->Project("rdg","tshmc._dist",reco+mcsel+active+gambig+primary);
  _tdiag->Project("rdb","tshmc._dist",reco+mcsel+active+bambig+primary);
  _tdiag->Project("rdn","tshmc._dist",reco+mcsel+active+nambig+primary);

  _tdiag->Project("rpullg","_resid/_residerr",reco+mcsel+active+gambig+primary);
  _tdiag->Project("rpullb","_resid/_residerr",reco+mcsel+active+bambig+primary);
  _tdiag->Project("rpulln","_resid/_residerr",reco+mcsel+active+nambig+primary);
  _tdiag->Project("rpulld","_resid/_residerr",reco+mcsel+active+delta);

  TCanvas* residcan = new TCanvas("residcan","Residuals",1200,800);
  residcan->Divide(2,1);
  residcan->cd(1);
  rdg->Draw();
  rdb->Draw("same");
  rdn->Draw("same");

  TLegend* leg = new TLegend(0.3,0.3,0.8,0.5);
  leg->AddEntry(rpullg,"Correct ambiguity","l");
  leg->AddEntry(rpullb,"Incorrect ambiguity","l");
  leg->AddEntry(rpulln,"No ambiguity assigned","l");
  leg->Draw();

  TPad* ppad = dynamic_cast<TPad*>(residcan->cd(2));
  ppad->Divide(1,4);
  ppad->cd(1);
  rpullg->Fit("gaus");
  ppad->cd(2);
  rpullb->Fit("gaus");
  ppad->cd(3);
  rpulln->Fit("gaus");
  ppad->cd(4);
  rpulld->Fit("gaus");

  residcan->cd(0);

}

void KalFit::Con() {
  TH1F* con1 = new TH1F("con1","#chi^{2} fit consistency",500,0.0,1.0);
  TH1F* con2 = new TH1F("con2","#chi^{2} fit consistency",500,0.0,1.0);
  TH1F* lcon1 = new TH1F("lcon1","log_{10}(#chi^{2}) fit consistency",100,-10,0);
  TH1F* lcon2 = new TH1F("lcon2","log_{10}(#chi^{2}) fit consistency",100,-10,0);
  con1->SetLineColor(kBlue);
  con2->SetLineColor(kRed);
  lcon1->SetLineColor(kBlue);
  lcon2->SetLineColor(kRed);
//  fcon1->SetStats(0);
//  fcon2->SetStats(0);

  _tdiag->Project("con1","fit.con",mcsel+"fit.status==1");
  _tdiag->Project("con2","fit.con",mcsel+"fit.status==2");
  _tdiag->Project("lcon1","log10(fit.con)",mcsel+"fit.status==1");
  _tdiag->Project("lcon2","log10(fit.con)",mcsel+"fit.status==2");

  TCanvas* fcan = new TCanvas("fcan","fit consistency",500,800);
  fcan->Clear();
  fcan->Divide(1,2);
  fcan->cd(1);
  con1->Draw();
  con2->Draw("same");
  fcan->cd(2);
  lcon1->Draw();
  lcon2->Draw("same");

  TLegend* leg = new TLegend(0.1,0.6,0.4,0.8);
  leg->AddEntry(con1,"Fully Converged Fit","l");
  leg->AddEntry(con2,"Unconverged Fit","l");
  leg->Draw();

}

void KalFit::Error(){
    TF1* sgau = new TF1("sgau",splitgaus,-1.,1.,7);
    sgau->SetParName(0,"Norm");
    sgau->SetParName(1,"Mean");
    sgau->SetParName(2,"SigH");
    sgau->SetParName(3,"SigL");
    sgau->SetParName(4,"TFH");
    sgau->SetParName(5,"TSigH");
    sgau->SetParName(6,"TSigL");

    TCut final = (reco+mcsel);
    TH1F* momres1 = new TH1F("momres1","momentum resolution at start of tracker;MeV/c",151,-2.5,2.5);
    TH1F* momres2 = new TH1F("momres2","momentum resolution at start of tracker;MeV/c",151,-2.5,2.5);

    _tdiag->Project("momres1","fit.mom-mcent.mom",final+"fit.momerr<0.15");
    _tdiag->Project("momres2","fit.mom-mcent.mom",final+"fit.momerr>0.2");

    double integral = momres1->GetEntries()*momres1->GetBinWidth(1);
    sgau->SetParameters(integral,0.0,momres1->GetRMS(),momres1->GetRMS(),0.01,2*momres1->GetRMS(),2*momres1->GetRMS());
    sgau->SetParLimits(5,1.0*momres1->GetRMS(),1.0);
    sgau->SetParLimits(6,1.0*momres1->GetRMS(),1.0);
    sgau->SetParLimits(4,0.0,0.8);
    momres1->Fit("sgau","L");

    integral = momres2->GetEntries()*momres1->GetBinWidth(1);
    sgau->SetParameters(integral,0.0,momres2->GetRMS(),momres2->GetRMS(),0.01,2*momres2->GetRMS(),2*momres2->GetRMS());
    sgau->SetParLimits(5,1.0*momres2->GetRMS(),1.0);
    sgau->SetParLimits(6,1.0*momres2->GetRMS(),1.0);
    sgau->SetParLimits(4,0.0,0.8);
    momres2->Fit("sgau","L");

   TCanvas* mecan = new TCanvas("mecan","Momentum Res",800,500);
   mecan->Divide(2,1);
   mecan->cd(1);
   momres1->Draw();
   mecan->cd(2);
   momres2->Draw();

}

void KalFit::Drift(){

  TH1F* rdg = new TH1F("rdg","Reco Hit Drift Radius;radius (mm);N hits",100,-0.05,2.55);
  TH1F* rdb = new TH1F("rdb","Reco Hit Drift Radius;radius (mm);N hits",100,-0.05,2.55);
  rdg->SetStats(0);
  rdb->SetStats(0);

  _tdiag->Project("rdg","_rdrift","fit.status==1&&fit.momerr<0.15&&_active");
  _tdiag->Project("rdb","_rdrift","fit.status==1&&fit.momerr>0.2&&_active");

  TCanvas* dcan = new TCanvas("dcan","Drift Radius",800,500);
  dcan->Divide(2,1);
  dcan->cd(1);
  rdg->Draw();
  dcan->cd(2);
  rdb->Draw();

}

void KalFit::NHits(){
  TH1F* nch = new TH1F("nch","Number of Conversion Electron Tracker Hits;N hits;N Conversion Tracks",100,-0.5,99.5);
  TH1F* ncha = new TH1F("ncha","Number of Conversion Electron Tracker Hits;N hits;N Conversion Tracks",100,-0.5,99.5);
  nch->SetStats(0);
  nch->SetStats(0);
  TCut goodmc("mcmid.t0%1695>550&&mcmid.mom>100");
  _tdiag->Project("nch","mc.ndigigood",tmom+tpitch+goodmc);
  _tdiag->Project("ncha","mc.ngood",tmom+tpitch+goodmc+"fit.status>0");
  nch->SetLineColor(kRed);
  ncha->SetLineColor(kBlue);
  TLegend* leg = new TLegend(0.6,0.7,0.9,0.9);
  leg->AddEntry(nch,"Generated","l");
  leg->AddEntry(ncha,"Reconstructed","l");
  TCanvas* ncan = new TCanvas("ncan","Number of hits",800,600);
  ncan->Clear();
//  nch->SetMaximum(2000);
  nch->Draw();
  ncha->Draw("same");
  leg->Draw();
}

void KalFit::Nactive() {
  gStyle->SetOptStat(111111);
  TCut pure("nactive-mc.nactive==0");
  TCut goodtrk =mcsel+reco+rpitch;

  TH1F* dna = new TH1F("dna","N non-conversion hits active in fit;N_{active}-N_{active,conversion}",10,-0.5,9.5);
  _tdiag->Project("dna","nactive-mc.nactive",goodtrk);

  TH1F* momres0 = new TH1F("momres0","Momentum resolution at start of tracker;p_{reco}-p_{true}(MeV/c)",101,-4,4);
  TH1F* momres1 = new TH1F("momres1","Momentum resolution at start of tracker;p_{reco}-p_{true}(MeV/c)",101,-4,4);
  TH1F* momres2 = new TH1F("momres2","Momentum resolution at start of tracker;p_{reco}-p_{true}(MeV/c)",101,-4,4);
  momres0->SetLineColor(kBlack);
  momres1->SetMarkerColor(kCyan);
  momres1->SetMarkerStyle(4);
  momres2->SetMarkerColor(kOrange);
  momres2->SetMarkerStyle(5);
 //  momres->SetStats(0);
  _tdiag->Project("momres0","fit.mom-mcent.mom",goodtrk);
  _tdiag->Project("momres1","fit.mom-mcent.mom",goodtrk&&pure);
  _tdiag->Project("momres2","fit.mom-mcent.mom",goodtrk&&!pure);

  TCanvas* ncan = new TCanvas("ncan","N hits",1000,800);
  ncan->Divide(1,2);
  ncan->cd(2);
  //gPad->SetLogy();
  dna->Draw();

  ncan->cd(1);
  TF1* degau = new TF1("degau",doublegausexp,-1.5,1.5,7);
  degau->SetParName(0,"Norm");
  degau->SetParName(1,"GTailFrac");
  degau->SetParName(2,"ETailFrac");
  degau->SetParName(3,"Mean");
  degau->SetParName(4,"CoreSig");
  degau->SetParName(5,"GTailSig");
  degau->SetParName(6,"ETailLambda");
  double integral = momres0->GetEntries()*momres0->GetBinWidth(1);
  degau->SetParameters(integral,0.1,0.2,0.0,0.3*momres0->GetRMS(),2*momres0->GetRMS(),0.5*momres0->GetRMS(),0.25);
  degau->SetParLimits(1,0.02,0.3);
  degau->SetParLimits(2,0.02,0.3);
  degau->SetParLimits(4,0.05,momres0->GetRMS());
  degau->SetParLimits(5,0.12,2*momres0->GetRMS());
  degau->SetParLimits(6,0.1,momres0->GetRMS());

  gPad->SetLogy();
  momres0->Fit(degau);
  momres1->Draw("psame");
  momres2->Draw("psame");
  TLegend* mleg = new TLegend(0.13,0.6,0.43,0.85);
  char title[80];
  snprintf(title,80,"%4.0f All Fits",momres0->GetEntries());
  mleg->AddEntry(momres0,title,"l");
  snprintf(title,80,"%4.0f Pure conversion hits",momres1->GetEntries());
  mleg->AddEntry(momres1,title,"p");
  snprintf(title,80,"%4.0f >=1 non-conversion hits",momres2->GetEntries());
  mleg->AddEntry(momres2,title,"p");
  mleg->Draw();

}

void KalFit::TrkQualRes() {
  TH1F* goodf = new TH1F("goodf","Momentum Resolution;Reco - True Momentum (MeV/c)",100,-4,4);
  TH1F* badf = new TH1F("badf","Momentum Resolution;Reco - True Momentum (MeV/c)",100,-4,4);
  goodf->SetLineColor(kGreen);
  goodf->SetLineWidth(2);
  goodf->SetFillColor(kGreen);
  goodf->SetFillStyle(3004);
  badf->SetLineColor(kMagenta);
  badf->SetFillColor(kMagenta);
  badf->SetLineWidth(2);
  badf->SetFillStyle(3005);
  goodf->SetStats(0);
  badf->SetStats(0);
  _tdiag->Project("goodf","fit.mom-mcent.mom",reco+rpitch+"trkqual>0.2");
  _tdiag->Project("badf","fit.mom-mcent.mom",reco+rpitch+"trkqual<0.2");
  badf->Scale(goodf->GetEntries()/badf->GetEntries());
  TCanvas* tqcan = new TCanvas("tqcan","TrkQual",1000,800);
  goodf->Draw();
  badf->Draw("same");

  double gmax = goodf->GetMaximum();
  int gbin1 = goodf->FindFirstBinAbove(0.5*gmax);
  int gbin2 = goodf->FindLastBinAbove(0.5*gmax);
  double gfwhm = goodf->GetBinCenter(gbin2) - goodf->GetBinCenter(gbin1);

  double bmax = badf->GetMaximum();
  int bbin1 = badf->FindFirstBinAbove(0.5*bmax);
  int bbin2 = badf->FindLastBinAbove(0.5*bmax);
  double bfwhm = badf->GetBinCenter(bbin2) - badf->GetBinCenter(bbin1);

  char ltitle[40];

  TLegend* leg = new TLegend(0.55,0.6,0.9,0.9);
  snprintf(ltitle,40,"TrkQual<0.2, FWHM =%3.0f KeV/c",1000.0*gfwhm);
  leg->AddEntry(goodf,ltitle,"F");
  snprintf(ltitle,40,"TrkQual>0.2, FWHM =%3.0f KeV/c",1000.0*bfwhm);
  leg->AddEntry(badf,ltitle,"F");
  leg->Draw();
}

void KalFit::Mom() {
  TH1F* cemomtp = new TH1F("cemomtp","Conversion Electron Momentum;p (MeV/c)",200,98,106);
  TH1F* cemomte = new TH1F("cemomte","Conversion Electron Momentum;p (MeV/c)",200,98,106);
  TH1F* cemomr = new TH1F("cemomr","Conversion Electron Momentum;p (MeV/c)",200,98,106);
  cemomtp->SetStats(0);
  cemomtp->SetLineColor(kGreen);
  cemomtp->SetFillColor(kGreen);
  cemomte->SetStats(0);
  cemomte->SetLineColor(kBlue);
  cemomte->SetFillColor(kBlue);
  cemomte->SetFillStyle(3005);
  cemomte->SetLineWidth(2.0);
  cemomr->SetStats(0);
  cemomr->SetLineColor(kRed);
  cemomr->SetFillColor(kRed);
  cemomr->SetFillStyle(3004);
  cemomr->SetLineWidth(2.0);
  _tdiag->Project("cemomtp","mcgen.mom",goodfit[icut]);
  _tdiag->Project("cemomte","mcent.mom",goodfit[icut]);
  _tdiag->Project("cemomr","fit.mom",goodfit[icut]);
  TCanvas* cemcan = new TCanvas("cemcan","Conversion Momentum",1000,800);
  cemcan->Divide(1,1);
  cemcan->cd(1);
//  cemomte->Draw();
//  cemomr->Draw("same");

  double tmax = cemomte->GetMaximum();
  int bin1 = cemomte->FindFirstBinAbove(0.5*tmax);
  int bin2 = cemomte->FindLastBinAbove(0.5*tmax);
  double tfwhm = cemomte->GetBinCenter(bin2) - cemomte->GetBinCenter(bin1);
  TArrow* tfwhma = new TArrow(cemomte->GetBinCenter(bin1),0.5*tmax,
  cemomte->GetBinCenter(bin2),0.5*tmax,0.02,"<|>");
  tfwhma->SetLineWidth(2);
  tfwhma->SetLineColor(kBlue);
  tfwhma->SetFillColor(kBlue);
  bin1 = cemomte->FindFirstBinAbove(0.9*tmax);
  bin2 = cemomte->FindLastBinAbove(0.9*tmax);
  cemomte->SetMaximum(1.2*tmax);
  cemomte->Fit("gaus","0","same",cemomte->GetBinCenter(bin1),cemomte->GetBinCenter(bin2));
  cemomte->Draw();
  tfwhma->Draw();


  double rmax = cemomr->GetMaximum();
  bin1 = cemomr->FindFirstBinAbove(0.5*rmax);
  bin2 = cemomr->FindLastBinAbove(0.5*rmax);
  double rfwhm = cemomr->GetBinCenter(bin2) - cemomr->GetBinCenter(bin1);
  TArrow* rfwhma = new TArrow(cemomr->GetBinCenter(bin1),0.5*rmax,
  cemomr->GetBinCenter(bin2),0.5*rmax,0.02,"<|>");
  rfwhma->SetLineWidth(2);
  rfwhma->SetLineColor(kRed);
  rfwhma->SetFillColor(kRed);
  bin1 = cemomr->FindFirstBinAbove(0.9*rmax);
  bin2 = cemomr->FindLastBinAbove(0.9*rmax);
  cemomr->Fit("gaus","0","same",cemomr->GetBinCenter(bin1),cemomr->GetBinCenter(bin2));
  cemomr->Draw("same");
  rfwhma->Draw();

  cemomtp->Draw("same");
  double emax = cemomtp->GetBinCenter(cemomtp->GetMaximumBin());
  double tmean = cemomte->GetFunction("gaus")->GetParameter(1);
  double tshift = tmean - emax;
  cout << "True shift = " << tshift << endl;

  TArrow* tfsa = new TArrow(tmean,tmax,emax,tmax,0.02,"<-|");
  tfsa->SetLineColor(kBlue);
  tfsa->SetLineWidth(2);
  tfsa->SetFillColor(kBlue);
  tfsa->Draw();

  double rmean = cemomr->GetFunction("gaus")->GetParameter(1);
  double rshift = rmean - emax;
  cout << "Reco shift = " << rshift << endl;
  TArrow* rfsa = new TArrow(rmean,rmax,emax,rmax,0.02,"<-|");
  rfsa->SetLineColor(kRed);
  rfsa->SetLineWidth(2);
  rfsa->SetFillColor(kRed);
  rfsa->Draw();

  TLegend* cemleg = new TLegend(0.1,0.6,0.6,0.9);
  cemleg->AddEntry(cemomtp,"CE Momentum at Production","F");
  char line[50];
  cemleg->AddEntry(cemomte,"True CE Momentum at Tracker Entrance","F");
  snprintf(line,50,"Shift = %2.0f KeV/c, FWHM = %3.0f KeV/c",1000.0*tshift, 1000.0*tfwhm);

  cemleg->AddEntry(tfsa,line,"L");
  cemleg->AddEntry(cemomr,"Reconstructed CE Momentum","F");
  snprintf(line,50,"Shift = %2.0f KeV/c, FWHM = %3.0f KeV/c",1000.0*rshift, 1000.0*rfwhm);
  cemleg->AddEntry(rfsa,line,"L");
cemleg->Draw();
}

void KalFit::Rad() {
  TH1F* rmax = new TH1F("rmax","Track R_{max}, nominal BField;d_{0}+2/#omega (mm)",100,300,900);
  TH1F* rmaxm = new TH1F("rmaxm","Track R_{max}, B *= 0.95;d_{0}+2/(0.95#times#omega) (mm)",100,300,900);
  rmax->SetStats(0);
  rmaxm->SetStats(0);

  _tdiag->Project("rmax","d0+2.0/om",reco+tnhits+tmom+goodfit[icut]+livegate);
  _tdiag->Project("rmaxm","d0+2.0/(0.95*om)",reco+tnhits+tmom+goodfit[icut]+livegate);


  TCanvas* radcan = new TCanvas("radcan","radcan",800,800);
  radcan->Divide(1,2);
  radcan->cd(1);
  rmax->Draw();
  radcan->cd(2);
  rmaxm->Draw();
  TLine* rmaxcut = new TLine(700,0.0,700,rmaxm->GetMaximum());
  rmaxcut->SetLineColor(kBlack);
  rmaxcut->SetLineStyle(2);
  rmaxcut->SetLineWidth(2);
  rmaxcut->Draw();
}

void KalFit::MomTails(int iwt) {
  TString weight;
  if(iwt==0)
    weight = "max(1.0,5.0*(min(fit.mom-mcent.mom,2.0)))";
  else if(iwt==1)
    weight="max(1.0,exp(2.0*min(fit.mom-mcent.mom,2.0)))";
  else if(iwt==2)
    weight ="max(1.0,pow(2.0*min(fit.mom-mcent.mom,2.0),5.0))";
  else
    weight="1";

  TF1* corewt = new TF1("sigwt","[0]",-0.15,0.1);
  corewt->SetLineColor(kBlue);
  TF1* tailwt_lin = new TF1("tailwt_lin","[0]*max(1.0,5.0*(min(x,2.0)))",0.5,5.0);
  tailwt_lin->SetLineColor(kRed);
  TF1* tailwt_exp = new TF1("tailwt_exp","[0]*max(1.0,exp(2.0*min(x,2.0)))",0.5,5.0);
  tailwt_exp->SetLineColor(kGreen);
  TF1* tailwt_pol = new TF1("tailwt_pol","[0]*max(1.0,pow(2.0*min(x,2.0),5.0))",0.5,5.0);
  tailwt_pol->SetLineColor(kCyan);

  double maxmom(5.0);
  TH1F* momres = new TH1F("momres","Momentum Resolution",200,-1,maxmom);

  _tdiag->Project("momres","fit.mom-mcent.mom",reco+tmom+rpitch);
  momres->SetStats(0);

  TCanvas* mrcan = new TCanvas("mrcan","Momentum Resolution",800,800);
  mrcan->SetLogy();
  momres->Draw();
  TBox* coresel= new TBox(-0.15,0.0,0.1,momres->GetMaximum());
  coresel->SetFillStyle(3609);
  coresel->SetFillColor(kBlue);
  coresel->SetLineColor(kBlue);
  coresel->SetLineStyle(1);
  coresel->SetLineWidth(2);
  coresel->Draw("1");

  TBox* tailsel= new TBox(0.5,0.0,maxmom,momres->GetBinContent(momres->FindBin(0.5)));
  tailsel->SetFillStyle(3644);
  tailsel->SetFillColor(kRed);
  tailsel->SetLineColor(kRed);
  tailsel->SetLineStyle(1);
  tailsel->SetLineWidth(2);
  tailsel->Draw("1");

  corewt->SetParameter(0,0.1*momres->GetMaximum());
  tailwt_lin->SetParameter(0,0.01*momres->GetBinContent(momres->FindBin(0.5)));
  tailwt_exp->SetParameter(0,0.01*momres->GetBinContent(momres->FindBin(0.5)));
  tailwt_pol->SetParameter(0,0.01*momres->GetBinContent(momres->FindBin(0.5)));

  corewt->Draw("same");
  tailwt_lin->Draw("same");
  tailwt_exp->Draw("same");
  tailwt_pol->Draw("same");

  TLegend* mrleg = new TLegend(0.7,0.7,0.9,0.9);;
  mrleg->AddEntry(coresel,"Core Selection","F");
  mrleg->AddEntry(tailsel,"Tail Selection","F");
//  mrleg->AddEntry(corewt,"Core Weight","L");
  mrleg->AddEntry(tailwt_lin,"Linear Tail Weight","L");
  mrleg->AddEntry(tailwt_exp,"Exponential Tail Weight","L");
  mrleg->AddEntry(tailwt_pol,"Polynomial Tail Weight","L");
  mrleg->Draw();

  gStyle->SetOptStat(0);
  TCut core = reco + tmom + rpitch + TCut("fit.mom-mcent.mom<0.1&&fit.mom-mcent.mom>-0.15");
  TCut tail = reco + tmom + rpitch + TCut("fit.mom-mcent.mom>0.5");
  TH1F* cnact = new TH1F("cnact","N Active hits",91,9.5,100.5);
  TH1F* tnact = new TH1F("tnact","N Active hits",91,9.5,100.5);
  cnact->SetLineColor(kBlue);
  tnact->SetLineColor(kRed);

  TH1F* cnactf = new TH1F("cnactf","Active hit fraction",100,0.5,1.01);
  TH1F* tnactf = new TH1F("tnactf","Active hit fraction",100,0.5,1.01);
  cnactf->SetLineColor(kBlue);
  tnactf->SetLineColor(kRed);

  TH1F* cndoubf = new TH1F("cndoubf","Double hit fraction",100,0.5,1.01);
  TH1F* tndoubf = new TH1F("tndoubf","Double hit fraction",100,0.5,1.01);
  cndoubf->SetLineColor(kBlue);
  tndoubf->SetLineColor(kRed);

  TH1F* cfitcon = new TH1F("cfitcon","Log_{10} Fit Consistency",100,-6,0.0);
  TH1F* tfitcon = new TH1F("tfitcon","Log_{10} Fit Consistency",100,-6,0.0);
  cfitcon->SetLineColor(kBlue);
  tfitcon->SetLineColor(kRed);

  TH1F* cmomerr = new TH1F("cmomerr","Mom Err;MeV/c",100,0.0,0.5);
  TH1F* tmomerr = new TH1F("tmomerr","Mom Err;MeV/c",100,0.0,0.5);
  cmomerr->SetLineColor(kBlue);
  tmomerr->SetLineColor(kRed);

  TH1F* ct0err = new TH1F("ct0err","t0 Err;nsec",100,0.0,1.5);
  TH1F* tt0err = new TH1F("tt0err","t0 Err;nsec",100,0.0,1.5);
  ct0err->SetLineColor(kBlue);
  tt0err->SetLineColor(kRed);

  TH1F* cd0 = new TH1F("cd0","d0;mm",100,-150,150);
  TH1F* td0 = new TH1F("td0","d0;mm",100,-150,150);
  cd0->SetLineColor(kBlue);
  td0->SetLineColor(kRed);

  TH1F* crmax = new TH1F("crmax","Rmax;mm",100,400,680);
  TH1F* trmax = new TH1F("trmax","Rmax;mm",100,400,680);
  crmax->SetLineColor(kBlue);
  trmax->SetLineColor(kRed);

  TH1F* ctandip = new TH1F("ctandip","Tan #lambda",100,0.4,1.2);
  TH1F* ttandip = new TH1F("ttandip","Tan #lambda",100,0.4,1.2);
  ctandip->SetLineColor(kBlue);
  ttandip->SetLineColor(kRed);

  TH1F* ctrkqual = new TH1F("ctrkqual","TrkQual",100,-0.1,1.2);
  TH1F* ttrkqual = new TH1F("ttrkqual","TrkQual",100,-0.1,1.2);
  ctrkqual->SetLineColor(kBlue);
  ttrkqual->SetLineColor(kRed);

// MC variables

  TH1F* cmcambigf = new TH1F("cmcambigf","Correct Ambiguity Fraction",100,0.5,1.01);
  TH1F* tmcambigf = new TH1F("tmcambigf","Correct Ambiguity Fraction",100,0.5,1.01);
  cmcambigf->SetLineColor(kBlue);
  tmcambigf->SetLineColor(kRed);

  TH1F* cmcdp = new TH1F("cmcdp","MC Tracker Momentum Change",100,0.0,3.0);
  TH1F* tmcdp = new TH1F("tmcdp","MC Tracker Momentum Change",100,0.0,3.0);
  cmcdp->SetLineColor(kBlue);
  tmcdp->SetLineColor(kRed);

  TH1F* cmcdt0 = new TH1F("cmcdt0","MC t0 Res;nsec",100,-3.0,6.0);
  TH1F* tmcdt0 = new TH1F("tmcdt0","MC t0 Res;nsec",100,-3.0,6.0);
  cmcdt0->SetLineColor(kBlue);
  tmcdt0->SetLineColor(kRed);

  TH1F* cmchitfrac = new TH1F("cmchitfrac","MC Hit Fraction",100,0.5,1.01);
  TH1F* tmchitfrac = new TH1F("tmchitfrac","MC Hit Fraction",100,0.5,1.01);
  cmchitfrac->SetLineColor(kBlue);
  tmchitfrac->SetLineColor(kRed);

  _tdiag->Project("cnact","nactive",core*weight);
  _tdiag->Project("tnact","nactive",tail*weight);

  _tdiag->Project("cfitcon","log10(fit.con)",core*weight);
  _tdiag->Project("tfitcon","log10(fit.con)",tail*weight);

  _tdiag->Project("cmomerr","fit.momerr",core*weight);
  _tdiag->Project("tmomerr","fit.momerr",tail*weight);

  _tdiag->Project("ct0err","t0err",core*weight);
  _tdiag->Project("tt0err","t0err",tail*weight);

  _tdiag->Project("cd0","d0",core*weight);
  _tdiag->Project("td0","d0",tail*weight);

  _tdiag->Project("crmax","d0+2/om",core*weight);
  _tdiag->Project("trmax","d0+2/om",tail*weight);

  _tdiag->Project("ctandip","td",core*weight);
  _tdiag->Project("ttandip","td",tail*weight);

  _tdiag->Project("cnactf","nactive/nhits",core*weight);
  _tdiag->Project("tnactf","nactive/nhits",tail*weight);

  _tdiag->Project("cndoubf","ndactive/nactive",core*weight);
  _tdiag->Project("tndoubf","ndactive/nactive",tail*weight);

  _tdiag->Project("ctrkqual","trkqual",core*weight);
  _tdiag->Project("ttrkqual","trkqual",tail*weight);

// MC variables

  _tdiag->Project("cmcambigf","mc.nambig/mc.nactive",core*weight);
  _tdiag->Project("tmcambigf","mc.nambig/mc.nactive",tail*weight);

  _tdiag->Project("cmcdp","mcent.mom-mcxit.mom",core*weight);
  _tdiag->Project("tmcdp","mcent.mom-mcxit.mom",tail*weight);

  _tdiag->Project("cmcdt0","t0-mcmid.t0%1695",core*weight);
  _tdiag->Project("tmcdt0","t0-mcmid.t0%1695",tail*weight);

  _tdiag->Project("cmchitfrac","mc.ngood/mc.ndigigood",core*weight);
  _tdiag->Project("tmchitfrac","mc.ngood/mc.ndigigood",tail*weight);


  double factor = cnact->Integral()/tnact->Integral();
  tnact->Scale(factor);
  tnactf->Scale(factor);
  tndoubf->Scale(factor);
  tfitcon->Scale(factor);
  tmomerr->Scale(factor);
  tt0err->Scale(factor);
  td0->Scale(factor);
  trmax->Scale(factor);
  ttandip->Scale(factor);
  ttrkqual->Scale(factor);
  tmcambigf->Scale(factor);
  tmcdp->Scale(factor);
  tmcdt0->Scale(factor);
  tmchitfrac->Scale(factor);

  TLegend* leg = new TLegend(0.5,0.7,0.9,0.9);
  leg->AddEntry(cnact,"Res. Core","L");
  leg->AddEntry(tnact,"Res. Tail (scaled)","L");

  TCanvas* mtcan1 = new TCanvas("mtcan1","Mom res tail",1000,800);
  mtcan1->Divide(3,3);
  unsigned ipad(1);
  mtcan1->cd(ipad++);
  tnact->Draw();
  cnact->Draw("same");
  leg->Draw();
  mtcan1->cd(ipad++);
  cfitcon->Draw();
  tfitcon->Draw("same");
  mtcan1->cd(ipad++);
  cmomerr->Draw();
  tmomerr->Draw("same");
  mtcan1->cd(ipad++);
  ct0err->Draw();
  tt0err->Draw("same");
  mtcan1->cd(ipad++);
  td0->Draw();
  cd0->Draw("same");
  mtcan1->cd(ipad++);
  trmax->Draw();
  crmax->Draw("same");
  mtcan1->cd(ipad++);
  cnactf->Draw();
  tnactf->Draw("same");
  mtcan1->cd(ipad++);
  tndoubf->Draw();
  cndoubf->Draw("same");
  mtcan1->cd(ipad++);
  ttrkqual->Draw();
  ctrkqual->Draw("same");

  TCanvas* mtcan2 = new TCanvas("mtcan2","Mom res tail",1000,800);
  mtcan2->Divide(2,2);
  ipad = 1;
  mtcan2->cd(ipad++);
  cmcambigf->Draw();
  tmcambigf->Draw("same");
  mtcan2->cd(ipad++);
  tmcdp->Draw();
  cmcdp->Draw("same");
  leg->Draw();
  mtcan2->cd(ipad++);
  cmcdt0->Draw();
  tmcdt0->Draw("same");
  mtcan2->cd(ipad++);
  cmchitfrac->Draw();
  tmchitfrac->Draw("same");

}

void KalFit::MomTailsHits() {
  TCut core = reco + TCut("abs(fit.mom-mcent.mom)<0.2");
  TCut tail = reco + TCut("fit.mom-mcent.mom>0.5");
  TH1F* cmchambig = new TH1F("cmchambig","Hit Ambiguity reco*true",3,-1.5,1.5);
  TH1F* tmchambig = new TH1F("tmchambig","Hit Ambiguity reco*true",3,-1.5,1.5);
  cmchambig->SetLineColor(kBlue);
  tmchambig->SetLineColor(kRed);

  TH1F* crdrift = new TH1F("crdrift","Reco Drift Radius",100,-0.001,2.5);
  TH1F* trdrift = new TH1F("trdrift","Reco Drift Radius",100,-0.001,2.5);
  crdrift->SetLineColor(kBlue);
  trdrift->SetLineColor(kRed);

  TH1F* cnpanel = new TH1F("cnpanel","N panel",5,-0.5,4.5);
  TH1F* tnpanel = new TH1F("tnpanel","N panel",5,-0.5,4.5);
  cnpanel->SetLineColor(kBlue);
  tnpanel->SetLineColor(kRed);

  TH1F* ctdrift = new TH1F("ctdrift","MC Drift Radius",100,-0.001,2.5);
  TH1F* ttdrift = new TH1F("ttdrift","MC Drift Radius",100,-0.001,2.5);
  ctdrift->SetLineColor(kBlue);
  ttdrift->SetLineColor(kRed);

  _tdiag->Project("cnpanel","_npanel",core+"_active");
  _tdiag->Project("tnpanel","_npanel",tail+"_active");

  _tdiag->Project("cmchambig","_ambig*tshmc._ambig",core+"_active");
  _tdiag->Project("tmchambig","_ambig*tshmc._ambig",tail+"_active");

  _tdiag->Project("crdrift","_rdrift",core+"_active");
  _tdiag->Project("trdrift","_rdrift",tail+"_active");

  _tdiag->Project("ctdrift","tshmc._dist",core+"_active");
  _tdiag->Project("ttdrift","tshmc._dist",tail+"_active");

  double factor=0.5*cnpanel->GetEntries()/tnpanel->GetEntries();
  trdrift->Scale(factor);
  tnpanel->Scale(factor);
  tmchambig->Scale(factor);
  ttdrift->Scale(factor);

  TCanvas* mthcan = new TCanvas("mthcan","Mom res tail",1000,800);
  mthcan->Divide(2,2);
  unsigned ipad = 1;
  mthcan->cd(ipad++);
  crdrift->Draw();
  trdrift->Draw("same");
  mthcan->cd(ipad++);
  cmchambig->Draw();
  tmchambig->Draw("same");
  mthcan->cd(ipad++);
  cnpanel->Draw();
  tnpanel->Draw("same");
  mthcan->cd(ipad++);
  ctdrift->Draw();
  ttdrift->Draw("same");

}

void
KalFit::StrawMat() {
  TH1F* naddmat = new TH1F("naddmat","N Added Straws",30,-0.5,29.5);
  TH2F* matvshit = new TH2F("matvshit","N Straw vs N Hits",100,-0.5,99.5,100,-0.5,99.5);
  TH1F* addmatfrac = new TH1F("addmatfrac","Fraction of Added Straws",100,-0.01,0.5);
  addmatfrac->SetStats(0);
  matvshit->SetStats(0);

  TH1F* lofracres = new TH1F("lofracres","Momentum Resolution;Reco - True Mom. (MeV/c)",100,-2,2);
  TH1F* hifracres = new TH1F("hifracres","Momentum Resolution;Reco - True Mom. (MeV/c)",100,-2,2);
  lofracres->SetStats(0);
  hifracres->SetStats(0);
  lofracres->SetLineColor(kRed);
  hifracres->SetLineColor(kBlack);

  TH1F* hitdoca = new TH1F("hitdoca","DOCA to Wire;DOCA (mm)",100,-0.05,2.65);
  TH1F* adddoca = new TH1F("adddoca","DOCA to Wire;DOCA (mm)",100,-0.05,2.65);
  hitdoca->SetStats(0);
  hitdoca->SetLineColor(kRed);
  adddoca->SetStats(0);
  adddoca->SetLineColor(kBlue);

  TH1F* hitstraw = new TH1F("hitstraw","Straw Number;straw #",100,-0.5,99.5);
  TH1F* addstraw = new TH1F("addstraw","Straw number;straw #",100,-0.5,99.5);
  hitstraw->SetStats(0);
  hitstraw->SetLineColor(kRed);
  addstraw->SetStats(0);
  addstraw->SetLineColor(kBlue);

  _tdiag->Project("addmatfrac","(fit.nmatactive-fit.nactive)/fit.nmatactive","fit.status>0");
  _tdiag->Project("matvshit","fit.nmatactive:fit.nactive","fit.status>0");
  _tdiag->Project("naddmat","fit.nmatactive-fit.nactive","fit.status>0");
  _tdiag->Project("adddoca","tm._doca","fit.status>0&&tm._active&&(!tm._thita)");
  _tdiag->Project("hitdoca","tm._doca","fit.status>0&&tm._thita");
  _tdiag->Project("addstraw","tm._straw","fit.status>0&&tm._active&&(!tm._thita)");
  _tdiag->Project("hitstraw","tm._straw","fit.status>0&&tm._thita");

  _tdiag->Project("lofracres","fit.mom-mcent.mom","fit.status>0&&(fit.nmatactive-fit.nactive)/fit.nmatactive<0.1");
  _tdiag->Project("hifracres","fit.mom-mcent.mom","fit.status>0&&(fit.nmatactive-fit.nactive)/fit.nmatactive>0.1");

  TLegend* leg = new TLegend(0.6,0.7,0.9,.9);
  leg->AddEntry(hitdoca,"Hit Straw","L");
  leg->AddEntry(adddoca,"Added Straw","L");

  TCanvas* mcan = new TCanvas("mcan","mcan",1200,800);
  mcan->Divide(3,2);
  mcan->cd(1);
  matvshit->Draw("colorz");
  mcan->cd(2);
  naddmat->Draw();
  mcan->cd(3);
  addmatfrac->Draw();
  mcan->cd(4);
  hitdoca->Draw();
  adddoca->Draw("same");
  leg->Draw();
  mcan->cd(5);
  hitstraw->Draw();
  addstraw->Draw("same");
  leg->Draw();
  mcan->cd(6);
  lofracres->Draw();
  hifracres->Draw("same");
  TLegend* mleg = new TLegend(0.6,0.7,0.9,0.9);
  mleg->AddEntry(lofracres,"N_{added}/N<0.1","L");
  mleg->AddEntry(hifracres,"N_{added}/N>0.1","L");
  mleg->Draw();
}
