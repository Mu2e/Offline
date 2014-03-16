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

// basic parameters

double tdlow(0.57735027);
double tdhigh(1.0);
double t0min(710);
double momlow(103.35);
double momhigh(104.75);
int minnhits(20);
size_t icut=2;
unsigned minnactive[4] = {20,22,25,30};
double maxt0err[4] = {1.5,0.95,0.9,0.8};
double maxmomerr[4] = {0.3,0.28,0.25,0.22};
double minfitcon[4] = {1e-6,1e-3,2e-3,1e-2};

TCut ncuts[4], t0cuts[4], momcuts[4], fitcuts[4];
TCut reco,goodfit,cosmic,rmom,rpitch,livegate;
TCut tpitch, tt0,tmom,nmch,mcsel;

bool donecuts(false);
void KalCuts() {
  for(size_t ic=0;ic<4;++ic){
    char cutstring[100];
    snprintf(cutstring,100,"nactive>=%i",minnactive[ic]);
    ncuts[ic] = TCut(cutstring);
    snprintf(cutstring,100,"t0err<%3.2f",maxt0err[ic]);
    t0cuts[ic] = TCut(cutstring);
    snprintf(cutstring,100,"fitmomerr<%4.3f",maxmomerr[ic]);
    momcuts[ic] = TCut(cutstring);
    snprintf(cutstring,100,"fitcon>%5.4f",minfitcon[ic]);
    fitcuts[ic] = TCut(cutstring);
  }
  char ctext[80];
  snprintf(ctext,80,"td>%4.3f&&td<%4.3f",tdlow,tdhigh);
  rpitch = TCut(ctext);
  snprintf(ctext,80,"t0>%f",t0min);
  livegate = TCut(ctext);
  snprintf(ctext,80,"mcenttd>%4.3f&&mcenttd<%4.3f",tdlow-0.02,tdhigh+0.02);
  tpitch = TCut(ctext);
  snprintf(ctext,80,"mct0>%f",t0min);
  tt0 = TCut(ctext);
  tmom = TCut("mcentmom>100");
  snprintf(ctext,80,"ncgood>=%i",minnhits);
  nmch = TCut(ctext);
  mcsel = nmch+tmom+tpitch;
//  mcsel = nmch+tmom;

  reco = TCut("fitstatus>0");
  goodfit = reco+ncuts[icut]+t0cuts[icut]+momcuts[icut]+fitcuts[icut];
  cosmic = TCut("abs(d0)<105 && d0+2/om>450 && d0+2/om<680");
  snprintf(ctext,80,"fitmom>%f&&fitmom<%f",momlow,momhigh);
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

void KalFitHit (TTree* hits ) {
  if(!donecuts)KalCuts();
  TCut goodhit = reco+ncuts[icut]+momcuts[icut]+"_active";

  TCanvas* dcan = new TCanvas("driftcan","driftcan",1200,800);
  TH1F* dres = new TH1F("dres","Drift radius resolution;mm",100,-1,1);
  TH1F* dpull = new TH1F("dpull","Drift radius pull",100,-10,10);
   TH1F* rpull = new TH1F("rpull","residual pull",100,-10,10);
  hits->Project("dres","_rdrift-_mcdist",goodhit);
  hits->Project("dpull","(_rdrift-_mcdist)/_rdrifterr",goodhit);
  hits->Project("rpull","_resid/_residerr",goodhit);
  dcan->Clear();
  dcan->Divide(2,2);
  dcan->cd(1);
  dres->Fit("gaus");
  dcan->cd(2);
  dpull->Fit("gaus");
  dcan->cd(3);
  rpull->Fit("gaus");


  TCanvas* t2dcan = new TCanvas("t2dcan","t2dcan",1200,800);
  TH2F* drad = new TH2F("drad","Drift radius;true drift radius (mm);reco drift radius (mm)",
      55,-0.1,2.6,55,-0.1,2.6);
  TProfile* pdresid = new TProfile("pdresid","Drift Residual vs drift radius;reco drift radius (mm);reco - true radius (mm)",25,0.0,2.5,0,2.5,"s");
  TProfile* pdresidt = new TProfile("pdresidt","Drift Residual vs true drift radius;MC true drift radius (mm);reco - true radius (mm)",25,0.0,2.5,0.,2.5,"s");
  TH2F* dresid = new TH2F("dresid","Drift Residual vs drift radius;reco drift radius (mm);reco - true radius (mm)",25,0.0,2.5,25,-0.5,0.8);
  TH2F* dresidt = new TH2F("dresidt","Drift Residual vs true drift radius;MC true drift radius (mm);reco - true radius (mm)",25,0.0,2.5,25,-0.5,0.8);
  hits->Project("drad","_rdrift:_mcdist",goodhit);
  hits->Project("pdresid","_rdrift-_mcdist:_rdrift",goodhit);
  hits->Project("pdresidt","_rdrift-_mcdist:_mcdist",goodhit);
  hits->Project("dresid","_rdrift-_mcdist:_rdrift",goodhit);
  hits->Project("dresidt","_rdrift-_mcdist:_mcdist",goodhit);
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
//  hits->Project("t0res","_ht-_mcht","_active");
//  hits->Project("t0pull","(_ht-_mcht)/_t0err","_active");
//  tcan->Clear();
//  tcan->Divide(2,2);
//  tcan->cd(1);
//  hits->Draw("_ht:_mcht>>dt0","_active");
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

  hits->Project("tdres","_tddist-_mclen","_active");
  hits->Project("tdpull","(_tddist-_mclen)/_tdderr","_active");
  tdcan->Clear();
  tdcan->Divide(2,2);
  tdcan->cd(1);
  hits->Draw("_tddist:_mclen>>dtd","_active");
  tdcan->cd(2);
  dtd->FitSlicesY(0,20,80);
  TH1D* dtd_1 = (TH1D*)gDirectory->Get("dtd_1");
  dtd_1->SetTitle("Average reco #Delta-t V vs true V;true V (mm);reco V (mm)");
  dtd_1->Fit("pol1");
  tdcan->cd(3);
  tdres->Fit("gaus");
  tdcan->cd(4);
  tdpull->Fit("gaus");

}

void KalFitTrk (TTree* trks ) {
  if(!donecuts)KalCuts();
  TCanvas* tcan = new TCanvas("tt0can","trk_t0can",1200,800);
  TH1F* t00res = new TH1F("t00res","Initial t0 resolution;nsec",100,-20,20);
  TH1F* t0res = new TH1F("t0res","Final t0 resolution;nsec",100,-10,10);
  TH1F* t0pull = new TH1F("t0pull","Track t0 pull",100,-10,10);
//  TH2F* dt0 = new TH2F("dt0","Track t0;true t0 (nsec);Initial t0 (nsec)",
//      100,500,4000,100,500,4000);
  trks->Project("t00res","t00-mct0","fitstatus>0");
  trks->Project("t0res","t0-mct0","fitstatus>0");
  trks->Project("t0pull","(t0-mct0)/t0err","fitstatus>0");
  tcan->Clear();
  tcan->Divide(2,2);
  tcan->cd(1);
  trks->Draw("t00:mct0>>dt0","fitstatus>0");
  tcan->cd(2);
  t00res->Fit("gaus");
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
  trks->Project("d0pull","(d0-mcd0)/d0err","fitstatus>0");
  trks->Project("p0pull","(p0-mcp0)/p0err","fitstatus>0");
  trks->Project("ompull","(om-mcom)/omerr","fitstatus>0");
  trks->Project("z0pull","(z0-mcz0)/z0err","fitstatus>0");
  trks->Project("tdpull","(td-mctd)/tderr","fitstatus>0");
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
//  TH2F* mom = new TH2F("mom","momentum at first hit;true mom (MeV);reco mom (MeV)",
//    100,80,110,100,80,110);
  TH1F* mres = new TH1F("mres","momentum resolution at first hit;MeV",100,-2,2);
  TH1F* mpull = new TH1F("mpull","momentum pull at first hit",100,-10,10);
  TH1F* chisq = new TH1F("chisq","Chisq/NDof",100,0,10);
  trks->Project("mres","fitmom-mcmom","fitstatus>0");
  trks->Project("mpull","(fitmom-mcmom)/fitmomerr","fitstatus>0");
  trks->Project("chisq","chisq/ndof","fitstatus>0");
  fcan->Clear();
  fcan->Divide(2,2);
  fcan->cd(1);
  trks->Draw("fitmom:mcmom>>mom","fitstatus>0");
  fcan->cd(2);
  mres->Draw();
  fcan->cd(3);
  mpull->Fit("gaus");
  fcan->cd(4);
  chisq->Draw();
}

void KalFitAccPlots(TTree* trks) {
  if(!donecuts)KalCuts();

  TH1F* nmc = new TH1F("nmc","N Straw Hits from CE;N straws",81,-0.5,80.5);
  TH1F* mcmom = new TH1F("mcmom","CE true momentum at tracker;CE momentum (MeV)",57,49,106);
  
  TH1F* fitcon = new TH1F("fitcon","log_{10} fit consistency",101,-8,0);
  TH1F* momerr = new TH1F("momerr","Fit momentum error;momentum error (MeV)",100,0,0.5);
  TH1F* t0err = new TH1F("t0err","Fit t_{0} error; t_{0} error (ns)",100,0,2.0);
  TH1F* na = new TH1F("na","Fit N active hits;N hits",71,-0.5,70.5);

  TH1F* t0 = new TH1F("t0","Track Fit t_{0};t_{0} (ns)",100,300,1800);
  TH1F* td = new TH1F("td","Track Fit tan(#lambda);tan(#lambda)",100,0.5,1.5);
  TH1F* d0 = new TH1F("d0","Track fit d_{0};d_{0} (mm)",100,-150,150);
  TH1F* rmax = new TH1F("rmax","Track fit rmax;d_{0}+2/#omega (mm)",100,300,900);

  TH1F* fitmom = new TH1F("fitmom","Track fit momentum;fit momentum (MeV)",100,98,107);

  trks->Project("nmc","ncgood");
  trks->Project("mcmom","mcentmom",nmch);
  
  trks->Project("fitcon","log10(fitcon)",reco+nmch+tmom);
  trks->Project("momerr","fitmomerr",reco+nmch+tmom);
  trks->Project("t0err","t0err",reco+nmch+tmom);
  trks->Project("na","nactive",reco+nmch+tmom);

  trks->Project("t0","t0",reco+nmch+tmom+goodfit);
  trks->Project("td","td",reco+nmch+tmom+goodfit+livegate);

  trks->Project("d0","d0",reco+nmch+tmom+goodfit+livegate);
  trks->Project("rmax","d0+2.0/om",reco+nmch+tmom+goodfit+livegate);

  trks->Project("fitmom","fitmom",reco+nmch+tmom+goodfit+livegate+cosmic);

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

void KalFitAcc(TTree* trks) {
  if(!donecuts)KalCuts();
  unsigned nbins(10);
  double bmax = nbins-0.5;
  TH1F* acc = new TH1F("acc","CE Acceptance;;cummulative acceptance",nbins,-0.5,bmax);
  TH1F* racc = new TH1F("racc","CE Acceptance;;relative acceptance",nbins,-0.5,bmax);
//  acc->Sumw2();
//  racc->Sumw2();
  unsigned ibin(1);
  acc->GetXaxis()->SetBinLabel(ibin++,"All CE");
  acc->GetXaxis()->SetBinLabel(ibin++,">=20 CE SH");
  acc->GetXaxis()->SetBinLabel(ibin++,"CE p>100 MeV/c");
  acc->GetXaxis()->SetBinLabel(ibin++,"CE pitch");
  acc->GetXaxis()->SetBinLabel(ibin++,"KF Track fit");
  acc->GetXaxis()->SetBinLabel(ibin++,"Fit Quality");
  acc->GetXaxis()->SetBinLabel(ibin++,"Livegate");
  acc->GetXaxis()->SetBinLabel(ibin++,"Reco pitch");
  acc->GetXaxis()->SetBinLabel(ibin++,"Cosmic Rejection");
  acc->GetXaxis()->SetBinLabel(ibin++,"Momentum window");


  ibin = 1;
  racc->GetXaxis()->SetBinLabel(ibin++,"All CE");
  racc->GetXaxis()->SetBinLabel(ibin++,">=20 CE SH");
  racc->GetXaxis()->SetBinLabel(ibin++,"CE p>100 MeV/c");
  racc->GetXaxis()->SetBinLabel(ibin++,"CE pitch");
  racc->GetXaxis()->SetBinLabel(ibin++,"KF Track fit");
  racc->GetXaxis()->SetBinLabel(ibin++,"Fit Quality");
  racc->GetXaxis()->SetBinLabel(ibin++,"Livegate");
  racc->GetXaxis()->SetBinLabel(ibin++,"Reco pitch");
  racc->GetXaxis()->SetBinLabel(ibin++,"Cosmic Rejection");
  racc->GetXaxis()->SetBinLabel(ibin++,"Momentum window");
  
  ibin = 0;
  const char* binnames[11] ={"0.0","1.0","2.0","3.0","4.0","5.0","6.0","7.0","8.0","9.0","10.0"};
  trks->Project("acc",binnames[ibin++]);
  trks->Project("+acc",binnames[ibin++],nmch);
  trks->Project("+acc",binnames[ibin++],nmch+tmom);
  trks->Project("+acc",binnames[ibin++],nmch+tmom+tpitch);
  trks->Project("+acc",binnames[ibin++],nmch+tmom+tpitch+reco);
  trks->Project("+acc",binnames[ibin++],nmch+tmom+tpitch+reco+goodfit);
  trks->Project("+acc",binnames[ibin++],nmch+tmom+tpitch+reco+goodfit+livegate);
  trks->Project("+acc",binnames[ibin++],nmch+tmom+tpitch+reco+goodfit+livegate+rpitch);
  trks->Project("+acc",binnames[ibin++],nmch+tmom+tpitch+reco+rpitch+livegate+cosmic+goodfit);
  trks->Project("+acc",binnames[ibin++],nmch+tmom+tpitch+reco+rpitch+livegate+cosmic+goodfit+rmom);

  double all = acc->GetBinContent(1);
  double prev = all;
  for(ibin=1;ibin<=nbins;ibin++){
    racc->SetBinContent(ibin,acc->GetBinContent(ibin)/prev);
    prev = acc->GetBinContent(ibin);
  }
  cout << "Found " << all << "Entries." << endl;
  racc->SetMaximum(1.1);
  acc->Scale(1.0/all);
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

void KalFitRes(TTree* trks,int mincut=0,int maxcut=3) {
  if(!donecuts)KalCuts();
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

  TF1* cball = new TF1("cball",crystalball,-2.0,1.5,7);
  cball->SetParName(0,"Norm");
  cball->SetParName(1,"x0");
  cball->SetParName(2,"sigma");
  cball->SetParName(3,"n");
  cball->SetParName(4,"alpha");
  cball->SetParName(5,"tailfrac");
  cball->SetParName(6,"taillambda");

  TH1F* momres[4];
  TF1*  fitmomres[4];
  TH1F* effnorm = new TH1F("effnorm","effnorm",100,0,150);
  trks->Project("effnorm","mcentmom",mcsel+livegate);
 
  TCanvas* rcan = new TCanvas("rcan","Momentum Resolution",1200,800);
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
    momres[ires] = new TH1F(mname,"momentum resolution at start of tracker;MeV",251,-4,4);
//  momres[ires]->SetStats(0);
    TCut quality = ncuts[ires]+t0cuts[ires]+momcuts[ires]+fitcuts[ires];
    TCut final = reco+quality+rpitch+cosmic+livegate;
    trks->Project(mname,"fitmom-mcentmom",final);
    double integral = momres[ires]->GetEntries()*momres[ires]->GetBinWidth(1);
    cout << "Integral = " << integral << " mean = " << momres[ires]->GetMean() << " rms = " << momres[ires]->GetRMS() << endl;
    cball->SetParameters(3*integral,momres[ires]->GetMean()+0.07,0.3*momres[ires]->GetRMS(),3.0,1.0,0.02,0.2);
    cball->SetParLimits(5,0.001,0.4);
    cball->SetParLimits(6,0.1,momres[ires]->GetRMS());

    momres[ires]->SetMinimum(0.5);
    momres[ires]->Fit("cball","LRQ");
    momres[ires]->Fit("cball","LRM");
    fitmomres[ires] = new TF1(*cball);
    fitmomres[ires]->SetName(fitname);
    gDirectory->Append(fitmomres[ires]);

    TLine* zero = new TLine(0.0,0.0,0.0,momres[ires]->GetBinContent(momres[ires]->GetMaximumBin()));
    zero->SetLineStyle(2);
    zero->Draw();
  
    double keff = momres[ires]->GetEntries()/effnorm->GetEntries();

//    TPaveText* ttext = new TPaveText(0.1,0.75,0.4,0.9,"NDC");  
//    ttext->AddText("Truth Cuts");
//    ttext->AddText(nmch.GetTitle());
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
    sprintf(line,"%s",ncuts[ires].GetTitle());
    rtext->AddText(line);
    sprintf(line,"%s",t0cuts[ires].GetTitle());
    rtext->AddText(line);
    sprintf(line,"%s",momcuts[ires].GetTitle());
    rtext->AddText(line);
    sprintf(line,"%s",fitcuts[ires].GetTitle());
    rtext->AddText(line);
    sprintf(line,"Eff=%4.4f",keff);
    rtext->AddText(line);
    rtext->Draw();
 
  }
  rcan->cd(0);
}

void KalFitRes2(TTree* trks,int ires) {
  if(!donecuts)KalCuts();
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
  TCut quality = ncuts[ires] && t0cuts[ires] && momcuts[ires] && fitcuts[ires];
  TCut final = reco+quality+mcsel;
  trks->Project("momres","fitmom-mcentmom",final);
  double integral = momres->GetEntries()*momres->GetBinWidth(1);
  cout << "Integral = " << integral << " mean = " << momres->GetMean() << " rms = " << momres->GetRMS() << endl;
  cball2->SetParameters(3*integral,momres->GetMean()+0.07,0.3*momres->GetRMS(),3.0,1.0,0.02,0.2);
  cball2->SetParLimits(5,0.001,0.4);
  cball2->SetParLimits(6,0.1,momres->GetRMS());

  momres->Fit("cball2","RNQ");
  momres->Fit("cball2","LRNQ");

  TH1F* effnorm = new TH1F("effnorm","effnorm",100,0,150);
  trks->Project("effnorm","mcentmom",mcsel);

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

void KalFitAmbig(TTree* t, int acut=0) {
  if(!donecuts)KalCuts();
  gStyle->SetOptStat(1111);

  TCut gambig("_mcambig==_ambig");
  TCut bambig("_mcambig!=_ambig&&_ambig!=0");
  TCut nambig("_ambig==0");
  TCut active("_active>0");
// apply requested cuts
  TCut quality;
  if(acut>=0)
    quality = ncuts[acut] && t0cuts[acut] && momcuts[acut] && fitcuts[acut];

  TCut goodtrk = (reco+quality+mcsel);

//  TCut goodtrk ="fitstatus>0";

  TH1F* rdg = new TH1F("rdg","Hit fraction vs drift radius;true radius (mm);hit fraction",100,0.0,2.7);
  TH1F* rdn = new TH1F("rdn","Hit fraction vs drift radius;true radius (mm);hit fraction",100,0.0,2.7);
  TH1F* rdb = new TH1F("rdb","Hit fraction vs drift radius;true radius (mm);hit fraction",100,0.0,2.7);
  TH1F* rda = new TH1F("rda","Drift radius;true radius (mm)",100,0.0,2.7);
  TH1F* rdi = new TH1F("rdi","Drift radius;true radius (mm)",100,0.0,2.7);
  rdg->SetLineColor(kGreen);
  rdn->SetLineColor(kBlue);
  rdb->SetLineColor(kRed);
  rda->SetLineColor(kBlack);
  rdi->SetLineColor(kCyan);
  rdg->SetStats(0);
  rdn->SetStats(0);
//  rdb->SetStats(0);
  rdi->SetStats(0);
  rda->SetStats(0);
  rdg->Sumw2();
  rdn->Sumw2();
  rdb->Sumw2();
  rda->Sumw2();

  t->Project("rdg","_mcdist",goodtrk+active+gambig);
  t->Project("rdn","_mcdist",goodtrk+active+nambig);
  t->Project("rdb","_mcdist",goodtrk+active+bambig);
  t->Project("rda","_mcdist",goodtrk+active);
  t->Project("rdi","_mcdist",goodtrk+!active);
  Double_t ntotal = rda->GetEntries();
  Double_t nright = rdg->GetEntries();
  Double_t nneutral = rdn->GetEntries();
  Double_t nwrong = rdb->GetEntries();
  rdg->Divide(rda);  
  rdn->Divide(rda);
  rdb->Divide(rda);

//  TH1F* frdg = new TH1F("frdg","True Drift radius, failed fits;radius (mm);N hits",100,-0.05,2.55);
//  TH1F* frdb = new TH1F("frdb","True Drift radius, failed fits;radius (mm);N hits",100,-0.05,2.55);
//  frdg->SetLineColor(kBlue);
//  frdb->SetLineColor(kRed);
//  frdg->SetStats(0);
//  frdb->SetStats(0);

//  t->Project("frdg","_mcdist",mcsel+active+ambig+(!goodfit));
//  t->Project("frdb","_mcdist",mcsel+active+(!ambig)+(!goodfit));

  TH1F* momres0 = new TH1F("momres0","Momentum resolution at start of tracker;p_{reco}-p_{true}(MeV/c)",151,-4,4);
  TH1F* momres1 = new TH1F("momres1","Momentum resolution at start of tracker;p_{reco}-p_{true}(MeV/c)",151,-4,4);
  TH1F* momres2 = new TH1F("momres2","Momentum resolution at start of tracker;p_{reco}-p_{true}(MeV/c)",151,-4,4);
  momres0->SetLineColor(kBlack);
  momres1->SetMarkerColor(kCyan);
  momres1->SetMarkerStyle(4);
  momres2->SetMarkerColor(kOrange);
  momres2->SetMarkerStyle(5);
 //  momres->SetStats(0);
  t->Project("momres0","fitmom-mcentmom",goodtrk);
  t->Project("momres1","fitmom-mcentmom",goodtrk&&"fitstatus==1");
  t->Project("momres2","fitmom-mcentmom",goodtrk&&"fitstatus==2");

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
  t->Project("afg","fitmom-mcentmom",goodtrk+active+gambig);
  t->Project("afn","fitmom-mcentmom",goodtrk+active+nambig);
  t->Project("afb","fitmom-mcentmom",goodtrk+active+bambig);
  t->Project("afa","fitmom-mcentmom",goodtrk+active);
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
  TLegend* drleg = new TLegend(0.15,0.35,0.55,0.5);
  char dtitle[100];
  snprintf(dtitle,100,"%4.0f Active hits",rda->GetEntries());
  drleg->AddEntry(rda,dtitle,"l");
  snprintf(dtitle,100,"%4.0f Inactive hits",rdi->GetEntries());
  drleg->AddEntry(rdi,dtitle,"l");
  drleg->Draw();

  ambigcan->cd(2);

  TF1* ex = new TF1("ex","[0]*exp(-x/[1])/[1]+[2]*x+[3]",4);
  ex->SetParameters(0.0,0.1,-0.01,0.04);
  ex->SetParName(0,"ExpNorm");
  ex->SetParName(1,"Lambda");
  ex->SetParName(2,"Slope");
  ex->SetParName(3,"Intercept");
  rdb->SetMaximum(1.1);
  rdb->SetMinimum(0.0);
  rdb->Fit(ex,"L");
  rdn->Draw("same");
  rdg->Draw("same");

  TLegend* leg = new TLegend(0.16,0.35,0.625,0.6);
  char title[80];
  snprintf(title,80,"Correct ambiguity %4.3f",nright/ntotal);
  leg->AddEntry(rdg,title,"l");
  snprintf(title,80,"0 ambiguity %4.3f",nneutral/ntotal);
  leg->AddEntry(rdn,title,"l");
  snprintf(title,80,"Incorrect ambiguity %4.3f",nwrong/ntotal);
  leg->AddEntry(rdb,title,"l");
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

void KalFitResid(TTree* t) {
  if(!donecuts)KalCuts();


  TCut delta("_mcproc==17");
  TCut primary("_mcproc==56");
  TCut gambig("_mcambig==_ambig&&_ambig!=0");
  TCut bambig("_mcambig!=_ambig&&_ambig!=0");
  TCut nambig("_ambig==0");
  TCut active("fitstatus==1 && _active>0");

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

  t->Project("rdg","_mcdist",reco+mcsel+active+gambig+primary);
  t->Project("rdb","_mcdist",reco+mcsel+active+bambig+primary);
  t->Project("rdn","_mcdist",reco+mcsel+active+nambig+primary);

  t->Project("rpullg","_resid/_residerr",reco+mcsel+active+gambig+primary);
  t->Project("rpullb","_resid/_residerr",reco+mcsel+active+bambig+primary);
  t->Project("rpulln","_resid/_residerr",reco+mcsel+active+nambig+primary);
  t->Project("rpulld","_resid/_residerr",reco+mcsel+active+delta);

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

void KalFitCon(TTree* t) {
  if(!donecuts)KalCuts();
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

  t->Project("con1","fitcon",mcsel+"fitstatus==1");
  t->Project("con2","fitcon",mcsel+"fitstatus==2");
  t->Project("lcon1","log10(fitcon)",mcsel+"fitstatus==1");
  t->Project("lcon2","log10(fitcon)",mcsel+"fitstatus==2");

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

void KalFitError(TTree* t){
  if(!donecuts)KalCuts();
    TF1* sgau = new TF1("sgau",splitgaus,-1.,1.,7);
    sgau->SetParName(0,"Norm");
    sgau->SetParName(1,"Mean");
    sgau->SetParName(2,"SigH");
    sgau->SetParName(3,"SigL");
    sgau->SetParName(4,"TFH");
    sgau->SetParName(5,"TSigH");
    sgau->SetParName(6,"TSigL");

//    TCut quality = ncuts[ires] && fitcuts[ires];
    TCut final = (reco+mcsel);
    TH1F* momres1 = new TH1F("momres1","momentum resolution at start of tracker;MeV",151,-2.5,2.5);
    TH1F* momres2 = new TH1F("momres2","momentum resolution at start of tracker;MeV",151,-2.5,2.5);
   
    t->Project("momres1","fitmom-mcentmom",final+"fitmomerr<0.15");
    t->Project("momres2","fitmom-mcentmom",final+"fitmomerr>0.2");
    
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

void KalFitDrift(TTree* t){
  if(!donecuts)KalCuts();

  TH1F* rdg = new TH1F("rdg","Reco Hit Drift Radius;radius (mm);N hits",100,-0.05,2.55);
  TH1F* rdb = new TH1F("rdb","Reco Hit Drift Radius;radius (mm);N hits",100,-0.05,2.55);
  rdg->SetStats(0);
  rdb->SetStats(0);

  t->Project("rdg","_rdrift","fitstatus==1&&fitmomerr<0.15&&_active");
  t->Project("rdb","_rdrift","fitstatus==1&&fitmomerr>0.2&&_active");

  TCanvas* dcan = new TCanvas("dcan","Drift Radius",800,500);
  dcan->Divide(2,1);
  dcan->cd(1);
  rdg->Draw();
  dcan->cd(2);
  rdb->Draw();

}

void KalFitNHits(TTree* t){
  if(!donecuts)KalCuts();
  TH1F* nch = new TH1F("nch","Number of Conversion Electron Tracker Hits;N hits;N Conversion Tracks",100,-0.5,99.5);
  TH1F* ncha = new TH1F("ncha","Number of Conversion Electron Tracker Hits;N hits;N Conversion Tracks",100,-0.5,99.5);
  nch->SetStats(0);
  nch->SetStats(0);
  t->Project("nch","ncgood",tmom+tpitch);
  t->Project("ncha","ncgood",tmom+tpitch+"fitstatus>0");
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

void KalFitNactive(TTree* t) {
  if(!donecuts)KalCuts();
  gStyle->SetOptStat(111111);
  TCut pure("nactive-ncactive==0");
  TCut goodtrk =mcsel+reco+rpitch;

  TH1F* dna = new TH1F("dna","N non-conversion hits active in fit;N_{active}-N_{active,conversion}",10,-0.5,9.5);
  t->Project("dna","nactive-ncactive",goodtrk);

  TH1F* momres0 = new TH1F("momres0","Momentum resolution at start of tracker;p_{reco}-p_{true}(MeV/c)",101,-4,4);
  TH1F* momres1 = new TH1F("momres1","Momentum resolution at start of tracker;p_{reco}-p_{true}(MeV/c)",101,-4,4);
  TH1F* momres2 = new TH1F("momres2","Momentum resolution at start of tracker;p_{reco}-p_{true}(MeV/c)",101,-4,4);
  momres0->SetLineColor(kBlack);
  momres1->SetMarkerColor(kCyan);
  momres1->SetMarkerStyle(4);
  momres2->SetMarkerColor(kOrange);
  momres2->SetMarkerStyle(5);
 //  momres->SetStats(0);
  t->Project("momres0","fitmom-mcentmom",goodtrk);
  t->Project("momres1","fitmom-mcentmom",goodtrk&&pure);
  t->Project("momres2","fitmom-mcentmom",goodtrk&&!pure);

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
