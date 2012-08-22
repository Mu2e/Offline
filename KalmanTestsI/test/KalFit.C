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

// basic parameters

double tdlow(0.57735027);
double tdhigh(1.0);
double t0min(720);
double momlow(103.3);
double momhigh(104.7);
int minnhits(25);
size_t icut=2;

TCut ncuts[4], t0cuts[4], momcuts[4], fitcuts[4], fitcutsText[4];
TCut reco,goodfit,cosmic,rmom,rpitch,livegate;
TCut tpitch, tt0,tmom,nmch,mcsel;

bool donecuts(false);
void KalCuts() {
  ncuts[0] = "fitinfo.nhits>=30";
  ncuts[1] = "fitinfo.nhits>=30";
  ncuts[2] = "fitinfo.nhits>=40";
  ncuts[3] = "fitinfo.nhits>=50";
  t0cuts[0] = "fitinfo.errt0<2.0";
  t0cuts[1] = "fitinfo.errt0<1.5";
  t0cuts[2] = "fitinfo.errt0<1.0";
  t0cuts[3] = "fitinfo.errt0<0.9";
  momcuts[0] = "fitinfo.fitmomerr<0.3";
  momcuts[1] = "fitinfo.fitmomerr<0.2";
  momcuts[2] = "fitinfo.fitmomerr<0.18";
  momcuts[3] = "fitinfo.fitmomerr<0.15";
  fitcutsText[0] = "fitcon>1e-6";
  fitcutsText[1] = "fitcon>1e-4";
  fitcutsText[2] = "fitcon>1e-3";
  fitcutsText[3] = "fitcon>1e-2";
  fitcuts[0] = "TMath::Prob(fitinfo.chi2,(fitinfo.nhits-5))>1e-6";
  fitcuts[1] = "TMath::Prob(fitinfo.chi2,(fitinfo.nhits-5))>1e-4";
  fitcuts[2] = "TMath::Prob(fitinfo.chi2,(fitinfo.nhits-5))>1e-3";
  fitcuts[3] = "TMath::Prob(fitinfo.chi2,(fitinfo.nhits-5))>1e-2";

  char ctext[80];
  snprintf(ctext,80,"recopar.paramVec.m[4]>%4.3f&&recopar.paramVec.m[4]<%4.3f",tdlow,tdhigh);
  rpitch = TCut(ctext);
  snprintf(ctext,80,"fitinfo.t0fit>%f",t0min);
  livegate = TCut(ctext);
  snprintf(ctext,80,"startpar.paramVec.m[4]>%4.3f&&startpar.paramVec.m[4]<%4.3f",tdlow-0.02,tdhigh+0.02);
  tpitch = TCut(ctext);
  snprintf(ctext,80,"fitinfo.t0>%f",t0min);
  tt0 = TCut(ctext);
  tmom = TCut("fitinfo.momin>100");
  snprintf(ctext,80,"fitinfo.nhitstot>=%i",minnhits);
  nmch = TCut(ctext);
  mcsel = nmch+tmom+tpitch;
//  mcsel = nmch+tmom;

  reco = TCut("fitinfo.fit>0");
  goodfit = reco+ncuts[icut]+t0cuts[icut]+momcuts[icut]+fitcuts[icut];
  cosmic = TCut("abs(recopar.paramVec.m[0])<105 && recopar.paramVec.m[0]+2.0/recopar.paramVec.m[2]>460 && recopar.paramVec.m[0]+2.0/recopar.paramVec.m[2]<660");
  snprintf(ctext,80,"fitinfo.fitmom>%f&&fitinfo.fitmom<%f",momlow,momhigh);
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
  // par[5] : fraction of 2nd Gaussian
  // par[6] : tail gaussian sigma

  if ( (x[0]- par[1])/fabs(par[2]) > -1.*par[4]) {
    double g = par[0]*TMath::Gaus(x[0], par[1], par[2]);
//    double g2 = par[5]*par[0]*TMath::Gaus(x[0], par[1], par[6]);
//    return g1+g2;
    double e = par[0]*par[5]*exp(-(x[0]-par[1])/par[6])/par[6];
    return g+e;
  }
  else {
    double A = pow(par[3]/fabs(par[4]), par[3])*exp(-0.5*par[4]*par[4]);
    double B = par[3]/fabs(par[4]) - fabs(par[4]);
    return par[0]*A*pow(B-(x[0]-par[1])/fabs(par[2]), -1.*par[3]);
  }
}

Double_t doubleCrystalball (Double_t *x, Double_t *par) {
	  // par[0] : norm
	  // par[1] : x0
	  // par[2] : sigma
	  // par[3] : n
	  // par[4] : alpha
	  // par[5] : fraction core
	  // par[6] : tail gaussian sigma
	  // par[7] : tail n
	  // par[8] : tail alpha

	  const double invSqrt2pi=0.398942280401432702863;

	  double DeltaX = (x[0]- par[1]);
	  double absSigma = fabs(par[2]);
	  double fval = par[0]*invSqrt2pi/par[2];
	  if ( DeltaX/absSigma > -1.*par[4]) {
	    fval *= TMath::Gaus(x[0], par[1], par[2]);
	  }
	  else {
	    double absAlpha = fabs(par[4]);
	    double A = pow(par[3]/absAlpha, par[3])*exp(-0.5*par[4]*par[4]);
	    double B = par[3]/absAlpha - absAlpha;
	    fval *= A*pow(B-DeltaX/absSigma, -1.*par[3]);
	  }

	  double tailFval=par[0]*invSqrt2pi/par[6];
	  double tailAbsSigma = fabs(par[6]);
	  DeltaX*=-1.0;
	  if ( DeltaX/tailAbsSigma > -1.*par[8]) {
	    tailFval *= TMath::Gaus(x[0], par[1], par[6]);
	  }
	  else {
	    double tailAbsAlpha = fabs(par[8]);
	    double tailA = pow(par[7]/tailAbsAlpha, par[7])*exp(-0.5*par[8]*par[8]);
	    double tailB = par[7]/tailAbsAlpha - tailAbsAlpha;
	    tailFval *= tailA*pow(tailB-DeltaX/tailAbsSigma, -1.*par[7]);
	  }

	  return fval*par[5] + (1.0-par[5])*tailFval;
}

/*
void KalFitHit (TTree* hits ) {
  TCanvas* dcan = new TCanvas("driftcan","driftcan",1200,800);
  TH1F* dres = new TH1F("dres","Drift radius resolution;mm",100,-1,1);
  TH1F* dpull = new TH1F("dpull","Drift radius pull",100,-10,10);
//  TH2F* drad = new TH2F("drad","Drift radius;true drift radius (mm);reco drift radius (mm)",
//      100,-0.3,2.8,100,-0.3,2.8);
  TH1F* rpull = new TH1F("rpull","residual pull",100,-10,10);
  hits->Project("dres","rdrift-mcrdrift","active");
  hits->Project("dpull","(rdrift-mcrdrift)/rdrifterr","active");
  hits->Project("rpull","resid/residerr","active");
  dcan->Clear();
  dcan->Divide(2,2);
  dcan->cd(1);
  hits->Draw("rdrift:mcrdrift>>drad","active");
  dcan->cd(2);
  dres->Fit("gaus");
  dcan->cd(3);
  dpull->Fit("gaus");
  dcan->cd(4);
  rpull->Fit("gaus");

  TCanvas* tcan = new TCanvas("ht0can","hit_t0can",1200,800);
  TH1F* t0res = new TH1F("t0res","hit t0 resolution;nsec",100,-10,10);
  TH1F* t0pull = new TH1F("t0pull","hit t0 pull",100,-10,10);
//  TH2F* dt0 = new TH2F("dt0","Hit t0;true t0 (nsec);reco t0 (nsec)",
//      100,500,4000,100,500,4000);
  hits->Project("t0res","hitt0-mchitt0","active");
  hits->Project("t0pull","(hitt0-mchitt0)/hitt0err","active");
  tcan->Clear();
  tcan->Divide(2,2);
  tcan->cd(1);
  hits->Draw("hitt0:mchitt0>>dt0","active");
  tcan->cd(2);
  t0res->Fit("gaus");
  tcan->cd(3);
  t0pull->Fit("gaus");

  TCanvas* tdcan = new TCanvas("tdcan","tdcan",1200,800);
  TH1F* tdres = new TH1F("tdres","#Deltat V resolution;mm",100,-200,200);
  TH1F* tdpull = new TH1F("tdpull","#Deltat V pull",100,-10,10);
//  TH2F* dtd = new TH2F("dtd","Hit V position;true V (mm);#Deltat V (mm)",
//      100,-600,600,100,-600,600);
//  TH2F* pocatd = new TH2F("pocatd","Hit POCA V;true V (mm);POCA V (mm)",
//      100,-600,600,100,-600,600);

  hits->Project("tdres","dmid-mcdmid","active");
  hits->Project("tdpull","(dmid-mcdmid)/dmiderr","active");
  tdcan->Clear();
  tdcan->Divide(2,2);
  tdcan->cd(1);
  hits->Draw("dmid:mcdmid>>dtd","active");
  tdcan->cd(2);
  hits->Draw("hflt:mcdmid>>pocatd","active");
  tdcan->cd(3);
  tdres->Fit("gaus");
  tdcan->cd(4);
  tdpull->Fit("gaus");

}

void KalFitTrk (TTree* trks ) {

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

  trks->Project("nmc","nchits");
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
  TLine* fitconcut = new TLine(-4,0.0,-4,fitcon->GetMaximum());
  fitconcut->SetLineColor(kBlack);
  fitconcut->SetLineStyle(2);
  fitconcut->SetLineWidth(2);
  fitconcut->Draw();

  fcan->cd(2);
  na->Draw();
  TLine* nacut = new TLine(20,0.0,20,na->GetMaximum());
  nacut->SetLineColor(kBlack);
  nacut->SetLineStyle(2);
  nacut->SetLineWidth(2);
  nacut->Draw();

  fcan->cd(3);
  t0err->Draw();
  TLine* t0errcut = new TLine(1.5,0.0,1.5,t0err->GetMaximum());
  t0errcut->SetLineColor(kBlack);
  t0errcut->SetLineStyle(2);
  t0errcut->SetLineWidth(2);
  t0errcut->Draw();
  
  fcan->cd(4);
  momerr->Draw();
  TLine* momerrcut = new TLine(0.2,0.0,0.2,momerr->GetMaximum());
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
*/

void KalFitRes(TTree* trks, int nGenEv=-1, TString fitOpt="L") {

	if (nGenEv<=0) nGenEv=((TH1F*)d->Get("g4run/totalMultiplicity"))->GetEntries();

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

//  TF1* cball = new TF1("cball",crystalball,-2.0,1.5,7);
//  cball->SetParName(0,"Norm");
//  cball->SetParName(1,"x0");
//  cball->SetParName(2,"sigma");
//  cball->SetParName(3,"n");
//  cball->SetParName(4,"alpha");
//  cball->SetParName(5,"tailfrac");
//  cball->SetParName(6,"taillambda");

  TF1* dcball = new TF1("dcball",doubleCrystalball,-2.5,2.5,9);
  dcball->SetParName(0,"Norm");
  dcball->SetParName(1,"x0");
  dcball->SetParName(2,"sigma");
  dcball->SetParName(3,"n");
  dcball->SetParName(4,"alpha");
  dcball->SetParName(5,"frac_{core}");
  dcball->SetParName(6,"sigma_{tail}");
  dcball->SetParName(7,"n_{tail}");
  dcball->SetParName(8,"alpha_{tail}");
  dcball->SetNpx(1000);

  TH1F* momres[4];
  TH1F* effnorm = new TH1F("effnorm","effnorm",100,0,150);
  trks->Project("effnorm","fitinfo.momin",mcsel);
 
  TCanvas* rcan = new TCanvas("rcan","Momentum Resolution",1200,800);
  rcan->Clear();
  rcan->Divide(2,2);
  gStyle->SetOptFit(111111);
  gStyle->SetOptStat("oumr");
  for(unsigned ires=0;ires<4;ires++){
    rcan->cd(ires+1);
    gPad->SetLogy();
    char mname[50];
    snprintf(mname,50,"momres%i",ires);
    momres[ires] = new TH1F(mname,"momentum resolution at start of tracker;MeV",251,-2.5,2.5);
//  momres[ires]->SetStats(0);
    TCut quality = ncuts[ires] && t0cuts[ires] /*&& momcuts[ires]*/ && fitcuts[ires];
    TCut final = (reco+quality+rpitch+mcsel);
    trks->Project(mname,"fitinfo.fitmom-fitinfo.momin",final);
    double integral = momres[ires]->GetEntries()*momres[ires]->GetBinWidth(1);
//    sgau->SetParameters(integral,0.0,0.5*momres[ires]->GetRMS(),0.5*momres[ires]->GetRMS(),0.01,2*momres[ires]->GetRMS(),2*momres[ires]->GetRMS());
//    sgau->SetParLimits(5,0.1,1.0);
//    sgau->SetParLimits(6,0.1,1.0);
//    sgau->SetParLimits(4,0.0,0.8);
//    momres[ires]->Fit("sgau","LIR");
//    degau->SetParameters(integral,0.1,0.2,0.0,0.3*momres[ires]->GetRMS(),2*momres[ires]->GetRMS(),0.5*momres[ires]->GetRMS(),0.25);
//    degau->SetParLimits(1,0.02,0.3);
//    degau->SetParLimits(2,0.02,0.3);
//    degau->SetParLimits(4,0.05,momres[ires]->GetRMS());
//    degau->SetParLimits(5,0.12,2*momres[ires]->GetRMS());
//    degau->SetParLimits(6,0.1,momres[ires]->GetRMS());
//    degau->SetParLimits(7,1,3.0);
  
//    cball->SetParameters(integral,0.0,0.1,1.0,1.0,0.05,0.5);
//    cball->SetParLimits(5,0.001,0.4);
//    cball->SetParLimits(6,0.1,momres[ires]->GetRMS());
//    momres[ires]->Fit("cball","LIR");

    dcball->SetParameters(integral,momres[ires]->GetMean(1)*0.5,momres[ires]->GetRMS(1)*0.25,1.0,1.0
  		  ,0.8,momres[ires]->GetRMS(1),1.0,1.0);
    dcball->SetParLimits(0,1,integral*2);
    dcball->SetParLimits(1,-0.2,0.4);
    dcball->SetParLimits(2,0.03,momres[ires]->GetRMS(1)*0.7);
    dcball->SetParLimits(3,0,20.);
    dcball->SetParLimits(4,0,5);
    dcball->SetParLimits(5,0,1.0);
    dcball->SetParLimits(6,0.1,momres[ires]->GetRMS(1)*2.0);
    dcball->SetParLimits(7,0,30.);
    dcball->SetParLimits(8,0,10);

    momres[ires]->Fit("dcball");
    momres[ires]->Fit("dcball");
    momres[ires]->Fit("dcball",fitOpt.Data());
    //momres[ires]->Fit("dcball","LIR");

    double covMat[9][9];
    gMinuit->mnemat(&covMat[0][0],9);
    double corrF = covMat[2][6]/(sqrt(covMat[2][2]*covMat[6][6]));
    double wSigmaCore = dcball->GetParameter(2)*dcball->GetParameter(5);
    double wSigmaTail = dcball->GetParameter(6)*(1.0-dcball->GetParameter(5));
    double meanSigma = sqrt(wSigmaCore*wSigmaCore + wSigmaTail*wSigmaTail + 2.0*corrF*wSigmaCore*wSigmaTail);

    int intFirstBin = 1 + (int)(((dcball->GetParameter(1)-3.0*meanSigma) - momres[ires]->GetXaxis()->GetXmin() )/momres[ires]->GetBinWidth(1));
    int intLastBin  = 1 + (int)(((dcball->GetParameter(1)+3.0*meanSigma) - momres[ires]->GetXaxis()->GetXmin() )/momres[ires]->GetBinWidth(1));

    float nEvtIn3sigma = momres[ires]->Integral(intFirstBin,intLastBin);


    TLine* zero = new TLine(0.0,0.0,0.0,momres[ires]->GetBinContent(momres[ires]->GetMaximumBin()));
    zero->SetLineStyle(2);
    zero->Draw();
  
    double keff = momres[ires]->GetEntries()/effnorm->GetEntries();
    //cout<<"norm "<<effnorm->GetEntries()<<" sel entries "<<momres[ires]->GetEntries()<<" eff "<<keff<<endl;

    TPaveText* ttext = new TPaveText(0.1,0.65,0.4,0.9,"NDC");
    ttext->SetTextSize(0.04);
    ttext->AddText("Truth Cuts");
    ttext->AddText(nmch.GetTitle());
    ttext->AddText(tmom.GetTitle());
    ttext->AddText(Form("%4.3f<tan(#lambda_{in})<%4.3f",tdlow-0.02,tdhigh+0.02));
    //ttext->AddText(tpitch.GetTitle());
    ttext->AddText(Form("accept = %3.2f %%",(effnorm->GetEntries()/((float)nGenEv))*100.0));
    ttext->Draw();
 
    TPaveText* rtext = new TPaveText(0.1,0.4,0.4,0.65,"NDC");
    rtext->SetTextSize(0.04);
    rtext->AddText("Reco Cuts");
    char line[40];
    snprintf(line,80,"%4.3f<tan(#lambda)<%4.3f",tdlow,tdhigh);
    rtext->AddText(line);
    //snprintf(line,80,"t0>%5.1f nsec",t0min);
    //rtext->AddText(line);
    sprintf(line,"%s",ncuts[ires].GetTitle());
    rtext->AddText(line);
    //sprintf(line,"%s",t0cuts[ires].GetTitle());
    //rtext->AddText(line);
    //sprintf(line,"%s",momcuts[ires].GetTitle());
    //rtext->AddText(line);
    //sprintf(line,"%s",fitcuts[ires].GetTitle());
    sprintf(line,"%s",fitcutsText[ires].GetTitle());
    rtext->AddText(line);
    sprintf(line,"Eff = %3.2f %%",keff*100.0);
    rtext->AddText(line);
    rtext->Draw();
    if (!fitOpt.Contains("L")) {
    TPaveText* vtext = new TPaveText(0.1,0.30,0.4,0.4,"NDC");
    vtext->SetTextSize(0.04);
    vtext->AddText(Form("#bar{#sigma} = %.3f",meanSigma));
    vtext->AddText(Form("Eff (in x0#pm 3#bar{#sigma}) = %3.2f %%",(nEvtIn3sigma/effnorm->GetEntries())*100.0));
    vtext->Draw();
    }
 
  }
  rcan->cd(0);
}

/*
void KalFitAmbig(TTree* t, int acut=0) {
  if(!donecuts)KalCuts();
  gStyle->SetOptStat(1111);

  TCut gambig("_mcambig==_ambig");
  TCut bambig("_mcambig!=_ambig&&_ambig!=0");
  TCut nambig("_ambig==0");
  TCut active("_active>0");
// apply requested cuts
  TCut quality;
  if(acut>0)
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
  TCut delta("_mcproc==17");
  TCut ambig("_mcambig==_ambig");
  TCut active("fitstatus==1 && _active>0");

  TH1F* rdg = new TH1F("rdg","True Drift radius;radius (mm);N hits",100,-0.05,2.55);
  TH1F* rdb = new TH1F("rdb","True Drift radius;radius (mm);N hits",100,-0.05,2.55);
  rdg->SetLineColor(kBlue);
  rdb->SetLineColor(kRed);
  rdg->SetStats(0);
  rdb->SetStats(0);

  TH1F* rpullg = new TH1F("rpullg","Residual Pull;Pull;N hits",100,-6,6);
  TH1F* rpullb = new TH1F("rpullb","Residual Pull;Pull;N hits",100,-6,6);
  TH1F* rpulld = new TH1F("rpulld","Residual Pull;Pull;N hits",100,-6,6);
  rpullg->SetLineColor(kBlue);
  rpullb->SetLineColor(kRed);
  rpulld->SetLineColor(kGreen);

  t->Project("rdg","_mcdist",mcsel+active+ambig+(!delta));
  t->Project("rdb","_mcdist",mcsel+active+(!ambig)+(!delta));

  t->Project("rpullg","_resid/_residerr",mcsel+active+ambig+(!delta));
  t->Project("rpullb","_resid/_residerr",mcsel+active+(!ambig)+(!delta));
  t->Project("rpulld","_resid/_residerr",mcsel+active+ambig+delta);

  TCanvas* residcan = new TCanvas("residcan","Residuals",1200,800);
  residcan->Divide(2,1);
  residcan->cd(1);
  rdg->Draw();
  rdb->Draw("same");

  TLegend* leg = new TLegend(0.3,0.3,0.8,0.5);
  leg->AddEntry(rpullg,"Correct ambiguity","l");
  leg->AddEntry(rpullb,"Incorrect ambiguity","l");
  leg->AddEntry(rpulld,"eIoni (#delta ray)","l");
  leg->Draw();

  residcan->cd(2);
  rpullg->Fit("gaus");
  rpullb->Draw("same");
  rpulld->Draw("same");

  residcan->cd(0);

}
*/

void KalFitCon(TTree* t) {
	  if(!donecuts)KalCuts();
  TH1F* fprob1 = new TH1F("fprob1","fit consistency ( Prob(#chi^{2},nDOF) )",250,0,1);
  TH1F* fprob2 = new TH1F("fprob2","fit consistency ( Prob(#chi^{2},nDOF) )",250,0,1);
  fprob1->SetLineColor(kBlue);
  fprob2->SetLineColor(kRed);
  fprob1->SetStats(0);
  fprob2->SetStats(0);

  t->Project("fprob1","TMath::Prob(fitinfo.chi2,(fitinfo.nhits-5))",mcsel+"fitinfo.fit==1");
  t->Project("fprob2","TMath::Prob(fitinfo.chi2,(fitinfo.nhits-5))",mcsel+"fitinfo.fit!=1");

  TH1F* fcon1 = new TH1F("fcon1","log_{10}(fit consistency)",100,-10,0);
  TH1F* fcon2 = new TH1F("fcon2","log_{10}(fit consistency)",100,-10,0);
  fcon1->SetLineColor(kBlue);
  fcon2->SetLineColor(kRed);
  fcon1->SetStats(0);
  fcon2->SetStats(0);

  t->Project("fcon1","log10(TMath::Prob(fitinfo.chi2,(fitinfo.nhits-5)))",mcsel+"fitinfo.fit==1");
  t->Project("fcon2","log10(TMath::Prob(fitinfo.chi2,(fitinfo.nhits-5)))",mcsel+"fitinfo.fit!=1");


  TLegend* leg = new TLegend(0.1,0.6,0.4,0.8);
  leg->AddEntry(fcon1,"Fully Converged Fit","l");
  leg->AddEntry(fcon2,"Unconverged Fit","l");
  TLegend* leg1 = new TLegend(0.1,0.6,0.4,0.8);
  leg1->AddEntry(fprob1,"Fully Converged Fit","l");
  leg1->AddEntry(fprob2,"Unconverged Fit","l");

  TCanvas* fccan = new TCanvas("fccan","fit consistency",1000,500);
  fccan->Clear();
  fccan->Divide(2,1);
  fccan->cd(1)->SetLogy(1);
  fprob1->Draw();
  fprob2->Draw("same");
  leg1->Draw();
  fccan->cd(2)->SetLogy(1);
  fcon1->Draw();
  fcon2->Draw("same");
  leg->Draw();

}

/*
void KalFitError(TTree* t){
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
  TH1F* nch = new TH1F("nch","Number of Conversion Electron Tracker Hits;N hits;N Conversion Tracks",100,-0.5,99.5);
  TH1F* ncha = new TH1F("ncha","Number of Conversion Electron Tracker Hits;N hits;N Conversion Tracks",100,-0.5,99.5);
  nch->SetStats(0);
  nch->SetStats(0);
  t->Project("nch","nchits",tmom+tpitch);
  t->Project("ncha","nchits",tmom+tpitch+"fitstatus>0");
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
*/
