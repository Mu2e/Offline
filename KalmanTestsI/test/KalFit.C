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
#include <fstream>
#include <iomanip>

// basic parameters

double tdlow(/*0.5*/0.57735027);
double tdhigh(1.0/*1.6*/);
double t0min(710);
double momlow(103.3/*103.3*/);
double momhigh(104.7);
int minnhits(50);
//size_t icut=0;
bool useMomErrCut(true);

TCut ncuts[4], t0cuts[4], momcuts[4], fitcuts[4];
TCut ncutsText[4], t0cutsText[4], momcutsText[4], fitcutsText[4];
TCut reco,goodfit,cosmic,rmom,rmomCorr,rpitch,livegate;
TCut tpitch, tt0,tmom,nmch,mcsel;
bool donecuts(false);
TCut onlyFrstSD("fitinfo.iseed==0");

Long64_t nentries = 1000000000;

TCut misrec("recopar.parvec.m[2]>0.0");

void KalCuts(size_t icut=0, int cutOption=0, double liveTime=-1, int useExtapMom=0, double momCorr=0.0, double momErrScale=1.0/*0.61*/) {
  if (liveTime>0.) t0min=liveTime;

  switch (cutOption) {
  case 1:
  {
          minnhits=40;
          ncutsText[0] = Form("nhits>=%i",minnhits+0);
          ncutsText[1] = Form("nhits>=%i",minnhits+0);
          ncutsText[2] = Form("nhits>=%i",minnhits+10);
          ncutsText[3] = Form("nhits>=%i",minnhits+20);//+30);
          ncuts[0] = Form("fitinfo.nhits>=%i",minnhits+0);
          ncuts[1] = Form("fitinfo.nhits>=%i",minnhits+0);
          ncuts[2] = Form("fitinfo.nhits>=%i",minnhits+10);
          ncuts[3] = Form("fitinfo.nhits>=%i",minnhits+20);//+30);
          t0cutsText[0] = "errt0<3";
          t0cutsText[1] = "errt0<2";
          t0cutsText[2] = "errt0<1.2";
          t0cutsText[3] = "errt0<1.2";
          t0cuts[0] = "fitinfo.errt0<3";
          t0cuts[1] = "fitinfo.errt0<2";
          t0cuts[2] = "fitinfo.errt0<1.2";
          t0cuts[3] = "fitinfo.errt0<1.2";
          momcutsText[0] = "fitmomerr<0.4";
          momcutsText[1] = "fitmomerr<0.3";
          momcutsText[2] = "fitmomerr<0.2";
          momcutsText[3] = "fitmomerr<0.2";
          momcuts[0] = Form("%f*fitinfo.fitmomerr<0.4",momErrScale);
          momcuts[1] = Form("%f*fitinfo.fitmomerr<0.3",momErrScale);
          momcuts[2] = Form("%f*fitinfo.fitmomerr<0.2",momErrScale);
          momcuts[3] = Form("%f*fitinfo.fitmomerr<0.2",momErrScale);
          fitcutsText[0] = "fitcon>1e-7";//1e-6;//1e-3
          fitcutsText[1] = "fitcon>1e-6";//1e-5;//1e-2
          fitcutsText[2] = "fitcon>1e-5";//1e-4;//1e-2
          fitcutsText[3] = "fitcon>1e-4";//5e-3;//2e-2
          //  fitcuts[0] = "TMath::Prob(fitinfo.chi2,(fitinfo.nhits-5))>1e-6";//1e-3
          //  fitcuts[1] = "TMath::Prob(fitinfo.chi2,(fitinfo.nhits-5))>1e-4";//1e-2
          //  fitcuts[2] = "TMath::Prob(fitinfo.chi2,(fitinfo.nhits-5))>1e-3";//1e-2
          //  fitcuts[3] = "TMath::Prob(fitinfo.chi2,(fitinfo.nhits-5))>1e-2";//2e-2
          fitcuts[0] = "fitinfo.fitcon>1e-7";//1e-6;//1e-3
          fitcuts[1] = "fitinfo.fitcon>1e-6";//1e-5;//1e-2
          fitcuts[2] = "fitinfo.fitcon>1e-5";//1e-4;//1e-2
          fitcuts[3] = "fitinfo.fitcon>1e-4";//5e-3;//2e-2

  }
          break;
  default:
  {
          ncutsText[0] = Form("nhits>=%i",minnhits+10);
          ncutsText[1] = Form("nhits>=%i",minnhits+10);
          ncutsText[2] = Form("nhits>=%i",minnhits+20);
          ncutsText[3] = Form("nhits>=%i",minnhits+20);//+30);
          ncuts[0] = Form("fitinfo.nhits>=%i",minnhits+10);
          ncuts[1] = Form("fitinfo.nhits>=%i",minnhits+10);
          ncuts[2] = Form("fitinfo.nhits>=%i",minnhits+20);
          ncuts[3] = Form("fitinfo.nhits>=%i",minnhits+20);//+30);
          t0cutsText[0] = "errt0<3";
          t0cutsText[1] = "errt0<2";
          t0cutsText[2] = "errt0<1.2";
          t0cutsText[3] = "errt0<1.2";
          t0cuts[0] = "fitinfo.errt0<3";
          t0cuts[1] = "fitinfo.errt0<2";
          t0cuts[2] = "fitinfo.errt0<1.2";
          t0cuts[3] = "fitinfo.errt0<1.2";
          momcutsText[0] = "fitmomerr<0.4";
          momcutsText[1] = "fitmomerr<0.3";
          momcutsText[2] = "fitmomerr<0.2";
          momcutsText[3] = "fitmomerr<0.2";
          momcuts[0] = Form("%f*fitinfo.fitmomerr<0.4",momErrScale);
          momcuts[1] = Form("%f*fitinfo.fitmomerr<0.3",momErrScale);
          momcuts[2] = Form("%f*fitinfo.fitmomerr<0.2",momErrScale);
          momcuts[3] = Form("%f*fitinfo.fitmomerr<0.2",momErrScale);
          fitcutsText[0] = "fitcon>1e-6";//1e-7;//1e-3
          fitcutsText[1] = "fitcon>1e-5";//1e-6;//1e-2
          fitcutsText[2] = "fitcon>1e-4";//1e-5;//1e-2
          fitcutsText[3] = "fitcon>5e-3";//1e-4;//2e-2
          //  fitcuts[0] = "TMath::Prob(fitinfo.chi2,(fitinfo.nhits-5))>1e-6";//1e-3
          //  fitcuts[1] = "TMath::Prob(fitinfo.chi2,(fitinfo.nhits-5))>1e-4";//1e-2
          //  fitcuts[2] = "TMath::Prob(fitinfo.chi2,(fitinfo.nhits-5))>1e-3";//1e-2
          //  fitcuts[3] = "TMath::Prob(fitinfo.chi2,(fitinfo.nhits-5))>1e-2";//2e-2
          fitcuts[0] = "fitinfo.fitcon>1e-6";//1e-7;//1e-3
          fitcuts[1] = "fitinfo.fitcon>1e-5";//1e-6;//1e-2
          fitcuts[2] = "fitinfo.fitcon>1e-4";//1e-5;//1e-2
          fitcuts[3] = "fitinfo.fitcon>5e-3";//1e-4;//2e-2

  }
  break;
  }

  char ctext[80];
  snprintf(ctext,80,"recopar.parvec.m[4]>%4.3f&&recopar.parvec.m[4]<%4.3f",tdlow,tdhigh);
  rpitch = TCut(ctext);
  snprintf(ctext,80,"fitinfo.t0fit>%f",t0min);
  livegate = TCut(ctext);
  snprintf(ctext,80,"startpar.parvec.m[4]>%4.3f&&startpar.parvec.m[4]<%4.3f",tdlow,tdhigh+0.02);
  tpitch = TCut(ctext);
  snprintf(ctext,80,"fitinfo.t0>%f",t0min);
  tt0 = TCut(ctext);
  tmom = TCut("fitinfo.momin>100");
  snprintf(ctext,80,"fitinfo.nhitstot>=%i",minnhits);
  nmch = TCut(ctext);
  mcsel = nmch+tmom+tpitch;
//  mcsel = nmch+tmom;

  //reco = TCut("fitinfo.fit>0");
  reco = TCut("fitinfo.fit==1");
  reco += misrec;
  goodfit = reco+ncuts[icut]+t0cuts[icut]+fitcuts[icut];
  if (useMomErrCut) goodfit += momcuts[icut];
  cosmic = TCut("abs(recopar.parvec.m[0])<105 && recopar.parvec.m[0]+2.0/recopar.parvec.m[2]>400 && recopar.parvec.m[0]+2.0/recopar.parvec.m[2]<660");

  snprintf(ctext,80,"fitinfo.fitmom>=%f&&fitinfo.fitmom<=%f",momlow,momhigh);
  rmom = TCut(ctext);

  if (useExtapMom>0) {
          snprintf(ctext,80,"fitinfo.fitmombeam>=%f&&fitinfo.fitmombeam<=%f",momlow,momhigh);
  } else if (momCorr>0) {
          snprintf(ctext,80,"(fitinfo.fitmom+%f)>=%f&&(fitinfo.fitmom+%f)<=%f",momCorr,momlow,momCorr,momhigh);
  }
  rmomCorr = TCut(ctext);

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


//Double_t crystalball (Double_t *x, Double_t *par) {
//  // par[0] : norm
//  // par[1] : x0
//  // par[2] : sigma
//  // par[3] : n
//  // par[4] : alpha
//  // par[5] : fraction of 2nd Gaussian
//  // par[6] : tail gaussian sigma
//
//  const double sqrttwPi = 2.50662827463100024;
//  double gausNorm = sqrttwPi*par[2];
//  //par[0]/=gausNorm;
//  //double tailLamb = par[6];
//  double tailLamb = 0.5*gausNorm;
//
//
//  if ( (x[0]- par[1])/fabs(par[2]) > -1.*par[4]) {
//    double g = par[0]*TMath::Gaus(x[0], par[1], par[2])/gausNorm;
////    double g2 = par[5]*par[0]*TMath::Gaus(x[0], par[1], tailLamb);
////    return g1+g2;
//    //double e = 0.0;
//    if ((x[0]- par[1])>0.0) e = par[0]*par[5]*exp(-(x[0]-par[1])/tailLamb)/tailLamb;
//    else return g;
//    //return g+e;
//    return ((1.0-par[5])*g+0.5*e);
//  }
//  else {
//    double A = pow(par[3]/fabs(par[4]), par[3])*exp(-0.5*par[4]*par[4]);
//    double B = par[3]/fabs(par[4]) - fabs(par[4]);
//    return par[0]*A*pow(B-(x[0]-par[1])/fabs(par[2]), -1.*par[3])/gausNorm;
//  }
//}

Double_t crystalball (Double_t *x, Double_t *par) {
  // par[0] : norm
  // par[1] : x0
  // par[2] : sigma
  // par[3] : n
  // par[4] : alpha
  // par[5] : fraction of exponential tail
  // par[6] : tail exponential lambda

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
        double fval = par[0]*invSqrt2pi/absSigma;
        if ( DeltaX/absSigma > -1.*par[4]) {
                fval *= TMath::Gaus(x[0], par[1], par[2]);
        }
        else {
                double absAlpha = fabs(par[4]);
                double A = pow(par[3]/absAlpha, par[3]);
                if (A<1.e+100/*1.8e+278*/) {
                        A*=exp(-0.5*par[4]*par[4]);
                        double B = par[3]/absAlpha - absAlpha;
                        double tmpVal = B-DeltaX/absSigma;
                        if (tmpVal>0) {
                                fval *= A*pow(B-DeltaX/absSigma, -1.*par[3]);
                        } else {
                                fval = 0.0;
                        }
                } else {
                        fval=0.0;
                }
        }

        double tailAbsSigma = fabs(par[6]);
        double tailFval = par[0]*invSqrt2pi/tailAbsSigma;
        DeltaX*=-1.0;
        if ( DeltaX/tailAbsSigma > -1.*par[8]) {
                tailFval *= TMath::Gaus(x[0], par[1], par[6]);
        }
        else {
                double tailAbsAlpha = fabs(par[8]);
                double tailA = pow(par[7]/tailAbsAlpha, par[7]);
               if (tailA<1.0e+100/*1.8e+278*/) {
                        tailA*=exp(-0.5*par[8]*par[8]);
                        double tailB = par[7]/tailAbsAlpha - tailAbsAlpha;
                        double tmpVal = tailB-DeltaX/tailAbsSigma;
                        if (tmpVal>0) {
                                tailFval *= tailA*pow(tailB-DeltaX/tailAbsSigma, -1.*par[7]);
                        } else {
                                tailFval = 0.0;
                        }
                } else {
                        tailFval = 0.0;
                }
        }

        fval*=par[5];
        fval+= (1.0-par[5])*tailFval;
        return fval;
}

Double_t twoNgaus(Double_t *x, Double_t *par) {
  Double_t retval;
  retval = par[0]*(par[3]*TMath::Gaus(x[0],par[1],par[2],true) + (1.0-par[3])*TMath::Gaus(x[0],par[4],par[5],true));
  return retval;
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
  trks->Project("t00res","t00-mct0","fitinfo.fit>0");
  trks->Project("t0res","t0-mct0","fitinfo.fit>0");
  trks->Project("t0pull","(t0-mct0)/t0err","fitinfo.fit>0");
  tcan->Clear();
  tcan->Divide(2,2);
  tcan->cd(1);
  trks->Draw("t00:mct0>>dt0","fitinfo.fit>0");
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
  trks->Project("d0pull","(d0-mcd0)/d0err","fitinfo.fit>0");
  trks->Project("p0pull","(p0-mcp0)/p0err","fitinfo.fit>0");
  trks->Project("ompull","(om-mcom)/omerr","fitinfo.fit>0");
  trks->Project("z0pull","(z0-mcz0)/z0err","fitinfo.fit>0");
  trks->Project("tdpull","(td-mctd)/tderr","fitinfo.fit>0");
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
  trks->Project("mres","fitmom-mcmom","fitinfo.fit>0");
  trks->Project("mpull","(fitmom-mcmom)/fitmomerr","fitinfo.fit>0");
  trks->Project("chisq","chisq/ndof","fitinfo.fit>0");
  fcan->Clear();
  fcan->Divide(2,2);
  fcan->cd(1);
  trks->Draw("fitmom:mcmom>>mom","fitinfo.fit>0");
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
*/
void KalFitAcc(TTree* trks, int nGenEv=-1, TString addCut="", size_t cutType=0, int cutOption=0, double liveTime=-1, int useExtapMom=0, double pCorr=0.0, Long64_t NEntries=-1, Long64_t skipentries=0) {

  if(!donecuts)KalCuts(cutType,cutOption,liveTime,useExtapMom,pCorr);
  TCut addcut(addCut.Data());
  if (NEntries>0 && NEntries<1000000000) {nentries = NEntries;}

  unsigned nbins(10);
  double bmax = nbins-0.5;
  TH1F* acc = new TH1F("acc","CE Acceptance;;cummulative acceptance",nbins,-0.5,bmax);
  TH1F* racc = new TH1F("racc","CE Acceptance;;relative acceptance",nbins,-0.5,bmax);
//  acc->Sumw2();
//  racc->Sumw2();
  unsigned ibin(1);
  acc->GetXaxis()->SetBinLabel(ibin++,"All CE");
  acc->GetXaxis()->SetBinLabel(ibin++,Form(">=%i CE Hit",minnhits));
  acc->GetXaxis()->SetBinLabel(ibin++,"CE p>100 MeV/c");
  acc->GetXaxis()->SetBinLabel(ibin++,"CE pitch");
  acc->GetXaxis()->SetBinLabel(ibin++,"KF Track fit");
  acc->GetXaxis()->SetBinLabel(ibin++,"Fit Quality");
  acc->GetXaxis()->SetBinLabel(ibin++,Form("Livegate (t_{0}>%3.0f)",t0min));
  acc->GetXaxis()->SetBinLabel(ibin++,"Reco pitch");
  acc->GetXaxis()->SetBinLabel(ibin++,"Cosmic Rejection");
  acc->GetXaxis()->SetBinLabel(ibin++,"Momentum window");


  ibin = 1;
  racc->GetXaxis()->SetBinLabel(ibin++,"All CE");
  racc->GetXaxis()->SetBinLabel(ibin++,Form(">=%i CE Hit",minnhits));
  racc->GetXaxis()->SetBinLabel(ibin++,"CE p>100 MeV/c");
  racc->GetXaxis()->SetBinLabel(ibin++,"CE pitch");
  racc->GetXaxis()->SetBinLabel(ibin++,"KF Track fit");
  racc->GetXaxis()->SetBinLabel(ibin++,"Fit Quality");
  racc->GetXaxis()->SetBinLabel(ibin++,Form("Livegate (t_{0}>%3.0f)",t0min));
  racc->GetXaxis()->SetBinLabel(ibin++,"Reco pitch");
  racc->GetXaxis()->SetBinLabel(ibin++,"Cosmic Rejection");
  racc->GetXaxis()->SetBinLabel(ibin++,"Momentum window");
  
  ibin = 0;
  const char* binnames[11] ={"0.0","1.0","2.0","3.0","4.0","5.0","6.0","7.0","8.0","9.0","10.0"};
  trks->Project("acc",binnames[ibin++],onlyFrstSD,"",nentries,skipentries);
  trks->Project("+acc",binnames[ibin++],onlyFrstSD+nmch,"",nentries,skipentries);
  trks->Project("+acc",binnames[ibin++],onlyFrstSD+nmch+tmom,"",nentries,skipentries);
  trks->Project("+acc",binnames[ibin++],onlyFrstSD+nmch+tmom+tpitch,"",nentries,skipentries);
  trks->Project("+acc",binnames[ibin++],onlyFrstSD+addcut+nmch+tmom+tpitch+reco,"",nentries,skipentries);
  trks->Project("+acc",binnames[ibin++],onlyFrstSD+addcut+nmch+tmom+tpitch+reco+goodfit,"",nentries,skipentries);
  trks->Project("+acc",binnames[ibin++],onlyFrstSD+addcut+nmch+tmom+tpitch+reco+goodfit+livegate,"",nentries,skipentries);
  trks->Project("+acc",binnames[ibin++],onlyFrstSD+addcut+nmch+tmom+tpitch+reco+goodfit+livegate+rpitch,"",nentries,skipentries);
  trks->Project("+acc",binnames[ibin++],onlyFrstSD+addcut+nmch+tmom+tpitch+reco+rpitch+livegate+cosmic+goodfit,"",nentries,skipentries);

  TCut goodMatCorr("(fitinfo.fitmombeam-fitinfo.fitmom)>0&&(fitinfo.fitmombeam-fitinfo.fitmom)<1.0");
  if (useExtapMom==2) {
          Long64_t nSelected = trks->Project("+acc",binnames[ibin++],onlyFrstSD+addcut+nmch+tmom+tpitch+reco+rpitch+livegate+cosmic+goodfit+goodMatCorr+rmomCorr,"",nentries,skipentries);
          cout<<"nCE in window (good mat corr) = "<<nSelected<<endl;
          nSelected += trks->Project("+acc",binnames[ibin],onlyFrstSD+addcut+nmch+tmom+tpitch+reco+rpitch+livegate+cosmic+goodfit+!goodMatCorr+rmom,"",nentries,skipentries);
          cout<<"nCE in window (plus not good mat corr) = "<<nSelected<<endl;
          acc->SetBinContent(10,nSelected);
  } else if (useExtapMom==1) {
          trks->Project("+acc",binnames[ibin++],onlyFrstSD+addcut+nmch+tmom+tpitch+reco+rpitch+livegate+cosmic+goodfit+rmomCorr,"",nentries,skipentries);
  } else {
          if (pCorr>0.0){
                  TCut ezcorr("(recopar.parvec.m[3]+(recopar.parvec[1]+380.*TMath::Pi()*recopar.parvec.m[2])/recopar.parvec.m[2]*recopar.parvec[4])>=(-550-1500+500)");//
                  TCut nezcorr("(recopar.parvec.m[3]+(recopar.parvec[1]+380.*TMath::Pi()*recopar.parvec.m[2])/recopar.parvec.m[2]*recopar.parvec[4])<(-550-1500+500)");//
                  Long64_t nSelected = trks->Project("+acc",binnames[ibin++],onlyFrstSD+addcut+nmch+tmom+tpitch+reco+rpitch+livegate+cosmic+goodfit+ezcorr+rmomCorr,"",nentries,skipentries);
                  nSelected += trks->Project("+acc",binnames[ibin],onlyFrstSD+addcut+nmch+tmom+tpitch+reco+rpitch+livegate+cosmic+goodfit+nezcorr+rmom,"",nentries,skipentries);
                  acc->SetBinContent(10,nSelected);
          } else {
                  trks->Project("+acc",binnames[ibin++],onlyFrstSD+addcut+nmch+tmom+tpitch+reco+rpitch+livegate+cosmic+goodfit+rmom,"",nentries,skipentries);
          }
  }
  cout<<"nCE in window = "<<acc->GetBinContent(10)<<endl;

  double all;
  if (nGenEv>0) {
    all=nGenEv;
    acc->SetBinContent(1,all);
    racc->SetBinContent(1,all);
  } else {
    all=acc->GetBinContent(1);
  }
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


void KalPlotMomSpec(TTree* trks, TString addCut="", size_t cutType=0, int cutOption=0, double scale=1.0, double liveTime=-1, bool useExtapMom=true, double pCorr=0.0, Long64_t NEntries=-1, Long64_t skipentries=0) {

  if(!donecuts)KalCuts(cutType,cutOption,liveTime,useExtapMom,pCorr);
  TCut addcut(addCut.Data());
  if (NEntries>0 && NEntries<1000000000) {nentries = NEntries;}

  unsigned nbins(100);
  //double bmax = nbins-0.5;
  TH1F* momspec = new TH1F("momspec","Reconstructed e spectrum;MeV/c;",nbins,101,106);

  if (useExtapMom==2) {
          TCut goodMatCorr("(fitinfo.fitmombeam-fitinfo.fitmom)>0&&(fitinfo.fitmombeam-fitinfo.fitmom)<1.0");
          trks->Project("momspec","fitinfo.fitmombeam",onlyFrstSD+addcut+reco+rpitch+livegate+cosmic+goodfit+goodMatCorr,"",nentries,skipentries);
          trks->Project("+momspec","fitinfo.fitmom",onlyFrstSD+addcut+reco+rpitch+livegate+cosmic+goodfit+!goodMatCorr,"",nentries,skipentries);
  } else if (useExtapMom==1) {
          trks->Project("momspec","fitinfo.fitmombeam",onlyFrstSD+addcut+reco+rpitch+livegate+cosmic+goodfit,"",nentries,skipentries);
  } else {
          if (pCorr>0.0){
                  TCut ezcorr("(recopar.parvec.m[3]+(recopar.parvec[1]+380.*TMath::Pi()*recopar.parvec.m[2])/recopar.parvec.m[2]*recopar.parvec[4])>=(-550-1500+500)");//
                  TCut nezcorr("(recopar.parvec.m[3]+(recopar.parvec[1]+380.*TMath::Pi()*recopar.parvec.m[2])/recopar.parvec.m[2]*recopar.parvec[4])<(-550-1500+500)");//
                  trks->Project("momspec",Form("fitinfo.fitmom+%f",pCorr),onlyFrstSD+addcut+reco+rpitch+livegate+cosmic+goodfit+ezcorr,"",nentries,skipentries);
                  trks->Project("+momspec","fitinfo.fitmom",onlyFrstSD+addcut+reco+rpitch+livegate+cosmic+goodfit+nezcorr,"",nentries,skipentries);
          } else {
                  trks->Project("momspec","fitinfo.fitmom",onlyFrstSD+addcut+reco+rpitch+livegate+cosmic+goodfit,"",nentries,skipentries);
          }
  }
  momspec->Sumw2();
  momspec->Scale(scale);

  TCanvas* acan = new TCanvas("pcan","Reco Mom",1200,800);
  pcan->Clear();
  pcan->cd();
  momspec->Draw("E");
}


void KalFitRes(TTree* trks, int nGenEv=-1, int nCuts=4, int cutOption=0, TString fitOpt="L", TString addCut="", int fType=2, int useExtapMom=0, Long64_t NEntries=-1, Long64_t skipentries=0) {
          if(!donecuts)KalCuts(0,cutOption,-1,useExtapMom);
          if (NEntries>0 && NEntries<1000000000) {nentries = NEntries;}

          TF1 *fitf;
          if (fType==0) {
                  fitf = new TF1("fitf",splitgaus,-1.5,1.5,7);
                  fitf->SetParName(0,"Norm");
                  fitf->SetParName(1,"Mean");
                  fitf->SetParName(2,"SigH");
                  fitf->SetParName(3,"SigL");
                  fitf->SetParName(4,"TFH");
                  fitf->SetParName(5,"TSigH");
                  fitf->SetParName(6,"TSigL");
          } else if (fType==1) {
                  fitf = new TF1("fitf",doublegausexp,-1.5,1.5,8);
                  fitf->SetParName(0,"Norm");
                  fitf->SetParName(1,"GTailFrac");
                  fitf->SetParName(2,"ETailFrac");
                  fitf->SetParName(3,"Mean");
                  fitf->SetParName(4,"CoreSig");
                  fitf->SetParName(5,"GTailSig");
                  fitf->SetParName(6,"ETailLambda");
                  fitf->SetParName(7,"ETailPower");
          } else if (fType==2) {
                  fitf = new TF1("fitf",crystalball,-2.0,1.5,7);
                  fitf->SetParName(0,"Norm");
                  fitf->SetParName(1,"x0");
                  fitf->SetParName(2,"sigma");
                  fitf->SetParName(3,"n");
                  fitf->SetParName(4,"alpha");
                  fitf->SetParName(5,"tailfrac");
                  fitf->SetParName(6,"taillambda");
          } else if (fType==3) {
                  fitf = new TF1("fitf",doubleCrystalball,-2.0,1.5,9);
                  fitf->SetParName(0,"Norm");
                  fitf->SetParName(1,"x0");
                  fitf->SetParName(2,"sigma");
                  fitf->SetParName(3,"n");
                  fitf->SetParName(4,"alpha");
                  fitf->SetParName(5,"frac_{core}");
                  fitf->SetParName(6,"sigma_{tail}");
                  fitf->SetParName(7,"n_{tail}");
                  fitf->SetParName(8,"alpha_{tail}");
                  fitf->SetNpx(1000);
          }

  TH1F* momres[4];
  TH1F* effnorm = new TH1F("effnorm","effnorm",100,0,150);
  trks->Project("effnorm","fitinfo.momin",onlyFrstSD+mcsel,"",nentries,skipentries);

  if (nCuts<1) nCuts=1;
  if (nCuts>4) nCuts=4;
  TCanvas* rcan;
  if (nCuts==1) {
    rcan = new TCanvas("rcan","Momentum Resolution");
  } else if (nCuts==2) {
    rcan = new TCanvas("rcan","Momentum Resolution",1200,400);
    rcan->Divide(1,2);
  } else {
    rcan = new TCanvas("rcan","Momentum Resolution",1200,800);
    rcan->Divide(2,2);
  }
  //rcan->Clear();
  cout<<"nCuts "<<nCuts<<endl;
  gStyle->SetOptFit(111111);
  gStyle->SetOptStat("oumr");
  fstream outfile("Parameters.txt",ios::out);
  for(unsigned ires=0;ires<nCuts;ires++){
    rcan->cd(ires+1);
    gPad->SetLogy();
    char mname[50];
    snprintf(mname,50,"momres%i",ires);
    momres[ires] = new TH1F(mname,"momentum resolution at start of tracker;MeV/c",251,-2.5,2.5);
//  momres[ires]->SetStats(0);
    TCut quality = ncuts[ires] && t0cuts[ires] && fitcuts[ires];
    if (useMomErrCut) quality+=momcuts[ires];
    TCut final = (onlyFrstSD+reco+quality+rpitch+mcsel+TCut(addCut.Data()));
    TCut goodMatCorr("(fitinfo.fitmombeam-fitinfo.fitmom)>0&&(fitinfo.fitmombeam-fitinfo.fitmom)<1.0");
    if (useExtapMom==2) {
            trks->Project(mname,"fitinfo.fitmombeam-fitinfo.momtrackerin",final+goodMatCorr,"",nentries,skipentries);
            trks->Project(Form("+%s",mname),"fitinfo.fitmom-fitinfo.momin",final+!goodMatCorr,"",nentries,skipentries);
    } else if (useExtapMom==1) {
            trks->Project(mname,"fitinfo.fitmombeam-fitinfo.momtrackerin",final,"",nentries,skipentries);
    } else {
            trks->Project(mname,"fitinfo.fitmom-fitinfo.momin",final,"",nentries,skipentries);
    }
    double integral = momres[ires]->GetEntries()*momres[ires]->GetBinWidth(1);
    if (fType==0) {
            fitf->SetParameters(integral,0.0,0.5*momres[ires]->GetRMS(),0.5*momres[ires]->GetRMS(),0.01,2*momres[ires]->GetRMS(),2*momres[ires]->GetRMS());
            fitf->SetParLimits(5,0.1,1.0);
            fitf->SetParLimits(6,0.1,1.0);
            fitf->SetParLimits(4,0.0,0.8);
            momres[ires]->Fit("fitf","LIR");
    } else if (fType==1) {
            fitf->SetParameters(integral,0.1,0.2,0.0,0.3*momres[ires]->GetRMS(),2*momres[ires]->GetRMS(),0.5*momres[ires]->GetRMS(),0.25);
            fitf->SetParLimits(1,0.02,0.3);
            fitf->SetParLimits(2,0.02,0.3);
            fitf->SetParLimits(4,0.05,momres[ires]->GetRMS());
            fitf->SetParLimits(5,0.12,2*momres[ires]->GetRMS());
            fitf->SetParLimits(6,0.1,momres[ires]->GetRMS());
            fitf->SetParLimits(7,1,3.0);
    } else if (fType==2) {
            fitf->SetParameters(integral,0.0,0.1,1.0,1.0,0.05,0.5);
            fitf->SetParLimits(0,1,integral*4);
            fitf->SetParLimits(1,-2,0.4);
            fitf->SetParLimits(2,0.001,momres[ires]->GetRMS(1)*0.7);
            fitf->SetParLimits(3,0,30.);
            fitf->SetParLimits(4,1e-6,5);
            fitf->SetParLimits(5,0.00001,0.9);
            fitf->SetParLimits(6,0.001,10.0*momres[ires]->GetRMS());
            //momres[ires]->Fit("fitf", "LIR");
            momres[ires]->Fit("fitf");
            if (!fitOpt.IsNull()) momres[ires]->Fit("fitf",fitOpt.Data());

            outfile << "\nCut" << ires+1 << "_1CBall_x0\t" << setprecision(8) << fitf->GetParameter(1);
            outfile << "\nCut" << ires+1 << "_1CBall_x0_err\t" << setprecision(8) << fitf->GetParError(1);
            outfile << "\nCut" << ires+1 << "_1CBall_sigma\t" << setprecision(8) << fitf->GetParameter(2);
            outfile << "\nCut" << ires+1 << "_1CBall_sigma_err\t" << setprecision(8) << fitf->GetParError(2);
            outfile << "\nCut" << ires+1 << "_1CBall_frac\t" << setprecision(8) << fitf->GetParameter(5);
            outfile << "\nCut" << ires+1 << "_1CBall_frac_err\t" << setprecision(8) << fitf->GetParError(5);
            outfile << "\nCut" << ires+1 << "_1CBall_taillambda\t" << setprecision(8) << fitf->GetParameter(6);
            outfile << "\nCut" << ires+1 << "_1CBall_taillambda_err\t" << setprecision(8) << fitf->GetParError(6) << endl;
   } else if (fType==3) {

            fitf->SetParameters(integral,momres[ires]->GetMean(1)*0.5,momres[ires]->GetRMS(1)*0.25,1.0,1.0
                            ,0.8,momres[ires]->GetRMS(1),1.0,1.0);
            fitf->SetParLimits(0,1,integral*2);
            fitf->SetParLimits(1,-0.2,0.4);
            fitf->SetParLimits(2,0.001,momres[ires]->GetRMS(1)*0.7);
            fitf->SetParLimits(3,0,30.);
            fitf->SetParLimits(4,1e-4,10);
            fitf->SetParLimits(5,0,1.0);
            fitf->SetParLimits(6,0.01,momres[ires]->GetRMS(1)*3.0);
            fitf->SetParLimits(7,0,150.);
            fitf->SetParLimits(8,1e-4,20);

            //momres[ires]->Fit("fitf");
            momres[ires]->Fit("fitf");
            if (!fitOpt.IsNull()) momres[ires]->Fit("fitf",fitOpt.Data());
            //momres[ires]->Fit("fitf","LIR");
            TF1* dcballL = new TF1("dcballL",doubleCrystalball,-2.5,2.5,9);
            dcballL->SetNpx(1000);
            TF1* dcballR = new TF1("dcballR",doubleCrystalball,-2.5,2.5,9);
            dcballR->SetNpx(1000);
            dcballL->SetParameters(fitf->GetParameters());
            dcballL->SetParameter(0,fitf->GetParameter(0)*fitf->GetParameter(5));
            dcballL->SetParameter(5,1);
            dcballL->SetLineColor(kBlack);
            dcballL->SetLineWidth(1);
            dcballL->Draw("same");
            dcballR->SetParameters(fitf->GetParameters());
            dcballR->SetParameter(0,fitf->GetParameter(0)*(1.0-fitf->GetParameter(5)));
            dcballR->SetParameter(5,0);
            dcballR->SetLineColor(kGreen);
            dcballR->SetLineWidth(1);
            dcballR->Draw("same");

            /*
    double covMat[9][9];
    gMinuit->mnemat(&covMat[0][0],9);
    double corrF = covMat[2][6]/(sqrt(covMat[2][2]*covMat[6][6]));
             */
            double wSigmaCore = fitf->GetParameter(2)*fitf->GetParameter(5);
            double wSigmaTail = fitf->GetParameter(6)*(1.0-fitf->GetParameter(5));
            //double meanSigma = sqrt(wSigmaCore*wSigmaCore + wSigmaTail*wSigmaTail + 2.0*corrF*wSigmaCore*wSigmaTail);
            double meanSigma = sqrt(wSigmaCore*fitf->GetParameter(2) + wSigmaTail*fitf->GetParameter(6));

            int intFirstBin = 1 + (int)(((fitf->GetParameter(1)-3.0*meanSigma) - momres[ires]->GetXaxis()->GetXmin() )/momres[ires]->GetBinWidth(1));
            int intLastBin  = 1 + (int)(((fitf->GetParameter(1)+3.0*meanSigma) - momres[ires]->GetXaxis()->GetXmin() )/momres[ires]->GetBinWidth(1));

            float nEvtIn3sigma = momres[ires]->Integral(intFirstBin,intLastBin);

            outfile << "\nCut" << ires+1 << "_2CBall_x0\t" << setprecision(8) << fitf->GetParameter(1);
            outfile << "\nCut" << ires+1 << "_2CBall_x0_err\t" << setprecision(8) << fitf->GetParError(1);
            outfile << "\nCut" << ires+1 << "_2CBall_sigma1\t" << setprecision(8) << fitf->GetParameter(2);
            outfile << "\nCut" << ires+1 << "_2CBall_sigma1_err\t" << setprecision(8) << fitf->GetParError(2);
            outfile << "\nCut" << ires+1 << "_2CBall_frac\t" << setprecision(8) << fitf->GetParameter(5);
            outfile << "\nCut" << ires+1 << "_2CBall_frac_err\t" << setprecision(8) << fitf->GetParError(5);
            outfile << "\nCut" << ires+1 << "_2CBall_sigmatail\t" << setprecision(8) << fitf->GetParameter(6);
            outfile << "\nCut" << ires+1 << "_2CBall_sigmatail_err\t" << setprecision(8) << fitf->GetParError(6) << endl;
            outfile << "\nCut" << ires+1 << "_2CBall_meansigma\t" << setprecision(8) << meanSigma << endl;

    }
    TLine* zero = new TLine(0.0,0.0,0.0,momres[ires]->GetBinContent(momres[ires]->GetMaximumBin()));
    zero->SetLineStyle(2);
    zero->Draw();

    double keff = momres[ires]->GetEntries()/effnorm->GetEntries();
    //cout<<"norm "<<effnorm->GetEntries()<<" sel entries "<<momres[ires]->GetEntries()<<" eff "<<keff<<endl;

    TPaveText* ttext = new TPaveText(0.1,0.65,0.38,0.9,"NDC");
    ttext->SetTextSize(0.04);
    ttext->AddText("Truth Cuts");
    ttext->AddText(nmch.GetTitle());
    ttext->AddText(tmom.GetTitle());
    ttext->AddText(Form("%4.3f<tan(#lambda_{in})<%4.3f",tdlow,tdhigh+0.02));
    //ttext->AddText(tpitch.GetTitle());
    //ttext->AddText(Form("accept = %3.2f %%",(effnorm->GetEntries()/((float)nGenEv))*100.0));
    ttext->Draw();

    TPaveText* rtext = new TPaveText(0.1,0.35,0.38,0.65,"NDC");
    rtext->SetTextSize(0.04);
    rtext->AddText("Reco Cuts");
    char line[40];
    snprintf(line,80,"%4.3f<tan(#lambda)<%4.3f",tdlow,tdhigh);
    rtext->AddText(line);
    //snprintf(line,80,"t0>%5.1f nsec",t0min);
    //rtext->AddText(line);
    sprintf(line,"%s",ncutsText[ires].GetTitle());
    //sprintf(line,"%s",ncuts[ires].GetTitle());
    rtext->AddText(line);
    sprintf(line,"%s",t0cutsText[ires].GetTitle());
    //sprintf(line,"%s",t0cuts[ires].GetTitle());
    rtext->AddText(line);
    sprintf(line,"%s",momcutsText[ires].GetTitle());
    //sprintf(line,"%s",momcuts[ires].GetTitle());
    if (useMomErrCut) rtext->AddText(line);
    sprintf(line,"%s",fitcutsText[ires].GetTitle());
    //sprintf(line,"%s",fitcuts[ires].GetTitle());
    rtext->AddText(line);
    if (!addCut.IsNull()) {
            TString addCutText = addCut;
            addCutText.ReplaceAll("fitinfo.","");
            addCutText.ReplaceAll(">","\>");
            addCutText.ReplaceAll("<","\<");
            int iSubs=0;
            TString subAddCut[10];
            while (addCutText.Contains("&") && iSubs<10) {
                    Ssiz_t andpos = addCutText.First("&");
                    TSubString splitAddCut = addCutText(0,andpos);
                    subAddCut[iSubs]=splitAddCut;
                    cout<<"andpos "<< andpos<<" "<<splitAddCut<<endl;
                    if (andpos>0) {
                            rtext->AddText(subAddCut[iSubs].Data());
                            addCutText=addCutText.Remove(0,andpos+2);
                    }
                    ++iSubs;
            }
            cout<<"addCutText "<<addCutText<<endl;
      rtext->AddText(addCutText.Data());
    }
    sprintf(line,"Eff = %3.2f %%",keff*100.0);
    rtext->AddText(line);
    rtext->Draw();
    if (fType==3) {
            TPaveText* vtext = new TPaveText(0.1,0.25,0.38,0.35,"NDC");
            vtext->SetTextSize(0.04);
            vtext->AddText(Form("#bar{#sigma} = %.3f",meanSigma));
            vtext->AddText(Form("Eff (in x0#pm 3#bar{#sigma}) = %3.2f %%",(nEvtIn3sigma/effnorm->GetEntries())*100.0));
            vtext->Draw();
    }


  }

  outfile.close();
  rcan->cd(0);
}


void KalFitMomPull(TTree* trks, int nGenEv=-1, int nCuts=4, int cutOption=0, TString fitOpt="L", TString addCut="", int fType=2, int useExtapMom=0, Long64_t NEntries=-1, Long64_t skipentries=0) {
          if(!donecuts)KalCuts(0,cutOption,-1,useExtapMom);
          if (NEntries>0 && NEntries<1000000000) {nentries = NEntries;}

          TF1 *fitf;
          if (fType==0) {
                  fitf = new TF1("fitfpl","gaus",-4,4);
                  fitf->SetParName(0,"Norm");
                  fitf->SetParName(1,"Mean");
                  fitf->SetParName(2,"Sigma");
          } else /*if (fType==1)*/ {
                  fitf = new TF1("fitfpl",twoNgaus,-4,4,6);
                  fitf->SetParName(0,"Norm");
                  fitf->SetParName(1,"Mean");
                  fitf->SetParName(2,"Sigma");
                  fitf->SetParName(3,"fract_{core}");
                  fitf->SetParName(4,"Mean_{tail}");
                  fitf->SetParName(5,"Sigma_{tail}");
          }

  TH1F* mompull[4];
  TH1F* effnorm = new TH1F("effnorm","effnorm",100,0,150);
  trks->Project("effnorm","fitinfo.momin",onlyFrstSD+mcsel,"",nentries,skipentries);

  if (nCuts<1) nCuts=1;
  if (nCuts>4) nCuts=4;
  TCanvas* pllcan;
  if (nCuts==1) {
    pllcan = new TCanvas("pllcan","Momentum Pull");
  } else if (nCuts==2) {
    pllcan = new TCanvas("pllcan","Momentum Pull",1200,400);
    pllcan->Divide(1,2);
  } else {
    pllcan = new TCanvas("pllcan","Momentum Pull",1200,800);
    pllcan->Divide(2,2);
  }
  //pllcan->Clear();
  cout<<"nCuts "<<nCuts<<endl;
  gStyle->SetOptFit(111111);
  gStyle->SetOptStat("oumr");
  for(unsigned ires=0;ires<nCuts;ires++){
    pllcan->cd(ires+1);
    gPad->SetLogy();
    char mname[50];
    snprintf(mname,50,"mompull%i",ires);
    mompull[ires] = new TH1F(mname,"Momentum Pull at start of tracker",200,-5,5);
//  mompull[ires]->SetStats(0);
    TCut quality = ncuts[ires] && t0cuts[ires] && fitcuts[ires];
    if (useMomErrCut) quality += momcuts[ires];
    TCut final = (onlyFrstSD+reco+quality+rpitch+mcsel+TCut(addCut.Data()));
    TCut goodMatCorr("(fitinfo.fitmombeam-fitinfo.fitmom)>0&&(fitinfo.fitmombeam-fitinfo.fitmom)<1.0");
    if (useExtapMom==2) {
            trks->Project(mname,"(fitinfo.fitmombeam-fitinfo.momtrackerin)/fitinfo.fitmomerr",final+goodMatCorr,"",nentries,skipentries);
            trks->Project(Form("+%s",mname),"(fitinfo.fitmom-fitinfo.momin)/fitinfo.fitmomerr",final+!goodMatCorr,"",nentries,skipentries);
    } else if (useExtapMom==1) {
            trks->Project(mname,"(fitinfo.fitmombeam-fitinfo.momtrackerin)/fitinfo.fitmomerr",final,"",nentries,skipentries);
    } else {
            trks->Project(mname,"(fitinfo.fitmom-fitinfo.momin)/fitinfo.fitmomerr",final,"",nentries,skipentries);
    }
    double integral = mompull[ires]->GetEntries()*mompull[ires]->GetBinWidth(1);
    if (fType==0) {
            fitf->SetParameters(integral,0.0,0.5*mompull[ires]->GetRMS());
            mompull[ires]->Fit("fitfpl","LIR");
    } else /*if (fType==1)*/ {
            fitf->SetParameters(integral,0.0,0.5*mompull[ires]->GetRMS(),0.9,0.0,mompull[ires]->GetRMS());
            fitf->SetParLimits(0,0.5*integral,2.*integral);
            fitf->SetParLimits(1,-0.8,0.8);
            fitf->SetParLimits(2,0.0,mompull[ires]->GetRMS());
            fitf->SetParLimits(3,0,1.0);
            fitf->SetParLimits(4,-5,5);
            fitf->SetParLimits(6,0,2.0*mompull[ires]->GetRMS());
            mompull[ires]->Fit("fitfpl","LIR");
    }

    TLine* zero = new TLine(0.0,0.0,0.0,mompull[ires]->GetBinContent(mompull[ires]->GetMaximumBin()));
    zero->SetLineStyle(2);
    zero->Draw();

    double keff = mompull[ires]->GetEntries()/effnorm->GetEntries();
    //cout<<"norm "<<effnorm->GetEntries()<<" sel entries "<<mompull[ires]->GetEntries()<<" eff "<<keff<<endl;

    TPaveText* ttext = new TPaveText(0.1,0.65,0.4,0.9,"NDC");
    ttext->SetTextSize(0.04);
    ttext->AddText("Truth Cuts");
    ttext->AddText(nmch.GetTitle());
    ttext->AddText(tmom.GetTitle());
    ttext->AddText(Form("%4.3f<tan(#lambda_{in})<%4.3f",tdlow,tdhigh+0.02));
    //ttext->AddText(tpitch.GetTitle());
    //ttext->AddText(Form("accept = %3.2f %%",(effnorm->GetEntries()/((float)nGenEv))*100.0));
    ttext->Draw();
 
    TPaveText* rtext = new TPaveText(0.1,0.4,0.4,0.65,"NDC");
    rtext->SetTextSize(0.04);
    rtext->AddText("Reco Cuts");
    char line[40];
    snprintf(line,80,"%4.3f<tan(#lambda)<%4.3f",tdlow,tdhigh);
    rtext->AddText(line);
    //snprintf(line,80,"t0>%5.1f nsec",t0min);
    //rtext->AddText(line);
    sprintf(line,"%s",ncutsText[ires].GetTitle());
    //sprintf(line,"%s",ncuts[ires].GetTitle());
    rtext->AddText(line);
    sprintf(line,"%s",t0cutsText[ires].GetTitle());
    //sprintf(line,"%s",t0cuts[ires].GetTitle());
    rtext->AddText(line);
    sprintf(line,"%s",momcutsText[ires].GetTitle());
    //sprintf(line,"%s",momcuts[ires].GetTitle());
    if(useMomErrCut) rtext->AddText(line);
    //sprintf(line,"%s",fitcuts[ires].GetTitle());
    sprintf(line,"%s",fitcutsText[ires].GetTitle());
    rtext->AddText(line);
    if (!addCut.IsNull()) {
            TString addCutText = addCut;
            addCutText.ReplaceAll("fitinfo.","");
            addCutText.ReplaceAll(">","\>");
            addCutText.ReplaceAll("<","\<");
            int iSubs=0;
            TString subAddCut[10];
            while (addCutText.Contains("&") && iSubs<10) {
                    Ssiz_t andpos = addCutText.First("&");
                    TSubString splitAddCut = addCutText(0,andpos);
                    subAddCut[iSubs]=splitAddCut;
                    cout<<"andpos "<< andpos<<" "<<splitAddCut<<endl;
                    if (andpos>0) {
                            rtext->AddText(subAddCut[iSubs].Data());
                            addCutText=addCutText.Remove(0,andpos+2);
                    }
                    ++iSubs;
            }
            cout<<"addCutText "<<addCutText<<endl;
      rtext->AddText(addCutText.Data());
    }
    sprintf(line,"Eff = %3.2f %%",keff*100.0);
    rtext->AddText(line);
    rtext->Draw();

  }

  pllcan->cd(0);
}

/*
void KalFitAmbig(TTree* t, int acut=0) {
  if(!donecuts)KalCuts();
  gStyle->SetOptStat(1111);

  TCut gambig("hitinfo._mcambig==hitinfo._ambig");
  TCut bambig("hitinfo._mcambig!=hitinfo._ambig&&hitinfo._ambig!=0");
  TCut nambig("hitinfo._ambig==0");
  TCut active("hitinfo._active>0");
// apply requested cuts
  TCut quality;
  if(acut>0)
    quality = ncuts[acut] && t0cuts[acut] && momcuts[acut] && fitcuts[acut];

  TCut goodtrk = (reco+quality+mcsel);

//  TCut goodtrk ="fitinfo.fit>0";

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

  t->Project("rdg","hitinfo._mcdist",goodtrk+active+gambig);
  t->Project("rdn","hitinfo._mcdist",goodtrk+active+nambig);
  t->Project("rdb","hitinfo._mcdist",goodtrk+active+bambig);
  t->Project("rda","hitinfo._mcdist",goodtrk+active);
  t->Project("rdi","hitinfo._mcdist",goodtrk+!active);
  Double_t ntotal = rda->GetEntries();
  Double_t nright = rdg->GetEntries();
  Double_t nneutral = rdn->GetEntries();
  Double_t nwrong = rdb->GetEntries();
  rdg->Divide(rda);  
  rdn->Divide(rda);
  rdb->Divide(rda);

//  TH1F* frdg = new TH1F("frdg","True Drift radius, failed fits;radius (mm);N hits",140,-0.05,6.95);
//  TH1F* frdb = new TH1F("frdb","True Drift radius, failed fits;radius (mm);N hits",140,-0.05,6.95);
//  frdg->SetLineColor(kBlue);
//  frdb->SetLineColor(kRed);
//  frdg->SetStats(0);
//  frdb->SetStats(0);

//  t->Project("frdg","hitinfo._mcdist",mcsel+active+ambig+(!goodfit));
//  t->Project("frdb","hitinfo._mcdist",mcsel+active+(!ambig)+(!goodfit));

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
  t->Project("momres1","fitmom-mcentmom",goodtrk&&"fitinfo.fit==1");
  t->Project("momres2","fitmom-mcentmom",goodtrk&&"fitinfo.fit==2");

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

}*/

void KalFitResid(TTree* t,int cutOption=0, Long64_t NEntries=-1, Long64_t skipentries=0) {
  if(!donecuts)KalCuts(0,cutOption,-1,false);
  if (NEntries>0 && NEntries<1000000000) {nentries = NEntries;}

  TCut delta("hitinfo._mcproc==17");
  TCut ambig("hitinfo._mcambig==hitinfo._ambig");
  TCut active("fitinfo.fit==1 && hitinfo._active>0");

  TH1F* rdg = new TH1F("rdg","True Drift radius;radius (mm);N hits",140,-0.05,6.95);
  TH1F* rdb = new TH1F("rdb","True Drift radius;radius (mm);N hits",140,-0.05,6.95);
  rdg->SetLineColor(kBlue);
  rdb->SetLineColor(kRed);
  rdg->SetStats(0);
  rdb->SetStats(0);

  TH1F* rdgNrm = new TH1F("rdgNrm","True Drift radius normalized per cell dim;norm radius (%);N hits",285,-0.005,1.42);
  TH1F* rdbNrm = new TH1F("rdbNrm","True Drift radius normalized per cell dim;norm radius (%);N hits",285,-0.005,1.42);
  rdgNrm->SetLineColor(kBlue);
  rdbNrm->SetLineColor(kRed);
  rdgNrm->SetStats(0);
  rdbNrm->SetStats(0);

  TH1F* rpullg = new TH1F("rpullg","Residual Pull;Pull;N hits",100,-6,6);
  TH1F* rpullb = new TH1F("rpullb","Residual Pull;Pull;N hits",100,-6,6);
  TH1F* rpulld = new TH1F("rpulld","Residual Pull;Pull;N hits",100,-6,6);
  rpullg->SetLineColor(kBlue);
  rpullb->SetLineColor(kRed);
  rpulld->SetLineColor(kGreen);

  t->Project("rdg","hitinfo._mcdist",onlyFrstSD+mcsel+active+ambig+(!delta),"",nentries,skipentries);
  t->Project("rdb","hitinfo._mcdist",onlyFrstSD+mcsel+active+(!ambig)+(!delta),"",nentries,skipentries);

  double invNcell = 1.0/384.0;
  double scaleFactor = (1.0+TMath::Pi()*invNcell)/(1.0-TMath::Pi()*invNcell);
  double frstCellRad = 6.27*0.5;

  t->Project("rdgNrm",Form("hitinfo._mcdist/(TMath::Power(%f,hitinfo._Slayer)*%f)",scaleFactor,frstCellRad),onlyFrstSD+mcsel+active+ambig+(!delta),"",nentries,skipentries);
  t->Project("rdbNrm",Form("hitinfo._mcdist/(TMath::Power(%f,hitinfo._Slayer)*%f)",scaleFactor,frstCellRad),onlyFrstSD+mcsel+active+(!ambig)+(!delta),"",nentries,skipentries);

  t->Project("rpullg","hitinfo._resid/hitinfo._residerr",onlyFrstSD+mcsel+active+ambig+(!delta),"",nentries,skipentries);
  t->Project("rpullb","hitinfo._resid/hitinfo._residerr",onlyFrstSD+mcsel+active+(!ambig)+(!delta),"",nentries,skipentries);
  t->Project("rpulld","hitinfo._resid/hitinfo._residerr",onlyFrstSD+mcsel+active+ambig+delta,"",nentries,skipentries);

  TCanvas* residcan = new TCanvas("residcan","Residuals",1200,800/*1200*/);
  residcan->Divide(2,1/*2*/);
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

//  residcan->cd(3);
//  rdgNrm->Draw();
//  rdbNrm->Draw("same");

  residcan->cd(0);

}


void KalFitCon(TTree* t, double chi2corrCoeff = 0, Long64_t NEntries=-1, Long64_t skipentries=0) {
  if(!donecuts)KalCuts();
  if (NEntries>0 && NEntries<1000000000) {nentries = NEntries;}

  TH1F* fprob1 = new TH1F("fprob1","fit consistency ( Prob(#chi^{2},nDOF) )",250,0,1);
  TH1F* fprob2 = new TH1F("fprob2","fit consistency ( Prob(#chi^{2},nDOF) )",250,0,1);
  fprob1->SetLineColor(kBlue);
  fprob2->SetLineColor(kRed);
  fprob1->SetStats(0);
  fprob2->SetStats(0);

//  t->Project("fprob1","TMath::Prob(fitinfo.chi2,(fitinfo.nhits-5))",mcsel+"fitinfo.fit==1");
//  t->Project("fprob2","TMath::Prob(fitinfo.chi2,(fitinfo.nhits-5))",mcsel+"fitinfo.fit!=1");
  if (chi2corrCoeff==0 || chi2corrCoeff<-1.0|| chi2corrCoeff>1.0) {
          t->Project("fprob1","fitinfo.fitcon",mcsel+"fitinfo.fit==1","",nentries,skipentries);
          t->Project("fprob2","fitinfo.fitcon",mcsel+"fitinfo.fit!=1","",nentries,skipentries);
  } else {
          t->Project("fprob1",Form("TMath::Prob(fitinfo.chi2+%f*fitinfo.ndof,fitinfo.ndof)",chi2corrCoeff),mcsel+"fitinfo.fit==1","",nentries,skipentries);
          t->Project("fprob2",Form("TMath::Prob(fitinfo.chi2+%f*fitinfo.ndof,fitinfo.ndof)",chi2corrCoeff),mcsel+"fitinfo.fit!=1","",nentries,skipentries);
  }

  TH1F* fcon1 = new TH1F("fcon1","log_{10}(fit consistency)",100,-10,0);
  TH1F* fcon2 = new TH1F("fcon2","log_{10}(fit consistency)",100,-10,0);
  fcon1->SetLineColor(kBlue);
  fcon2->SetLineColor(kRed);
  fcon1->SetStats(0);
  fcon2->SetStats(0);

//  t->Project("fcon1","log10(TMath::Prob(fitinfo.chi2,(fitinfo.nhits-5)))",mcsel+"fitinfo.fit==1");
//  t->Project("fcon2","log10(TMath::Prob(fitinfo.chi2,(fitinfo.nhits-5)))",mcsel+"fitinfo.fit!=1");
  if (chi2corrCoeff==0 || chi2corrCoeff<-1.0|| chi2corrCoeff>1.0) {
          t->Project("fcon1","log10(fitinfo.fitcon)",mcsel+"fitinfo.fit==1","",nentries,skipentries);
          t->Project("fcon2","log10(fitinfo.fitcon)",mcsel+"fitinfo.fit!=1","",nentries,skipentries);
  } else {
          t->Project("fcon1",Form("log10(TMath::Prob(fitinfo.chi2+%f*fitinfo.ndof,fitinfo.ndof))",chi2corrCoeff),mcsel+"fitinfo.fit==1","",nentries,skipentries);
          t->Project("fcon2",Form("log10(TMath::Prob(fitinfo.chi2+%f*fitinfo.ndof,fitinfo.ndof))",chi2corrCoeff),mcsel+"fitinfo.fit!=1","",nentries,skipentries);
 }


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

  TH1F* rdg = new TH1F("rdg","Reco Hit Drift Radius;radius (mm);N hits",140,-0.05,6.95);
  TH1F* rdb = new TH1F("rdb","Reco Hit Drift Radius;radius (mm);N hits",140,-0.05,6.95);
  rdg->SetStats(0);
  rdb->SetStats(0);

  t->Project("rdg","hitinfo._rdrift","fitinfo.fit==1&&fitmomerr<0.15&&hitinfo._active");
  t->Project("rdb","hitinfo._rdrift","fitinfo.fit==1&&fitmomerr>0.2&&hitinfo._active");

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
  t->Project("ncha","nchits",tmom+tpitch+"fitinfo.fit>0");
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
