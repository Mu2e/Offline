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


TCut utd,dtd,umom,dmom, umomerr, dmomerr,ut0err,dt0err,unact,dnact,ureco,dreco,ugood,dgood,goodpair,goodmc;
TCut ud0,dd0;
TF1* cball;
TF1* diffcball;

double momlo(-4.0);
double momhi(4.0);

double difflo(-10.0);
double diffhi(5.0);

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

void init() {
  double tdlow(0.57735027);
  double tdhigh(1.0);
  char ctext[80];
  snprintf(ctext,80,"utd<%4.3f&&utd>%4.3f",-tdlow,-tdhigh);
  utd = TCut(ctext);
  snprintf(ctext,80,"dtd>%4.3f&&dtd<%4.3f",tdlow,tdhigh);
  dtd = TCut(ctext);

  double momlow(80);
  double momhigh(120);
  snprintf(ctext,80,"umom>%4.3f&&umom<%4.3f",momlow,momhigh);
  umom = TCut(ctext);
  snprintf(ctext,80,"dmom>%4.3f&&dmom<%4.3f",momlow,momhigh);
  dmom = TCut(ctext);

  ut0err = TCut("ut0err<0.9");
  dt0err = TCut("dt0err<0.9");

  umomerr = TCut("umomerr<0.18");
  dmomerr = TCut("dmomerr<0.18");

  unact = TCut("unactive>=25");
  dnact = TCut("dnactive>=25");

  ureco = TCut("ufitstat>0");
  dreco = TCut("dfitstat>0");

  ugood = ureco+unact+umomerr+ut0err+umom+utd;
  dgood = dreco+dnact+dmomerr+dt0err+dmom+dtd;

  goodpair = ugood+dgood;
  goodmc = TCut("umcmom>0.0&&dmcmom>0.0");

  ud0 = TCut("ud0>100");
  dd0 = TCut("dd0>100");

  cball = new TF1("cball",crystalball,-2.0,2.0,7);
  cball->SetParName(0,"Norm");
  cball->SetParName(1,"x0");
  cball->SetParName(2,"sigma");
  cball->SetParName(3,"n");
  cball->SetParName(4,"alpha");
  cball->SetParName(5,"tailfrac");
  cball->SetParName(6,"taillambda");

  diffcball = new TF1("diffcball",crystalball,-5.0,1.,7);
  diffcball->SetParName(0,"Norm");
  diffcball->SetParName(1,"x0");
  diffcball->SetParName(2,"sigma");
  diffcball->SetParName(3,"n");
  diffcball->SetParName(4,"alpha");
  diffcball->SetParName(5,"tailfrac");
  diffcball->SetParName(6,"taillambda");

}


void momfit(TTree* ref) {
  TCanvas* rcan = new TCanvas("rcan","Momentum Resolution",1200,800);
  rcan->Clear();
  rcan->Divide(2,2);
  gStyle->SetOptFit(111111);
  gStyle->SetOptStat("oumr");

  rcan->cd(1);
  gPad->SetLogy();
  TH1F* umomres = new TH1F("umomres","Upstream momentum resolution;MeV",251,momlo,momhi);
  ref->Project("umomres","umcmom-umom",ugood+goodmc);
  double uint = umomres->GetEntries()*umomres->GetBinWidth(1);
  cball->SetParameters(uint,0.0,0.1,3.0,0.8,0.02,0.3);
  cball->SetParLimits(5,0.001,0.4);
  cball->SetParLimits(6,0.1,umomres->GetRMS());
  umomres->Fit("cball","LIR");

  rcan->cd(2);
  gPad->SetLogy();
  TH1F* dmomres = new TH1F("dmomres","Downstream momentum resolution;MeV",251,momlo,momhi);
  ref->Project("dmomres","dmom-dmcmom",dgood+goodmc);
  double dint = dmomres->GetEntries()*dmomres->GetBinWidth(1);
  cball->SetParameters(dint,0.0,0.1,3.0,0.8,0.02,0.3);
  cball->SetParLimits(5,0.001,0.4);
  cball->SetParLimits(6,0.1,dmomres->GetRMS());
  dmomres->Fit("cball","LIR");

  rcan->cd(3);
  gPad->SetLogy();
  TH1F* umomresd = new TH1F("umomresd","Upstream momentum resolution, d0 cut;MeV",251,momlo,momhi);
  ref->Project("umomresd","umcmom-umom",ugood+goodmc+ud0);
  double uintd = umomresd->GetEntries()*umomresd->GetBinWidth(1);
  cball->SetParameters(uintd,0.0,0.1,3.0,0.8,0.02,0.3);
  cball->SetParLimits(5,0.001,0.4);
  cball->SetParLimits(6,0.1,umomresd->GetRMS());
  umomresd->Fit("cball","LIR");

  rcan->cd(4);
  gPad->SetLogy();
  TH1F* dmomresd = new TH1F("dmomresd","Downstream momentum resolution, d0 cut;MeV",251,momlo,momhi);
  ref->Project("dmomresd","dmom-dmcmom",dgood+goodmc+dd0);
  double dintd = dmomresd->GetEntries()*dmomresd->GetBinWidth(1);
  cball->SetParameters(dintd,0.0,0.1,3.0,0.8,0.02,0.3);
  cball->SetParLimits(5,0.001,0.4);
  cball->SetParLimits(6,0.1,dmomresd->GetRMS());
  dmomresd->Fit("cball","LIR");

}

void diffres(TTree* ref) {
  TH1F* momdiff = new TH1F("momdiff","Reco Downstream - Upstream momentum",251,difflo,diffhi);
  TH1F* mcmomdiff = new TH1F("mcmomdiff","True Downstream - Upstream momentum",251,difflo,diffhi);
  TH1F* dmomdiff = new TH1F("dmomdiff","Reco Downstream - Upstream momentum, d0 cut",251,difflo,diffhi);
  TH1F* dmcmomdiff = new TH1F("dmcmomdiff","True Downstream - Upstream momentum, d0 cut",251,difflo,diffhi);

  ref->Project("momdiff","dmom-umom",ugood+dgood+goodmc);
  ref->Project("mcmomdiff","dmcmom-umcmom",ugood+dgood+goodmc);

  ref->Project("dmomdiff","dmom-umom",ugood+dgood+goodmc+ud0+dd0);
  ref->Project("dmcmomdiff","dmcmom-umcmom",ugood+dgood+goodmc+ud0+dd0);

  TCanvas* dcan = new TCanvas("dcan","Momentum Difference",1200,800);
  dcan->Clear();
  dcan->Divide(2,2);

  dcan->cd(1);
  gPad->SetLogy();
  double difint = momdiff->GetEntries()*momdiff->GetBinWidth(1);
  diffcball->SetParameters(difint,0.0,0.1,10.0,1.0,0.05,0.5);
  diffcball->SetParLimits(5,0.000,0.4);
  diffcball->SetParLimits(6,0.01,momdiff->GetRMS());
  momdiff->Fit("diffcball","LIR");

  dcan->cd(2);
  gPad->SetLogy();
  double mcdifint = mcmomdiff->GetEntries()*mcmomdiff->GetBinWidth(1);
  diffcball->SetParameters(mcdifint,0.0,0.1,10.0,1.0,0.05,0.5);
  diffcball->SetParLimits(5,0.00,0.4);
  diffcball->SetParLimits(6,0.01,mcmomdiff->GetRMS());
  mcmomdiff->Fit("diffcball","LIR");

  dcan->cd(3);
  gPad->SetLogy();
  double ddifint = dmomdiff->GetEntries()*dmomdiff->GetBinWidth(1);
  diffcball->SetParameters(ddifint,0.0,0.1,10.0,1.0,0.05,0.5);
  diffcball->SetParLimits(5,0.00,0.4);
  diffcball->SetParLimits(6,0.01,dmomdiff->GetRMS());
  dmomdiff->Fit("diffcball","LIR");

  dcan->cd(4);
  gPad->SetLogy();
  double dmcdifint = dmcmomdiff->GetEntries()*dmcmomdiff->GetBinWidth(1);
  diffcball->SetParameters(dmcdifint,0.0,0.1,10.0,1.0,0.05,0.5);
  diffcball->SetParLimits(5,0.00,0.4);
  diffcball->SetParLimits(6,0.01,dmcmomdiff->GetRMS());
  dmcmomdiff->Fit("diffcball","LIR");

}
