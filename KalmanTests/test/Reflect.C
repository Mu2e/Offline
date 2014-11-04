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
class Reflect {
public:
  void init();
  void momres();
  void diffres();
  Reflect(TTree* ref,double momlow, double momhigh);
private:
  TTree* _ref;
  TCut utd,dtd,ud0, dd0, umom,dmom, umomerr, dmomerr,ut0err,dt0err,unact,dnact,ureco,dreco,ugood,dgood,goodpair,goodmc,ufitc,dfitc, utq, dtq;
  TCut ud0low,dd0low,ud0hi,dd0hi;
  TF1* cball;
  TF1* diffcball;
  TF1* truecball;

  bool isinit;
  double momlo;
  double momhi;
  double difflo;
  double diffhi;

  double _momlow, _momhigh;
};

Reflect::Reflect(TTree* ref,double momlow, double momhigh) : _ref(ref), _momlow(momlow), _momhigh(momhigh) 
{ init(); }

void Reflect::init() {
  momlo = -4.0;
  momhi = 4.0;
  difflo = -10.0;
  diffhi = 5.0;

  double tdlow(0.57735027);
  double tdhigh(1.0);
  double d0low(-80.0);
  double d0high(80.0);
  char ctext[80];
  snprintf(ctext,80,"utrk.td<%4.3f&&utrk.td>%4.3f",-tdlow,-tdhigh);
  utd = TCut(ctext);
  snprintf(ctext,80,"dtrk.td>%4.3f&&dtrk.td<%4.3f",tdlow,tdhigh);
  dtd = TCut(ctext);
  snprintf(ctext,80,"utrk.d0>%4.3f&&utrk.d0<%4.3f",d0low,d0high);
  ud0 = TCut(ctext);
  snprintf(ctext,80,"dtrk.d0>%4.3f&&dtrk.d0<%4.3f",d0low,d0high);
  dd0 = TCut(ctext);
// note upstream selection is on MC, downstream on fit, due to inverted direction
  snprintf(ctext,80,"umc.mom>%4.3f&&umc.mom<%4.3f",_momlow,_momhigh);
  umom = TCut(ctext);
  snprintf(ctext,80,"dtrk.fitmom>%4.3f&&dtrk.fitmom<%4.3f",_momlow,_momhigh);
  dmom = TCut(ctext);

  ut0err = TCut("utrk.t0err<0.9");
  dt0err = TCut("dtrk.t0err<0.9");

  umomerr = TCut("utrk.fitmomerr<0.25");
  dmomerr = TCut("dtrk.fitmomerr<0.25");

  unact = TCut("utrk.nactive>=25");
  dnact = TCut("dtrk.nactive>=25");

  ureco = TCut("utrk.fitstat>0");
  dreco = TCut("dtrk.fitstat>0");

  ufitc = TCut("utrk.fitcon>2e-3");
  dfitc = TCut("dtrk.fitcon>2e-3");

  utq = TCut("utrk.trkqual>0.4");
  dtq = TCut("dtrk.trkqual>0.4");

  ugood = umom+utq+utd+ud0;
  dgood = dmom+dtq+dtd+ud0;

  goodpair = ugood+dgood;
  goodmc = TCut("umc.mom>0.0&&dmc.mom>0.0&&mc.pdg==11&&mc.pdg==utrk.fitpart");

  ud0low = TCut("utrk.d0<100");
  dd0low = TCut("dtrk.d0<100");
  ud0hi = TCut("utrk.d0>100");
  dd0hi = TCut("dtrk.d0>100");

  cball = new TF1("cball",crystalball,-2.0,2.0,7);
  cball->SetParName(0,"Norm");
  cball->SetParName(1,"x0");
  cball->SetParName(2,"sigma");
  cball->SetParName(3,"n");
  cball->SetParName(4,"alpha");
  cball->SetParName(5,"tailfrac");
  cball->SetParName(6,"taillambda");

  diffcball = new TF1("diffcball",crystalball,-5.0,2.,7);
  diffcball->SetParName(0,"Norm");
  diffcball->SetParName(1,"x0");
  diffcball->SetParName(2,"sigma");
  diffcball->SetParName(3,"n");
  diffcball->SetParName(4,"alpha");
  diffcball->SetParName(5,"tailfrac");
  diffcball->SetParName(6,"taillambda");
 
  truecball = new TF1("truecball",crystalball,-5.0,2.,7);
  truecball->SetParName(0,"Norm");
  truecball->SetParName(1,"x0");
  truecball->SetParName(2,"sigma");
  truecball->SetParName(3,"n");
  truecball->SetParName(4,"alpha");
  truecball->SetParName(5,"tailfrac");
  truecball->SetParName(6,"taillambda");
  truecball->FixParameter(5,0.0);
  truecball->FixParameter(6,0.2);

  isinit = true;
}


void Reflect::momres() {
  if(!isinit)init();
  TCanvas* rcan = new TCanvas("rcan","Momentum Resolution",1200,500);
  rcan->Clear();
  rcan->Divide(2,1);
  gStyle->SetOptFit(111111);
  gStyle->SetOptStat("oumr");

  gPad->SetLogy();
  TH1F* dmomres = new TH1F("dmomres","Downstream Single Track Momentum Resolution;fitmom-mcmom (MeV/c)",201,momlo,momhi);
  TH1F* umomres = new TH1F("umomres","Upstream Single Track Momentum Resolution;mcmom-fitmom (MeV/c)",201,momlo,momhi);
  _ref->Project("umomres","umc.mom-utrk.fitmom",ugood+goodmc);
  cout << umomres->GetEntries() <<endl;
  _ref->Project("dmomres","dtrk.fitmom-dmc.mom",dgood+goodmc);
  cout << dmomres->GetEntries() <<endl;
  double upint = 2*umomres->GetEntries()*umomres->GetBinWidth(1);
  cball->SetParameters(upint,0.0,0.1,3.0,0.8,0.02,0.3);
  cball->SetParLimits(5,0.001,0.4);
  cball->SetParLimits(6,0.1,umomres->GetRMS());
  rcan->cd(1);
  gPad->SetLogy();
  umomres->Fit(cball,"RQN");
  umomres->Fit(cball,"LRQN");
  umomres->Fit(cball,"LRQ");
  umomres->Fit(cball,"LRM");
  rcan->cd(2);
  gPad->SetLogy();
  dmomres->Fit(cball,"RQN");
  dmomres->Fit(cball,"LRQN");
  dmomres->Fit(cball,"LRQ");
  dmomres->Fit(cball,"LRM");

//  rcan->cd(2);
//  gPad->SetLogy();
//  TH1F* momresd = new TH1F("momresd","Single Track Momentum Resolution, d0 cut;MeV/c",201,momlo,momhi);
//  _ref->Project("momresd","umc.mom-utrk.fitmom",ugood+goodmc+ud0low);
//  _ref->Project("+dmomresd","dtrk.fitmom-dmc.mom",dgood+goodmc+dd0low);
//  double upintd = 2*momresd->GetEntries()*momresd->GetBinWidth(1);
//  cball->SetParameters(upintd,0.0,0.1,3.0,0.8,0.02,0.3);
//  cball->SetParLimits(5,0.001,0.4);
//  cball->SetParLimits(6,0.1,momresd->GetRMS());
//  momresd->Fit("cball","LIR");
}

void Reflect::diffres() {
  if(!isinit)init();
  TH1F* momdiff = new TH1F("momdiff","Reco Downstream - Upstream momentum",201,difflo,diffhi);
  TH1F* mcmomdiff = new TH1F("mcmomdiff","True Downstream - Upstream momentum",201,difflo,diffhi);
  TH1F* dmomdiff = new TH1F("dmomdiff","Reco Downstream - Upstream momentum, d0 cut",201,difflo,diffhi);
  TH1F* dmcmomdiff = new TH1F("dmcmomdiff","True Downstream - Upstream momentum, d0 cut",201,difflo,diffhi);

  TCut notarget("abs(utrk.d0>100)||(abs(uz0>200)&&utrk.td>-.7)");

  _ref->Project("momdiff","dtrk.fitmom-utrk.fitmom",ugood+dgood+goodmc);
  _ref->Project("mcmomdiff","dmc.mom-umc.mom",ugood+dgood+goodmc);

  _ref->Project("dmomdiff","dtrk.fitmom-utrk.fitmom",ugood+dgood+goodmc+dd0low+ud0low);
  _ref->Project("dmcmomdiff","dmc.mom-umc.mom",ugood+dgood+goodmc+dd0low+ud0low);

  TCanvas* dcan = new TCanvas("dcan","Momentum Difference",1200,800);
  dcan->Clear();
  dcan->Divide(2,2);

  dcan->cd(1);
  gPad->SetLogy();
  double difint = momdiff->GetEntries()*momdiff->GetBinWidth(1);
  double mean = momdiff->GetMean();
  double rms = momdiff->GetRMS();
  diffcball->SetParameters(difint,mean,0.4,10.0,1.0,0.01,0.5);
  diffcball->SetParLimits(5,0.000,0.4);
  diffcball->SetParLimits(6,0.01,momdiff->GetRMS());
  momdiff->Fit("diffcball","LIR");

  dcan->cd(2);
  gPad->SetLogy();
  double mcdifint = mcmomdiff->GetEntries()*mcmomdiff->GetBinWidth(1);
  mean = mcmomdiff->GetMean();
  truecball->SetParameters(mcdifint,mean,0.4,10.0,1.0,0.0,0.5);
  mcmomdiff->Fit("truecball","LIR");

  dcan->cd(3);
  gPad->SetLogy();
  double ddifint = dmomdiff->GetEntries()*dmomdiff->GetBinWidth(1);
  mean = dmomdiff->GetMean();
  diffcball->SetParameters(ddifint,mean,0.28,3.5,0.65,0.01,0.5);
  diffcball->SetParLimits(5,0.00,0.4);
  diffcball->SetParLimits(6,0.01,dmomdiff->GetRMS());
  dmomdiff->Fit("diffcball","LIR");

  dcan->cd(4);
  gPad->SetLogy();
  double dmcdifint = dmcmomdiff->GetEntries()*dmcmomdiff->GetBinWidth(1);
  mean = dmomdiff->GetMean();
  truecball->SetParameters(dmcdifint,mean,0.25,3.5,0.65,0.01,0.5);
  dmcmomdiff->Fit("truecball","LIR");

}
