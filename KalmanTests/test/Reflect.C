#include "TH1F.h"
#include "TF1.h"
#include "TTree.h"
#include "TLine.h"
#include "TArrow.h"
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
  void params(TTree* ce=0);
  Reflect(TTree* ref,double momlow, double momhigh,int pdg=11);
private:
  TTree* _ref;
  TCut utd,dtd,ud0, dd0, umom,dmom, umomerr, dmomerr,ut0err,dt0err,unact,dnact,ureco,dreco,ugood,dgood,goodpair,goodmc,ufitc,dfitc, utq, dtq, updg, dpdg;
  TCut ud0low,dd0low,ud0hi,dd0hi;
  TF1* cball;
  TF1* diffcball;
  TF1* truecball;

  bool isinit;
  double difflow, diffhigh;
  double d0low,d0high;
  double tdlow,tdhigh;
  double tqlow,tqhigh;
  double _momlow, _momhigh;
  int _pdg;
};

Reflect::Reflect(TTree* ref,double momlow, double momhigh,int pdg) : _ref(ref), _momlow(momlow), _momhigh(momhigh), _pdg(pdg) 
{ init(); }

void Reflect::init() {
  difflow = -10.0;
  diffhigh = 5.0;

  tdlow = 0.57735027;
  tdhigh = 1.0;
  d0low = -80.0;
  d0high =80.0;
  tqlow = 0.4;
  tqhigh =1.5;
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
  snprintf(ctext,80,"umcent.mom>%4.3f&&umcent.mom<%4.3f",_momlow,_momhigh);
  umom = TCut(ctext);
  snprintf(ctext,80,"dtrk.mom>%4.3f&&dtrk.mom<%4.3f",_momlow,_momhigh);
  dmom = TCut(ctext);

  snprintf(ctext,80,"utrk.trkqual> %4.3f & &utrk.trkqual < %4.3f",tqlow,tqhigh);
  utq = TCut(ctext);
  snprintf(ctext,80,"dtrk.trkqual> %4.3f & &dtrk.trkqual < %4.3f",tqlow,tqhigh);
  dtq = TCut(ctext);

  snprintf(ctext,80,"utrk.pdg==%i",_pdg);
  updg = TCut(ctext);
  snprintf(ctext,80,"dtrk.pdg==%i",_pdg);
  dpdg = TCut(ctext);

  ut0err = TCut("utrk.t0err<0.9");
  dt0err = TCut("dtrk.t0err<0.9");

  umomerr = TCut("utrk.momerr<0.25");
  dmomerr = TCut("dtrk.momerr<0.25");

  unact = TCut("utrk.nactive>=25");
  dnact = TCut("dtrk.nactive>=25");

  ureco = TCut("utrk.status>0");
  dreco = TCut("dtrk.status>0");

  ufitc = TCut("utrk.con>2e-3");
  dfitc = TCut("dtrk.con>2e-3");


  ugood = umom+utq+utd+ud0+unact+updg;
  dgood = dmom+dtq+dtd+dd0+dnact+dpdg;

  goodpair = ugood+dgood;
  goodmc = TCut("umcent.mom>0.0&&dmcent.mom>0.0&&mc.pdg==utrk.pdg");

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
  TCanvas* rcan = new TCanvas("rcan","Momentum Resolution",1200,1200);
  rcan->Clear();
  rcan->Divide(2,2);
  gStyle->SetOptFit(111111);
  gStyle->SetOptStat("oumr");

  double momlo(-4.0);
  double momhigh(4.0);
  gPad->SetLogy();
  TH1F* dmomres = new TH1F("dmomres","Downstream Single Track Momentum Resolution, Tracker Entrance;Trk-MC mom (MeV/c)",201,momlo,momhigh);
  TH1F* umomres = new TH1F("umomres","Upstream Single Track Momentum Resolution, Tracker Entrance;MC-Trk mom (MeV/c)",201,momlo,momhigh);
  TH1F* dmomresxit = new TH1F("dmomresxit","Downstream Single Track Momentum Resolution, Tracker Exit;MC - Trkmom (MeV/c)",201,momlo,momhigh);
  TH1F* umomresxit = new TH1F("umomresxit","Upstream Single Track Momentum Resolution, Tracker Exit;Trk -MC mom (MeV/c)",201,momlo,momhigh);
  _ref->Project("umomres","umcent.mom-utrk.mom",goodpair+goodmc);
  cout << umomres->GetEntries() <<endl;
  _ref->Project("dmomres","dtrk.mom-dmcent.mom",goodpair+goodmc);
  cout << dmomres->GetEntries() <<endl;
  _ref->Project("umomresxit","utrkxit.mom-umcxit.mom",goodpair+goodmc);
  _ref->Project("dmomresxit","dmcxit.mom-dtrkxit.mom",goodpair+goodmc);
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

  rcan->cd(3);
  gPad->SetLogy();
  umomresxit->Fit(cball,"RQN");
  umomresxit->Fit(cball,"LRQN");
  umomresxit->Fit(cball,"LRQ");
  umomresxit->Fit(cball,"LRM");
  rcan->cd(4);
  gPad->SetLogy();
  dmomresxit->Fit(cball,"RQN");
  dmomresxit->Fit(cball,"LRQN");
  dmomresxit->Fit(cball,"LRQ");
  dmomresxit->Fit(cball,"LRM");

//  rcan->cd(2);
//  gPad->SetLogy();
//  TH1F* momresd = new TH1F("momresd","Single Track Momentum Resolution, d0 cut;MeV/c",201,momlo,momhigh);
//  _ref->Project("momresd","umcent.mom-utrk.mom",ugood+goodmc+ud0low);
//  _ref->Project("+dmomresd","dtrk.mom-dmcent.mom",dgood+goodmc+dd0low);
//  double upintd = 2*momresd->GetEntries()*momresd->GetBinWidth(1);
//  cball->SetParameters(upintd,0.0,0.1,3.0,0.8,0.02,0.3);
//  cball->SetParLimits(5,0.001,0.4);
//  cball->SetParLimits(6,0.1,momresd->GetRMS());
//  momresd->Fit("cball","LIR");
}

void Reflect::diffres() {
  if(!isinit)init();
  TH1F* momdiff = new TH1F("momdiff","Reco Downstream - Upstream momentum",201,difflow,diffhigh);
  TH1F* mcmomdiff = new TH1F("mcmomdiff","True Downstream - Upstream momentum",201,difflow,diffhigh);
  TH1F* dmomdiff = new TH1F("dmomdiff","Reco Downstream - Upstream momentum, d0 cut",201,difflow,diffhigh);
  TH1F* dmcmomdiff = new TH1F("dmcmomdiff","True Downstream - Upstream momentum, d0 cut",201,difflow,diffhigh);

  TCut notarget("abs(utrk.d0>100)||(abs(uz0>200)&&utrk.td>-.7)");

  _ref->Project("momdiff","dtrk.mom-utrk.mom",ugood+dgood+goodmc);
  _ref->Project("mcmomdiff","dmcent.mom-umcent.mom",ugood+dgood+goodmc);

  _ref->Project("dmomdiff","dtrk.mom-utrk.mom",ugood+dgood+goodmc+dd0low+ud0low);
  _ref->Project("dmcmomdiff","dmcent.mom-umcent.mom",ugood+dgood+goodmc+dd0low+ud0low);

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

void Reflect::params( TTree* ce) {
  TH1F* crdtd = new TH1F("crdtd","Downstream tan(#lambda)",100,0.2,1.2);
  TH1F* crdmom = new TH1F("crdmom","Downstream Momentum",100,30,150);
  TH1F* crdd0 = new TH1F("crdd0","Downstream d0;mm",100,-500,500);
  TH1F* crdp0 = new TH1F("crdp0","Downstream #phi0",100,-3.142,3.142);
  TH1F* crdtq = new TH1F("crdtq","Downstream TrkQual",100,-0.2,1.5);
  TH1F* crdna = new TH1F("crdna","Downstream N Active",100,-0.5,99.5);
  crdtd->SetLineColor(kBlue);
  crdmom->SetLineColor(kBlue);
  crdd0->SetLineColor(kBlue);
  crdp0->SetLineColor(kBlue);
  crdtq->SetLineColor(kBlue);
  crdna->SetLineColor(kBlue);
  crdtd->SetStats(0);
  crdmom->SetStats(0);
  crdd0->SetStats(0);
  crdp0->SetStats(0);
  crdtq->SetStats(0);
  crdna->SetStats(0);
  _ref->Project("crdtd","dtrk.td",dreco+dpdg+goodmc);
  _ref->Project("crdmom","dtrk.mom",dreco+dpdg+goodmc);
  _ref->Project("crdd0","dtrk.d0",dreco+dpdg+goodmc);
  _ref->Project("crdp0","dtrk.p0",dreco+dpdg+goodmc);
  _ref->Project("crdtq","dtrk.trkqual",dreco+dpdg+goodmc);
  _ref->Project("crdna","dtrk.nactive",dreco+dpdg+goodmc);
  
  unsigned ncr = crdtd->GetEntries();

  TH1F* cetd(0);
  TH1F* cemom(0);
  TH1F* ced0(0);
  TH1F* cep0(0);
  TH1F* cetq(0);
  TH1F* cena(0);
  TLegend* leg = new TLegend(0.6,0.7,0.9,0.9);
  leg->AddEntry(crdtd,"CR Branch","L");
  if(ce != 0){
    cetd = new TH1F("cetd","Downstream tan(#lambda)",100,0.2,1.2);
    cemom = new TH1F("cemom","Downstream Momentum",100,30,150);
    ced0 = new TH1F("ced0","Downstream d0;mm",100,-500,500);
    cep0 = new TH1F("cep0","Downstream #phi0",100,-3.142,3.142);
    cetq = new TH1F("cetq","Downstream TrkQual",100,-0.2,1.5);
    cena = new TH1F("cena","Downstream N Active",100,-0.5,99.5);
    cetd->SetLineColor(kRed);
    cemom->SetLineColor(kRed);
    ced0->SetLineColor(kRed);
    cep0->SetLineColor(kRed);
    cetq->SetLineColor(kRed);
    cena->SetLineColor(kRed);
    cetd->SetStats(0);
    cemom->SetStats(0);
    ced0->SetStats(0);
    cep0->SetStats(0);
    cetq->SetStats(0);
    cena->SetStats(0);
    ce->Project("cetd","fit.td","fit.status>0");
    ce->Project("cemom","fit.mom","fit.status>0");
    ce->Project("ced0","fit.d0","fit.status>0");
    ce->Project("cep0","fit.p0","fit.status>0");
    ce->Project("cetq","fit.trkqual","fit.status>0");
    ce->Project("cena","fit.nactive","fit.status>0");

    unsigned nce = cetd->GetEntries();
    float scale(ncr/nce);
    cetd->Scale(scale);
    cemom->Scale(scale);
    ced0->Scale(scale);
    cep0->Scale(scale);
    cetq->Scale(scale);
    cena->Scale(scale);
    leg->AddEntry(cetd,"Conversion","L");
  }

  TLine *ltdlow , *ltdhigh, *lmomlow, *lmomhigh, *ld0low, *ld0high, *ltqlow, *ltqhigh;
  if(ce != 0){
    ltdlow = new TLine(tdlow,0,tdlow,cetd->GetMaximum());
    ltdhigh= new TLine(tdhigh,0,tdhigh,cetd->GetMaximum());
    lmomlow = new TLine(_momlow,0,_momlow,cemom->GetMaximum());
    lmomhigh= new TLine(_momhigh,0,_momhigh,cemom->GetMaximum());
    ld0low = new TLine(d0low,0,d0low,ced0->GetMaximum());
    ld0high= new TLine(d0high,0,d0high,ced0->GetMaximum());
    ltqlow = new TLine(tqlow,0,tqlow,cetq->GetMaximum());
    ltqhigh= new TLine(tqhigh,0,tqhigh,cetq->GetMaximum());
  } else {
    ltdlow = new TLine(tdlow,0,tdlow,crdtd->GetMaximum());
    ltdhigh= new TLine(tdhigh,0,tdhigh,crdtd->GetMaximum());
    lmomlow = new TLine(_momlow,0,_momlow,crdmom->GetMaximum());
    lmomhigh= new TLine(_momhigh,0,_momhigh,crdmom->GetMaximum());
    ld0low = new TLine(d0low,0,d0low,crdd0->GetMaximum());
    ld0high= new TLine(d0high,0,d0high,crdd0->GetMaximum());
    ltqlow = new TLine(tqlow,0,tqlow,crdtq->GetMaximum());
    ltqhigh= new TLine(tqhigh,0,tqhigh,crdtq->GetMaximum());
  }

  TCanvas* pcan = new TCanvas("pcan","Parameters",1200,800);
  pcan->Divide(3,2);
  pcan->cd(1);
  if(cetd)cetd->Draw("same");
  crdtd->Draw("same");
  leg->Draw("same");
  ltdlow->Draw();
  ltdhigh->Draw();
  pcan->cd(2);
  if(cemom)cemom->Draw("same");
  crdmom->Draw("same");
  lmomlow->Draw();
  lmomhigh->Draw();
  pcan->cd(3);
  if(ced0)ced0->Draw("same");
  crdd0->Draw("same");
  ld0low->Draw();
  ld0high->Draw();
  pcan->cd(4);
  crdtq->Draw("same");
  if(cetq)cetq->Draw("same");
  ltqlow->Draw();
  ltqhigh->Draw();
  pcan->cd(5);
  if(cena)cena->Draw("same");
  crdna->Draw("same");
  pcan->cd(6);
  crdp0->Draw("same");
  if(cep0)cep0->Draw("same");

}
