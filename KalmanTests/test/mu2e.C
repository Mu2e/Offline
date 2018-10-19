#include "TTree.h"
#include "TStyle.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TCut.h"
#include "TPaveText.h"
#include "TLine.h"
#include "TArrow.h"
#include "TRandom3.h"
#include "TMath.h"
#include <iostream>
#include <math.h>
#include <vector>
#include "KalmanTests/test/DIOCZ.h"

using namespace std;

Double_t DIOCZ_R(Double_t *x, Double_t *par) {
  double norm = par[0];
  double eloss = par[1];
  double acc0 = par[2];
  double acc1 = par[3];
  double acc2 = par[4];
  static const double mal(25133);
  //    double mmu(105.654);
  static const double emu(105.194);
//  static const double emue(104.973);
  //    double me(0.511);
  static const double a4(1e-17);
  static const double a5(8.6434e-17);
  static const double a6(1.16874e-17);
  static const double a7(-1.87828e-19);
  static const double a8(9.16327e-20);
  double ee = x[0];
  double delta = emu + eloss - ee - ee*ee/(2*mal);
  if(delta>0.0)
    return norm*(acc0+acc1*(ee-100.0) + a4*acc2*pow(delta,3))*(a5*pow(delta,5) + a6*pow(delta,6) + a7*pow(delta,7) + a8*pow(delta,8));
  else
    return 0.0;
}

Double_t RECODIO(Double_t *x, Double_t *par) {
  double ee = x[0];
  double mal(25133);
  double emu = 105.194 - par[0];
  double delta = emu - ee - ee*ee/(2*mal);
  if(delta > 0.0)
    return par[1]*(pow(delta,5)+par[2]*pow(delta,6)+par[3]*pow(delta,7)+par[4]*pow(delta,8));
  else
    return 0.0;
}

Double_t RECOFITDIO(Double_t *x, Double_t *par) {
  double ee = x[0];
  double mal(25133);
  double emu = 105.194 - par[0] - par[5];
  double delta = emu - ee - ee*ee/(2*mal);
  if(delta > 0.0)
    return par[6]*par[1]*(pow(delta,5)+par[2]*pow(delta,6)+par[3]*pow(delta,7)+par[4]*pow(delta,8));
  else
    return 0.0;
}

// reconstruction acceptance, as a function of the true DIO momentum
Double_t recoacc(Double_t *x, Double_t *par) {
  double accmid = par[0];
  double accslope = par[1];
  double mom = x[0];
  static const double midmom(100.0);
  return accmid + (mom-midmom)*accslope;
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


class mu2e {
  public:
    mu2e(TTree* d, TTree* c, double dgenrange, double nd, double nc,bool weightd=true,double np=3.6e20,double mustopfrac=1.87e-3) : dio(d), con(c),diogenrange(dgenrange),
    ndio(nd),ncon(nc),weightdio(weightd),nproton(np),nstopped(np*mustopfrac),capfrac(0.609),rmue(1e-16),trueconvmom(104.973),
    tdlow(0.57735027),tdhigh(1.0),t0min(700),t0max(1695),rpc(0.025), ap(0.083333),cmu(0.041666),mu2ecut(2),
    reco("dem.status>0"),_init(false),nactive("dem.nactive>=20")
  {
  }
    void init();

    void fillmu2e(unsigned nbins=151,double mmin=101.0,double mmax=106.0);
    void drawmu2e(double momlow,double momhigh,bool log,unsigned ilow=0,unsigned ihi=3,const char* suffix=".png");
    void drawdio(double momlow,double momhigh,const char* suffix=".png");
    void smearDIO(unsigned ntrials=1e5,unsigned nres=1e5);
    void doExperiments(double momlow, double momhigh,double cprob,unsigned ispec,unsigned nexp=18, unsigned npave=0);
    void fitReco(unsigned icut);
    TTree *dio, *con;
    double diogenrange;
    double ndio, ncon;
    bool weightdio;
    unsigned _nbins;
    double _mmin, _mmax;
    double nproton,nstopped,capfrac,rmue;
    double decayfrac,ndecay,ncap,mevperbin,conscale;
    double trueconvmom,tdlow,tdhigh,t0min,t0max;
    double rpc, ap, cmu,flat; // rates per MeV
    double dioint,dioscale;
    unsigned mu2ecut;
    TCut reco, pitch, livegate, cosmic;
    TCut pid, upstream;
    TCut quality[4], final[4];
    TCut mcdio, mccon;
    TF1* _diocz_f;
    TF1* _cball;
    TF1* _dscb;
    TF1* _racc;
    TF1* _reco_f;
//    TF1* _flat_f[4];
    TH1F* _diospec[4];
    TH1F* _conspec[4];
    TH1D* _recodio;
    TH1F* _momres;
    TLegend* _leg;
    TPaveText* _info;
    bool _init;
    double _conint[4];
    double _conint_err[4];
    TCut nactive;
};

void mu2e::init(){
  if(_init)return;
  using namespace std;
  decayfrac = 1.0 - capfrac;
  ndecay = nstopped*decayfrac;
  ncap = nstopped*capfrac;
  conscale = ncap*rmue/ncon;
  cout << "Conversion scale factor =" << conscale << endl;
  // dio spectrum
  _diocz_f = new TF1("_diocz_f",DIOCZ,95.0,trueconvmom,1);
  _diocz_f->SetLineColor(kGreen);
  _diocz_f->SetParameter(0,1.0);
  // integrate the DIO spectrum over the range specified.  This is relative to the free decay rate
  dioint = _diocz_f->Integral(trueconvmom-diogenrange,trueconvmom);
  dioscale = 0.0;
  if(ndio>0){
    if(weightdio){ 
      cout << "Weighting DIO" << endl;
      dioscale =ndecay*diogenrange/ndio;
    } else {
      cout << "Not weighting DIO" << endl;
      dioscale = dioint*ndecay/ndio;
    }
  }
  cout << "DIO scale factor = " << dioscale << endl;

//  flat = rpc+ap+cmu;
//  cout << "Flat rate = " << flat << " counts/MeV/c" << endl;
  // basic cuts
  char ctext[80];
  snprintf(ctext,80,"dem.td>%f&&dem.td<%f",tdlow,tdhigh);
  pitch = TCut(ctext);
  snprintf(ctext,80,"dem.t0>%f&&dem.t0<%f",t0min,t0max);
  livegate = TCut(ctext);
  cosmic = TCut("dem.d0<105&&dem.d0>-80 && (dem.d0+2/dem.om)>450 && (dem.d0+2/dem.om)<680");
  pid = TCut("demc.eclust>10.0&&demc.eclust<120&&demc.uvchisq<100.0&&max(demc.dtllr,0.0)+max(demc.epllr,0.0)>1.5");
  upstream = TCut("tcnt.nuem<=0");
  // insure this is the primary generated DIO particle, not something from the background frame
  mcdio = TCut("demmc.gen==28");
  mccon = TCut("demmc.gen==2");
  double trkqualcut[4] = {0.3,0.4,0.5,0.6};
  // cuts for different tightness of selection
  for(unsigned icut=0;icut<4;icut++){
    snprintf(ctext,80,"dem.trkqual>%f",trkqualcut[icut]);
    quality[icut] = TCut(ctext);
    final[icut] = (reco+pitch+livegate+quality[icut]+cosmic+pid+upstream);
    cout << "final cut  " << icut << " = " <<  final[icut].GetTitle() << endl;
  }
  _init = true;
}

void mu2e::fillmu2e(unsigned nbins,double mmin,double mmax) {
  init();
  _nbins = nbins;
  _mmin = mmin;
  _mmax = mmax;
  mevperbin = (_mmax-_mmin)/_nbins;
  for(unsigned icut=0;icut<4;icut++){
    char dioname[50];
    snprintf(dioname,50,"diospec%i",icut);
    char conname[50];
    snprintf(conname,50,"conspec%i",icut);
//    char flatname[50];
//    snprintf(flatname,50,"flat_f%i",icut);

    char dtitle[100];
    snprintf(dtitle,100,"Reconstructed e^{-} Momentum;p (MeV/c);N Events/%2.2f MeV/c",mevperbin);
    _diospec[icut] = new TH1F(dioname,dtitle,_nbins,_mmin,_mmax);
    _diospec[icut]->SetStats(0);
    _diospec[icut]->SetLineColor(kBlue);
    _diospec[icut]->Sumw2();

    _conspec[icut] = new TH1F(conname,dtitle,_nbins,_mmin,_mmax);
    _conspec[icut]->SetStats(0);
    _conspec[icut]->SetLineColor(kRed);
    _conspec[icut]->Sumw2();

    if(dio){
      dio->Project(dioname,"dem.mom","evtinfo.evtwt"*final[icut]);
      _diospec[icut]->Scale(dioscale);
      _diospec[icut]->SetMinimum(0.08*_diocz_f->Eval(trueconvmom-0.1)*ndecay*mevperbin);
      _diospec[icut]->SetMaximum(0.08*_diocz_f->Eval(_mmin)*ndecay*mevperbin);
    }
    con->Project(conname,"dem.mom","evtinfo.evtwt"*final[icut]);
    _conspec[icut]->Scale(conscale);
    
//    _flat_f[icut] = new TF1(flatname,"[0]",_mmin,_mmax);
//    _flat_f[icut]->SetLineColor(kGreen);
//    double acc(1.0);
//    if(icut>0)acc = _conspec[icut]->Integral()/_conspec[0]->Integral();
//    _flat_f[icut]->SetParameter(0,flat*mevperbin*acc);
  }

  _leg = new TLegend(0.6,0.7,0.9,0.9);
  _leg->AddEntry(_diospec[0],"DIO","L");
  _leg->AddEntry(_conspec[0],"Conversion","L");
//  _leg->AddEntry(_flat_f[0],"RPC+AP+cosmic","L");

  _info = new TPaveText(0.4,0.8,0.7,0.9,"NDC");
  char text[80];
  snprintf(text,80,"%5.2e stopped #mu^{-}",nstopped);
  TString snstop(text);
  _info->AddText(snstop);
  snprintf(text,80,"%e Conversion Rate",rmue);
  TString srmue(text);
  _info->AddText(srmue);
  _info->SetBorderSize(0);
}

void mu2e::drawmu2e(double momlow, double momhigh,bool logy,unsigned ilow,unsigned ihi,const char* suffix) {
  char ctext[80];
  snprintf(ctext,80,"dem.mom>%f&&dem.mom<%f",momlow,momhigh);
  TCut momwin(ctext);
  // plot results
  TCanvas* allcan = new TCanvas("mu2eall","mu2e results",1200,800);
  allcan->Clear();
  int ncan=ihi-ilow;
  if(ncan >=3)
    allcan->Divide(2,2);
  else if(ncan>1)
    allcan->Divide(2,1);
  else
    allcan->Divide(1,1);
// setup a function to fit the DIO spectrum
  TF1* diofit = new TF1("diofit",DIOCZ_R,_mmin,_mmax,5);
  diofit->SetNumberFitPoints(1000);
  double rawscale_f = (_mmax-_mmin)*ndecay/_nbins;
  diofit->FixParameter(0,rawscale_f);
  diofit->SetParameter(1,-1.0);
  diofit->SetParameter(2,0.12);
  diofit->SetParameter(3,0.0);
  diofit->SetParameter(4,0.0);
  diofit->SetRange(102.5,105.5);
  diofit->SetLineColor(kBlue);

  for(unsigned icut=ilow;icut<ihi+1;icut++){
    double conmax = 1.50*_conspec[icut]->GetBinContent(_conspec[icut]->GetMaximumBin());
    allcan->cd(icut+1);
    if(logy){
      gPad->SetLogy();
    } else {
      _conspec[icut]->SetMinimum(-0.01);
      _conspec[icut]->SetMaximum(conmax);
    }
    _conspec[icut]->Draw();
    TF1* diof = diofit;
    if(dio){
      _diospec[icut]->Fit("diofit","R0");
      _diospec[icut]->Fit("diofit","RM0");
      diof = (TF1*)_diospec[icut]->FindObject("diofit");
      _diospec[icut]->Draw("same");
      if(diof !=0)diof->Draw("same");
    }
    
    int istart = _conspec[icut]->FindFixBin(momlow+0.5*mevperbin);
    int istop = _conspec[icut]->FindFixBin(momhigh-0.5*mevperbin);
        cout << "Integration low edge " << _diospec[icut]->GetBinLowEdge(istart)
        << " high edge " << _diospec[icut]->GetBinLowEdge(istop)+mevperbin << " for cut " << icut << endl;
    double dint_err, cint_err;
    double dint = _diospec[icut]->IntegralAndError(istart,istop,dint_err);
    double cint = _conspec[icut]->IntegralAndError(istart,istop,cint_err);
//    double fint = _flat_f[icut]->Integral(momlow,momhigh)/mevperbin;
    _conint[icut] = cint;
    _conint_err[icut] = cint_err;

    TPaveText* inttext = new TPaveText(0.15,0.4,0.4,0.9,"NDC");
    char itext[50];
  
    snprintf(itext,50,"%3.1e Protons On Target",nproton);
    TText* l = inttext->AddText(itext);
    snprintf(itext,50,"%3.1e Stopped #mu^{-}",nstopped);
    l = inttext->AddText(itext);

    snprintf(itext,50,"R_{#mue} = %2.2g",rmue);
    l = inttext->AddText(itext);
    l->SetTextColor(kRed);
//    inttext->AddLine();
//    snprintf(itext,50,"%4.2f MeV/c < P < %4.2f MeV/c",momlow,momhigh);
//    inttext->AddText(itext);
    snprintf(itext,50,"#Sigma CE = %4.3f #pm %2.2f",cint,cint_err);
    l = inttext->AddText(itext);
    l->SetTextColor(kRed);

    double ses = rmue/cint;
    double ses_err = rmue*cint_err/cint;
    snprintf(itext,50,"CE SES= %3.2e #pm %2.2e",ses,ses_err);
    l = inttext->AddText(itext);
    l->SetTextColor(kRed);

    snprintf(itext,50,"#Sigma DIO = %3.3f #pm %2.3f",dint,dint_err);
    l = inttext->AddText(itext);
//    l->SetTextColor(kBlue);
//    snprintf(itext,50,"#int RPC+AP+Cosmic = %2.2f",fint);
//    l = inttext->AddText(itext);
    l->SetTextColor(kBlue);
    inttext->Draw();

    double intfactor = _nbins/(_mmax-_mmin);
    double diofitint = diof->Integral(momlow,momhigh)*intfactor;
    double diofitint_err(0.0);
    if(dio)diofitint_err = diof->IntegralError(momlow,momhigh)*intfactor;
    snprintf(itext,50,"#int DIO fit = %3.2f #pm %2.2f",diofitint,diofitint_err);
    l = inttext->AddText(itext);
    l->SetTextColor(kBlue);
    inttext->Draw();

    TPaveText* cuttext = new TPaveText(0.1,0.2,0.4,0.4,"NDC");  
    char line[80];
    snprintf(line,80,"%s",pitch.GetTitle());
    cuttext->AddText(line);
    snprintf(line,80,"%s",livegate.GetTitle());
    cuttext->AddText(line);
    snprintf(line,80,"%s",quality[icut].GetTitle());
    cuttext->AddText(line);
    snprintf(line,80,"%s",cosmic.GetTitle());
    cuttext->AddText(line);
    //cuttext->Draw();
    cout <<  final[icut].GetTitle() << endl;

    TLine* momlowl = new TLine(momlow,0.0,momlow,0.8*conmax);
    momlowl->SetLineColor(kBlack);
    momlowl->SetLineStyle(2);
    momlowl->SetLineWidth(2);
    momlowl->Draw();

    TLine* momhighl = new TLine(momhigh,0.0,momhigh,0.8*conmax);
    momhighl->SetLineColor(kBlack);
    momhighl->SetLineStyle(2);
    momhighl->SetLineWidth(2);
    momhighl->Draw();

    TArrow* sigline = new TArrow(momlow,0.8*conmax,momhigh,0.8*conmax,0.01,"<>");
    sigline->SetLineWidth(2);
    sigline->Draw();

    TText* sigwin = new TText(0.5*(momlow+momhigh),0.90*conmax,"Signal Window");
    sigwin->SetTextAlign(21);
    sigwin->Draw();
    snprintf(line,80,"%4.2f < p < %4.2f MeV/c",momlow,momhigh); 
    TText* sigwin2 = new TText(0.5*(momlow+momhigh),0.84*conmax,line);
    sigwin2->SetTextAlign(21);
    sigwin2->SetTextSize(0.025);
    sigwin2->Draw();
  }
  allcan->cd(0);
  std::string ssuf(suffix);
  allcan->SaveAs((std::string("mu2e")+ssuf).c_str());
}

void mu2e::drawdio(double momlow,double momhigh,const char* suffix) {
  char ctext[80];
  snprintf(ctext,80,"dem.mom>%f&&dem.mom<%f",momlow,momhigh);
  TCut momwin(ctext);
  std::string ssuf(suffix);
  TCanvas* dioc = new TCanvas("dioc","dio",1200,800);
  dioc->Divide(2,2);
  Double_t dmhi = trueconvmom;
  Double_t dmlow = trueconvmom - diogenrange;
  TH1F* diogen = new TH1F("diogen","True DIO momentum;MeV/c",_nbins,dmlow,dmhi);
  TH1F* evtwt = new TH1F("evtwt","True DIO momentum;MeV/c",_nbins,dmlow,dmhi);
  //  evtwt->Sumw2();
  if(dio)dio->Project("diogen","demmcgen.mom");
  if(dio)dio->Project("evtwt","demmcgen.mom","evtinfo.evtwt");
  evtwt->Scale(dioscale);
  evtwt->SetLineColor(kBlue);
  diogen->SetLineColor(kRed);
  evtwt->SetStats(0);
  diogen->SetStats(0);


  Int_t colors[4] = {kRed,kBlue,kGreen,kBlack};
  TH1F* diogenwin[4] = {0,0,0,0};
  TH1F* diodiffwin[4] = {0,0,0,0};
  const char* dopt[4] = {"","same","same","same"};
  const char* cutset[4] = {"All Reconstructed Tracks","Loose Track Selection","Default Track Selection","Tight Track Selection"};
  TLegend* dgenwinleg = new TLegend(.15,.6,.45,.9);
  for(unsigned icut=0;icut<4;icut++){
    char diogenname[50], diodiffname[50];
    snprintf(diogenname,50,"diogenwin%i",icut);
    diogenwin[icut] = new TH1F(diogenname,"Generated Momentum of DIO in Signal Box;Generated DIO Momentum (MeV/c)",100,dmlow,dmhi);
    diogenwin[icut]->SetStats(0);

    snprintf(diodiffname,50,"diodiffwin%i",icut);
    diodiffwin[icut] = new TH1F(diodiffname,"Reco - True Momentum of DIO in Signal Box;#Delta Momentum (MeV/c)",100,-1,2.5);
    diodiffwin[icut]->SetStats(0);

    if(dio)dio->Project(diogenname,"mcent.mom","evtinfo.evtwt"*(final[icut]+momwin));
    diogenwin[icut]->SetFillColor(colors[icut]);
    if(dio)dio->Project(diodiffname,"dem.mom-mcent.mom","evtinfo.evtwt"*(final[icut]+momwin));
    diodiffwin[icut]->SetFillColor(colors[icut]);
    dgenwinleg->AddEntry(diogenwin[icut],cutset[icut],"f");
  }

  dioc->cd(1);
  gPad->SetLogy();
  // dead-reconing on spectrum, accounting for bins
  double diofscale = ndecay*diogenrange*diogen->GetEntries()/(_nbins*ndio);
  cout << "dio function scale = " << diofscale << endl;
  _diocz_f->SetParameter(0,diofscale);
  evtwt->Draw();

  _diocz_f->Draw("same");
  cout << "dio function @103 = " << _diocz_f->Eval(103.0) << endl;
  diogen->Draw("same");
  TLegend* dioleg = new TLegend(.2,.4,.6,.6);
  dioleg->AddEntry(diogen,"Generated","l");
  dioleg->AddEntry(evtwt,"Weighted","l");
  dioleg->AddEntry(_diocz_f,"Czarnecki etal","l");
  dioleg->Draw();

  dioc->cd(2);
  for(unsigned icut=0;icut<4;icut++){
    _diospec[icut]->SetFillColor(colors[icut]);
    if(icut==0)
      _diospec[icut]->Draw("Hist");
    else
      _diospec[icut]->Draw("Histsame");
  }
  dgenwinleg->Draw();
  TLine* momlowl = new TLine(momlow,0.0,momlow,1.5*_diospec[0]->GetBinContent(_diospec[0]->GetMaximumBin()));
  momlowl->SetLineColor(kBlack);
  momlowl->SetLineStyle(2);
  momlowl->SetLineWidth(2);
  momlowl->Draw();

  TLine* momhighl = new TLine(momhigh,0.0,momhigh,1.5*_diospec[0]->GetBinContent(_diospec[0]->GetMaximumBin()));
  momhighl->SetLineColor(kBlack);
  momhighl->SetLineStyle(2);
  momhighl->SetLineWidth(2);
  momhighl->Draw();

  dioc->cd(3);
  gPad->SetLogy();
  for(unsigned icut=0;icut<4;icut++){
    diogenwin[icut]->Draw(dopt[icut]);
  }
  dgenwinleg->Draw();

  dioc->cd(4);
  gPad->SetLogy();
  for(unsigned icut=0;icut<4;icut++){
    diodiffwin[icut]->Draw(dopt[icut]);
  }

  dioc->SaveAs((std::string("diocan")+ssuf).c_str());

  TCanvas* diores = new TCanvas("diores","DIO result",800,600);
  gPad->SetLogy();
  diodiffwin[mu2ecut]->Draw();
  double split = 0.720; // define tail as above 4 sigma
  TLine* td = new TLine(split,0.0,split,diodiffwin[mu2ecut]->GetMaximum());
  td->SetLineColor(kBlack);
  td->SetLineStyle(2);
  td->SetLineWidth(2);
  td->Draw();

  int istart = diodiffwin[mu2ecut]->FindFixBin(split);
  double core = diodiffwin[mu2ecut]->Integral(0,istart);
  double tail = diodiffwin[mu2ecut]->Integral(istart+1,100);
  double total = core+tail;
  core /= total;
  tail /= total;
  core *= 100;
  tail *= 100;
  cout <<"core = " << core << " tail = " << tail << endl;
  char ccore[30], ctail[30];
  snprintf(ccore,30,"Core =%3.1f%%",core);
  snprintf(ctail,30,"Tail =%3.1f%%",tail);
  TText* tcore = new TText(0.2,0.4,ccore);
  TText* ttail = new TText(0.6,0.4,ctail);
  tcore->SetNDC();
  ttail->SetNDC();
  tcore->Draw();
  ttail->Draw();
  diores->SaveAs((std::string("diores")+ssuf).c_str());
}

void mu2e::fitReco(unsigned icut) {
// fit the acceptance
  gStyle->SetOptFit(1);
  gStyle->SetOptStat(0);
  TH1F* amomd = new TH1F("amomd","Reco Acceptance vs Momentum;MeV/c",_nbins,_mmin,_mmax);
  TH1F* amom = new TH1F("amom","Reco Acceptance vs Momentum;MeV/c",_nbins,_mmin,_mmax);
  amom->Sumw2();
  amomd->Sumw2();
//  char fmomcut[80];
//  snprintf(fmomcut,80,"dem.mom>%5.3f",demmcgen.mom);
  if(dio)dio->Project("amom","demmcgen.mom",final[icut]+TCut("dem.mom>demmcgen.mom-4.0"));
  if(dio)dio->Project("amomd","demmcgen.mom");
  amom->Divide(amomd);
  amom->SetMinimum(0.08);
  amom->SetMaximum(0.18);
  _racc = new TF1("racc",recoacc,_mmin,_mmax,2);
  _racc->SetParName(0,"acc100");
  _racc->SetParName(1,"slope");

  _cball = new TF1("cball",crystalball,-10.0,5,7);
  _cball->SetParName(0,"Norm");
  _cball->SetParName(1,"x0");
  _cball->SetParName(2,"sigma");
  _cball->SetParName(3,"n");
  _cball->SetParName(4,"alpha");
  _cball->SetParName(5,"tailfrac");
  _cball->SetParName(6,"taillambda");

  _dscb = new TF1("dscb",fnc_dscb,-2.0,2.5,7);
  _dscb->SetParName(0,"Norm");
  _dscb->SetParName(1,"x0");
  _dscb->SetParName(2,"sigma");
  _dscb->SetParName(3,"ANeg");
  _dscb->SetParName(4,"PNeg");
  _dscb->SetParName(5,"APos");
  _dscb->SetParName(6,"PPos");

// set parameters according to cutset 'C'
  _momres = new TH1F("momres","Track Momentum Resolution;p_{eco}-p_{True} (MeV/c)",_nbins,-5,2.0);
//  _momres->Sumw2();
  dio->Project("momres","dem.mom-demmcgen.mom",final[icut]);
//  _momres->Scale(conscale);
  TCanvas* fcan = new TCanvas("fcan","Fits",1000,800);
  fcan->Clear();
  fcan->Divide(2,1);
  fcan->cd(1);
  amom->Fit("racc");
  fcan->cd(2);
  gPad->SetLogy();
  double integral = _momres->GetEntries()*_momres->GetBinWidth(1);
//  _cball->SetParameters(integral,_momres->GetMean()+0.07,0.3*_momres->GetRMS(),3.0,1.0,0.001,0.2);
//  _cball->SetParLimits(5,0.00001,0.01);
//  _cball->SetParLimits(6,0.1,_momres->GetRMS());
//  _cball->FixParameter(5,0.0);
//  _cball->FixParameter(6,0.1);

  _dscb->SetParameters(3*integral,_momres->GetMean()+0.07,0.3*_momres->GetRMS(),0.9,3.5,1.5,6.0);

  _momres->Fit(_dscb,"0");
  _momres->Fit(_dscb,"L0");
  _momres->Fit(_dscb,"L");
}

void mu2e::smearDIO(unsigned ntrials,unsigned nres) {
// raw dio spectrum
  TF1* diocz_f = new TF1("diocz_f",DIOCZ,95.0,trueconvmom,1);
  diocz_f->SetLineColor(kCyan);
  diocz_f->SetParameter(0,1.0);
  double buffer(4.0);
// resolution function
// override GRandom
  gRandom = new TRandom3(324913);
  double mevpb = (_mmax-_mmin)/_nbins;
  char dtitle[100];
  snprintf(dtitle,100,"DIO e^{-} Momentum;P_{TRUE} (MeV/c);N Events/%2.2f MeV/c",mevpb);
  TH1F* rawdio = new TH1F("rawdio",dtitle,_nbins,_mmin,_mmax);
  rawdio->Sumw2();
  snprintf(dtitle,100,"All Reco Effects Applied DIO e^{-} Momentum;P_{RECO} (MeV/c);N Events/%2.2f MeV/c",mevpb);
  _recodio = new TH1D("recodio",dtitle,_nbins,_mmin,_mmax);
  _recodio->Sumw2();
  TH1F* rawcon = new TH1F("rawcon","Raw Conversion e^{-} Momentum;P_{true} (MeV/c);N Events",_nbins,_mmin,_mmax);
  TH1F* acccon = new TH1F("acccon","Acceptance Applied Conversion e^{-} Momentum;P_{reco} (MeV/c);N Events",_nbins,_mmin,_mmax);
  TH1F* econ = new TH1F("econ","Acc, #Delta E Applied Conversion e^{-} Momentum;P_{reco} (MeV/c);N Events",_nbins,_mmin,_mmax);
  rawcon->SetLineColor(kRed);
  acccon->SetLineColor(kRed);
  econ->SetLineColor(kRed);
  rawcon->Fill(trueconvmom,rmue*ncap);
  acccon->Fill(trueconvmom,rmue*ncap*_racc->Eval(trueconvmom));
  econ->Fill(trueconvmom+_dscb->GetParameter(1),rmue*ncap*_racc->Eval(trueconvmom));
  double genrange = trueconvmom-_mmin+buffer;
  for(unsigned itrial=0;itrial<ntrials;++itrial){
// make the range bigger than expected, to account for smearing
    double diomom = gRandom->Uniform(_mmin-buffer,trueconvmom);
    double evtwt = diocz_f->Eval(diomom);
    double acc = _racc->Eval(diomom);
    double recowt = evtwt*acc;
    rawdio->Fill(diomom,evtwt);
    for(unsigned ires=0;ires<nres;++ires){
      double mres = _dscb->GetRandom();
      double recomom = diomom+mres;
      _recodio->Fill(recomom,recowt); 
    }
  }
// scale these to the # of expected events for mu2e
  double rawscale = genrange*ndecay/ntrials;
  cout << "raw scale = " << rawscale << endl;
  rawdio->Scale(rawscale);
  double recoscale = rawscale/nres;
  cout << "reco scale = " << recoscale << endl;
  _recodio->Scale(recoscale);
  TCanvas* diocan = new TCanvas("diocan","dio spectra",800,600);
  diocan->Clear();
  diocan->Divide(2,1);
  diocan->cd(1);
  gPad->SetLogy();
  rawdio->SetMinimum(0.01);
  TF1* raw_f = new TF1("raw_f",DIOCZ,_mmin,_mmax,1);
  double rawscale_f = (_mmax-_mmin)*ndecay/_nbins;
  raw_f->SetParameter(0,rawscale_f);
  raw_f->SetLineColor(kBlue);
  double a5(8.6434e-17);
  double a6(1.16874e-17);
  double a7(-1.87828e-19);
  double a8(9.16327e-20);
  rawdio->SetMinimum(1e-3);
  rawdio->Draw();
  raw_f->Draw("same");
  diocan->cd(2);
  gPad->SetLogy();
  _recodio->SetMinimum(1e-3);
  _reco_f = new TF1("_reco_f",RECODIO,_mmin,trueconvmom,5);
  _reco_f->SetLineColor(kBlue);
  double recoscale_f  = genrange*ndecay*_racc->Eval(103.0)/_nbins;
  _reco_f->SetParameter(0,0.0);
  _reco_f->SetParameter(1,recoscale_f);
  _reco_f->SetParameter(2,a6/a5);
  _reco_f->SetParameter(3,a7/a5);
  _reco_f->SetParameter(4,a8/a5);
  _recodio->Fit("_reco_f","N");
  _recodio->Fit("_reco_f","MN");
  _recodio->Fit("_reco_f","M");

  TF1* dspeca_f = new TF1("dspeca_f",DIOCZ_R,_mmin,_mmax,6);
  dspeca_f->FixParameter(0,rawscale_f);
  dspeca_f->SetParameter(1,-1.0);
  dspeca_f->SetParameter(2,_racc->GetParameter(0));
  dspeca_f->SetParameter(3,_racc->GetParameter(1));
  dspeca_f->SetParameter(4,0.0);
  dspeca_f->SetLineColor(kBlue);

  TF1* dspece_f = new TF1("dspece_f",DIOCZ_R,_mmin,_mmax,6);
  dspece_f->FixParameter(0,rawscale_f);
  dspece_f->SetParameter(1,_dscb->GetParameter(1)); // parameter 1 is the offset of the mean in the resolution function = average energy loss
  dspece_f->SetParameter(2,_racc->GetParameter(0));
  dspece_f->SetParameter(3,_racc->GetParameter(1));
  dspece_f->SetParameter(4,0.0);
  dspece_f->SetLineColor(kBlue);

  char ttitle[100];
  char atitle[100];
  char etitle[100];
  char rtitle[100];
  char title[100];
  snprintf(title,100,"e^{-} Momentum;P (MeV/c);N Events/%2.2f MeV/c",mevpb);
  snprintf(ttitle,100,"Theory Predictions;e^{-} Momentum (MeV/c);N Events/%2.2f MeV/c",mevpb);
  snprintf(atitle,100,"After Reco Acceptance;e^{-} Momentum (MeV/c);N Events/%2.2f MeV/c",mevpb);
  snprintf(etitle,100,"After Reco Acceptance+#DeltaE;e^{-} Momentum (MeV/c);N Events/%2.2f MeV/c",mevpb);
  snprintf(rtitle,100,"After Reco Acceptance+#DeltaE+Resolution;e^{-} Momentum (MeV/c);N Events/%2.2f MeV/c",mevpb);
  TH1F* dspec = new TH1F("dspec",title,_nbins,_mmin,_mmax);
  TH1F* dspect = new TH1F("dspect",ttitle,_nbins,_mmin,_mmax);
  TH1F* dspeca = new TH1F("dspeca",atitle,_nbins,_mmin,_mmax);
  TH1F* dspece = new TH1F("dspece",etitle,_nbins,_mmin,_mmax);
  TH1F* dspecs = new TH1F("dspecs",rtitle,_nbins,_mmin,_mmax);
  dspec->SetMinimum(1e-3);
  dspec->SetMaximum(raw_f->Eval(_mmin)*3);
  dspect->SetMinimum(1e-3);
  dspect->SetMaximum(raw_f->Eval(_mmin)*3);
  dspeca->SetMinimum(1e-3);
  dspeca->SetMaximum(raw_f->Eval(_mmin)*3);
  dspece->SetMinimum(1e-3);
  dspece->SetMaximum(raw_f->Eval(_mmin)*3);
  dspecs->SetMinimum(1e-3);
  dspecs->SetMaximum(raw_f->Eval(_mmin)*3);
  TCanvas* diocan2 = new TCanvas("diocan2","Momentum  Effects",1000,1000);
  diocan2->Divide(2,2);
  diocan2->cd(1);
  gPad->SetLogy();
  dspect->Draw();
  raw_f->Draw("same");
  TLegend* tleg = new TLegend(0.3,0.75,0.9,0.9);
  tleg->AddEntry(raw_f,"DIO Prediction","L");
  tleg->AddEntry(rawcon,"Conversion, R_{#mue}=10^{-16}","L");
  tleg->Draw();
  rawcon->Draw("same");
  diocan2->cd(2);
  gPad->SetLogy();
  dspeca->Draw();
  acccon->Draw("same");
  dspeca_f->Draw("same");
  TLegend* aleg = new TLegend(0.3,0.75,0.9,0.9);
  aleg->AddEntry(dspeca_f,"DIO Prediction","L");
  aleg->AddEntry(acccon,"Conversion, R_{#mue}=10^{-16}","L");
  aleg->Draw();
  diocan2->cd(3);
  gPad->SetLogy();
  dspece->Draw();
  econ->Draw("same");
  dspece_f->Draw("same");
  TLegend* eleg = new TLegend(0.3,0.75,0.9,0.9);
  eleg->AddEntry(dspece_f,"DIO Prediction","L");
  eleg->AddEntry(econ,"Conversion, R_{#mue}=10^{-16}","L");
  eleg->Draw();
  diocan2->cd(4);
  gPad->SetLogy();
  dspecs->Draw();
  _reco_f->Draw("same");
  TF1* reccon = new TF1(*_dscb);
  reccon->SetLineColor(kRed);
  reccon->SetName("reccon");
  reccon->SetParameter(1,_dscb->GetParameter(1)+trueconvmom);
// must scale functions
  double csum = _dscb->Integral(-5.0,1.0)/_momres->GetBinWidth(1);
  cout << "CBall sum = " << csum;
  double norm = conscale*_dscb->GetParameter(0)*_momres->GetBinWidth(1)/dspec->GetBinWidth(1);
  reccon->SetParameter(0,norm);
  reccon->SetRange(_mmin,_mmax);
  cout << "reccon name " << reccon->GetName() << " norm = " << reccon->GetParameter(0) << " shift = " << reccon->GetParameter(1) << endl;
//  reccon->Draw();
  reccon->Draw("same");
  TLegend* rleg = new TLegend(0.3,0.75,0.9,0.9);
  rleg->AddEntry(_reco_f,"DIO Prediction","L");
  rleg->AddEntry(reccon,"Conversion, R_{#mue}=10^{-16}","L");
  rleg->Draw();
 
  TCanvas* diocan3 = new TCanvas("diocan3","DIO Spectrum Effects",600,600);
  diocan3->Divide(1,1);
  diocan3->cd(1);
  gPad->SetLogy();
  dspec->Draw();
  raw_f->Draw("same");
  TLegend* leg = new TLegend(0.5,0.7,0.9,0.9);
  leg->AddEntry(raw_f,"Czarnecki etal","L");
  dspeca_f->Draw("same");
  leg->AddEntry(dspeca_f,"#times Acc_{Detector}","L");
  dspece_f->Draw("same");
  leg->AddEntry(dspece_f,"#times Acc_- #DeltaE_{Detector}","L");
  _reco_f->Draw("same");
  leg->AddEntry(_reco_f,"(#times Acc_- #DeltaE)#otimes #sigma_{Detector}","L");
  leg->Draw();

}

void mu2e::doExperiments(double momlow, double momhigh,double cprob,unsigned ispec, unsigned nexp, unsigned npave) {
  unsigned ncans(0);
  if(npave > 0)
    ncans = ceil(float(nexp)/npave);
  std::vector<TCanvas*> cans(ncans,0);
  double avgbkg(0.4);
  double plevel(3e-7); // 5 sigma

// override GRandom
  gRandom = new TRandom3(nexp*nexp*7+13);

  double conmean = _conspec[ispec]->Integral()*cprob/rmue;
  double diomean = _diospec[ispec]->Integral();
//  double flatmean = _flat_f[ispec]->Integral(_mmin,_mmax)/mevperbin;
  cout << "conv mean = " << conmean << " dio mean = " << diomean 
//  << " flat mean = " << flatmean i
  << endl;

  int istart = _diospec[ispec]->FindFixBin(momlow+0.5*mevperbin);
  int istop = _diospec[ispec]->FindFixBin(momhigh-0.5*mevperbin);
  TLine* momlowl = new TLine(momlow,0.0,momlow,_diospec[ispec]->GetBinContent(_diospec[ispec]->GetMaximumBin()));
  momlowl->SetLineColor(kBlack);
  momlowl->SetLineStyle(2);
  momlowl->SetLineWidth(2);
  TLine* momhighl = new TLine(momhigh,0.0,momhigh,_diospec[ispec]->GetBinContent(_diospec[ispec]->GetMaximumBin()));
  momhighl->SetLineColor(kBlack);
  momhighl->SetLineStyle(2);
  momhighl->SetLineWidth(2);

  unsigned ican=0;
  unsigned ipave=1;
  char dioname[50];
  char conname[50];
  //char flatname[50];
  TLegend* leg(0);
// fit each experiment DIO spectrum OUTSIDE THE SIGNAL BOX
  TF1* dioexp_f = new TF1("dioexp_f",RECOFITDIO,_mmin,momlow-1.0,7);
  dioexp_f->SetLineColor(kBlue);
  for(unsigned ip=0;ip<5;++ip)
    dioexp_f->FixParameter(ip,_reco_f->GetParameter(ip));
  dioexp_f->SetParameter(5,0.0);
  dioexp_f->SetParameter(6,1.0);

  TH1F* ndio_h = new TH1F("ndio_h","N Events in Signal Window;N Events;N Toys",15,-0.5,14.5);
  TH1F* nobk_h = new TH1F("nobk_h","N Events in Signal Window;N Events;N Toys",15,-0.5,14.5);
  TH1F* ncon_h = new TH1F("ncon_h","N Events in Signal Window;N Events;N Toys",15,-0.5,14.5);
  ndio_h->SetLineColor(kBlue);
  ncon_h->SetLineColor(kRed);
  nobk_h->SetLineColor(kGreen);
  TH1F* diofscale = new TH1F("diofscale","DIO Fit Scale;Scale (dimensionles);N Toys",100,0.85,1.15);
  TH1F* diofshift = new TH1F("diofshift","DIO Fit Shift;#Delta P (MeV/c);N Toys",100,-0.15,0.15);
  TH2F* diofparams = new TH2F("diofparams","DIO Fit Parameters;Scale (dimensionless);Shift (MeV/c)",100,0.85,1.15,100,-0.15,0.15);
  TH1F* diofint = new TH1F("diofint","DIO Fit Signal Window Integral;N Events;N Toys",100,0,0.6);
  TH1F* diofint_err = new TH1F("diofint_err","DIO Fit Signal Window Integral Error;N Events;N Toys",100,0,0.2);


  char etitle[100];
  snprintf(etitle,100,"Toy Mu2e Experiment;Momentum (MeV/c);Events/%3.3f MeV/c",mevperbin);
  int npass(0);
  TArrow* sigline = new TArrow(momlow,4.0,momhigh,4.0,0.01,"<>");
  sigline->SetLineWidth(2);
  TText* sigwin = new TText(0.5*(momlow+momhigh),10.0,"Signal Window");
  sigwin->SetTextAlign(21);
  char line[80];
  snprintf(line,80,"%4.1f < p < %4.1f MeV/c",momlow,momhigh); 
  TText* sigwin2 = new TText(0.5*(momlow+momhigh),6.0,line);
  sigwin2->SetTextAlign(21);
  sigwin2->SetTextSize(0.025);
  for(unsigned iexp=0;iexp<nexp;++iexp){
    snprintf(conname,50,"conexp%i",iexp);
    TH1F* conexp_h = new TH1F(conname,etitle,_nbins,_mmin,_mmax);
    conexp_h->SetStats(0);
    conexp_h->SetLineColor(kRed);
    conexp_h->SetFillColor(kRed);
//    conexp_h->SetMarkerStyle(2);
//    conexp_h->SetMarkerColor(kRed);
    
    snprintf(dioname,50,"dioexp%i",iexp);
    TH1F* dioexp_h = new TH1F(dioname,etitle,_nbins,_mmin,_mmax);
    dioexp_h->SetStats(0);
    dioexp_h->SetLineColor(kBlue);
    dioexp_h->SetMarkerColor(kBlue);
    dioexp_h->SetMinimum(1e-1);
    dioexp_h->SetMaximum(2.0*dioexp_f->Eval(_mmin));

//    snprintf(flatname,50,"flatexp%i",iexp);
//    TH1F* flatexp_h = new TH1F(flatname,etitle,_nbins,_mmin,_mmax);
//    flatexp_h->SetStats(0);
//    flatexp_h->SetLineColor(kGreen);
//    flatexp_h->SetMarkerStyle(7);
//    flatexp_h->SetMarkerColor(kGreen);
    
    unsigned nconexp = gRandom->Poisson(conmean);
    unsigned ndioexp = gRandom->Poisson(diomean);
//    unsigned nflatexp = gRandom->Poisson(flatmean);

    for(unsigned idio=0;idio<ndioexp;++idio){
      dioexp_h->Fill(_diospec[ispec]->GetRandom());
    }
    for(unsigned icon=0;icon<nconexp;++icon){
      conexp_h->Fill(_conspec[ispec]->GetRandom());
    }
//    for(unsigned iflat=0;iflat<nflatexp;++iflat){
//      flatexp_h->Fill(_flat_f[ispec]->GetRandom());
//    }
    
    dioexp_h->Fit("dioexp_f","qLR0","P");
    TF1* ff = dioexp_h->GetFunction("dioexp_f");
    if(ncans > 0){
      if(cans[ican]==0){
	char ctext[80];
	snprintf(ctext,80,"expcan_%i",ican);
	cans[ican] = new TCanvas(ctext,"Simulated events",900,900);
	cans[ican]->Divide(npave,npave);
	ipave=1;
      }
      cans[ican]->cd(ipave);
      gPad->SetLogy();
      dioexp_h->Draw("H");
      ff->Draw("same");
//      flatexp_h->Draw("sameP");
      conexp_h->Draw("same");
      momlowl->Draw();
      momhighl->Draw();
    } 

    double dint = dioexp_h->Integral(istart,istop);
    double cint = conexp_h->Integral(istart,istop);
//    double fint = flatexp_h->Integral(istart,istop);
//    double tint = dint+cint+fint;
    double tint = dint+cint;
    double bkgprob = TMath::PoissonI(tint,avgbkg);
    if(bkgprob < plevel)++npass;

    double shift = ff->GetParameter(5);
    double scale = ff->GetParameter(6);
    double dfint = ff->Integral(momlow,momhigh)/mevperbin;
    double dfint_err = ff->IntegralError(momlow,momhigh)/mevperbin;
    ndio_h->Fill(dint);
//    nobk_h->Fill(fint);
    ncon_h->Fill(cint);
    diofscale->Fill(scale);
    diofshift->Fill(shift);
    diofparams->Fill(scale,shift);
    diofint->Fill(dfint);
    diofint_err->Fill(dfint_err);

    if(ncans > 0){

      cout << "DIO fit shift = " << shift << " scale = " << scale << " inetgral " << dfint << " +- " << dfint_err << endl;

      TPaveText* inttext = new TPaveText(0.6,0.7,0.9,0.9,"NDC");
      char text[80];
      snprintf(text,80,"%5.2e stopped muons",nstopped);
      TString snstop(text);
      inttext->AddText(snstop);
      snprintf(text,80,"R_{#mu e} = %e",cprob);
      TString srmue(text);
      inttext->AddText(srmue);
      inttext->SetBorderSize(1);
      char itext[50];
      snprintf(itext,50,"%3.1f < P < %3.1f MeV/c",momlow,momhigh);
      inttext->AddText(itext);
      snprintf(itext,50,"N DIO = %2.0f",dint);
      inttext->AddText(itext);
//      snprintf(itext,50,"N RPC+AP+Cosmic = %2.0f",fint);
//      inttext->AddText(itext);
      snprintf(itext,50,"N Conversion = %2.0f",cint);
      inttext->AddText(itext);
//      snprintf(itext,50,"DIO Fit Shift = %3.3f #pm %3.3f",ff->GetParameter(5),ff->GetParError(5));
//      inttext->AddText(itext);
//      snprintf(itext,50,"DIO Fit Scale = %3.3f #pm %3.3f",ff->GetParameter(6),ff->GetParError(6));
//      inttext->AddText(itext);
//      snprintf(itext,50,"DIO Fit Integral = %3.3f #pm %3.3f",dfint,dfint_err);
//      inttext->AddText(itext);
      inttext->Draw();
      if(leg==0){
	leg = new TLegend(0.15,0.8,0.4,0.9);
	leg->AddEntry(dioexp_h,"DIO","L");
	leg->AddEntry(conexp_h,"Conversions","F");
//	leg->AddEntry(flatexp_h,"RPC+AP+Cosmics","P");
      }
      leg->Draw();
      sigline->Draw();
      sigwin->Draw();
      sigwin2->Draw();

      ipave++;
      if(ipave>npave*npave){
	char cfile[50];
	snprintf(cfile,50,"mu2e_exp_%i.png",ican);
	cans[ican]->SaveAs(cfile);
	++ican;
      }
    } else {
      delete conexp_h;
      delete dioexp_h;
//      delete flatexp_h;
    }
  }
  TCanvas* expcan1 = new TCanvas("expcan1","Toy Experiments",800,800);
  expcan1->Clear();
  expcan1->Divide(1,2);
//  expcan1->cd(1);
//  nobk_h->Draw();
//  ncon_h->Draw("same");
//  ndio_h->Draw("same");

//  TLegend* nleg = new TLegend(0.2,0.6,0.9,0.9);
//  char title[80];
//  snprintf(title,80,"Conversions, <N>=%3.3f",ncon_h->GetMean());
//  nleg->AddEntry(ncon_h,title,"L");
//  snprintf(title,80,"DIO, <N>=%3.3f",ndio_h->GetMean());
//  nleg->AddEntry(ndio_h,title,"L");
//  snprintf(title,80,"RPC+AP+Cosmics, <N>=%3.3f",nobk_h->GetMean());
//  nleg->AddEntry(nobk_h,title,"L");
//  nleg->Draw();
  expcan1->cd(1);
  diofint->SetStats(1);
  diofint->Draw();
//  diofint->Fit("gaus");
  expcan1->cd(2);
  diofint_err->SetStats(1);
  diofint_err->Draw();
//  diofint_err->Fit("gaus");

  TCanvas* expcan2 = new TCanvas("expcan2","Toy Experiments",800,800);
  expcan2->Clear();
  expcan2->Divide(2,2);
  expcan2->cd(1);

  expcan2->cd(2);
  diofscale->Fit("gaus");
  expcan2->cd(3);
  diofshift->Fit("gaus");
  expcan2->cd(4);
  diofparams->Draw();

  double passfrac = float(npass)/float(nexp);
  cout << "Pass prob cut " << passfrac << endl;
}


