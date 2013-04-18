#include "TTree.h"
#include "TH1F.h"
#include "TF1.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TCut.h"
#include "TPaveText.h"
#include "TLine.h"
#include "TRandom3.h"
#include <iostream>
#include <math.h>
#include <vector>

// the following approximation is from Czarnecki etal, 'Muon decay in orbit:spectrum of high-energy electrons',
// for E>85 MeV
Double_t DIOCZ(Double_t *x, Double_t *par) {
  double ee = x[0];
  double norm = par[0];
  double mal(25133);
  //    double mmu(105.654);
  double emu(105.194);
  //    double emue(104.973);
  //    double me(0.511);
  double a5(8.6434e-17);
  double a6(1.16874e-17);
  double a7(-1.87828e-19);
  double a8(9.16327e-20);
  double delta = emu - ee - ee*ee/(2*mal);
  return norm*(a5*pow(delta,5) + a6*pow(delta,6) + a7*pow(delta,7) + a8*pow(delta,8));
}

class mu2e {
  public:
    mu2e(TTree* d, TTree* c, double dgenrange, double nd, double nc,bool weightd=true) : dio(d), con(c),diogenrange(dgenrange),
    ndio(nd),ncon(nc),weightdio(weightd),nbins(151),nstopped(7.56e17),capfrac(0.609),conprob(1e-16),trueconvmom(104.973),
    tdlow(0.57735027),tdhigh(1.0),mmin(101),mmax(106),t0min(710),rpc(0.025), ap(.083333),cmu(0.041666),reco("fitstatus>0")
  {
    init();
  }
    void init();

    void drawmu2e(double momlow,double momhigh,const char* suffix=".png");
    void drawdio(double momlow,double momhigh,const char* suffix=".png");
    void doExperiments(double momlow, double momhigh,double cprob,unsigned ispec,unsigned nexp=18, unsigned npave=3);
    TTree *dio, *con;
    double diogenrange;
    double ndio, ncon;
    bool weightdio;
    unsigned nbins;
    double nstopped,capfrac,conprob;
    double decayfrac,ndecay,ncap,mevperbin,conscale;
    double trueconvmom,tdlow,tdhigh,mmin,mmax,t0min;
    double rpc, ap, cmu,flat; // rates per MeV
    double dioint,dioscale;
    unsigned mu2ecut;
    TCut reco, pitch, livegate, cosmic;
    TCut ncuts[4], t0cuts[4], momcuts[4], fitcuts[4], quality[4], final[4];
    TF1* diocz_f;
    TF1* flat_f[4];
    TH1F* diospec[4];
    TH1F* conspec[4];
};

void mu2e::init(){
  decayfrac = 1.0 - capfrac;
  ndecay = nstopped*decayfrac;
  ncap = nstopped*capfrac;
  mevperbin = (mmax-mmin)/nbins;
  conscale = ncap*conprob/ncon;
  cout << "Conversion scale factor =" << conscale << endl;
  // dio spectrum
  diocz_f = new TF1("diocz_f",DIOCZ,85.0,105,1);
  diocz_f->SetLineColor(kGreen);
  diocz_f->SetParameter(0,1.0);
  // integrate the DIO spectrum over the range specified.  This is relative to the free decay rate
  dioint = diocz_f->Integral(trueconvmom-diogenrange,trueconvmom);
  if(weightdio){ 
    dioscale =ndecay*diogenrange/ndio;
  } else {
    dioscale = dioint*ndecay/ndio;
  }
  cout << "DIO scale factor = " << dioscale << endl;
  flat = rpc+ap+cmu;
  cout << "Flat rate = " << flat << " counts/MeV/c" << endl;
  // basic cuts
  char ctext[80];
  snprintf(ctext,80,"td>%f&&td<%f",tdlow,tdhigh);
  pitch = TCut(ctext);
  snprintf(ctext,80,"t0>%f",t0min);
  livegate = TCut(ctext);
  cosmic = TCut("d0<105&&d0>-80 && d0+2/om>450 && d0+2/om<680");

  // cuts for different tightness of selection
  ncuts[0] = "nactive>=20";
  ncuts[1] = "nactive>=22";
  ncuts[2] = "nactive>=25";
  ncuts[3] = "nactive>=30";
  t0cuts[0] = "t0err<1.5";
  t0cuts[1] = "t0err<0.95";
  t0cuts[2] = "t0err<0.9";
  t0cuts[3] = "t0err<0.8";
  momcuts[0] = "fitmomerr<0.3";
  momcuts[1] = "fitmomerr<0.28";
  momcuts[2] = "fitmomerr<0.25";
  momcuts[3] = "fitmomerr<0.22";
  fitcuts[0] = "fitcon>1e-6";
  fitcuts[1] = "fitcon>1e-3";
  fitcuts[2] = "fitcon>2e-3";
  fitcuts[3] = "fitcon>1e-2";

  for(unsigned icut=0;icut<4;icut++){
    quality[icut] = ncuts[icut] && t0cuts[icut] && momcuts[icut] && fitcuts[icut];
    final[icut] = (reco+pitch+livegate+quality[icut]+cosmic);
  } 
  mu2ecut=2;
}

void mu2e::drawmu2e(double momlow, double momhigh,const char* suffix) {
  std::string ssuf(suffix);
  for(unsigned icut=0;icut<4;icut++){
    char dioname[50];
    snprintf(dioname,50,"diospec%i",icut);
    char conname[50];
    snprintf(conname,50,"conspec%i",icut);
    char flatname[50];
    snprintf(flatname,50,"flat_f%i",icut);

    diospec[icut] = new TH1F(dioname,"Reconstructed momentum;MeV/c;Events per experiment",nbins,mmin,mmax);
    diospec[icut]->SetStats(0);
    diospec[icut]->SetLineColor(kBlue);
    diospec[icut]->Sumw2();

    conspec[icut] = new TH1F(conname,"Reconstructed momentum;MeV/c;Events per experiment",nbins,mmin,mmax);
    conspec[icut]->SetStats(0);
    conspec[icut]->SetLineColor(kRed);
    conspec[icut]->Sumw2();

    dio->Project(dioname,"fitmom","diowt"*final[icut]);
    diospec[icut]->Scale(dioscale);

    con->Project(conname,"fitmom",final[icut]);
    conspec[icut]->Scale(conscale);

    
    flat_f[icut] = new TF1(flatname,"[0]",mmin,mmax);
    flat_f[icut]->SetLineColor(kGreen);
    double acc(1.0);
    if(icut>0)acc = conspec[icut]->Integral()/conspec[0]->Integral();
    flat_f[icut]->SetParameter(0,flat*mevperbin*acc);
  }

  TLegend* leg = new TLegend(0.7,0.8,0.9,0.9);
  leg->AddEntry(diospec[0],"DIO","L");
  leg->AddEntry(conspec[0],"Conversion","L");
  leg->AddEntry(flat_f[0],"RPC+AP+cosmic","L");

  TPaveText* info = new TPaveText(0.4,0.8,0.7,0.9,"NDC");
  char text[80];
  snprintf(text,80,"%g stopped muons",nstopped);
  TString snstop(text);
  info->AddText(snstop);
  snprintf(text,80,"%g Conversion Rate",conprob);
  TString sconprob(text);
  info->AddText(sconprob);
  info->SetBorderSize(0);

  // plot results
  TCanvas* mu2ecan = new TCanvas("mu2e","mu2e result",900,600);
  mu2ecan->Clear();
  mu2ecan->Divide(1,1);
  TCanvas* allcan = new TCanvas("mu2eall","mu2e results",1200,800);
  allcan->Clear();
  allcan->Divide(2,2);
  for(unsigned icut=0;icut<4;icut++){
    allcan->cd(icut+1);
    TH1* diocopy = diospec[icut]->DrawCopy();
    diocopy->SetMinimum(-1.0/nbins);
    diocopy->SetMaximum(60.0/nbins);
    conspec[icut]->Draw("same");
    flat_f[icut]->Draw("same");

    int istart = diospec[icut]->FindFixBin(momlow+0.5*mevperbin);
    int istop = diospec[icut]->FindFixBin(momhigh-0.5*mevperbin);
    //    cout << "Integration low edge " << diospec[icut]->GetBinLowEdge(istart) << " for cut at " << momlow << endl;
    //    cout << "Integration high edge " << diospec[icut]->GetBinLowEdge(istop)+mevperbin << " for cut at " << momhigh << endl;
    double dint_err, cint_err;
    double dint = diospec[icut]->IntegralAndError(istart,istop,dint_err);
    double cint = conspec[icut]->IntegralAndError(istart,istop,cint_err);
    double fint = flat_f[icut]->Integral(momlow,momhigh)/mevperbin;

    TPaveText* inttext = new TPaveText(0.5,0.6,0.9,0.8,"NDC");
    char itext[50];
    snprintf(itext,50,"%4.2f MeV/c < P < %4.2f MeV/c",momlow,momhigh);
    inttext->AddText(itext);
    snprintf(itext,50,"DIO integral = %4.2f #pm %4.2f",dint,dint_err);
    inttext->AddText(itext);
    snprintf(itext,50,"Conv. integral = %4.2f #pm %4.2f",cint,cint_err);
    inttext->AddText(itext);
    snprintf(itext,50,"RPC+AP+cosmic integral = %4.2f",fint);
    inttext->AddText(itext);
    inttext->Draw();

    TPaveText* cuttext = new TPaveText(0.1,0.5,0.4,0.8,"NDC");  
    char line[40];
    snprintf(line,80,"%4.3f<tan(#lambda)<%4.3f",tdlow,tdhigh);
    cuttext->AddText(line);
    snprintf(line,80,"t0>%5.1f nsec",t0min);
    cuttext->AddText(line);
    sprintf(line,"%s",ncuts[icut].GetTitle());
    cuttext->AddText(line);
    sprintf(line,"%s",t0cuts[icut].GetTitle());
    cuttext->AddText(line);
    sprintf(line,"%s",momcuts[icut].GetTitle());
    cuttext->AddText(line);
    sprintf(line,"%s",fitcuts[icut].GetTitle());
    cuttext->AddText(line);
    cuttext->Draw();

    TLine* momlowl = new TLine(momlow,0.0,momlow,1.5*conspec[icut]->GetBinContent(conspec[icut]->GetMaximumBin()));
    momlowl->SetLineColor(kBlack);
    momlowl->SetLineStyle(2);
    momlowl->SetLineWidth(2);
    momlowl->Draw();

    TLine* momhighl = new TLine(momhigh,0.0,momhigh,1.5*conspec[icut]->GetBinContent(conspec[icut]->GetMaximumBin()));
    momhighl->SetLineColor(kBlack);
    momhighl->SetLineStyle(2);
    momhighl->SetLineWidth(2);
    momhighl->Draw();
    leg->Draw();
    info->Draw();
    if(icut == mu2ecut){
      mu2ecan->cd(0);
      diocopy->Draw();
      conspec[icut]->Draw("same");
      flat_f[icut]->Draw("same");
      inttext->Draw();
      cuttext->Draw();
      momlowl->Draw();
      momhighl->Draw();
      leg->Draw();
      info->Draw();
    }
  }
  allcan->cd(0);
  allcan->SaveAs((std::string("mu2e_all")+ssuf).c_str());
  mu2ecan->SaveAs((std::string("mu2e")+ssuf).c_str());
}


void mu2e::drawdio(double momlow,double momhigh,const char* suffix) {
  std::string ssuf(suffix);
  TCanvas* dioc = new TCanvas("dioc","dio",1200,800);
  dioc->Divide(2,2);
  Double_t dmhi = trueconvmom;
  Double_t dmlow = trueconvmom - diogenrange;
  TH1F* diogen = new TH1F("diogen","True DIO momentum;MeV",nbins,dmlow,dmhi);
  TH1F* diowt = new TH1F("diowt","True DIO momentum;MeV",nbins,dmlow,dmhi);
  //  diowt->Sumw2();
  dio->Project("diogen","mcmom");
  dio->Project("diowt","mcmom","diowt");
  diowt->Scale(dioscale);
  diowt->SetLineColor(kBlue);
  diogen->SetLineColor(kRed);
  diowt->SetStats(0);
  diogen->SetStats(0);

  char ctext[80];
  snprintf(ctext,80,"fitmom>%f&&fitmom<%f",momlow,momhigh);
  TCut momwin(ctext);

  Int_t colors[4] = {kRed,kBlue,kGreen,kBlack};
  TH1F* diogenwin[4] = {0,0,0,0};
  TH1F* diodiffwin[4] = {0,0,0,0};
  const char* dopt[4] = {"","same","same","same"};
  const char* cutset[4] = {"Cutset A","Cutset B","Cutset C","Cutset D"};
  TLegend* dgenwinleg = new TLegend(.5,.6,.7,.9);
  for(unsigned icut=0;icut<4;icut++){
    char diogenname[50], diodiffname[50];
    snprintf(diogenname,50,"diogenwin%i",icut);
    diogenwin[icut] = new TH1F(diogenname,"True momentum of DIO in signal box;MeV",100,dmlow,dmhi);
    diogenwin[icut]->SetStats(0);

    snprintf(diodiffname,50,"diodiffwin%i",icut);
    diodiffwin[icut] = new TH1F(diodiffname,"Reco - true momentum of DIO in signal box;MeV",100,-1,2);
    diodiffwin[icut]->SetStats(0);

    dio->Project(diogenname,"mcentmom","diowt"*(final[icut]+momwin));
    diogenwin[icut]->SetFillColor(colors[icut]);
    dio->Project(diodiffname,"fitmom-mcentmom","diowt"*(final[icut]+momwin));
    diodiffwin[icut]->SetFillColor(colors[icut]);
    dgenwinleg->AddEntry(diogenwin[icut],cutset[icut],"f");

  }

  dioc->cd(1);
  gPad->SetLogy();
  // dead-reconing on spectrum, accounting for bins
  double diofscale = ndecay*(dmhi-dmlow)/nbins;
  diocz_f->SetParameter(0,diofscale);
  diowt->Draw();
  diocz_f->Draw("same");
  diogen->Draw("same");
  TLegend* dioleg = new TLegend(.2,.4,.6,.6);
  dioleg->AddEntry(diogen,"Generated","l");
  dioleg->AddEntry(diowt,"Weighted","l");
  dioleg->AddEntry(diocz_f,"Czarnecki etal","l");
  dioleg->Draw();

  dioc->cd(2);
  for(unsigned icut=0;icut<4;icut++){
    diospec[icut]->SetFillColor(colors[icut]);
    if(icut==0)
      diospec[icut]->Draw("Hist");
    else
      diospec[icut]->Draw("Histsame");
  }
  dgenwinleg->Draw();
  TLine* momlowl = new TLine(momlow,0.0,momlow,1.5*diospec[0]->GetBinContent(diospec[0]->GetMaximumBin()));
  momlowl->SetLineColor(kBlack);
  momlowl->SetLineStyle(2);
  momlowl->SetLineWidth(2);
  momlowl->Draw();

  TLine* momhighl = new TLine(momhigh,0.0,momhigh,1.5*diospec[0]->GetBinContent(diospec[0]->GetMaximumBin()));
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
  double split = 0.35;
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
  diores->SaveAs("diores.png");
  char ccore[30], ctail[30];
  snprintf(ccore,30,"%4.1f%%",core);
  snprintf(ctail,30,"%4.1f%%",tail);
  TText* tcore = new TText(0.2,0.4,ccore);
  TText* ttail = new TText(0.6,0.4,ctail);
  tcore->SetNDC();
  ttail->SetNDC();
  tcore->Draw();
  ttail->Draw();
}

void mu2e::doExperiments(double momlow, double momhigh,double cprob,unsigned ispec, unsigned nexp, unsigned npave) {
  unsigned ncans = ceil(float(nexp)/npave);
  std::vector<TCanvas*> cans(ncans,0);
  std::vector<TH1F*> dioexp(nexp,0);
  std::vector<TH1F*> conexp(nexp,0);
  std::vector<TH1F*> flatexp(nexp,0);

  double conmean = conspec[ispec]->Integral()*cprob/conprob;
  double diomean = diospec[ispec]->Integral();
  double flatmean = flat_f[ispec]->Integral(mmin,mmax)/mevperbin;
  cout << "conv mean = " << conmean << " dio mean = " << diomean << " flat mean = " << flatmean << endl;

  TRandom3 rand(ispec*nexp*npave+11);
  // 1st value seems corrupt?
  for(unsigned idum=0;idum<nexp*npave;++idum){
    gRandom->Rndm();
  }

  int istart = diospec[ispec]->FindFixBin(momlow+0.5*mevperbin);
  int istop = diospec[ispec]->FindFixBin(momhigh-0.5*mevperbin);
  TLine* momlowl = new TLine(momlow,0.0,momlow,diospec[ispec]->GetBinContent(diospec[ispec]->GetMaximumBin()));
  momlowl->SetLineColor(kBlack);
  momlowl->SetLineStyle(2);
  momlowl->SetLineWidth(2);
  TLine* momhighl = new TLine(momhigh,0.0,momhigh,diospec[ispec]->GetBinContent(diospec[ispec]->GetMaximumBin()));
  momhighl->SetLineColor(kBlack);
  momhighl->SetLineStyle(2);
  momhighl->SetLineWidth(2);

  unsigned ican=0;
  unsigned ipave=1;
  char dioname[50];
  char conname[50];
  char flatname[50];
  TLegend* leg(0);
  for(unsigned iexp=0;iexp<nexp;++iexp){
    snprintf(conname,50,"conexp%i",iexp);
    conexp[iexp] = new TH1F(conname,"Sample Mu2e Experiment;Momentum (MeV/c);Events",nbins,mmin,mmax);
    conexp[iexp]->SetStats(0);
    conexp[iexp]->SetLineColor(kRed);
    conexp[iexp]->SetFillColor(kRed);
    conexp[iexp]->SetMarkerStyle(2);
    conexp[iexp]->SetMarkerColor(kRed);
    
    snprintf(dioname,50,"dioexp%i",iexp);
    dioexp[iexp] = new TH1F(dioname,"Sample Mu2e Experiment;Momentum (MeV/c);Events",nbins,mmin,mmax);
    dioexp[iexp]->SetStats(0);
    dioexp[iexp]->SetLineColor(kBlue);
    dioexp[iexp]->SetMarkerStyle(5);
    dioexp[iexp]->SetMarkerColor(kBlue);

    snprintf(flatname,50,"flatexp%i",iexp);
    flatexp[iexp] = new TH1F(flatname,"Sample Mu2e Experiment;Momentum (MeV/c);Events",nbins,mmin,mmax);
    flatexp[iexp]->SetStats(0);
    flatexp[iexp]->SetLineColor(kGreen);
    flatexp[iexp]->SetMarkerStyle(7);
    flatexp[iexp]->SetMarkerColor(kGreen);
    
    unsigned nconexp = rand.Poisson(conmean);
    unsigned ndioexp = rand.Poisson(diomean);
    unsigned nflatexp = rand.Poisson(flatmean);

    for(unsigned idio=0;idio<ndioexp;++idio){
      dioexp[iexp]->Fill(diospec[ispec]->GetRandom());
    }
    for(unsigned icon=0;icon<nconexp;++icon){
      conexp[iexp]->Fill(conspec[ispec]->GetRandom());
    }
    for(unsigned iflat=0;iflat<nflatexp;++iflat){
      flatexp[iexp]->Fill(flat_f[ispec]->GetRandom());
    }

    if(cans[ican]==0){
      char ctext[80];
      snprintf(ctext,80,"expcan%i",ican);
      cans[ican] = new TCanvas(ctext,"Simulated events",900,900);
      cans[ican]->Divide(npave,npave);
      ipave=1;
    }
    cans[ican]->cd(ipave);
    dioexp[iexp]->Draw("HP");
    flatexp[iexp]->Draw("sameHP");
    conexp[iexp]->Draw("sameHP");
    momlowl->Draw();
    momhighl->Draw();
    
    double dint = dioexp[iexp]->Integral(istart,istop);
    double cint = conexp[iexp]->Integral(istart,istop);
    double fint = flatexp[iexp]->Integral(istart,istop);

    TPaveText* inttext = new TPaveText(0.5,0.5,0.9,0.9,"NDC");
    char text[80];
    snprintf(text,80,"%g stopped muons",nstopped);
    TString snstop(text);
    inttext->AddText(snstop);
    snprintf(text,80,"R_{#mu e} = %g",cprob);
    TString sconprob(text);
    inttext->AddText(sconprob);
    inttext->SetBorderSize(1);
    char itext[50];
    snprintf(itext,50,"%3.1f < P < %3.1f MeV/c",momlow,momhigh);
    inttext->AddText(itext);
    snprintf(itext,50,"DIO integral = %2.0f",dint);
    inttext->AddText(itext);
    snprintf(itext,50,"RPC+AP+cosmic integral = %2.0f",fint);
    inttext->AddText(itext);
    snprintf(itext,50,"Conv. integral = %2.0f",cint);
    inttext->AddText(itext);
    inttext->Draw();

    if(leg==0){
      leg = new TLegend(0.15,0.7,0.45,0.9);
      leg->AddEntry(dioexp[0],"DIO","LP");
      leg->AddEntry(conexp[0],"Conversion","LP");
      leg->AddEntry(flatexp[0],"RPC+AP+cosmic","LP");
    }
    leg->Draw();

    ipave++;
    if(ipave>npave*npave)++ican;
 }
}
