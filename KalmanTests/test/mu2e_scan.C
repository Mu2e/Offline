#include "TTree.h"
#include "TH1F.h"
#include "TF1.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TCut.h"
#include "TPaveText.h"
#include "TLine.h"
#include "TGraphErrors.h"
#include "TVectorF.h"
#include <iostream>
#include "math.h"
#include "KalmanTests/test/DIOCZ.h"
using namespace std;

void mu2e_scan(TTree* dio, TTree* con, double diogenrange, double ndio, double ncon,bool weightdio=true,const char* suffix=".png") {
  // diogenrange is the momentum range over which the DIO events were generated
  double nstopped(5.76e17);
  double capfrac(0.609); 
  double decayfrac = 1.0 - capfrac;
  double ndecay = nstopped*decayfrac;
  double ncap = nstopped*capfrac;
  double conprob(1e-16);
  double trueconvmom(104.973);
  double momlow(103.38);
  double momhigh(104.8);

  unsigned nbins(1000);
  double mmin(101);
  double mmax(106);

  double conscale = ncap*conprob/ncon;
  cout << "Conversion scale factor =" << conscale << endl;
  // dio spectrum
  TF1* diocz_f = new TF1("diocz_f",DIOCZ,85.0,105,1);
  diocz_f->SetLineColor(kGreen);
  diocz_f->SetParameter(0,1.0);
  // integrate the DIO spectrum over the range specified.  This is relative to the free decay rate
  double dioint = diocz_f->Integral(trueconvmom-diogenrange,trueconvmom);
  double dioscale(1.0);
  if(weightdio){ 
    dioscale =ndecay*diogenrange/ndio;
  } else {
    dioscale = dioint*ndecay/ndio;
  }
  cout << "DIO scale factor = " << dioscale << endl;

  TH1F* diospec[4];
  TH1F* conspec[4];
  // basic cuts
  TCut reco("fitstatus>0");
  double tdlow(0.57735027);
  double tdhigh(1.0);
  double t0min(710);
  char ctext[80];
  snprintf(ctext,80,"td>%f&&td<%f",tdlow,tdhigh);
  TCut pitch(ctext);
  snprintf(ctext,80,"t0>%f",t0min);
  TCut livegate(ctext);
  TCut cosmic = TCut("d0<105&&d0>-80 && d0+2/om>450 && d0+2/om<680");

  // cuts for different tightness of selection
  TCut ncuts[4], t0cuts[4], momcuts[4], fitcuts[4];

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

    unsigned mu2ecut=2;
  unsigned icut=mu2ecut;
  char dioname[50];
  snprintf(dioname,50,"diospec%i",icut);
  char conname[50];
  snprintf(conname,50,"conspec%i",icut);

  diospec[icut] = new TH1F(dioname,"Reconstructed momentum;MeV/c;Events per experiment",nbins,mmin,mmax);
  diospec[icut]->SetStats(0);
  diospec[icut]->SetLineColor(kBlue);
  diospec[icut]->Sumw2();

  conspec[icut] = new TH1F(conname,"Reconstructed momentum;MeV/c;Events per experiment",nbins,mmin,mmax);
  conspec[icut]->SetStats(0);
  conspec[icut]->SetLineColor(kRed);
  conspec[icut]->Sumw2();

  TCut quality = ncuts[icut] && t0cuts[icut] && momcuts[icut] && fitcuts[icut];
  TCut final = (reco+pitch+livegate+quality+cosmic);

  dio->Project(dioname,"fitmom","evtwt"*final);
  diospec[icut]->Scale(dioscale);

  con->Project(conname,"fitmom",final);
  conspec[icut]->Scale(conscale);


  TLegend* leg = new TLegend(0.6,0.7,0.9,0.9);
  leg->AddEntry(diospec[icut],"DIO","L");
  leg->AddEntry(conspec[icut],"Conversion","L");

  TPaveText* info = new TPaveText(0.3,0.7,0.6,0.9,"NDC");
  char text[80];
  snprintf(text,80,"%g#times10^{17} stopped muons",nstopped*1e-17);
  TString snstop(text);
  info->AddText(snstop);
  snprintf(text,80,"R_{#mue}=%g#times10^{-16}",conprob*1e16);
  TString sconprob(text);
  info->AddText(sconprob);
  info->SetBorderSize(0);

  double mevperbin = (mmax-mmin)/nbins;
  int istart = diospec[icut]->FindFixBin(momlow+0.5*mevperbin);
  int istop = diospec[icut]->FindFixBin(momhigh-0.5*mevperbin);
  double dint_err, cint_err;
  double dint = diospec[icut]->IntegralAndError(istart,istop,dint_err);
  double cint = conspec[icut]->IntegralAndError(istart,istop,cint_err);

  TPaveText* inttext = new TPaveText(0.5,0.65,0.9,0.8,"NDC");
  char itext[50];
  snprintf(itext,50,"%4.2f MeV/c < P < %4.2f MeV/c",momlow,momhigh);
  inttext->AddText(itext);
  snprintf(itext,50,"DIO integral = %5.3f #pm %4.3f",dint,dint_err);
  inttext->AddText(itext);
  snprintf(itext,50,"Conv. integral = %5.3f #pm %4.3f",cint,cint_err);
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

  // plot results
  TCanvas* mu2ecan = new TCanvas("mu2e_result","mu2e result",900,600);
  mu2ecan->Clear();
  mu2ecan->Divide(1,1);
  mu2ecan->cd(0);
  diospec[icut]->SetMinimum(-5.0/nbins);
  diospec[icut]->SetMaximum(60.0/nbins);
  diospec[icut]->Draw();
  conspec[icut]->Draw("same");
  inttext->Draw();
  cuttext->Draw();
  momlowl->Draw();
  momhighl->Draw();
  leg->Draw();
  info->Draw();
  string ssuf(suffix);
  mu2ecan->SaveAs((string("mu2e_finebins")+ssuf).c_str());


  int nscan(25);
  vector<Double_t> diox,dioy,dioxerr,dioyerr;
  vector<Double_t> conx,cony,conxerr,conyerr;
  double startx = diospec[icut]->GetBinCenter(istart);
  for(int ibin=-nscan;ibin<=nscan;++ibin){
    double dints_err, cints_err;
    double dints = diospec[icut]->IntegralAndError(istart+ibin,istop,dints_err);
    double cints = conspec[icut]->IntegralAndError(istart+ibin,istop,cints_err);
    dioy.push_back(dints-dint);
    diox.push_back(1000.0*(diospec[icut]->GetBinCenter(istart+ibin)-startx));
    dioyerr.push_back(dints_err);
    dioxerr.push_back(0.0);
    cony.push_back(cints-cint);
    conx.push_back(1000.0*(conspec[icut]->GetBinCenter(istart+ibin)-startx));
    conyerr.push_back(cints_err);
    conxerr.push_back(0.0);
  }
  TGraphErrors* dioscan = new TGraphErrors(diox.size(),&diox[0],&dioy[0],&dioxerr[0],&dioyerr[0]);
  TGraphErrors* conscan = new TGraphErrors(conx.size(),&conx[0],&cony[0],&conxerr[0],&conyerr[0]);
  dioscan->SetLineColor(kBlue);
  conscan->SetLineColor(kRed);
  dioscan->SetMarkerColor(kBlue);
  conscan->SetMarkerColor(kRed);
  conscan->SetMaximum(0.5);
  conscan->SetMinimum(-0.5);
  dioscan->SetTitle("Rate Sensitivity to lower Momentum Cut;#Delta P (KeV/c);Integral Change");
  conscan->SetTitle("Rate Sensitivity to lower Momentum Cut;#Delta P (KeV/c);Integral Change");
  TCanvas* scancan = new TCanvas("scan","Scan of lower threshold",800,800);
  scancan->Divide(1,1);
  scancan->cd(1);
  conscan->Draw("ALP");
  dioscan->Draw("same");
  info->Draw();
  leg->Draw();

  int ilow = nscan-15;
  int ihi = nscan+15;

  TLine* klow = new TLine(diox[ilow],-0.5,diox[ilow],dioy[ilow]);
  TLine* khi = new TLine(diox[ihi],-0.5,diox[ihi],dioy[ihi]);
  TLine* klowx = new TLine(-150,dioy[ilow],diox[ilow],dioy[ilow]);
  TLine* khix = new TLine(-150,dioy[ihi],diox[ihi],dioy[ihi]);
  klow->SetLineStyle(2);
  khi->SetLineStyle(2);
  klowx->SetLineStyle(2);
  khix->SetLineStyle(2);
  klow->Draw("same");
  khi->Draw("same");
  klowx->Draw("same");
  khix->Draw("same");
  scancan->SaveAs((string("mu2e_scan")+ssuf).c_str());




// some numerical values
  vector<double> dioshift,conshift,dioshift_err,conshift_err,shift,shift_err;
  vector<double> fdioshift,fconshift,fdioshift_err,fconshift_err;
  double diocent = diospec[icut]->Integral(istart,istop);
  double concent = conspec[icut]->Integral(istart,istop);
  for(int ires=-nscan;ires<=nscan;++ires){
    double ddioint_err,dconint_err;
    double ddioint,dconint;
    if(istart+ires<istart){
      ddioint = diospec[icut]->IntegralAndError(istart+ires,istart,ddioint_err);
      dconint = conspec[icut]->IntegralAndError(istart+ires,istart,dconint_err);
    } else {
      ddioint = -diospec[icut]->IntegralAndError(istart,istart+ires,ddioint_err);
      dconint = -conspec[icut]->IntegralAndError(istart,istart+ires,dconint_err);
    }
    shift.push_back(ires*mevperbin);
    shift_err.push_back(0.0);
    dioshift.push_back(ddioint);
    conshift.push_back(dconint);
    dioshift_err.push_back(ddioint_err);
    conshift_err.push_back(dconint_err);
    fdioshift.push_back(ddioint/diocent);
    fconshift.push_back(dconint/concent);
    fdioshift_err.push_back(ddioint_err/diocent);
    fconshift_err.push_back(dconint_err/concent);
//    cout << "For lower-bound shift of " << ires*mevperbin << " MeV/c DIO integral change = " << ddioint
//      << " , conversion integral change = " << dconint << endl;
  }
  TGraphErrors* dioscan2 = new TGraphErrors(shift.size(),&shift[0],&dioshift[0],&shift_err[0],&dioshift_err[0]);
  TGraphErrors* conscan2 = new TGraphErrors(shift.size(),&shift[0],&conshift[0],&shift_err[0],&conshift_err[0]);
  TGraphErrors* fdioscan2 = new TGraphErrors(shift.size(),&shift[0],&fdioshift[0],&shift_err[0],&fdioshift_err[0]);
  TGraphErrors* fconscan2 = new TGraphErrors(shift.size(),&shift[0],&fconshift[0],&shift_err[0],&fconshift_err[0]);
  dioscan2->SetLineColor(kBlue);
  conscan2->SetLineColor(kRed);
  dioscan2->SetMarkerColor(kBlue);
  conscan2->SetMarkerColor(kRed);
  conscan2->SetMaximum(0.6);
  conscan2->SetMinimum(-0.8);
  fdioscan2->SetLineColor(kBlue);
  fconscan2->SetLineColor(kRed);
  fdioscan2->SetMarkerColor(kBlue);
  fconscan2->SetMarkerColor(kRed);
  fconscan2->SetMaximum(3.0);
  fconscan2->SetMinimum(-1.0);
  dioscan2->SetTitle("Rate Sensitivity to Lower Momentum Cut;#Delta P (MeV/c);#Delta I");
  conscan2->SetTitle("Rate Sensitivity to Lower Momentum Cut;#Delta P (MeV/c);#Delta I");
  fdioscan2->SetTitle("Fractional Rate Sensitivity to Lower Momentum Cut;#Delta P (MeV/c);#Delta I/I");
  fconscan2->SetTitle("Fractional Rate Sensitivity to Lower Momentum Cut;#Delta P (MeV/c);#Delta I/I");
  TCanvas* scancan2 = new TCanvas("scan2","Scan of lower threshold",800,800);
  scancan2->Divide(1,2);
  scancan2->cd(1);
  conscan2->Draw("ALP");
  dioscan2->Draw("same");
  info->Draw();
  leg->Draw();
  scancan2->cd(2);
  fconscan2->Draw("ALP");
  fdioscan2->Draw("same");
  info->Draw();
  leg->Draw();
  scancan2->SaveAs((string("mu2e_scan2")+ssuf).c_str());
}

