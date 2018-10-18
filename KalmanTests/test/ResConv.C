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

// the following approximation is from Czarnecki etal, 'Muon decay in orbit:spectrum of high-energy electrons',
// for E>85 MeV
double capfrac(0.609);
//const double nmu(6.7e17);
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

// The following is from Alexx Perloff, JetMetaAnalysis
double fnc_cemomf(double*xx,double*pp) {
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

void DrawCE(double nmu, double rmue, int nsteps,double resfac)  {
  double lowlim(100.0);
  double hilim(106.0);
  TH1F* cemom = new TH1F("cemom","Momentum;MeV/c;Events/MeV/c",200,lowlim,hilim);
  cemom->SetStats(0);
  TF1* cemomf = new TF1("cemomf",fnc_cemomf,lowlim,hilim,7);
  cemomf->SetParName(0,"Norm");
  cemomf->SetParName(1,"x0");
  cemomf->SetParName(2,"sigma");
  cemomf->SetParName(3,"ANeg");
  cemomf->SetParName(4,"PNeg");
  cemomf->SetParName(5,"APos");
  cemomf->SetParName(6,"PPos");
  double norm = nmu*capfrac*rmue*0.01107*31.375/5.0/resfac;
  cout << "setting norm to " << norm << endl;
  double mean = 105-0.608;
  cemomf->SetParameters(norm,mean,0.3235*resfac,0.6929,2.799,2.58,2.612+1.0);
  cemomf->SetNpx(1000);
  cemom->SetMaximum(1.5*cemomf->Eval(mean));
  TCanvas* cecan = new TCanvas("cecan","cecan",800,800);
  cecan->cd(1);
//  gPad->SetLogy();
  cemom->Draw();
  cemomf->Draw("same");
// now DIO
  TF1* diot = new TF1("diot",DIOCZ,lowlim,hilim,1);
  double dionorm = nmu*(1.0-capfrac);
  diot->SetParameter(0,dionorm);
  diot->SetLineColor(kBlue);
//  diot->Draw("same");
  TH1F* diomom = new TH1F("diomom","Momentum;MeV/c;Events/MeV/c",200,lowlim,hilim);
  diomom->SetStats(0);
  diomom->Sumw2();
  TF1* resf = new TF1("resf",fnc_cemomf,-5,5.0,7);
  resf->SetParName(0,"Norm");
  resf->SetParName(1,"x0");
  resf->SetParName(2,"sigma");
  resf->SetParName(3,"ANeg");
  resf->SetParName(4,"PNeg");
  resf->SetParName(5,"APos");
  resf->SetParName(6,"PPos");
  double rnorm = 0.01107*31.375/resfac;
  double rmean = -0.608;
  resf->SetParameters(rnorm,rmean,0.3235*resfac,0.6929,2.799,2.58,2.612);
  resf->SetNpx(1000);
  double step = (hilim-lowlim)/nsteps;
  double rstep = 5.0/nsteps;
  for(unsigned istep=0;istep<nsteps;++istep){
    double mom = lowlim + istep*step;
    double diov = diot->Eval(mom);
    for(int jstep=-nsteps;jstep<nsteps;++jstep){
      double dmom = jstep*rstep;
      double rfac = resf->Eval(dmom);
//      cout << "mom = " << mom << " diov = " << diov << " dmom = " << dmom << " rfac = " << rfac << endl;
      diomom->Fill(mom+dmom,diov*rfac);
    }
  }
  diomom->Scale(1000*0.5/(nsteps*nsteps)/(hilim-lowlim));
  diomom->SetLineColor(kBlue);
  diomom->Draw("same");
  // labels
  TPaveText* rtext = new TPaveText(0.1,0.7,0.4,0.9,"NDC");
  char line[80];
  sprintf(line,"N stopped #mu = %G",nmu);
  rtext->AddText(line);
  sprintf(line,"R_{#mue} = %G",rmue);
  rtext->AddText(line);
  rtext->Draw();
}
