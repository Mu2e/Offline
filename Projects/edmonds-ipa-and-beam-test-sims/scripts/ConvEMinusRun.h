#ifndef ConvEMinusRun_h_
#define ConvEMinusRun_h_

#include "BaseRun.h"

#include "TH1.h"
#include "TH3.h"
#include "TFile.h"
#include "TTree.h"
#include "TFitResult.h"
#include "TLine.h"
#include "TF1.h"
#include "TLatex.h"
#include "TMath.h"

#include <iostream>
#include <sstream>

class ConvEMinusRun : public BaseRun {

public:
  ConvEMinusRun(std::string filename, double min_mom = 100, double max_mom = 106, double bin_width = 0.1);

private:
  TH1F* fRecoHist;

  double fPeakEnergy;
  double fPeakEnergyError;
  double fPeakEnergyCount;
  double fPeakEnergyCountError;
  TF1* fPeakFn;

  double fLowEdge; // the low momentum edge of the FWHM
  double fHighEdge; // the high momentum edge of the FWHM
  double fFWHM;
  double fFWHMError;
  TF1* fLowEdgeFn;
  TF1* fHighEdgeFn;

  double fSignalEnergy;
  double fSignalRate;

  double fMinMomentum;
  double fMaxMomentum;
  double fBinWidth;

  void CalculateFWHMAndError();

  static double fnc_dscb(double*xx,double*pp); // crystal ball function
  TF1* fCrystalBallFit;

public:
  TH1F* GetRecoHist() { if(fRecoHist) {return fRecoHist;} else {std::cout << "No histogram for this run" << std::endl;} }
  double GetPeakEnergy() { return fPeakEnergy; }
  double GetPeakEnergyError() { return fPeakEnergyError; }
  double GetFWHM() { return fFWHM; }
  double GetFWHMError() { return fFWHMError; }
  double GetDeltaE() { return fSignalEnergy - fPeakEnergy; }

  void FitCrystalBall();
  void Draw(std::string drawopt  = "HIST E");
  double GetSignalRate() { return fSignalRate; }
};

ConvEMinusRun::ConvEMinusRun(std::string filename, double min_mom, double max_mom, double bin_width) 
                     : BaseRun(filename, "ConvEMinus"), fMinMomentum(min_mom), fMaxMomentum(max_mom), fBinWidth(bin_width) {

  int n_bins = (fMaxMomentum - fMinMomentum) / fBinWidth;
  fRecoHist = new TH1F("fRecoHist", "", n_bins,fMinMomentum,fMaxMomentum);
  fRecoHist->SetXTitle("Reco Momentum [MeV]");
  std::stringstream axistitle; axistitle << "Counts per " << fBinWidth << " MeV";
  fRecoHist->SetYTitle(axistitle.str().c_str());
  //    fRecoHist->SetStats(false);
  
  GetTrkDiagChain()->Draw("fit.mom>>fRecoHist", "fit.status>0", "goff");
  fRecoHist->SetDirectory(0);

  // Set some initial estimates of the peak energy and fwhm
  fPeakEnergy = fRecoHist->GetBinCenter(fRecoHist->GetMaximumBin());
  fPeakEnergyError = 0;
  fPeakEnergyCount = fRecoHist->GetBinContent(fRecoHist->GetMaximumBin());
  fPeakEnergyCountError = 0;
  
  fCrystalBallFit = NULL;
  FitCrystalBall();

  fSignalEnergy = 104.97;
  fSignalRate = 3e-17;
  SetNParticlesPerMicrobunch(n_POT_per_microbunch * n_stopped_muons_per_POT * n_decayed_muons_per_stopped_muon * fSignalRate);
}

void ConvEMinusRun::CalculateFWHMAndError() {
  double half_max = fPeakEnergyCount / 2;
  int bin1 = fRecoHist->FindFirstBinAbove(half_max);
  int bin2 = fRecoHist->FindLastBinAbove(half_max);

  // Get rough estimates of the high and low edge
  fLowEdge = fRecoHist->GetBinLowEdge(bin1);
  fHighEdge = fRecoHist->GetBinLowEdge(bin2) + fRecoHist->GetBinWidth(bin2);
  fFWHM = fHighEdge - fLowEdge;
  //  std::cout << "AE: At Start: Low Edge = " << fLowEdge << ", High Edge = " << fHighEdge << ", FWHM = " << fFWHM << std::endl;

  fLowEdge = fCrystalBallFit->GetX(half_max, fLowEdge-0.5, fLowEdge+0.5); // find the low edge
  fHighEdge = fCrystalBallFit->GetX(half_max, fHighEdge-0.5, fHighEdge+0.5); // find the low edge
  fFWHM = fHighEdge - fLowEdge;
  //    std::cout << "AE: Crystal Ball: Low Edge = " << fLowEdge << ", High Edge = " << fHighEdge << ", FWHM = " << fFWHM << std::endl;
  //    std::cout << half_max << ", " << fCrystalBallFit->Eval(fLowEdge) << ", " << fCrystalBallFit->Eval(fHighEdge) << std::endl;
  
  double higher_half_max = (fPeakEnergyCount + fPeakEnergyCountError) / 2;
  double higher_low_edge = fCrystalBallFit->GetX(higher_half_max, fLowEdge-0.5, fLowEdge+0.5); // find the low edge
  double higher_high_edge = fCrystalBallFit->GetX(higher_half_max, fHighEdge-0.5, fHighEdge+0.5); // find the high edge
  double higher_fwhm = higher_high_edge - higher_low_edge;
  //    std::cout << "AE: Crystal Ball (Higher): Low Edge = " << higher_low_edge << ", High Edge = " << higher_high_edge << ", FWHM = " << higher_fwhm << std::endl;
  //    std::cout << higher_half_max << ", " << fCrystalBallFit->Eval(higher_low_edge) << ", " << fCrystalBallFit->Eval(higher_high_edge) << std::endl;
  
  double lower_half_max = (fPeakEnergyCount - fPeakEnergyCountError) / 2;
  double lower_low_edge = fCrystalBallFit->GetX(lower_half_max, fLowEdge-0.5, fLowEdge+0.5); // find the low edge
  double lower_high_edge = fCrystalBallFit->GetX(lower_half_max, fHighEdge-0.5, fHighEdge+0.5); // find the high edge
  double lower_fwhm = lower_high_edge - lower_low_edge;
  //    std::cout << "AE: Crystal Ball (Lower): Low Edge = " << lower_low_edge << ", High Edge = " << lower_high_edge << ", FWHM = " << lower_fwhm << std::endl;
  //    std::cout << lower_half_max << ", " << fCrystalBallFit->Eval(lower_low_edge) << ", " << fCrystalBallFit->Eval(lower_high_edge) << std::endl;
  
  //    std::cout << "AE: Differences: " << std::fabs(fFWHM - higher_fwhm) << ", " << std::fabs(fFWHM - lower_fwhm) << std::endl;
  fFWHMError = ( std::fabs(fFWHM - higher_fwhm)+std::fabs(fFWHM - lower_fwhm) ) / 2;
}

void ConvEMinusRun::Draw(std::string drawopt) {
  fRecoHist->SetLineWidth(2);
  fRecoHist->SetLineColor(kBlue);
  fRecoHist->SetMaximum(fPeakEnergyCount+fPeakEnergyCountError+10);

  TLine* peak = new TLine(fPeakEnergy, 0, fPeakEnergy, fPeakEnergyCount);
  peak->SetLineColor(kBlue);
  peak->SetLineWidth(2);
  peak->SetLineStyle(2);

  TLine* peak_count_error = new TLine(fPeakEnergy, fPeakEnergyCount-fPeakEnergyCountError, fPeakEnergy, fPeakEnergyCount+fPeakEnergyCountError);
  peak_count_error->SetLineColor(kRed);
  peak_count_error->SetLineWidth(2);

  TLine* fwhm_line = new TLine(fLowEdge, fPeakEnergyCount/2, fHighEdge, fPeakEnergyCount/2);
  fwhm_line->SetLineColor(kBlue);
  fwhm_line->SetLineWidth(2);
  fwhm_line->SetLineStyle(2);

  fRecoHist->Draw(drawopt.c_str());

  std::stringstream latex;
  latex.precision(3);
  latex << "#splitline{#DeltaE = " << GetDeltaE() << " #pm ";
  latex.precision(1);
  latex << fPeakEnergyError << " MeV}{FWHM = ";
  latex.precision(2);
  latex << fFWHM << " #pm ";
  latex.precision(1);
  latex << fFWHMError << " MeV}";

  TLatex* text = new TLatex(101, 3*(fPeakEnergyCount/4), latex.str().c_str());
  text->SetTextSize(0.05);
  text->Draw("SAME");

  if (fCrystalBallFit) {
    fCrystalBallFit->SetLineWidth(3);
    fCrystalBallFit->SetLineColor(kRed);
    fCrystalBallFit->Draw("LSAME");
  }
  else if (fPeakFn && fLowEdgeFn && fHighEdgeFn) {
    fPeakFn->SetLineWidth(3);
    fLowEdgeFn->SetLineWidth(3);
    fHighEdgeFn->SetLineWidth(3);
    
    fPeakFn->SetLineColor(kBlack);
    fLowEdgeFn->SetLineColor(kRed);
    fHighEdgeFn->SetLineColor(kRed);

    fPeakFn->Draw("LSAME");
    fLowEdgeFn->Draw("LSAME");
    fHighEdgeFn->Draw("LSAME");
  }
  peak->Draw("LSAME");
  peak_count_error->Draw("LSAME");
  fwhm_line->Draw("LSAME");
}

// CrystalBall function taken from Offline/KalmanTests/test/KalFit.C
double ConvEMinusRun::fnc_dscb(double*xx,double*pp) {
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

void ConvEMinusRun::FitCrystalBall() {
  fCrystalBallFit= new TF1("dscb",fnc_dscb,fRecoHist->GetBinLowEdge(1), fRecoHist->GetBinLowEdge(fRecoHist->GetNbinsX())+fBinWidth, 7);
  fCrystalBallFit->SetParName(0,"Norm");
  fCrystalBallFit->SetParName(1,"x0");
  fCrystalBallFit->SetParName(2,"sigma");
  fCrystalBallFit->SetParName(3,"ANeg");
  fCrystalBallFit->SetParName(4,"PNeg");
  fCrystalBallFit->SetParName(5,"APos");
  fCrystalBallFit->SetParName(6,"PPos");

  double integral = fRecoHist->GetEntries() * fBinWidth;
  fCrystalBallFit->SetParameters(3*integral,fRecoHist->GetMean()+0.07,0.3*fRecoHist->GetRMS(),1.0,4.0,1.0,5.0);

  TFitResultPtr fit = fRecoHist->Fit("dscb", "SQR0");

  int fit_status = fit;
  if (fit_status != 0) {
    std::cout << "Crystal ball fit failed" << std::endl;
  }
  else {
    fPeakEnergy = fit->Parameter(1);
    fPeakEnergyError = fit->Error(1);
    fPeakEnergyCount = fit->Parameter(0);
    fPeakEnergyCountError = fit->Error(0);
    CalculateFWHMAndError(); // now recalculate the FWHM
  }

}

#endif
