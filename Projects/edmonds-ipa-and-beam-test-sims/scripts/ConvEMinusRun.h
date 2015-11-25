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
  ConvEMinusRun(std::string filename);

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
  void FitPeak();
  void FitCrystalBall();
  void Draw(std::string drawopt  = "HIST E");
  double GetSignalRate() { return fSignalRate; }
};

ConvEMinusRun::ConvEMinusRun(std::string filename) : BaseRun(filename, "eMinus") {

  double min_mom = 100; double max_mom = 106; fBinWidth = 0.1;
  int n_bins = (max_mom - min_mom) / fBinWidth;
  fRecoHist = new TH1F("fRecoHist", "", n_bins,min_mom,max_mom);
  fRecoHist->SetXTitle("Reco Momentum [MeV]");
  fRecoHist->SetYTitle("Counts per 0.1 MeV");
  //    fRecoHist->SetStats(false);
  
  fTrkDiagChain->Draw("fit.mom>>fRecoHist", "fit.status>0", "goff");
  fRecoHist->SetDirectory(0);

  // Set some initial estimates of the peak energy and fwhm
  fPeakEnergy = fRecoHist->GetBinCenter(fRecoHist->GetMaximumBin());
  fPeakEnergyError = 0;
  fPeakEnergyCount = fRecoHist->GetBinContent(fRecoHist->GetMaximumBin());
  fPeakEnergyCountError = 0;
  
  fCrystalBallFit = NULL;
  CalculateFWHMAndError();

  fSignalEnergy = 104.97;
  fSignalRate = 3e-17;
  fNParticlesPerMicrobunch = n_POT_per_microbunch * n_stopped_muons_per_POT * n_decayed_muons_per_stopped_muon * fSignalRate;
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

  if (fCrystalBallFit) {
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
    //    std::cout << "FWHM Error = " << fFWHMError << std::endl;
  }
  else {
    fLowEdgeFn = new TF1("low_edge", "pol1", fLowEdge-fBinWidth, fLowEdge+fBinWidth);
    fHighEdgeFn = new TF1("high_edge", "pol1", fHighEdge-fBinWidth, fHighEdge+fBinWidth);

    TFitResultPtr grad_fit_low = fRecoHist->Fit(fLowEdgeFn, "QSR0");
    TFitResultPtr grad_fit_high = fRecoHist->Fit(fHighEdgeFn, "QSR0");
    
    if ( (int) grad_fit_low == 0 && (int) grad_fit_high == 0) {
      
      // Calculate the central value
      double gradient_low = grad_fit_low->Parameter(1);
      double offset_low = grad_fit_low->Parameter(0);
      
      double gradient_high = grad_fit_high->Parameter(1);
      double offset_high = grad_fit_high->Parameter(0);
      
      fLowEdge = (half_max - offset_low)/gradient_low;
      fHighEdge = (half_max - offset_high)/gradient_high;
      
      fFWHM = fHighEdge - fLowEdge;
      
      // Calculate the extremes
      double half_max_plus_error = (fPeakEnergyCount + fPeakEnergyCountError)/2;
      double new_low_edge = (half_max_plus_error - offset_low)/gradient_low;
      double new_high_edge = (half_max_plus_error - offset_high)/gradient_high;
      double fwhm_plus = new_high_edge - new_low_edge;
      double fwhm_plus_error = std::fabs(fFWHM - fwhm_plus);
      
      double half_max_minus_error = (fPeakEnergyCount - fPeakEnergyCountError)/2;
      new_low_edge = (half_max_minus_error - offset_low)/gradient_low;
      new_high_edge = (half_max_minus_error - offset_high)/gradient_high;
      double fwhm_minus = new_high_edge - new_low_edge;
      double fwhm_minus_error = std::fabs(fFWHM - fwhm_minus);
      
      fFWHMError = (fwhm_plus_error + fwhm_minus_error)/2;
      //    std::cout << "AE: " << fPeakEnergyCount/2 << ", " << half_max_plus_error << ", " << half_max_minus_error << std::endl;
      //    std::cout << "AE: " << fFWHM << ", " << fwhm_plus << ", " << fwhm_minus << std::endl;
      //    std::cout << "AE: " << fFWHMError << ", " << fwhm_plus_error << ", " << fwhm_minus_error << std::endl;
    }
    else {
      std::cout << "FWHM gradient fits failed" << std::endl;
    }
  }
}

void ConvEMinusRun::FitPeak() {

  // Fit the current peak energy to a gaussian
  fPeakFn = new TF1("gaus", "gaus", fPeakEnergy-3*fBinWidth, fPeakEnergy+5*fBinWidth);
  TFitResultPtr fit = fRecoHist->Fit(fPeakFn, "SQR0");

  int fit_status = fit;
  if (fit_status != 0) {
    std::cout << "Peak energy fit failed" << std::endl;
  }
  else {
    fPeakEnergy = fit->Parameter(1);
    fPeakEnergyError = fit->Error(1);
    fPeakEnergyCount = fit->Parameter(0);
    fPeakEnergyCountError = fit->Error(0);
    CalculateFWHMAndError(); // now recalculate the FWHM
  }
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

/*TH1F* ConvEMinusRun::GetHitsPerLengthPlot(int device, int sector, int layer, int straw) {

  fFile->cd(); // make sure we're on the correct file so that we can draw onto the correct histogram and not create a new one
  fHitsPerLengthPlot->Reset();
  std::stringstream cut;
  if (device >= 0) {
    cut << "tsh._device==" << device << " && ";
  }
  if (sector >= 0) {
    cut << "tsh._sector==" << sector << " && ";
  }
  if (layer >= 0) {
    cut << "tsh._layer==" << layer << " && ";
  }
  if (straw >= 0) {
    cut << "tsh._straw==" << straw << " && ";
  }
  cut << "tshmc._pdg==11";

  fTrkDiagChain->Draw("tshmc._len>>fHitsPerLengthPlot", cut.str().c_str(), "goff");
  fHitsPerLengthPlot->Scale(1.0/fTotalCEStrawHits);

  return fHitsPerLengthPlot;
}
*/

#endif
