#include "Mu2eUtilities/inc/CaloPulseShape.hh"
#include "ConditionsService/inc/CalorimeterCalibrations.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"

#include "TFile.h"
#include "TH2F.h"
#include "TSpline.h"

#include <string>
#include <vector>
#include <iostream>
#include <algorithm>


namespace mu2e {


   CaloPulseShape::CaloPulseShape(double digiSampling, int pulseIntegralSteps, bool doIntegral) :
     digiSampling_(digiSampling), pulseIntegralSteps_(pulseIntegralSteps), doIntegral_(doIntegral), 
     nBinShape_(0.), integralVal_(), deltaT_(0.), digiStep_(digiSampling/pulseIntegralSteps), digitizedPulse_()
   {}

   //----------------------------------------------------------------------------------------------------------------------
   void CaloPulseShape::buildShapes()
   {
       ConditionsHandle<CalorimeterCalibrations> calorimeterCalibrations("ignored");
       std::string fileName = calorimeterCalibrations->pulseFileName();
       std::string histName = calorimeterCalibrations->pulseHistName();

       TH1F *pshape(0);
       TFile pulseFile(fileName.c_str());
       if (pulseFile.IsOpen()) pshape = (TH1F*) pulseFile.Get(histName.c_str());

       if (!pshape) throw cet::exception("CATEGORY")<<"CaloPulseShape:: Hitsogram "<<histName.c_str()
                                                   <<" from file "<<fileName.c_str()<<" does not exist";

       // smooth the histogram, normalize and resmaple with the desired binning
       int nbins = int((pshape->GetXaxis()->GetXmax()-pshape->GetXaxis()->GetXmin())/digiSampling_*pulseIntegralSteps_);
       TH1F pulseShape("ps_temp","ps_temp", nbins, pshape->GetXaxis()->GetXmin(), pshape->GetXaxis()->GetXmax());       
       TSpline3 spline(pshape);
       for (int i=1;i<pulseShape.GetNbinsX();++i) pulseShape.SetBinContent(i,spline.Eval(pulseShape.GetBinCenter(i)));

       // Integrated waveform over digitization period for each starting point if required
       if (doIntegral_) 
       {
	   TH1F pulseShapeIntegral(pulseShape);
	   for (int i=1;i<pulseShapeIntegral.GetNbinsX();++i) 
	      pulseShapeIntegral.SetBinContent(i,pulseShape.Integral(i, std::min(i+pulseIntegralSteps_-1,nbins-1)));
	   pulseShape = pulseShapeIntegral;
       }

       // Copy histogram into vector, including time padding      
       for (int j=1;j<=nbins;++j)
       {
	  int ibin = j-pulseIntegralSteps_;
	  integralVal_.push_back((ibin>0) ? pulseShape.GetBinContent(ibin) : 0.0);
       }

       // Normalize waveform such that the max vector is at 1
       double integralValMax = *max_element(integralVal_.begin(),integralVal_.end());
       for (auto& v : integralVal_) v/=integralValMax;
       
       //calculate the number of bins for the digitized waveform
       nBinShape_      = int(nbins/pulseIntegralSteps_);
       digitizedPulse_ = std::vector<double>(nBinShape_,0);
              
       // find difference between peak time and t0 for digitized waveform, so that we can shift the result of 
       // the fitted waveform back to t0. 
       int imax(1);
       for (;imax<nBinShape_;++imax) if (integralVal_[(imax+1)*pulseIntegralSteps_] < integralVal_[imax*pulseIntegralSteps_]) break;
       deltaT_ = ((imax-1)*pulseIntegralSteps_)*digiStep_;

       pulseFile.Close();
   }

   //----------------------------------------------------------------------------
   // forward shift in waveform = backward shift in time origin 
   const std::vector<double>& CaloPulseShape::digitizedPulse(double hitTime) const
   {
       int shiftBin = pulseIntegralSteps_ - int(hitTime/digiStep_)%pulseIntegralSteps_;
       for (int i=0;i<nBinShape_;++i) digitizedPulse_[i] = integralVal_[shiftBin+i*pulseIntegralSteps_];      
       return digitizedPulse_;
   }

   //----------------------------------------------------------------------------
   double CaloPulseShape::evaluateFromPeak(double tDifference) const
   {
       double t = tDifference+deltaT_;          
       int ibin = pulseIntegralSteps_ + int(t*pulseIntegralSteps_/digiSampling_);

       if (ibin < 0 || ibin >= int(integralVal_.size()-1)) return 0.0;
       
       double t0bin = (ibin-pulseIntegralSteps_)*digiStep_; //t0 is located at pulseIntegralSteps_            
       return (integralVal_[ibin+1]-integralVal_[ibin])/digiStep_*(t-t0bin)+integralVal_[ibin];                  
   }
   
   //----------------------------------------------------------------------------
   double CaloPulseShape::fromPeakToT0(double timePeak) const
   {
       return timePeak-deltaT_-0.5*digiSampling_;
   }

   //--------------------------------
   void CaloPulseShape::diag(bool fullDiag) const
   {
       std::cout<<"Number of digi bins "<<nBinShape_<<std::endl;
       
       std::cout<<"Cache content "<<std::endl;
       for (auto& h : digitizedPulse_) std::cout<<h<<" ";
       std::cout<<std::endl;
       
       if (!fullDiag) return;
       
       std::cout<<"Integral val "<<std::endl;
       for (auto& h : integralVal_) std::cout<<h<<" ";
       std::cout<<std::endl;
   }
}


//
// Additional notes to avoid scractching your head too long
//
// Detailed explanation of waveform construction and offsets
//
// Shifting the waveform forward by a time dt is equivalent to shifting the time origin backward by dt
// with original waveform: 
//  - value at pulseIntegralSteps_-nbin correspond to waveform shifted by time pulseIntegralSteps_+nbin, 
//    and we need pulseIntegralSteps_ bins to store the different values. 
//  - t0 value is located at bin pulseIntegralSteps_
//  - if the waveform "start" is located at a time smaller than digiSampling_, we extend the waveform 
//    with zeroes on the left side (equivalently, pad the integralVal_ vector)
//
// The t0 value correspond to the content of the digitized bin whose LOW EDGE is at time t0
//

