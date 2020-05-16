#include "CaloMC/inc/CaloPulseShape.hh"
#include "ConditionsService/inc/CalorimeterCalibrations.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"

#include <string>
#include <vector>
#include <iostream>
#include <algorithm>

#include "TFile.h"
#include "TH2F.h"
#include "TSpline.h"


namespace mu2e {

   CaloPulseShape::CaloPulseShape(double digiSampling, int pulseIntegralSteps) :
     digiSampling_(digiSampling), pulseIntegralSteps_(pulseIntegralSteps), nBinShape_(0), integralVal_(), digitizedPulse_()
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
       pulseShape.Scale(1.0/pulseShape.Integral());


       // find the start / end of the waveform
       double  pulseAmplitudeMax = pulseShape.GetBinContent(pulseShape.GetMaximumBin());
       int ifirst(pulseIntegralSteps_),ilast(nbins-pulseIntegralSteps_-1); 
       for (;ifirst<nbins;++ifirst) if (pulseShape.GetBinContent(ifirst)> pulseAmplitudeMax/1000) break;
       for (;ilast >0    ; --ilast) if (pulseShape.GetBinContent(ilast) > pulseAmplitudeMax/1000) break;
       
       if (ifirst==pulseIntegralSteps_) 
          std::cout<<"[CaloPulseShape] Warning, pulse histogram starts too early and will be clipped "<<std::endl;

       // calculate the integral over the sampling period for each start time
       integralVal_.clear();
       for (int j=0;j<ilast-ifirst+pulseIntegralSteps_;++j)
       {
	  int binStart = ifirst-pulseIntegralSteps_+j;
	  int binEnd   = binStart+pulseIntegralSteps_-1; 
	  integralVal_.push_back(pulseShape.Integral(binStart, binEnd));
       }
       
       nBinShape_      = (ilast-ifirst)/pulseIntegralSteps_+2;
       digitizedPulse_ = std::vector<double>(nBinShape_,0);

       pulseFile.Close();
   }
   

   //----------------------------------------------------------------------------
   // forward shift in waveform = backward shift in time origin 
   const std::vector<double>& CaloPulseShape::digitizedPulse(double hitTime) const
   {
       int shiftBin = pulseIntegralSteps_ - int(hitTime*pulseIntegralSteps_/digiSampling_)%pulseIntegralSteps_;
       for (int i=0;i<nBinShape_;++i) digitizedPulse_[i] = integralVal_[shiftBin+i*pulseIntegralSteps_];       
       return digitizedPulse_;
   }
   

   //--------------------------------
   void CaloPulseShape::diag() const
   {
       std::cout<<"Integral val "<<std::endl;
       for (auto& h : integralVal_) std::cout<<h<<" ";
       std::cout<<std::endl;
       
       std::cout<<"Number of digi bins "<<nBinShape_<<std::endl;
       
       std::cout<<"Cache content "<<std::endl;
       for (auto& h : digitizedPulse_) std::cout<<h<<" ";
       std::cout<<std::endl;
   }

}
