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

   CaloPulseShape::CaloPulseShape(double digiSampling, int pulseIntegralSteps, bool doIntegral) :
     digiSampling_(digiSampling), pulseIntegralSteps_(pulseIntegralSteps), doIntegral_(doIntegral), 
     nBinShape_(0), integralVal_(), digitizedPulse_()
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

       // Find the start of the waveform - this is arbitrary but we want to avoid taking too many empty bins before the waveform
       int ifirst(0),ilast(nbins-1); 
       double pulseAmplitudeMin = pulseShape.GetMaximum()/1000;
       for (;ifirst<nbins;++ifirst) if (pulseShape.GetBinContent(ifirst)> pulseAmplitudeMin) break;
       for (;ilast >0    ; --ilast) if (pulseShape.GetBinContent(ilast) > pulseAmplitudeMin) break;       

       // Tabulate the integrated values - need padding if waveform starts before histogram t0 + digitization
       nbins = ilast-ifirst+pulseIntegralSteps_;
       std::vector<double> integralVal(nbins,0.0);
       for (int j=0;j<nbins;++j)
       {
	  int ibin = ifirst-pulseIntegralSteps_+j;
	  integralVal_.push_back((ibin>0) ? pulseShape.GetBinContent(ibin) : 0.0);
       }

       // Normalize waveform such that the max vector is at 1
       double integralValMax = *max_element(integralVal_.begin(),integralVal_.end());
       for (auto& v : integralVal_) v/=integralValMax;
       
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
