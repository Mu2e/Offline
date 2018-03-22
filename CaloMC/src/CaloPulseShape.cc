//
// Utility class to hold Calorimeter readout pulse shape
//

#include "CaloMC/inc/CaloPulseShape.hh"
#include "ConditionsService/inc/CalorimeterCalibrations.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"

#include <string>
#include <vector>
#include <iostream>
#include <algorithm>

#include "TFile.h"
#include "TH2F.h"


namespace mu2e {


   CaloPulseShape::CaloPulseShape(double digiSampling, int pulseIntegralSteps) :
     digiSampling_(digiSampling), pulseIntegralSteps_(pulseIntegralSteps),pulseDigitized_()
   {
   }


   //----------------------------------------------------------------------------------------------------------------------
   void CaloPulseShape::buildShapes()
   {
       pulseDigitized_.clear();
       ConditionsHandle<CalorimeterCalibrations> calorimeterCalibrations("ignored");
       std::string fileName = calorimeterCalibrations->pulseFileName();
       std::string histName = calorimeterCalibrations->pulseHistName();

       TH1F *pshape(0);
       TFile pulseFile(fileName.c_str());
       if (pulseFile.IsOpen()) pshape = (TH1F*) pulseFile.Get(histName.c_str());

         if (!pshape) throw cet::exception("CATEGORY")<<"CaloPulseShape:: Hitsogram "<<histName.c_str()
	                                              <<" from file "<<fileName.c_str()<<" does not exist";

          double  pulseBinWidth = pshape->GetBinWidth(2);
          int     nBinTimeStamp = digiSampling_ / pulseBinWidth;
          int     integralStep  = digiSampling_ / pulseBinWidth / pulseIntegralSteps_;

          int ifirst(1),ilast(pshape->GetNbinsX()); //set bins prior to T0 to zero.
          for (; ifirst<pshape->GetNbinsX();++ifirst) {if (pshape->GetBinContent(ifirst) > 1e-3) break; pshape->SetBinContent(ifirst,0);}
          for (; ilast > 0                 ; --ilast) {if (pshape->GetBinContent(ilast)  > 1e-3) break; pshape->SetBinContent(ilast,0);}

          int  nTimeStamps = int((ilast-ifirst)*pulseBinWidth/digiSampling_) + 1;

          for (int i=0; i<pulseIntegralSteps_; ++i)
          {
               double sum(0);
               std::vector<double> pulseDigi;

               int binInit = ifirst - i*integralStep; //shift start time backward
               for (int j=0; j<nTimeStamps; ++j)
               {
                  int binMin = binInit+j*nBinTimeStamp;
                  double pulseIntegral = pshape->Integral(binMin, binMin+nBinTimeStamp-1);
                  sum += pulseIntegral;
                  pulseDigi.push_back(pulseIntegral);
               }
               for (auto &pulse : pulseDigi) pulse /= sum;
               pulseDigitized_.push_back(pulseDigi);
          }

	pulseFile.Close();
   }

   //----------------------------------------------------------------------------------------------------------------------
   const void CaloPulseShape::printShape() const
   {
       int ip(0);
       std::cout<<"[CaloPulseShape] Pulse "<<std::endl;
       for (const auto& pulseDigi : pulseDigitized_)
       {
	  std::cout<<"Pulse "<<ip<<std::endl;
          for (const auto& pulse : pulseDigi) {std::cout<<pulse<<" ";} std::cout<<std::endl;
          ++ip;
       }
   }



}
