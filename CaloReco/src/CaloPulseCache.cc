//
// Utility class to hold cached calorimeter readout pulse shape
//
#include "CaloReco/inc/CaloPulseCache.hh"

#include "cetlib_except/exception.h"
#include "ConditionsService/inc/CalorimeterCalibrations.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
 
#include <string> 
#include <vector>
#include <iostream>
#include <algorithm> 

#include "TH1F.h"
#include "TFile.h"


namespace mu2e {


   CaloPulseCache::CaloPulseCache() : cache_(), cacheSize_(0), deltaT_(0), factor_(0), step_(0)
   {
   }       


   //----------------------------------------------------------------------------------------------------------------------
   void CaloPulseCache::initialize()
   {   
      
       cache_.clear();
       
       ConditionsHandle<CalorimeterCalibrations> calorimeterCalibrations("ignored");
       std::string fileName = calorimeterCalibrations->pulseFileName();
       std::string histName = calorimeterCalibrations->pulseHistName();

       TFile pulseFile(fileName.c_str());
       TH1F *pshape(0);
       if (pulseFile.IsOpen()) pshape = (TH1F*) pulseFile.Get(histName.c_str());
       
         if (!pshape) throw cet::exception("CATEGORY")<<"CaloPulseCache:: Hitsogram "<<histName.c_str()
	                                              <<" from file "<<fileName.c_str()<<" does not exist";        
        
         
	 double sumNorm(0), contentMax(0), posMax(0), posMin(0);
         for (int i=0; i<pshape->GetNbinsX();++i)
         {
            double content = pshape->GetBinContent(i);
            if (content < 1e-5) continue;
            
            cache_.push_back(content);
            sumNorm += content*pshape->GetBinWidth(i);

            if (content > contentMax) {contentMax=content; posMax = pshape->GetBinCenter(i);}
            if (posMin < 1e-5) posMin =  pshape->GetBinCenter(i);
         }
               
         for (auto& val : cache_) val /= sumNorm;   
         double funcMax = *std::max_element(cache_.begin(),cache_.end());

         cacheSize_      = cache_.size();
         deltaT_         = posMax - posMin; 
         step_           = pshape->GetBinWidth(1);
         factor_         = 1.0/funcMax;
          
       pulseFile.Close();   

   }
   
   double CaloPulseCache::evaluate(double x)
   {
       int idx = int( (x+deltaT_)/step_ );
       if (idx < 0 || idx > cacheSize_-2) return 0;     
       return (cache_[idx+1]-cache_[idx])/step_*(x+deltaT_ - idx*step_) + cache_[idx];        
   }

   
   

}
