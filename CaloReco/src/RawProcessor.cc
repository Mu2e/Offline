// case checks: just one energy deposit, or two bins. 
// refit the rising peak

#include "CaloReco/inc/RawProcessor.hh"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"

#include "TH1.h"
#include "TCanvas.h"

#include <algorithm>
#include <string> 
#include <iostream>
#include <vector>



namespace mu2e {
        
  //-----------------------------------------------------------------------------
  RawProcessor::RawProcessor(const Config& config) :
     WaveformProcessor(),
     windowPeak_      (config.windowPeak()),
     minPeakAmplitude_(config.minPeakAmplitude()),
     shiftTime_       (config.shiftTime()),
     scaleFactor_     (config.scaleFactor()),
     diagLevel_       (config.diagLevel()),           
     nPeaks_(0),
     chi2_(999),
     res_(),
     resAmp_(),
     resAmpErr_(),
     resTime_(),
     resTimeErr_(),
     xvec_(),
     yvec_()
   {}       


   //---------------------------
   void RawProcessor::initialize() {}


   //---------------------------
   void RawProcessor::reset()
   {
      xvec_.clear();
      yvec_.clear();
      res_.clear();
      resAmp_.clear();
      resAmpErr_.clear();
      resTime_.clear();
      resTimeErr_.clear();

      nPeaks_  = 0;
      chi2_    = 999;
      ndf_     = 0;
   }



   //------------------------------------------------------------------------------------------
   void RawProcessor::extract(const std::vector<double> &xInput, const std::vector<double> &yInput)
   {

      reset();

      xvec_ = xInput;
      yvec_ = yInput;


      //find location of potential peaks
      std::vector<int> peakLocation;

      for (unsigned int iu=windowPeak_;iu<xvec_.size()-windowPeak_;++iu)
      {
	 auto maxp = std::max_element(&yvec_[iu-windowPeak_],&yvec_[iu+windowPeak_+1]);
	 if (maxp == &yvec_[iu] && yvec_[iu] > minPeakAmplitude_) peakLocation.push_back(iu);
      }
      int nPeak = peakLocation.size();

      if (diagLevel_ > 1) std::cout<<"[RawProcessor] Peaks found : "<<nPeak<<std::endl;


       nPeaks_ = nPeak;
       chi2_   = 0;
       ndf_    = 1;
       for (int i=0;i<nPeak;++i)
       { 
           resAmp_.push_back(scaleFactor_*yvec_[peakLocation[i]]);
	   resAmpErr_.push_back(0);
           resTime_.push_back(xvec_[peakLocation[i]] - shiftTime_);
           resTimeErr_.push_back(0);	  
       }

   }

   
   //---------------------------------------
   void RawProcessor::plot(const std::string& pname) const
   {
      TH1F h("test","Amplitude vs time",xvec_.size(),xvec_.front()-2.5,xvec_.back()+2.5);
      h.GetXaxis()->SetTitle("Time (ns)");
      h.GetYaxis()->SetTitle("Amplitude");
      for (unsigned int i=0;i<xvec_.size();++i) h.SetBinContent(i+1,yvec_[i]);

      TCanvas c1("c1","c1");
      h.Draw();
      c1.SaveAs(pname.c_str());
   }
}
