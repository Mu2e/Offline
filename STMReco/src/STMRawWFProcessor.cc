// case checks: just one energy deposit, or two bins.
// refit the rising peak

#include "Offline/STMReco/inc/STMRawWFProcessor.hh"
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
  STMRawWFProcessor::STMRawWFProcessor(const Config& config) :
     STMWaveformProcessor(),
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
   void STMRawWFProcessor::initialize() {}


   //---------------------------
   void STMRawWFProcessor::reset()
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
   void STMRawWFProcessor::extract(const std::vector<double> &xInput, const std::vector<double> &yInput)
   {

       reset();

       xvec_ = xInput;
       yvec_ = yInput;


       //find location of potential peaks
       std::vector<unsigned> peakLocation;
       for (unsigned iu=windowPeak_;iu+windowPeak_<xvec_.size();++iu)
       {
          if (yvec_[iu] < minPeakAmplitude_) continue;
          if (std::max_element(yvec_.begin()+iu-windowPeak_,yvec_.begin()+iu+windowPeak_+1) == yvec_.begin()+iu) peakLocation.push_back(iu);
       }
       if (diagLevel_ > 1) std::cout<<"[STMRawWFProcessor] Peaks found : "<<peakLocation.size()<<std::endl;


       nPeaks_ = peakLocation.size();
       chi2_   = 0;
       ndf_    = 1;
       for (size_t i=0;i<peakLocation.size();++i)
       {
           resAmp_.push_back(scaleFactor_*yvec_[peakLocation[i]]);
           resAmpErr_.push_back(0);
           resTime_.push_back(xvec_[peakLocation[i]] - shiftTime_);
           resTimeErr_.push_back(0);
       }

   }


   //---------------------------------------
   void STMRawWFProcessor::plot(const std::string& pname) const
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
