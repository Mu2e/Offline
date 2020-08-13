#include "CaloReco/inc/TemplateProcessor.hh"
#include "CaloReco/inc/TemplateUtil.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"

#include "TFile.h"
#include "TH2.h"

#include <algorithm>
#include <string> 
#include <iostream>
#include <vector>

namespace mu2e {

   TemplateProcessor::TemplateProcessor(const Config& config) :
      WaveformProcessor(),
      windowPeak_      (config.windowPeak()),
      minPeakAmplitude_(config.minPeakAmplitude()),
      numNoiseBins_    (config.numNoiseBins()),
      minDeltaPeakBin_ (config.minDeltaPeakBin()),
      doSecondaryPeak_ (config.doSecondaryPeak()),
      psdThreshold_    (config.psdThreshold()),
      refitLeadingEdge_(config.refitLeadingEdge()),
      timeCorr_        (config.timeCorr()),
      diagLevel_       (config.diagLevel()),
      fmutil_          (minPeakAmplitude_,config.digiSampling(),config.minDTPeaks(), config.fitPrintLevel()),
      chi2_            (999.),
      ndf_             (-1),
      shiftPar_        (),
      resAmp_          (),
      resAmpErr_       (),
      resTime_         (),
      resTimeErr_      ()
   {	                 
       if (diagLevel_ > 1)
       {
	  art::ServiceHandle<art::TFileService> tfs;
	  art::TFileDirectory tfdir = tfs->mkdir("FastFixedDiag");
	  _hTime    = tfdir.make<TH1F>("hTime",    "time",                  100, 0., 2000);
	  _hTimeErr = tfdir.make<TH1F>("hTimeErr", "time error",            100, 0.,   10);
	  _hEner    = tfdir.make<TH1F>("hAmp",     "Amplitude",             200, 0.,  200);
	  _hEnerErr = tfdir.make<TH1F>("hAmpErr",  "Amplitude error",       100, 0.,  100);
	  _hNpeak   = tfdir.make<TH1F>("hNpeak",   "Number of peak fitted",  10, 0.,   10);
	  _hChi2    = tfdir.make<TH1F>("hChi2",    "Chi2/ndf",              100, 0.,   20);
	  _hchi2Amp = tfdir.make<TH2F>("hchi2Amp", "Amp vs chi2/ndf",        50, 0.,   20, 100, 0, 2000);
       }
   }       


   
   
   void TemplateProcessor::initialize()
   {
       fmutil_.initialize();
       fmutil_.calcTimeCorr(shiftPar_);
   }


   void TemplateProcessor::extract(const std::vector<double>& xInput, const std::vector<double>& yInput)
   {       
       if (diagLevel_>2) std::cout<<"TemplateProcessor start"<<std::endl;
       reset();
       fmutil_.setXYVector(xInput, yInput);
       
       std::vector<unsigned> peakLocation;
       setPrimaryPeakPar(xInput, yInput, peakLocation);
       auto nFirstPeaks = peakLocation.size();
       if (nFirstPeaks==0) return;
       if (diagLevel_>2) std::cout<<"Found "<<nFirstPeaks<<" peaks"<<std::endl;
       fmutil_.setXYVector(xInput, yInput);
       fmutil_.fit();
      
       setSecondaryPeakPar(xInput, yInput, peakLocation);
       if (doSecondaryPeak_ && peakLocation.size() > nFirstPeaks)
       {
          if (diagLevel_>2) std::cout<<"Found "<<peakLocation.size()-nFirstPeaks<<" secondary peaks"<<std::endl;
          fmutil_.setXYVector(xInput, yInput);
          fmutil_.fit();
       }

       if (diagLevel_>2) std::cout<<"Fit done"<<std::endl;
       
       if (refitLeadingEdge_) fmutil_.refitEdge();

       for (unsigned i=fmutil_.nParBkg(); i<fmutil_.par().size(); i += fmutil_.nParFcn())
       {            
	   resAmp_.push_back(fmutil_.par()[i]);                          
	   resAmpErr_.push_back(fmutil_.parErr()[i]);
	   resTime_.push_back(fmutil_.fromPeakToT0(fmutil_.par()[i+1])); 
	   resTimeErr_.push_back(fmutil_.parErr()[i+1]);
       }
       chi2_ = fmutil_.chi2();
       ndf_  = yInput.size() - fmutil_.par().size();
              
       if (diagLevel_ > 1)
       {
          double chindf = (ndf_ > 0) ? chi2_/float(ndf_) : 999;
          _hNpeak->Fill(resAmp_.size());  
          _hChi2->Fill(chindf); 
          for (auto val : resAmp_)     _hchi2Amp->Fill(chindf,val);
          for (auto val : resAmp_)     _hEner->Fill(val);
          for (auto val : resAmpErr_)  _hEnerErr->Fill(val);
          for (auto val : resTime_)    _hTime->Fill(val);
          for (auto val : resTimeErr_) _hTimeErr->Fill(val);
       }          
   }


      
   //----------------------------------------------------------------------------------------------------------------------------------
   void TemplateProcessor::reset() {resAmp_.clear(); resAmpErr_.clear(); resTime_.clear(); resTimeErr_.clear();ndf_ = 0;chi2_ = 999; fmutil_.reset();}

   
   //---------------------------------------------------------------------------------------------------------------------------------------------
   void TemplateProcessor::setPrimaryPeakPar(const std::vector<double>& xvec, const std::vector<double>& yvec, std::vector<unsigned>& peakLocation)
   {        
        if (windowPeak_ > xvec.size()) return;
                
        std::vector<double> parInit,ywork(yvec);
        
 	//estimate noise, subtract baseline and get peaks
	float noise= std::accumulate(yvec.begin(),yvec.begin()+numNoiseBins_,0)/float(numNoiseBins_);
        parInit.push_back(noise);
        for (auto& val : ywork) val -= noise;        
       
	//first find "large" peaks, apply crude noise filtering (adjacent bins should be above mininum value) and remove small peaks
        for (auto i=windowPeak_;i<int(xvec.size())-windowPeak_;++i)
        {
            if (std::max_element(&ywork[i-windowPeak_],&ywork[i+windowPeak_+1]) != &ywork[i]) continue;
            if (ywork[i-1] < minPeakAmplitude_ || ywork[i] < minPeakAmplitude_ || ywork[i+1]< minPeakAmplitude_) continue;
            if (checkPeakDist(i,peakLocation)) continue;
            //if (!peakLocation.empty() && ywork[i] < 2*minPeakAmplitude_) continue;

            double xmaxEstimate = estimatePeakTime(xvec[i-1],xvec[i],xvec[i+1],ywork[i-1],ywork[i],ywork[i+1]);
	    parInit.push_back(ywork[i]); 
            parInit.push_back(xmaxEstimate);	     
	    peakLocation.push_back(i);
            
            if (diagLevel_ > 2) std::cout<<"[TemplateProcessor] Found primary peak "<<xmaxEstimate<<std::endl;

	    fmutil_.setPar(parInit);
	    for (unsigned j=0;j<xvec.size();++j) ywork[j] -= fmutil_.eval_logn(xvec[j],parInit.size()-fmutil_.nParFcn());		     
        }        
        
        if (peakLocation.empty()) 
        {
            unsigned ipeak(0);
            while (ipeak < ywork.size()) {if (ywork[ipeak] >= minPeakAmplitude_) break; ++ipeak;}
            parInit.push_back(ywork[ipeak]); 
            parInit.push_back(xvec[ipeak]);	     
	    peakLocation.push_back(ipeak);
            fmutil_.setPar(parInit);
        }
    }


    //--------------------------------------------------------------------------------------------------
    void TemplateProcessor::setSecondaryPeakPar(const std::vector<double>& xvec, const std::vector<double>& yvec, std::vector<unsigned>& peakLocation)
    {               
        if (windowPeak_ > xvec.size()) return;
	
	std::vector<double> parInit = fmutil_.par();

	std::vector<double> ywork(yvec);
	for (unsigned j=0;j<xvec.size();++j) ywork[j] -= fmutil_.eval_fcn(xvec[j]);	
	
        for (auto i=windowPeak_;i<xvec.size()-windowPeak_;++i)
        {
            if (std::max_element(&ywork[i-windowPeak_],&ywork[i+windowPeak_+1]) != &ywork[i]) continue;
            if (ywork[i-1] < minPeakAmplitude_ || ywork[i] < minPeakAmplitude_ || ywork[i+1] < minPeakAmplitude_) continue;
            if (checkPeakDist(i,peakLocation)) continue;            
            if (yvec[i]>0 && ywork[i]/yvec[i] < psdThreshold_) continue;  // pulseShape Discrimintion  

            double xmaxEstimate = estimatePeakTime(xvec[i-1],xvec[i],xvec[i+1],ywork[i-1],ywork[i],ywork[i+1]);
            parInit.push_back(ywork[i]); 
            parInit.push_back(xmaxEstimate);	     
	    peakLocation.push_back(i);

	    fmutil_.setPar(parInit);
	    for (unsigned j=0;j<xvec.size();++j) ywork[j] -= fmutil_.eval_logn(xvec[j],parInit.size()-fmutil_.nParFcn());	
        }
   }


   //------------------------------------------------------------
   double TemplateProcessor::estimatePeakTime(double x1, double x2, double x3, double y1, double y2, double y3)
   {       
       return (x1*y1+x2*y2+x3*y3)/(y1+y2+y3)+0*(shiftPar_[0]+shiftPar_[1]*y3/y2)*(x2-x1);
   }

   //------------------------------------------------------------------------------------------
   bool TemplateProcessor::checkPeakDist(unsigned i, const std::vector<unsigned>& peakLocation)
   {
      for (auto pl : peakLocation) 
	 if ( ((i>pl) ? i-pl : pl-i) < minDeltaPeakBin_) return true;
      
      return false;	 
   }
   

   //---------------------------------------------------------------------------------------------------------------------------------------
   void TemplateProcessor::plot(const std::string& name) const {fmutil_.plotFit(name);}
  
   //---------------------------------------------------------------------------------------------------------------------------------------
   void TemplateProcessor::dump(const std::string& name, const std::vector<double>& val) const 
   {
      std::cout<<name<<std::endl;
      for (auto h : val) std::cout<<h<<" ";
      std::cout<<std::endl;     
   }
}
