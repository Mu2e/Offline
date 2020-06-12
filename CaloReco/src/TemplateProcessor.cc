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
      pulseLowBuffer_  (config.pulseLowBuffer()),
      pulseHighBuffer_ (config.pulseHighBuffer()),
      minDeltaPeakBin_ (config.minDeltaPeakBin()),
      doSecondaryPeak_ (config.doSecondaryPeak()),
      refitLeadingEdge_(config.refitLeadingEdge()),
      timeCorr_        (config.timeCorr()),
      diagLevel_       (config.diagLevel()),
      fmutil_          (minPeakAmplitude_,config.digiSampling()),
      chi2_            (999.),
      ndf_             (-1),
      shiftPar_        (),
      resAmp_          (),
      resAmpErr_       (),
      resTime_         (),
      resTimeErr_      ()
   {	                 
       if (diagLevel_ > 2)
       {
	  art::ServiceHandle<art::TFileService> tfs;
	  art::TFileDirectory tfdir = tfs->mkdir("FastFixedDiag");
	  _hTime    = tfdir.make<TH1F>("hTime",    "time",                  100, 0., 2000);
	  _hTimeErr = tfdir.make<TH1F>("hTimeErr", "time error",            100, 0.,   10);
	  _hEner    = tfdir.make<TH1F>("hEner",    "Amplitude",             100, 0., 2000);
	  _hEnerErr = tfdir.make<TH1F>("hEnerErr", "Amplitude error",       100, 0.,  100);
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

       std::vector<unsigned> peakLocation;
       std::vector<double> xfit, yfit;

       setPrimaryPeakPar(xInput, yInput, peakLocation);
       auto nFirstPeaks = peakLocation.size();
       if (nFirstPeaks==0) return;
       buildXRange(peakLocation, xInput, yInput, xfit,yfit);              
       fmutil_.setXYVector(xfit,yfit);
       fmutil_.fit();
      
       setSecondaryPeakPar(xInput, yInput, peakLocation);
       if (doSecondaryPeak_ && peakLocation.size() > nFirstPeaks)
       {
	  buildXRange(peakLocation, xInput, yInput, xfit, yfit);
          fmutil_.setXYVector(xfit,yfit);
          fmutil_.fit();
       }

       if (diagLevel_>2) std::cout<<"Fit done"<<std::endl;
       
       if (refitLeadingEdge_ && fmutil_.par().size()>1) refitLeadingEdge(xInput, yInput);

       if (fmutil_.par().size())
       {
           chi2_ = fmutil_.chi2();
           ndf_  = yfit.size() - fmutil_.par().size();          
           for (unsigned i=0;i<fmutil_.par().size();++i)
           {            
	       if (i%2==0) {resAmp_.push_back(fmutil_.par()[i]);                        resAmpErr_.push_back(fmutil_.parErr()[i]);}
	       else        {resTime_.push_back(fmutil_.fromPeakToT0(fmutil_.par()[i])); resTimeErr_.push_back(fmutil_.parErr()[i]);}
           }
       } 
       
       
       if (diagLevel_ > 2)
       {
          double chindf = (ndf_ > 0) ? chi2_/float(ndf_) : 999;
          _hNpeak->Fill(fmutil_.nPeaks());  
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
        for (auto i=windowPeak_;i<int(xvec.size())-windowPeak_;++i)
        {
            if (std::max_element(&ywork[i-windowPeak_],&ywork[i+windowPeak_+1]) != &ywork[i]) continue;
            if (ywork[i-1] < minPeakAmplitude_ || ywork[i] < minPeakAmplitude_ || ywork[i+1] < minPeakAmplitude_) continue;
            if (checkPeakDist(i,peakLocation)) continue;

            double xmaxEstimate = estimatePeakTime(xvec[i-1],xvec[i],xvec[i+1],ywork[i-1],ywork[i],ywork[i+1]);
	    parInit.push_back(ywork[i]); 
            parInit.push_back(xmaxEstimate);	     
	    peakLocation.push_back(i);
           
            if (diagLevel_>2) std::cout<<"TemplateProcessor found primary peak "<<i<<" "<<xmaxEstimate<<std::endl;

	    fmutil_.setPar(parInit);
	    for (unsigned j=0;j<xvec.size();++j) ywork[j] -= fmutil_.eval_logn(xvec[j],parInit.size()-fmutil_.nParFcn());		     
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
             //if (/yvec[i]>0 && ywork[i]/yvec[i] < 0.2) continue;   //cheap PSD

             double xmaxEstimate = estimatePeakTime(xvec[i-1],xvec[i],xvec[i+1],ywork[i-1],ywork[i],ywork[i+1]);
             parInit.push_back(ywork[i]); 
             parInit.push_back(xmaxEstimate);	     
	     peakLocation.push_back(i);
	     
             if (diagLevel_>2) std::cout<<"TemplateProcessor found secondary peak "<<i<<" "<<xmaxEstimate<<std::endl;
	     
             fmutil_.setPar(parInit);
	     for (unsigned j=0;j<xvec.size();++j) ywork[j] -= fmutil_.eval_logn(xvec[j],parInit.size()-fmutil_.nParFcn());	
        }
   }


   //--------------------------------------------------------------------------------------------------
   void TemplateProcessor::setPeakPar(const std::vector<double>& xvec, const std::vector<double>& yvec, std::vector<unsigned>& peakLocation)
   {        
       if (windowPeak_ > xvec.size()) return;

       std::vector<double> parInit, ywork(yvec);
       for (auto i=windowPeak_;i<xvec.size()-windowPeak_;++i)
       {
            if (std::max_element(&ywork[i-windowPeak_],&ywork[i+windowPeak_+1]) != &ywork[i]) continue;
            if (ywork[i-1] < minPeakAmplitude_ || ywork[i] < minPeakAmplitude_ || ywork[i+1] < minPeakAmplitude_) continue;
            if (checkPeakDist(i,peakLocation)) continue;

            double xmaxEstimate = estimatePeakTime(xvec[i-1],xvec[i],xvec[i+1],ywork[i-1],ywork[i],ywork[i+1]);
            parInit.push_back(ywork[i]); 
            parInit.push_back(xmaxEstimate);	     
	    peakLocation.push_back(i);

	    fmutil_.setPar(parInit);
	    for (unsigned j=0;j<xvec.size();++j) ywork[j] -= fmutil_.eval_logn(xvec[j],parInit.size()-fmutil_.nParFcn());		     
       }

       //second pass on the residuals to catch secondary peaks
       if (!doSecondaryPeak_) return;

       for (auto i=windowPeak_;i<xvec.size()-windowPeak_;++i)
       {
            if (std::max_element(&ywork[i-windowPeak_],&ywork[i+windowPeak_+1]) != &ywork[i]) continue;
            if (ywork[i-1] < minPeakAmplitude_ || ywork[i] < minPeakAmplitude_ || ywork[i+1] < minPeakAmplitude_) continue;
            if (checkPeakDist(i,peakLocation)) continue;
            //if (/yvec[i]>0 && ywork[i]/yvec[i] < 0.2) continue;   //cheap PSD

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
   
   //-------------------------------------------------------------
   void TemplateProcessor::buildXRange(const std::vector<unsigned>& peakLoc, const std::vector<double>& xvec, const std::vector<double>& yvec, 
                                       std::vector<double>& xfit,std::vector<double>& yfit)
   {      
       std::set<unsigned> tempX;
       for (unsigned ipeak : peakLoc)
       {               
	    unsigned is(ipeak),ie(ipeak);
	    while (is>0 && is+pulseLowBuffer_ > ipeak)               {if (yvec[is]<minPeakAmplitude_/2) break; --is;}
	    while (ie<xvec.size() && ie < ipeak+pulseHighBuffer_ )  {if (yvec[ie]<minPeakAmplitude_/2) break; ++ie;}	     
	    for (unsigned ip=is; ip<ie; ++ip) tempX.insert(ip); 
       }

       xfit.clear(); yfit.clear();
       for (auto ix : tempX) {xfit.push_back(xvec[ix]); yfit.push_back(yvec[ix]);}  
   }
   
   //---------------------------------------------------------------------------------------------------------------------------------------
   void TemplateProcessor::refitLeadingEdge(const std::vector<double>& xvec, const std::vector<double>& yvec)
   {      
       int imax = int((fmutil_.par().at(1)-xvec[0])/(xvec[1]-xvec[0]));
       int ibLow(imax),ibUp(imax);
       for (int i=imax;i>0;--i)
       {
	  if (yvec[i]/yvec[imax] >0.9) ibUp  = i;
	  if (yvec[i]/yvec[imax] >0.1) ibLow = i;
       }

       if (ibUp-ibLow < 4) return;

       std::vector<double> xfit, yfit;
       for (int i = ibLow; i <= ibUp; ++i) {xfit.push_back(xvec[i]); yfit.push_back(yvec[i]);}

       fmutil_.setXYVector(xfit,yfit);
       fmutil_.refitMin();
       fmutil_.setXYVector(xvec,yvec); //reset original vectors
   }


   //---------------------------------------------------------------------------------------------------------------------------------------
   void TemplateProcessor::plot(std::string pname) const {fmutil_.plotFit(pname);}

}
