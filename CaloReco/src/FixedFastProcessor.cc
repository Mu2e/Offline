// This is an hybrid signal extraction algorithm with pre-calculated shape function

// If there is one peak, the normalization for a given time can be found analytically. This analytical 
// solution is used as imput to find the starting time with a quasi-Netwon method. Since the Netwon 
// method can oscillate around the minimum, an additional logarithmic scaling is used to make sure 
// the minimum is found  

// If there are more than one peak, we use a generic gradient descent method, namely minuit.

#include "CaloReco/inc/FixedFastProcessor.hh"
#include "CaloReco/inc/FastFixedUtil.hh"
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
        
   FixedFastProcessor::FixedFastProcessor(const Config& config) :
      WaveformProcessor(),
      windowPeak_      (config.windowPeak()),
      minPeakAmplitude_(config.minPeakAmplitude()),
      psdThreshold_    (config.psdThreshold()),
      pulseLowBuffer_  (config.pulseLowBuffer()),
      pulseHighBuffer_ (config.pulseHighBuffer()),
      minDiffTime_     (config.minDiffTime()),
      refitLeadingEdge_(config.refitLeadingEdge()),
      shiftTime_       (config.shiftTime()),
      printLevel_      (config.fitPrintLevel()),
      fitStrategy_     (config.fitStrategy()),
      diagLevel_       (config.diagLevel()),
      nParFcn_(2),
      fmutil_(nParFcn_, minPeakAmplitude_, minDiffTime_, printLevel_, fitStrategy_, diagLevel_),
      nPeaks_(0),
      chi2_(999),
      resAmp_(),
      resAmpErr_(),
      resTime_(),
      resTimeErr_()
   {	          
       if (diagLevel_ > 2)
       {
	  art::ServiceHandle<art::TFileService> tfs;
	  art::TFileDirectory tfdir = tfs->mkdir("FastFixedDiag");
	  _hTime    = tfdir.make<TH1F>("hTime",    "time",                  100, 0., 2000);
	  _hTimeErr = tfdir.make<TH1F>("hTimeErr", "time error",            100, 0.,   10);
	  _hEner    = tfdir.make<TH1F>("hEner",    "Amplitude",             100, 0., 5000);
	  _hEnerErr = tfdir.make<TH1F>("hEnerErr", "Amplitude error",       100, 0.,  100);
	  _hChi2    = tfdir.make<TH1F>("hChi2",    "Chi2/ndf",              100, 0.,   20);
	  _hNpeak   = tfdir.make<TH1F>("hNpeak",   "Number of peak fitted",  10, 0.,   10);
	  _hchi2Amp = tfdir.make<TH2F>("hchi2Amp", "Amp vs chi2/ndf",        50, 0.,   20, 100, 0, 5000);
	  _hpsd     = tfdir.make<TH1F>("hpsd",     "PSD",                    50, 0.,   1);
       }
   }       


   void FixedFastProcessor::initialize()
   {
       fmutil_.initialize();
   }


   void FixedFastProcessor::extract(const std::vector<double>& xInput, const std::vector<double>& yInput)
   {       
       reset();

       std::vector<double> parInit;
       std::vector<unsigned> xindices;
       findPeak(xInput, yInput, parInit, xindices);      

       if (diagLevel_ > 2) _hNpeak->Fill(parInit.size()/nParFcn_);       
       if (parInit.size() > 99) return;         
                   
       std::vector<double> xfit, yfit;
       for (auto i : xindices) {xfit.push_back(xInput[i]); yfit.push_back(yInput[i]);}
       fmutil_.setXYVector(xfit,yfit);
        
       if (parInit.size()==nParFcn_) fmutil_.fitNewton(parInit[1]);
       else                          fmutil_.fitMinuit(parInit);      
       
       std::vector<double> parFinal    = fmutil_.sfpar();
       std::vector<double> parFinalErr = fmutil_.esfpar();
       
       if (refitLeadingEdge_) refitTime(xInput, yInput, parFinal);

       chi2_   = fmutil_.chi2();
       ndf_    = xindices.size() - parFinal.size();
       nPeaks_ = parFinal.size()/ nParFcn_;
     
       for (unsigned i=0;i<parFinal.size();++i)
       {            
           resAmp_.push_back(parFinal[nParFcn_*i]);
           resAmpErr_.push_back(parFinalErr[nParFcn_*i]);
           resTime_.push_back(parFinal[nParFcn_*i+1] - shiftTime_ );
           resTimeErr_.push_back(parFinalErr[nParFcn_*i+1]);          
       
           if (diagLevel_ > 2)
           {
               _hTime->Fill(parFinal[nParFcn_*i+1]);
               _hTimeErr->Fill(parFinalErr[nParFcn_*i+1]);
               _hEner->Fill(parFinal[nParFcn_*i]);
               _hEnerErr->Fill(parFinalErr[nParFcn_*i]);
           }           
       }

       if (diagLevel_ > 2)
       {
          _hChi2->Fill(chi2_/float(ndf_)); 
          for (auto amp : resAmp_) _hchi2Amp->Fill(chi2_/float(ndf_),amp);
       }          
       
       return;
   }


      
   //---------------------------
   void FixedFastProcessor::reset()
   {
       resAmp_.clear();
       resAmpErr_.clear();
       resTime_.clear();
       resTimeErr_.clear();
       
       nPeaks_  = 0;
       chi2_    = 999;
       ndf_     = 0;
   }
     
   
   //----------------------------------------------------------------------------------------------------------------------
   void FixedFastProcessor::findPeak(const std::vector<double>& xvec, const std::vector<double>& yvec,std::vector<double>& parInit, 
                                     std::vector<unsigned>& xindices)
   {        
        if (windowPeak_ < 1 || windowPeak_ > xvec.size()) return;
                
        std::vector<unsigned> peakLocation;
        for (auto i=windowPeak_;i<xvec.size()-windowPeak_;++i)
        {
             if (std::max_element(&yvec[i-windowPeak_],&yvec[i+windowPeak_+1]) != &yvec[i]) continue;
             if (yvec[i-1] < minPeakAmplitude_ || yvec[i] < minPeakAmplitude_ || yvec[i+1] < minPeakAmplitude_) continue;

             double currentAmplitudeX = fmutil_.eval_fcn(xvec[i],parInit);
             double xmaxEstimate = meanParabol(xvec[i],xvec[i-1],xvec[i+1],yvec[i],yvec[i-1],yvec[i+1]);
            
             parInit.push_back((yvec[i] - currentAmplitudeX)); 
             parInit.push_back(xmaxEstimate);
	     peakLocation.push_back(i);
        }

        unsigned nPeakPrimary = peakLocation.size();
        if (diagLevel_ > 1) std::cout<<"[FixedFastProcessor] Primary peaks found : "<<peakLocation.size()<<std::endl;   
	if (peakLocation.empty()) return; 

        std::vector<double> residual;
        for (unsigned i=0;i<xvec.size();++i) residual.push_back( yvec[i] - fmutil_.eval_fcn(xvec[i],parInit) ); 

        for (unsigned int i=windowPeak_;i<xvec.size()-windowPeak_;++i)
        {
	     if (std::max_element(&residual[i-windowPeak_],&residual[i+windowPeak_+1]) != &residual[i]) continue;             
             if (residual[i] < minPeakAmplitude_) continue;
             
             double psd = (yvec[i] > 1e-3) ? residual[i]/yvec[i] : 0;	     
             if (diagLevel_ > 2) _hpsd->Fill(psd);
             if (psd < psdThreshold_) continue;
             
             //int dmin = *std::min_element(peakLocationInit.begin(),peakLocationInit.end(),
             //                            [&i](const unsigned& lhs,const unsigned& rhs){return(std::abs(lhs-i) < std::abs(rhs-i));});
                          
             double currentAmplitudeX = fmutil_.eval_fcn(xvec[i],parInit);
             parInit.push_back(yvec[i] - currentAmplitudeX); 
             parInit.push_back(xvec[i]); 
             peakLocation.push_back(i);
        }
                
        if (diagLevel_ > 1) std::cout<<"[FixedFastProcessor] Secondary peaks found : "<<peakLocation.size()-nPeakPrimary<<std::endl;   
	
        buildXRange(peakLocation, xvec, xindices);        
        return;           
   }

   //------------------------------------------------------------
   double FixedFastProcessor::meanParabol(double x1, double x2, double x3, double y1, double y2, double y3)
   {
       double a = ((y1-y2)/(x1-x2)-(y1-y3)/(x1-x3))/(x2-x3);
       double b = (y1-y2)/(x1-x2) - a*(x1+x2);
       if (std::abs(a) < 1e-6) return (x1+x2+x3)/3.0;       
       return -b/2.0/a ;
   }

  
   //-------------------------------------------------------------
   void FixedFastProcessor::buildXRange(const std::vector<unsigned>& peakLoc, const std::vector<double>& xvec, std::vector<unsigned>& xindices)
   {      
	std::set<unsigned> tempX;
	for (unsigned ipeak : peakLoc)
	{             
             unsigned is = std::max(ipeak-pulseLowBuffer_,0u);
	     unsigned ie = std::min(ipeak+pulseHighBuffer_,unsigned(xvec.size()));	     	     
	     for (unsigned ip=is; ip<ie; ++ip) tempX.insert(ip); 
	}

	for (auto i : tempX) xindices.push_back(i);  
   }
     


   
   void FixedFastProcessor::refitTime(const std::vector<double>& xvec, const std::vector<double>& yvec, std::vector<double>& parFinal)
   {      
       std::vector<double> ywork(yvec);
       for (unsigned ip=0; ip<parFinal.size(); ip+=nParFcn_)
       {
            //find rising edge of peak and refit it
            int imax = int((parFinal[ip+1]-xvec[0])/(xvec[1]-xvec[0]));
            int imin(imax);
            while (imin>0){ if (yvec[imin]/yvec[imax] < 0.1) break; --imin; }

            std::vector<double> xfit, yfit;
            for (int i = imin; i <= imax; ++i) {xfit.push_back(xvec[i]); yfit.push_back(yvec[i]);}
            fmutil_.setXYVector(xfit,yfit);

            parFinal[ip+1] = fmutil_.refineMin(parFinal[ip+1],0.01);

            //remove the peak contribution         
            for (unsigned i=0;i<xvec.size();++i) ywork[i] -= fmutil_.eval_logn(xvec[i],ip,parFinal);        
       }
   }
    
   void FixedFastProcessor::plot(std::string pname)
   {
       std::vector<double> param;
       for (unsigned i=0;i<resAmp_.size();++i)
       {
          param.push_back(resAmp_[i]);
          param.push_back(resTime_[i]);
       }
       
       fmutil_.plot(pname,param);
   }
     
}
