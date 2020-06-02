// This is an signal extraction method based on the waveform template. 
//
// Each peak in the waveform is described by two parameters: amplitide and peak time
// For a single peak, the amplitude can be found analytically for a given start time, and a 
// quasi-Netwon method can be used to fit the waveform. 
// If there are more than one peak, we use a generic gradient descent method, namely minuit.
//
// There is an additional option to refit the leding edge of the first peak to improve 
// timing accuracy
//
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
      refitLeadingEdge_(config.refitLeadingEdge()),
      timeCorr_        (config.timeCorr()),
      diagLevel_       (config.diagLevel()),
      fmutil_          (minPeakAmplitude_,config.digiSampling(),config.pulseIntegralSteps()),
      chi2_            (999.),
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
	  _hEner    = tfdir.make<TH1F>("hEner",    "Amplitude",             100, 0., 5000);
	  _hEnerErr = tfdir.make<TH1F>("hEnerErr", "Amplitude error",       100, 0.,  100);
	  _hNpeak   = tfdir.make<TH1F>("hNpeak",   "Number of peak fitted",  10, 0.,   10);
	  _hChi2    = tfdir.make<TH1F>("hChi2",    "Chi2/ndf",              100, 0.,   20);
	  _hchi2Amp = tfdir.make<TH2F>("hchi2Amp", "Amp vs chi2/ndf",        50, 0.,   20, 100, 0, 5000);
       }
   }       


   
   
   void TemplateProcessor::initialize()
   {
       fmutil_.initialize();
   }


   void TemplateProcessor::extract(const std::vector<double>& xInput, const std::vector<double>& yInput)
   {       
       reset();

       std::vector<double> parInit;
       std::vector<unsigned> xindices;
       findPeak(xInput, yInput, parInit, xindices);      

       if (parInit.size() > 99) return;         
                   
       std::vector<double> xfit, yfit;
       for (auto i : xindices) {xfit.push_back(xInput[i]); yfit.push_back(yInput[i]);}
       fmutil_.setXYVector(xfit,yfit);
       fmutil_.setPar(parInit);
       unsigned nParFcn(fmutil_.nParFcn());
        
       if (parInit.size()==nParFcn) fmutil_.fitNewton();
       else                         fmutil_.fitMinuit();      
              
       if (refitLeadingEdge_) refitLeadingEdge(xInput, yInput);

       std::vector<double> parFinal    = fmutil_.par();
       std::vector<double> parFinalErr = fmutil_.parErr();

       chi2_   = fmutil_.chi2();
       ndf_    = xindices.size() - parFinal.size();
     
       for (unsigned i=0;i<parFinal.size();++i)
       {            
           resAmp_.push_back(parFinal[nParFcn*i]);
           resAmpErr_.push_back(parFinalErr[nParFcn*i]);
           resTime_.push_back(fmutil_.fromPeakToT0(parFinal[nParFcn*i+1]) - timeCorr_ );
           resTimeErr_.push_back(parFinalErr[nParFcn*i+1]);          
           
           if (diagLevel_ > 2)
           {
               _hTime->Fill(parFinal[nParFcn*i+1]);
               _hTimeErr->Fill(parFinalErr[nParFcn*i+1]);
               _hEner->Fill(parFinal[nParFcn*i]);
               _hEnerErr->Fill(parFinalErr[nParFcn*i]);
           }           
       }
       
       if (diagLevel_ > 2)
       {
          _hNpeak->Fill(fmutil_.nPeaks());  
          _hChi2->Fill(chi2_/float(ndf_)); 
          for (auto amp : resAmp_) _hchi2Amp->Fill(chi2_/float(ndf_),amp);
       }          
   }


      
   //----------------------------------------------------------------------------------------------------------------------------------
   void TemplateProcessor::reset() {resAmp_.clear(); resAmpErr_.clear(); resTime_.clear(); resTimeErr_.clear();ndf_ = 0;chi2_ = 999;}


   //--------------------------------------------------------------------------------------------------
   void TemplateProcessor::findPeak(const std::vector<double>& xvec, const std::vector<double>& yvec,
                                    std::vector<double>& parInit, std::vector<unsigned>& xindices)
   {        
        if (windowPeak_ < 1 || windowPeak_ > xvec.size()) return;
	fmutil_.setPar(parInit);
                
        std::vector<double> ywork(yvec);
	std::vector<unsigned> peakLocation;
        for (auto i=windowPeak_;i<xvec.size()-windowPeak_;++i)
        {
             if (std::max_element(&ywork[i-windowPeak_],&ywork[i+windowPeak_+1]) != &ywork[i]) continue;
             if (ywork[i-1] < minPeakAmplitude_ || ywork[i] < minPeakAmplitude_ || ywork[i+1] < minPeakAmplitude_) continue;

             double xmaxEstimate = meanParabol(xvec[i],xvec[i-1],xvec[i+1],ywork[i],ywork[i-1],ywork[i+1]);
            
             parInit.push_back(ywork[i]); 
             parInit.push_back(xmaxEstimate);	     
	     peakLocation.push_back(i);
	     
	     fmutil_.setPar(parInit);
	     for (unsigned j=0;j<xvec.size();++j) ywork[j] -= fmutil_.eval_logn(xvec[j],parInit.size()-fmutil_.nParFcn());	
	     
	     i += minDeltaPeakBin_; // make sure we won;t find any other peak too close to this one.
        }
	
	//second pass on the residuals to catch secondary peaks
        for (auto i=windowPeak_;i<xvec.size()-windowPeak_;++i)
        {
             if (std::max_element(&ywork[i-windowPeak_],&ywork[i+windowPeak_+1]) != &ywork[i]) continue;
             if (ywork[i-1] < minPeakAmplitude_ || ywork[i] < minPeakAmplitude_ || ywork[i+1] < minPeakAmplitude_) continue;
             
	     double psd = (yvec[i] > 1e-3) ? ywork[i]/yvec[i] : 0;	     
             if (psd < 0.2) continue;

             double xmaxEstimate = meanParabol(xvec[i],xvec[i-1],xvec[i+1],ywork[i],ywork[i-1],ywork[i+1]);
            
             parInit.push_back(ywork[i]); 
             parInit.push_back(xmaxEstimate);	     
	     peakLocation.push_back(i);
	     
	     fmutil_.setPar(parInit);
	     for (unsigned j=0;j<xvec.size();++j) ywork[j] -= fmutil_.eval_logn(xvec[j],parInit.size()-fmutil_.nParFcn());	
	     
	     i += minDeltaPeakBin_; // make sure we won;t find any other peak too close to this one.
        }
	
        buildXRange(peakLocation, xvec, xindices);        
        return;           
   }


   //------------------------------------------------------------
   double TemplateProcessor::meanParabol(double x1, double x2, double x3, double y1, double y2, double y3)
   {
       double a = ((y1-y2)/(x1-x2)-(y1-y3)/(x1-x3))/(x2-x3);
       double b = (y1-y2)/(x1-x2) - a*(x1+x2);
       if (std::abs(a) < 1e-6) return (x1+x2+x3)/3.0;       
       return -b/2.0/a ;
   }

   
   //-------------------------------------------------------------
   void TemplateProcessor::buildXRange(const std::vector<unsigned>& peakLoc, const std::vector<double>& xvec, std::vector<unsigned>& xindices)
   {      
	std::set<unsigned> tempX;
	for (unsigned ipeak : peakLoc)
	{             
             unsigned is = (ipeak > pulseLowBuffer_) ? ipeak-pulseLowBuffer_ : 0u;
	     unsigned ie = std::min(ipeak+pulseHighBuffer_,unsigned(xvec.size()));	     	     
	     for (unsigned ip=is; ip<ie; ++ip) tempX.insert(ip); 
	}

	for (auto ix : tempX) xindices.push_back(ix);  
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
	fmutil_.refineMin(0.01);
        fmutil_.setXYVector(xvec,yvec); //reset original vectors
   }
      
   
   //---------------------------------------------------------------------------------------------------------------------------------------
   void TemplateProcessor::plot(std::string pname) const {fmutil_.plotFit(pname);}

}
