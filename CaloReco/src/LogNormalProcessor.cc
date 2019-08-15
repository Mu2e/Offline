
#include "CaloReco/inc/LogNormalProcessor.hh"
#include "art_root_io/TFileDirectory.h" 
#include "art_root_io/TFileService.h"
#include "ConditionsService/inc/ConditionsHandle.hh"


#include "TMinuit.h"
#include "TH1.h"
#include "TF1.h"
#include "TCanvas.h"

#include <algorithm>
#include <string> 
#include <iostream>
#include <vector>




//anonymous namespace containing the data structures and functions required by Minuit
namespace {

   std::vector<double> xvec_,yvec_;
   int                 nparTot_, nparFcn_, nparFix_;
   unsigned int        istart_,iend_;

   
   //-----------------------------------------------
   double logn(double x, double *par)
   {
	double logterms0 = 1.175*par[2]+sqrt(1.0+1.357225*par[2]*par[2]);
	double s0        = log(logterms0)/1.175;
	double Aterm     = par[2]/(2.506628274631*par[3]*s0);
	double logterm   = 1.0-(par[2]/par[3])*(x-par[1]);  

	if (logterm<0) logterm = 1e-6;

	double expterm = log(logterm)/s0;
	return par[0]*Aterm*exp(-0.5*expterm*expterm)/1.07755263; //1.07755263 for normalization      
   }

   
   //------------------------------------------------------
   double fitfunction(double x, double *par)
   {   
       double result(0);
       for (int i=0;i<nparTot_;i+=nparFcn_) result += logn(x,&par[i]);
       return result;
   }
   double fitfunction2(double* x, double *par) {return fitfunction(x[0],par);}
   
   //------------------------------------------------------------------------
   void myfcn(int& npar, double* , double &f, double *par, int)
   {   
       f=0;
       for (unsigned int i=istart_;i<iend_;++i)
       {     
           double x = xvec_[i];
           double y = yvec_[i];
	   double val = fitfunction(x, par);
           if (val > 1e-5) {f += (y-val)*(y-val)/val;}        
       }
   }
      
}





namespace mu2e {

  //-----------------------------------------------------------------------------
   LogNormalProcessor::LogNormalProcessor(fhicl::ParameterSet const& PSet) :

      WaveformProcessor(PSet),
      windowPeak_        (PSet.get<int>   ("windowPeak")),
      minPeakAmplitude_  (PSet.get<double>("minPeakAmplitude")),
      fixShapeSig_       (PSet.get<bool>  ("fixShapeSig")),
      psdThreshold_      (PSet.get<double>("psdThreshold")),
      pulseHighBuffer_   (PSet.get<int>   ("pulseHighBuffer")),
      timeFraction_      (PSet.get<double>("timeFraction")),
      shiftTime_         (PSet.get<double>("shiftTime")),
      printLevel_        (PSet.get<int>   ("fitPrintLevel",-1)),
      fitStrategy_       (PSet.get<int>   ("fitStrategy",1)),
      diagLevel_         (PSet.get<int>   ("diagLevel",0)),
      nPeaks_(0),
      chi2_(999),
      res_(),
      resAmp_(),
      resAmpErr_(),
      resTime_(),
      resTimeErr_()
   {
       nparTot_ = 0;
       nparFcn_ = 4;
       nparFix_ = fixShapeSig_ ? 2 : 0;

       double dummy[4]={1.0,1.0,-0.384473,11.1171};
       peakFactor_ = 1.0/logn(dummy[1],dummy);


       if (diagLevel_ > 2)
       {
          art::ServiceHandle<art::TFileService> tfs;
          art::TFileDirectory tfdir = tfs->mkdir("LogNormalDiag");
          _hTime     = tfdir.make<TH1F>("hTime",    "time",                  100, 0., 2000);
          _hTimeErr  = tfdir.make<TH1F>("hTimeErr", "time error",            100, 0.,   10);
          _hEner     = tfdir.make<TH1F>("hEner",    "Amplitude",             100, 0., 5000);
          _hEnerErr  = tfdir.make<TH1F>("hEnerErr", "Amplitude error",       100, 0.,  100);
          _hChi2     = tfdir.make<TH1F>("hChi2",    "Chi2/ndf",              100, 0.,   20);
          _hNpeak    = tfdir.make<TH1F>("hNpeak",   "Number of peak fitted",  10, 0.,   10);
          _hDelta    = tfdir.make<TH1F>("hDelta",   "Delta t",               100, -50.,   50);
       }
   }       

   //---------------------------
   void LogNormalProcessor::initialize()
   {
   }

   //---------------------------
   void LogNormalProcessor::reset()
   {
      xvec_.clear();
      yvec_.clear();
      res_.clear();
      resAmp_.clear();
      resAmpErr_.clear();
      resTime_.clear();
      resTimeErr_.clear();
      
      nparTot_ = 0;
      nPeaks_  = 0;
      istart_  = 0;
      iend_    = 0;
      chi2_    = 999;
      ndf_     = 0;
   }



   //------------------------------------------------------------------------------------------
   void LogNormalProcessor::extract(std::vector<double> &xInput, std::vector<double> &yInput)
   {

       reset();
       xvec_ = xInput;
       yvec_ = yInput;
       iend_ = xvec_.size();

       if (xvec_.size() < 2) return;


       double parInit[99]={0};
       findPeak(parInit);
       unsigned int nPeak = nparTot_ /nparFcn_;


       if (nparTot_>99) return;

       double chi2(999),fpar[99]={0},errfpar[99]={0};
       doFit(parInit, fpar, errfpar, chi2);



       //final results
       for (int i=0;i<nparTot_; ++i) res_.push_back(fpar[i]);

       chi2_   = chi2;
       nPeaks_ = 0;

       for (unsigned int i=0;i<nPeak;++i)
       { 
	   if (fpar[nparFcn_*i] < 1e-5) continue; 
	   ++nPeaks_;
	   	   
	   resAmp_.push_back(fpar[nparFcn_*i]);
	   resAmpErr_.push_back(errfpar[nparFcn_*i]);
           resTime_.push_back(fpar[nparFcn_*i+1] - shiftTime_);
           resTimeErr_.push_back(errfpar[nparFcn_*i+1]);	        

           if (diagLevel_ > 2)
           {
              _hTime->Fill(fpar[nparFcn_*i+1]);
              _hTimeErr->Fill(errfpar[nparFcn_*i+1]);
              _hEner->Fill(fpar[nparFcn_*i]);
              _hEnerErr->Fill(errfpar[nparFcn_*i]);
           }           
      }
      
      //finally, recalculate chi2
      ndf_ = (iend_- istart_) - (nparFcn_ - nparFix_)*nPeaks_;

   }
      
      
      
   //----------------------------------------------------------------------------------------------------------------------
   void LogNormalProcessor::findPeak(double* parInit)
   {

        std::vector<int> peakLocationInit,peakLocationRes;
        
        //find location of potential peaks: max element in the range i-window; i+window
        for (unsigned int i=windowPeak_;i<xvec_.size()-windowPeak_;++i)
        {
             auto maxp = std::max_element(&yvec_[i-windowPeak_],&yvec_[i+windowPeak_+1]);
             if (maxp == &yvec_[i] && yvec_[i] > minPeakAmplitude_) peakLocationInit.push_back(i);
        }

        for (unsigned int ipeak : peakLocationInit)
        {
	     double currentAmplitudeX = fitfunction(xvec_[ipeak],parInit);

             double loc(xvec_[ipeak]);
             if (ipeak>0 && ipeak<xvec_.size()-1) loc = meanParabol(ipeak,ipeak-1,ipeak+1);

	     parInit[nparTot_++] = peakFactor_*(yvec_[ipeak] - currentAmplitudeX); 
	     parInit[nparTot_++] = loc; 
	     parInit[nparTot_++] = -0.384473;
	     parInit[nparTot_++] = 11.1171;
        }

        
	if (diagLevel_ > 1) std::cout<<"[FixedFastProcessor] Peaks init found : "<<peakLocationInit.size()<<std::endl;        
	if (peakLocationInit.empty()) return; 
        
        unsigned int istart = peakLocationInit.front();
        unsigned int iend   = peakLocationInit.back();
        
        


        // find location of secondary peaks: calculate residualsals, remove regions close to primary peaks, scan and add components
        std::vector<double> residual;
        for (unsigned int i=0;i<xvec_.size();++i) 
        {
             double val = (yvec_[i] > 0 ) ? yvec_[i] - fitfunction(xvec_[i],parInit) : 0 ; 
             residual.push_back(val); 
        }
	for (auto ipeak : peakLocationInit) std::fill(&residual[ipeak-windowPeak_],&residual[ipeak+windowPeak_],0);


        for (unsigned int i=windowPeak_;i<xvec_.size()-windowPeak_;++i)
        {
             auto maxp = std::max_element(&residual[i-windowPeak_],&residual[i+windowPeak_+1]);
             if (maxp != &residual[i] || residual[i]/yvec_[i] < psdThreshold_ || residual[i] < minPeakAmplitude_) continue;
	     peakLocationRes.push_back(i);
        }


        for (unsigned int ipeak : peakLocationRes)
        {
	     double currentAmplitudeAtx = fitfunction(xvec_[ipeak],parInit);

	     parInit[nparTot_++] = peakFactor_*(yvec_[ipeak] - currentAmplitudeAtx); 
	     parInit[nparTot_++] = xvec_[ipeak]; 
	     parInit[nparTot_++] = -0.384473;
	     parInit[nparTot_++] = 11.1171;
             
             if (ipeak < istart) istart = ipeak;
             if (ipeak > iend  ) iend = ipeak; 
        }
	
        istart_ = (istart > 3) ? istart-3 : 0;
        iend_   = (iend+pulseHighBuffer_ < xvec_.size()) ? iend+pulseHighBuffer_ :  xvec_.size(); 

        return;           
   
   }
   
   
   
   //----------------------------------------------------------------------------------------------------------------------
   void LogNormalProcessor::doFit(double* parInit, double *sfpar, double *errsfpar, double& chi2)
   {

        chi2 = 999.0;
        int ierr(0),nvpar(999), nparx(999), istat(999);
        double arglist[2]={0,0},edm(999), errdef(999), temp(0);
        bool refit(false);


	
	TMinuit minuit(nparTot_); 
	minuit.SetFCN(myfcn);

	arglist[0] = printLevel_;
	minuit.mnexcm("SET PRI", arglist ,1,ierr);  
	arglist[0] = 1;
	minuit.mnexcm("SET NOW", arglist, 1, ierr);
	arglist[0] = fitStrategy_;
	minuit.mnexcm("SET STR", arglist ,1,ierr);  

	arglist[0] = 2000;
	arglist[1] = 1e-6;
	
		
	unsigned int nPeak = nparTot_/nparFcn_;
	for (unsigned int ip=0;ip<nPeak;++ip)
	{
	     double p0 = parInit[nparFcn_*ip];
	     double p1 = parInit[nparFcn_*ip+1];
	     minuit.mnparm(nparFcn_*ip+0, "par 0",    p0,      0.01,       0,    1e6, ierr);
	     minuit.mnparm(nparFcn_*ip+1, "par 1",    p1,      0.01,   p1-10,  p1+10, ierr);
	     minuit.mnparm(nparFcn_*ip+2, "par 2", -0.384473,   0.01,    -0.5,   -0.3, ierr);
	     minuit.mnparm(nparFcn_*ip+3, "par 3",  11.1171,   0.01,      10,     12, ierr);
	
  	     if (fixShapeSig_ || xvec_.size() < 4) 
	     {
	         minuit.FixParameter(nparFcn_*ip+2);
  	         minuit.FixParameter(nparFcn_*ip+3);
	     }
	}

	minuit.mnexcm("MIGRAD", arglist ,2,ierr);  
        if (!minuit.fCstatu.Contains("CONVERGED")) minuit.mnexcm("MIGRAD", arglist ,2,ierr);  

	
	//get list of fitted parameters at this stage
	std::vector<double> tempPar(nparTot_,0);
	for (int i=0;i<nparTot_;++i) minuit.GetParameter(i,tempPar[i],temp);    


        //remove small components and those too close to each other
	for (unsigned int ip=0;ip<nPeak;++ip)
        {    
	    
	    double minDTime(999);
	    for (unsigned int j=0;j<ip;++j) 
	       minDTime = std::min(minDTime,std::abs(tempPar[nparFcn_*ip+1] - tempPar[nparFcn_*j+1])); 

	    if (diagLevel_ > 2) _hDelta->Fill(minDTime);
	    if (tempPar[nparFcn_*ip] > minPeakAmplitude_ && minDTime > 6) continue;
	    
            refit = true;
	    minuit.mnparm(nparFcn_*ip, "par 0", 0,  0.01, 0, 1e6, ierr);
	    minuit.FixParameter(nparFcn_*ip);
	    minuit.FixParameter(nparFcn_*ip+1);
	    minuit.FixParameter(nparFcn_*ip+2);
	    minuit.FixParameter(nparFcn_*ip+3);
        }


	if (nPeak > 1 && refit)
	{  
	    minuit.mnexcm("MIGRAD", arglist ,2,ierr);  
	   if (!minuit.fCstatu.Contains("CONVERGED")) minuit.mnexcm("MIGRAD", arglist ,2,ierr);   
	}    


        minuit.mnstat(chi2,edm,errdef,nvpar,nparx,istat);        
        for (int i=0;i<nparTot_;++i) minuit.GetParameter(i,sfpar[i],errsfpar[i]);    
	
        if (diagLevel_ > 2) _hChi2->Fill(chi2/float(ndf_)); 
	
	return;	
   }
   
   
   
   //-----------------------------------
   double LogNormalProcessor::meanParabol(int i1, int i2, int i3)
   {
       double x1 = xvec_[i1];
       double x2 = xvec_[i2];
       double x3 = xvec_[i3];
       double y1 = yvec_[i1];
       double y2 = yvec_[i2];
       double y3 = yvec_[i3];

       double a = ((y1-y2)/(x1-x2)-(y1-y3)/(x1-x3))/(x2-x3);
       double b = (y1-y2)/(x1-x2) - a*(x1+x2);
       if (std::abs(a) < 1e-6) return (x1+x2+x3)/3.0;
       return -b/2.0/a;
   }


   //-----------------------------------
   double LogNormalProcessor::findTime(int ipeak)
   {

      double par[200]={0};
      for (int i=0;i<nparFcn_;++i) par[i] = res_[ipeak*nparFcn_+i];

      double xup  = res_[ipeak*nparFcn_+1];
      double xlow = res_[ipeak*nparFcn_+1] - 50;
      double f0 = timeFraction_*logn(xup,par);
      double df = 0.9;

      while(df > 0.01)
      {
	  double xmed = 0.5*(xup+xlow);
	  double f    = logn(xmed,par);
	  
	  if (f > f0) xup  = xmed;  
	  else        xlow = xmed;

	  df = std::abs(f-f0);
      }

      return 0.5*(xup+xlow);
   }



   //---------------------------------------
   void LogNormalProcessor::plot(std::string pname)
   {

      TH1F h("test","Amplitude vs time",xvec_.size(),xvec_.front()-2.5,xvec_.back()+2.5);
      h.GetXaxis()->SetTitle("Time (ns)");
      h.GetYaxis()->SetTitle("Amplitude");
      for (unsigned int i=0;i<xvec_.size();++i) h.SetBinContent(i+1,yvec_[i]);

      TF1 f("f",fitfunction2,xvec_[istart_],xvec_[iend_],nparTot_);
      for (unsigned int i=0;i<res_.size();++i)f.SetParameter(i,res_[i]);

      TCanvas c1("c1","c1");
      h.Draw();
      f.Draw("same");

      c1.SaveAs(pname.c_str());
   }


}
