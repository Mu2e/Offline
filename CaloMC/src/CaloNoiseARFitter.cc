#include "CaloMC/inc/CaloNoiseARFitter.hh"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "SeedService/inc/SeedService.hh"

#include "TMinuit.h"
#include <string> 
#include <vector>


namespace 
{
    unsigned npFit_ = 0;
    std::vector<double> wf_{};
    
    void myfcn(int& npar, double* , double &f, double *par, int)
    {   
        f=0;       
        for (unsigned i=npFit_+1;i<wf_.size();++i)
        {    
	    double xpred(0);
	    for (unsigned j=0;j<npFit_;++j) xpred += par[j]*wf_[i-j-1];        
	    f += (wf_[i]-xpred)*(wf_[i]-xpred);	
        }    
    }      
}



namespace mu2e {

   CaloNoiseARFitter::CaloNoiseARFitter(CLHEP::HepRandomEngine& engine, unsigned nParFit, int diagLevel) :
     nparFit_       (nParFit),
     param_         (),
     sigmaAR_       (0.0),
     status_        (0),
     randGauss_     (engine),
     diagLevel_     (diagLevel)
   {	                              
       npFit_ = nParFit;
   }       


   //------------------------------------------------------------------------------------------------------------------
   void CaloNoiseARFitter::setWaveform(const std::vector<double>& wf) {wf_=wf;}

   //------------------------------------------------------------------------------------------------------------------
   void CaloNoiseARFitter::fitARCoeff()
   {
       if (wf_.size()<2*nparFit_) 
          throw cet::exception("CATEGORY")<<"[CaloNoiseARFitter] Waveform size too small for AR fit, abort";
       
       param_.clear();

       TMinuit minuit(nparFit_); 
       minuit.SetFCN(myfcn);

       int ierr(0);
       double arglist[2]={0,0};
       arglist[0] = (diagLevel_>0) ? 0 : -1;
       minuit.mnexcm("SET PRI", arglist ,1,ierr);  
       arglist[0] = 1;
       minuit.mnexcm("SET NOW", arglist, 1, ierr);
       arglist[0] = 2;
       minuit.mnexcm("SET STR", arglist ,1,ierr);  

       for (unsigned ip=0;ip<nparFit_;++ip) minuit.mnparm(ip,std::to_string(ip).c_str(),0, 0.001,-100,100,ierr);

       arglist[0] = 2000;
       arglist[1] = 0.1;
       minuit.mnexcm("MIGRAD", arglist ,2,ierr);  
       if (!minuit.fCstatu.Contains("CONVERGED")) return;
       status_=1;       

       double val(0),err(0);
       for (unsigned i=0;i<nparFit_;++i) {minuit.GetParameter(i,val,err); param_.push_back(val);}    

       //calculate standard deviation (sigma) of waveform
       double sm_orig = calcSigma(wf_);

       //adjust the sigma parameter of the random Gaussian to match the original waveform variance
       unsigned niter(0);
       double sigma_low(0.0),sigma_up(sm_orig);
       double sl = ARsigma(sigma_low);
       double su = ARsigma(sigma_up);

       while (niter < 100)
       {
          double sigma = 0.5*(sigma_up+sigma_low);
          double sm    = ARsigma(sigma);
          if (abs(sm_orig-sm) < 0.1) break; 

          if (sm > sm_orig) {su = sm; sl=sl; sigma_up  = sigma;}
          else              {sl = sm; su=su; sigma_low = sigma;}
          ++niter;
       }
       sigmaAR_ = 0.5*(sigma_up+sigma_low);
   }

   //------------------------------------------------------------------------------------------------------------------
   double CaloNoiseARFitter::ARsigma(double sigma0)
   {
       std::vector<double> tnoise(2*wf_.size(),0);
       generateWF(tnoise, sigma0);
       return calcSigma(tnoise);
   }

   //------------------------------------------------------------------------------------------------------------------
   double CaloNoiseARFitter::calcSigma(const std::vector<double>& values)
   {
       double xm(0),xm2(0);       
       for (const auto& val : values) {xm += val; xm2 += val*val;}
       xm  /= values.size();
       xm2 /= values.size();
       return sqrt(xm2-xm*xm);      
   }


   //------------------------------------------------------------------------------------------------------------------
   void CaloNoiseARFitter::generateWF(std::vector<double>& wf)
   {
       generateWF(wf,sigmaAR_);
   }

   //------------------------------------------------------------------------------------------------------------------
   void CaloNoiseARFitter::generateWF(std::vector<double>& wf, double sigma0)
   {       
       if (status_==0) 
            throw cet::exception("CATEGORY")<<"[CaloNoiseARFitter] AR fit failed, noise can't be generated";
 
       //start with pure Gaussian RV, then run the AR sequence for a burn in period
       std::vector<double> burn_in(5*nparFit_,0.0);
       for (unsigned i=0;i<burn_in.size();++i)
       {
           burn_in[i] = randGauss_.fire(0.0,sigma0);
           if (i <nparFit_+1) continue; 
           for (unsigned j=0;j<param_.size();++j) {burn_in[i]+= param_[j]*burn_in[i-j-1];}
       }

       //produce the waveform, first copy the last elements of the burn_in chain and then continue AR generation 
       for (unsigned i=0;i<wf.size();++i)
       {    
           if (i <nparFit_+1) {wf[i] = burn_in[i]; continue;} 
           wf[i] = randGauss_.fire(0.0,sigma0);
           for (unsigned j=0;j<param_.size();++j) {wf[i] += param_[j]*wf[i-j-1];}
       }
   }
}
 

