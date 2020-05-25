#include "CaloReco/inc/FastFixedUtil.hh"
#include "TMinuit.h"
#include "TF1.h"
#include "TH2.h"
#include "TCanvas.h"
#include <algorithm>
#include <vector>



//An anonymous namespace to use Minuit
namespace 
{
   unsigned                  nParTot_;
   unsigned                  nParFcn_;
   std::vector<double>       xvec_,yvec_;
   mu2e::CaloPulseCache*     pulseCachePtr_;
   

   double logn(double x, double *par) {return par[0]*pulseCachePtr_->evaluate(x-par[1]); }

   double fitfunction(double x, double *par)
   {   
       double result(0);
       for (unsigned i=0;i<nParTot_;i+=nParFcn_) result += logn(x,&par[i]);
       return result;
   }      
   double fitfunction2(double* x, double *par) {return fitfunction(x[0],par);}
   
   void myfcn(int& npar, double* , double &f, double *par, int)
   {   
       f=0;       
       for (unsigned i=0;i<xvec_.size();++i)
       {    
           double x = xvec_[i];
           double y = yvec_[i];
           double val = fitfunction(x, par);
           if (y>1e-5) f += (y-val)*(y-val)/y;
       }
   }      
}






namespace mu2e {
        
   
   FastFixedUtil::FastFixedUtil(unsigned nParFcn, double minPeakAmplitude, double minDiffTime,
                                int printLevel, int fitStrategy, int diagLevel) :
      pulseCache_(CaloPulseCache()),
      diagLevel_(diagLevel),
      printLevel_(printLevel),
      fitStrategy_(fitStrategy),
      minPeakAmplitude_(minPeakAmplitude),
      minDiffTime_(minDiffTime),
      isFitDone_(false),
      sfpar_(),
      esfpar_(),
      chi2_(999.0)
   {
       nParTot_ = nParFcn;
       nParFcn_ = nParFcn;
       pulseCachePtr_ = &pulseCache_;
   }       

   void FastFixedUtil::initialize()
   {
       pulseCache_.initialize();
   }

   void FastFixedUtil::setXYVector(const std::vector<double>& xvec, const std::vector<double>& yvec)
   {
       xvec_ = xvec;
       yvec_ = yvec;      
   }

   void FastFixedUtil::fitMinuit(const std::vector<double>& parInit)
   {      
       sfpar_.clear();
       esfpar_.clear();
       nParTot_ = parInit.size();
       isFitDone_ = false;
       if (xvec_.empty() || parInit.size() > 99) return;

       int ierr(0),nvpar(999), nparx(999), istat(999);
       double arglist[2]={0,0}, edm(999), errdef(999);

       TMinuit minuit(nParTot_); 
       minuit.SetFCN(myfcn);

       arglist[0] = printLevel_;
       minuit.mnexcm("SET PRI", arglist ,1,ierr);  
       arglist[0] = 1;
       minuit.mnexcm("SET NOW", arglist, 1, ierr);
       arglist[0] = fitStrategy_;
       minuit.mnexcm("SET STR", arglist ,1,ierr);  

       unsigned nPeak = nParTot_/nParFcn_;
       for (unsigned ip=0;ip<nPeak;++ip)
       {
            double p0 = parInit[nParFcn_*ip];
            double p1 = parInit[nParFcn_*ip+1];
            minuit.mnparm(nParFcn_*ip+0, "par 0",  p0,  0.001,      0,    1e6, ierr);
            minuit.mnparm(nParFcn_*ip+1, "par 1",  p1,  0.001,  p1-15,  p1+15, ierr);
       }

       arglist[0] = 2000;
       arglist[1] = 0.1;
       minuit.mnexcm("MIGRAD", arglist ,2,ierr);  
       if (!minuit.fCstatu.Contains("CONVERGED")) minuit.mnexcm("MIGRAD", arglist ,2,ierr);          

       //get list of fitted parameters  and remove components too small or too close to each other
       std::vector<double> tempPar(nParTot_,0),tempErr(nParTot_,1);
       for (unsigned i=0;i<nParTot_;++i) minuit.GetParameter(i,tempPar[i],tempErr[i]);    

       bool refit(false);
       for (unsigned ip=0;ip<nPeak;++ip)
       {    
	   double minDTime = meanDistance(ip,tempPar);
           if (tempPar[nParFcn_*ip] < minPeakAmplitude_  || minDTime < minDiffTime_)
           { 
              refit = true;
              tempPar[nParFcn_*ip]=0;              
	      minuit.mnparm(nParFcn_*ip,   "par 0", 0,                      0.01, -1e6, 1e6, ierr);
	      minuit.mnparm(nParFcn_*ip+1, "par 1", tempPar[nParFcn_*ip+1], 0.01, 0, 1e6, ierr);
              minuit.FixParameter(nParFcn_*ip);
              minuit.FixParameter(nParFcn_*ip+1);
           }  	    	    
       }

       minuit.mnstat(chi2_,edm,errdef,nvpar,nparx,istat);        

       if (refit) minuit.mnexcm("MIGRAD", arglist ,2,ierr);  
       isFitDone_ = true;

       //only save the non-zero peaks       
       minuit.mnstat(chi2_,edm,errdef,nvpar,nparx,istat);        
       for (unsigned i=0;i<nParTot_; i+=nParFcn_)
       {
           double val(0),err(0);
           minuit.GetParameter(i,val,err);
 
           if (val<1e-3) continue;
           sfpar_.push_back(val);
           esfpar_.push_back(err);
  
           for (unsigned j=1;j<nParFcn_;++j)
           {
              minuit.GetParameter(i+j,val,err);
              sfpar_.push_back(val);
              esfpar_.push_back(err);
           } 
       }
    
    }

    //---------------------------------------------------------------
    //Newtonian method to find the minimum - usable for a single peak
    void FastFixedUtil::fitNewton(double tminInit)
    {    
        sfpar_.clear();
        esfpar_.clear();
        nParTot_ = nParFcn_;
        isFitDone_ = false;
        if (xvec_.empty()) return;

        int nTryMax(10);
        double EDM_p1(1e-3);
        double delta(0.01*EDM_p1); 

        int nTry(0);
        double tmin(tminInit);
        while (nTry<nTryMax)
        {
            double c1 = calcChi2(tmin+delta);
            double c2 = calcChi2(tmin);
            double c3 = calcChi2(tmin-delta);           
            double p1 = (c1-c3)/2.0/delta;
            double p2 = (c3-2.0*c2+c1)/delta/delta;                      

            if (fabs(p2) < 1e-8){nTry=nTryMax; break;}
            
            tmin -= p1/p2;                       
            if (p1 < EDM_p1) break;  
            ++nTry;          
        }
        isFitDone_ = true;

        if (nTry<nTryMax) 
        {
            tmin = refineMin(tmin,0.0011);
        } 
        else 
        {
            if (diagLevel_ > 1) std::cout<<"FastFixedUtil::doNewton non convergence"<<std::endl;        
            tmin = refineMin(tminInit,0.1);
        }          
 
        double alpha = calcAlpha(tmin);
        double v11   = 0.5*(calcChi2(tmin+delta,alpha)-2*calcChi2(tmin,alpha)+calcChi2(tmin-delta,alpha))/delta/delta;
        double v22   = 0.5*(calcChi2(tmin,alpha+delta)-2*calcChi2(tmin,alpha)+calcChi2(tmin,alpha-delta))/delta/delta;
        double v12   = 0.5*(calcChi2(tmin+delta,alpha+delta)-calcChi2(tmin,alpha+delta)-calcChi2(tmin+delta,alpha)+calcChi2(tmin,alpha))/delta/delta;
        double det   = v11*v22-v12*v12;

        sfpar_.push_back(alpha);
        sfpar_.push_back(tmin);
        esfpar_.push_back(sqrt(v11/det));
        esfpar_.push_back(sqrt(v22/det));
        chi2_ = calcChi2(tmin,alpha);

        return;           
    }  


    //------------------------------------------------------------
    double FastFixedUtil::refineMin(double t0init, double stepInit)
    {
        if (!isFitDone_) return -999.9;
        
        double t0(t0init);
        double step(stepInit);

        unsigned nSteps(0);
        while (step > 1e-3)
        {
           if (nSteps > 100) break;
           double current  = calcChi2(t0);
           double left     = calcChi2(t0-step);              
           double right    = calcChi2(t0+step);

           if      (left  < current) {t0-= step;}
           else if (right < current) {t0+= step;}
           else                      {step /= 2;}
           ++nSteps;
        }

        if (diagLevel_ > 2) std::cout<<"FastFixedUtil::refineMin  NSteps="<<nSteps<<std::endl;   
        return t0;      
    }

    //--------------------------------------------
    double FastFixedUtil::calcChi2(double testTime, double alpha)
    {       
        if (alpha < 1e-5) alpha = calcAlpha(testTime); //minimize formula for difference below

        double difference(0);
        for (unsigned i=0;i<xvec_.size();++i)
        {
	    double y = pulseCachePtr_->evaluate(xvec_[i]-testTime);
            if (yvec_[i] > 1e-5 ) difference += (yvec_[i]-alpha*y)*(yvec_[i]-alpha*y)/yvec_[i];     
        }

        return difference;
    }

    //--------------------------------------------
    double FastFixedUtil::calcAlpha(double testTime)
    {       
        double ytot(0),x2tot(0);
        for (unsigned i=0;i<xvec_.size();++i)
        {
	    double y = pulseCachePtr_->evaluate(xvec_[i]-testTime);
            if ( yvec_[i]> 0 ) {ytot += y; x2tot += y*y/yvec_[i];}       
        }

        return (x2tot > 0) ? ytot/x2tot : 0;
    }

    //------------------------------------------------------------
    double FastFixedUtil::meanDistance(unsigned ip, const std::vector<double>& tempPar)
    {
       double minDTime(999);
       unsigned int nPeak = nParTot_/nParFcn_;
       for (unsigned int j=0;j<nPeak;++j) 
	  if (j!=ip && tempPar[nParFcn_*j] > tempPar[nParFcn_*ip]) 
	      minDTime = std::min( minDTime,std::abs(tempPar[nParFcn_*ip+1] - tempPar[nParFcn_*j+1]) ); 
       return minDTime;       
    }

    //------------------------------------------------------------
    double FastFixedUtil::eval_logn(double x, int ioffset, const std::vector<double>& par)
    {
       if (par.empty() || par.size() > 98) return 0;
       double param[99]={0};
       for (unsigned j=ioffset;j<ioffset+nParFcn_;++j) param[j] = par[j]; 
       return logn(x,param);       
    }
    //------------------------------------------------------------
    double FastFixedUtil::eval_fcn(double x, const std::vector<double>& par)
    {
       
       if (par.empty() || par.size() > 98) return 0;
       double param[99]={0};
       for (unsigned j=0;j<par.size();++j) param[j] = par[j]; 
       return fitfunction(x,param);       
    }


    //---------------------------------------
    void FastFixedUtil::plot(std::string pname, const std::vector<double>& param)
    {
        double dx = xvec_[1]-xvec_[0];

        TH1F h("test","Amplitude vs time",xvec_.size(),xvec_.front()-0.5*dx,xvec_.back()+0.5*dx);
        h.GetXaxis()->SetTitle("Time (ns)");
        h.GetYaxis()->SetTitle("Amplitude");
        for (unsigned int i=0;i<xvec_.size();++i) h.SetBinContent(i+1,yvec_[i]);

        TF1 f("f",fitfunction2,xvec_.front(),xvec_.back(),param.size());
        for (unsigned i=0;i<param.size();++i) f.SetParameter(i,param[i]);

        TCanvas c1("c1","c1");
        h.Draw();
        f.Draw("same");
        std::cout<<"Save file as "<<pname<<std::endl;
        c1.SaveAs(pname.c_str());

        return;       
    }


}
