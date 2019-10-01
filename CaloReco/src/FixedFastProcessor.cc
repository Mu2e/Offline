// This is an hybrid signal extraction algorithm with pre-calculated shape function

// If there is one peak, the normalization for a given time can be found analytically. This analytical 
// solution is used as imput to find the starting time with a quasi-Netwon method. Since the Netwon 
// method can oscillate around the minimum, an additional logarithmic scaling is used to make sure 
// the minimum is found  

// If there are more than one peak, we use a generic gradient descent methos, namely minuit.


#include "CaloReco/inc/FixedFastProcessor.hh"
#include "CaloReco/inc/CaloPulseCache.hh"
#include "art_root_io/TFileDirectory.h" 
#include "art_root_io/TFileService.h"
#include "ConditionsService/inc/ConditionsHandle.hh"


#include "TFile.h"
#include "TMinuit.h"
#include "TH2.h"
#include "TF1.h"
#include "TCanvas.h"

#include <algorithm>
#include <string> 
#include <iostream>
#include <vector>




//anonymous namespace containing the data structures and functions required by Minuit

namespace {

   int                       nparTot_, nparFcn_;
   std::vector<unsigned int> xindices_;
   std::vector<double>       xvec_,yvec_;
   mu2e::CaloPulseCache*     pulseCachePtr_;
   
   //-----------------------------------------------
   double logn(double x, double *par) { return par[0]*pulseCachePtr_->evaluate(x-par[1]); }

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
       
       for (unsigned int i : xindices_)
       {    
           double x = xvec_[i];
           double y = yvec_[i];
           double val = fitfunction(x, par);
           if (y>1e-5) f += (y-val)*(y-val)/y;
       }
   }
      
}









namespace mu2e {

   //-----------------------------------------------------------------------------
   FixedFastProcessor::FixedFastProcessor(fhicl::ParameterSet const& PSet) :

      WaveformProcessor(PSet),
      windowPeak_         (PSet.get<int>         ("windowPeak")),
      minPeakAmplitude_   (PSet.get<double>      ("minPeakAmplitude")),
      psdThreshold_       (PSet.get<double>      ("psdThreshold")),
      pulseLowBuffer_     (PSet.get<unsigned int>("pulseLowBuffer")),
      pulseHighBuffer_    (PSet.get<unsigned int>("pulseHighBuffer")),
      minDiffTime_        (PSet.get<unsigned int>("minDiffTime")),
      shiftTime_          (PSet.get<double>      ("shiftTime")),
      printLevel_         (PSet.get<int>         ("fitPrintLevel",-1)),
      fitStrategy_        (PSet.get<int>         ("fitStrategy",1)),
      diagLevel_          (PSet.get<int>         ("diagLevel",0)),
      pulseCache_(CaloPulseCache()),
      nPeaks_(0),
      chi2_(999),
      res_(),
      resAmp_(),
      resAmpErr_(),
      resTime_(),
      resTimeErr_()
   {
       nparTot_ = 0;
       nparFcn_ = 2;
       pulseCachePtr_ = &pulseCache_;
	          
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
	  _hDelta   = tfdir.make<TH1F>("hDelta",   "Delta t",               100, 0.,   20);
       }
   }       


   //------------------------------------------------------------------------------------------    
   void FixedFastProcessor::initialize()
   {
       pulseCache_.initialize();
   }



   //------------------------------------------------------------------------------------------
   void FixedFastProcessor::extract(std::vector<double> &xInput, std::vector<double> &yInput)
   {       

       reset();
       xvec_ = xInput;
       yvec_ = yInput;       
       for (unsigned int i=0; i<xvec_.size(); ++i) xindices_.push_back(i);
       
       if (xInput.size() < 2) return;
      
       
       double parInit[99]={0};
       findPeak(parInit);
       if (nparTot_ > 99) return;       
       unsigned int nPeak = nparTot_/nparFcn_;
       
       
       if (diagLevel_ > 2)_hNpeak->Fill(nPeak);       


       double chi2(999),fpar[99]={0},errfpar[99]={0}; 
       if (nPeak==1) doNewton(parInit, fpar, errfpar, chi2);
       else          doFit(   parInit, fpar, errfpar, chi2);           
       
       
       
       //final results, keep only the good peaks
       for (int i=0;i<nparTot_; ++i) res_.push_back(fpar[i]);

       chi2_   = chi2;
       nPeaks_ = 0;
     
       for (unsigned int i=0;i<nPeak;++i)
       { 
           if (fpar[nparFcn_*i] < 1e-5) continue; 
           ++nPeaks_;
           
           resAmp_.push_back(fpar[nparFcn_*i]);
           resAmpErr_.push_back(errfpar[nparFcn_*i]);
           resTime_.push_back(fpar[nparFcn_*i+1] - shiftTime_ );
           resTimeErr_.push_back(errfpar[nparFcn_*i+1]);          
       
           if (diagLevel_ > 2)
           {
              _hTime->Fill(fpar[nparFcn_*i+1]);
              _hTimeErr->Fill(errfpar[nparFcn_*i+1]);
              _hEner->Fill(fpar[nparFcn_*i]);
              _hEnerErr->Fill(errfpar[nparFcn_*i]);
           }           
       }

       //finally, recalculate ndf = number of bins active in the fit - number of parameters       
       ndf_ = xindices_.size() - nparFcn_*nPeaks_;
       

       if (diagLevel_ > 2)
       {
          _hChi2->Fill(chi2/float(ndf_)); 
          for (auto amp : resAmp_) _hchi2Amp->Fill(chi2/float(ndf_),amp);
       }          
       
       return;
   }


      
   //---------------------------
   void FixedFastProcessor::reset()
   {
       xvec_.clear();
       yvec_.clear();
       xindices_.clear();
       res_.clear();
       resAmp_.clear();
       resAmpErr_.clear();
       resTime_.clear();
       resTimeErr_.clear();
       
       nparTot_ = 0;
       nPeaks_  = 0;
       chi2_    = 999;
       ndf_     = 0;
   }
      
      
   //----------------------------------------------------------------------------------------------------------------------
   void FixedFastProcessor::findPeak(double* parInit)
   {

        std::vector<unsigned int> peakLocationInit,peakLocationRes,peakLocation;
        	
        //find location of potential peaks: max element in the range i-window; i+window
        for (unsigned int i=windowPeak_;i<xvec_.size()-windowPeak_;++i)
        {
             if (std::max_element(&yvec_[i-windowPeak_],&yvec_[i+windowPeak_+1]) != &yvec_[i]) continue;
	     int imin = std::min(std::min(yvec_[i-1],yvec_[i+1]),yvec_[i]);

             if (imin < minPeakAmplitude_) continue;
             //if (yvec_[i] < minPeakAmplitude_) continue;
	     peakLocationInit.push_back(i);
	     peakLocation.push_back(i);
        }

        //fill initial paremeters and record fit range
        for (unsigned int ipeak : peakLocationInit)
        {
             double currentAmplitudeX = fitfunction(xvec_[ipeak],parInit);
             double loc               = meanParabol(ipeak,ipeak-1,ipeak+1);

             parInit[nparTot_++] = pulseCache_.factor()*(yvec_[ipeak] - currentAmplitudeX); 
             parInit[nparTot_++] = loc; 
        }

        


	if (diagLevel_ > 1) std::cout<<"[FixedFastProcessor] Peaks init found : "<<peakLocationInit.size()<<std::endl;        
	if (peakLocationInit.empty()) return; 
                



        // find location of secondary peaks: calculate residuals
	std::vector<double> residual;
        for (unsigned int i=0;i<xvec_.size();++i) 
        {
             double val = (yvec_[i] > 0 ) ? yvec_[i] - fitfunction(xvec_[i],parInit) : 0 ; 
             residual.push_back(val); 
        }


        //find secondat peaks (my best guess so far...)
        for (unsigned int i=windowPeak_;i<xvec_.size()-windowPeak_;++i)
        {
             if (std::max_element(&residual[i-windowPeak_],&residual[i+windowPeak_+1]) != &residual[i]) continue;
             
	     //int imin   = std::min(std::min(residual[i-1],residual[i+1]),residual[i]);
	     double psd = residual[i]/yvec_[i];
	     
             if (residual[i] < minPeakAmplitude_ || psd < psdThreshold_) continue;
	     peakLocationRes.push_back(i);
	     peakLocation.push_back(i);
        }

        //fill initial paremeters and record fit range        
	for (unsigned int ipeak : peakLocationRes)
        {
             double currentAMplitudeAtx = fitfunction(xvec_[ipeak],parInit);

             parInit[nparTot_++] = pulseCache_.factor()*(yvec_[ipeak] - currentAMplitudeAtx); 
             parInit[nparTot_++] = xvec_[ipeak]; 
        }
	
	
	
	// build vectors with x bins used in fit, more flexible than simple start / end range, can omit intermediate pts
	buildXRange(peakLocation);
	
        return;           
   }
   
   

   
   
   
   //----------------------------------------------------------------------------------------------------------------------
   void FixedFastProcessor::doFit(double* parInit, double *sfpar, double *errsfpar, double& chi2)
   {

        chi2 = 999.0;
        int ierr(0),nvpar(999), nparx(999), istat(999);
        double arglist[2]={0,0},edm(999), errdef(999);
	
	
	TMinuit minuit(nparTot_); 
        minuit.SetFCN(myfcn);

        arglist[0] = printLevel_;
        minuit.mnexcm("SET PRI", arglist ,1,ierr);  
        arglist[0] = 1;
        minuit.mnexcm("SET NOW", arglist, 1, ierr);
        arglist[0] = fitStrategy_;
        minuit.mnexcm("SET STR", arglist ,1,ierr);  
                
        
        unsigned int nPeak = nparTot_/nparFcn_;
	for (unsigned int ip=0;ip<nPeak;++ip)
        {
             double p0 = parInit[nparFcn_*ip];
             double p1 = parInit[nparFcn_*ip+1];
             minuit.mnparm(nparFcn_*ip+0, "par 0",  p0,  0.001,      0,    1e6, ierr);
             minuit.mnparm(nparFcn_*ip+1, "par 1",  p1,  0.001,  p1-15,  p1+15, ierr);
        }

        
	arglist[0] = 2000;
        arglist[1] = 0.1;
        minuit.mnexcm("MIGRAD", arglist ,2,ierr);  
        if (!minuit.fCstatu.Contains("CONVERGED")) minuit.mnexcm("MIGRAD", arglist ,2,ierr);          

        
	//get list of fitted parameters at this stage
	std::vector<double> tempPar(nparTot_,0),tempErr(nparTot_,1);
	for (int i=0;i<nparTot_;++i) minuit.GetParameter(i,tempPar[i],tempErr[i]);    


        //remove too small components or those too close to each other (in that case, remove the low peak)
        bool refit(false);
	std::vector<unsigned int> peakLoc;
	
	for (unsigned int ip=0;ip<nPeak;++ip)
        {    
	    double minDTime(999);
	    for (unsigned int j=0;j<nPeak;++j) 
	       if (j!=ip && tempPar[nparFcn_*j] > tempPar[nparFcn_*ip]) 
	           minDTime = std::min( minDTime,std::abs(tempPar[nparFcn_*ip+1] - tempPar[nparFcn_*j+1]) ); 

	    if (diagLevel_ > 2) _hDelta->Fill(minDTime);	     	    
	    if (tempPar[nparFcn_*ip] > minPeakAmplitude_  && minDTime > minDiffTime_) 
	    {
	       peakLoc.push_back((tempPar[nparFcn_*ip+1]-xvec_[0])/(xvec_[1]-xvec_[0]));	       
	       continue;
	    }

	    refit = true;
            tempPar[nparFcn_*ip] = 0;
	    minuit.mnparm(nparFcn_*ip, "par 0", 0,  0.01, 0, 1e6, ierr);
            minuit.FixParameter(nparFcn_*ip);
            minuit.FixParameter(nparFcn_*ip+1);	    	    
        }

        
	if (refit)
	{  
	    buildXRange(peakLoc);
	    minuit.mnexcm("MIGRAD", arglist ,2,ierr);  
	    if (!minuit.fCstatu.Contains("CONVERGED")) minuit.mnexcm("MIGRAD", arglist ,2,ierr);   
	}    


        minuit.mnstat(chi2,edm,errdef,nvpar,nparx,istat);        
        for (int i=0;i<nparTot_;++i) minuit.GetParameter(i,sfpar[i],errsfpar[i]);
        for (int i=0;i<nparTot_;++i) errsfpar[i] = std::abs(  errsfpar[i]); 

        return;
   }
   

   
   
   
   //------------------------------------------------------------------------------
   //Quasi-Newtonian method to find the minimum - please improve when you have time
   void FixedFastProcessor::doNewton(double* parInit, double *sfpar, double *errsfpar, double& chi2)
   {
    
       double tmin(parInit[1]);
       int nTry(0);

       for (;nTry<10;++nTry)
       {
           double p1 = (calcChi2(tmin+1e-4)-calcChi2(tmin-1e-4))/2e-4;
           double p2 = (calcChi2(tmin+1e-4)-2*calcChi2(tmin)+calcChi2(tmin-1e-4))/1e-4/1e-4;
           
           if (fabs(p2) < 1e-8) p2 = (p2>0) ? 1e-8 : -1e-8;
           tmin = tmin - p1/p2;            
           
           if (fabs(p1) < 1e-3) break;          
       }
       
 
       if (nTry==10 || calcChi2(tmin,calcAlpha(tmin)) > 2)
       {
           tmin = parInit[1]; 
	   tmin = refineMin(nTry, tmin, 1);
           tmin = refineMin(nTry, tmin, 0.1);
           tmin = refineMin(nTry, tmin, 0.01);         
       }
       
       double alpha = calcAlpha(tmin);
       double v11   = 0.5*(calcChi2(tmin+1e-5,alpha)-2*calcChi2(tmin,alpha)+calcChi2(tmin-1e-5,alpha))/1e-5/1e-5;
       double v22   = 0.5*(calcChi2(tmin,alpha+1e-5)-2*calcChi2(tmin,alpha)+calcChi2(tmin,alpha-1e-5))/1e-5/1e-5;
       double v12   = 0.5*(calcChi2(tmin+1e-5,alpha+1e-5)-calcChi2(tmin,alpha+1e-5)-calcChi2(tmin+1e-5,alpha)+calcChi2(tmin,alpha))/1e-5/1e-5;
       double det   = v11*v22-v12*v12;

       sfpar[0]    = alpha;
       sfpar[1]    = tmin;
       errsfpar[0] = v11/det > 1e-6 ? sqrt(v11/det) : sfpar[0];
       errsfpar[1] = v22/det > 1e-6 ? sqrt(v22/det) : sfpar[1];

       chi2 = calcChi2(tmin,alpha);
       
       //uncomment to refit time a
       //sfpar[1] = refitTime(tmin,alpha);
      
       return;           
   }  

   
   //------------------------------------------------------------
   double FixedFastProcessor::refineMin(int& nTry, double tmin, double step)
   {
       double t0(tmin);
       double min0  = calcChi2(t0);
       double min1  = calcChi2(t0+0.1*step);              
       double delta = min1<min0 ? step : -step;

       while (nTry < 30)
       {          
           double newmin = calcChi2(t0+delta);
           if (newmin > min0) break;          
           min0 = newmin;  
           t0  += delta;
           ++nTry;
       }
       
       return t0;      
   }
   
   //--------------------------------------------
   double FixedFastProcessor::calcChi2(double testTime, double alpha)
   {
       
       if (alpha<1e-5) alpha = calcAlpha(testTime);
       
       double difference(0);
       for (unsigned int i : xindices_)
       {
	   double y = pulseCachePtr_->evaluate(xvec_[i]-testTime);
           if (yvec_[i] > 1e-5 ) difference += (yvec_[i]-alpha*y)*(yvec_[i]-alpha*y)/yvec_[i];     
       }

       return difference;
   }

   //--------------------------------------------
   double FixedFastProcessor::calcAlpha(double testTime)
   {       
       double ytot(0),x2tot(0);
       for (unsigned int i : xindices_)
       {
	   double y = pulseCachePtr_->evaluate(xvec_[i]-testTime);
           if ( yvec_[i]> 0 ) {ytot += y; x2tot += y*y/yvec_[i];}       
       }

       return (x2tot > 0) ? ytot/x2tot : 0;
   }
    
   //--------------------------------------------
   double FixedFastProcessor::refitTime(double tmin, double alpha)
   {

       //this one refit the leading edge with the shape
       int imax = (tmin-xvec_[0])/(xvec_[1]-xvec_[0]);
       int imin(imax);
              
       double yimax = yvec_[imax];
       while( yvec_[imax] > 0.8*yimax && imax > 0) --imax;
       while( yvec_[imin] > minPeakAmplitude_) --imin;
       
       
       xindices_.clear();
       for (int i = imin; i <= imax; ++i) xindices_.push_back(i);
       
       int dummy(0);
       double tmin2 = refineMin(dummy, tmin, 0.002); 
       return tmin2;
       
       
       
       /*
       #include "CLHEP/Matrix/Vector.h"
       #include "CLHEP/Matrix/Matrix.h"

       //this one does a second order polynomial fit, pick your poison!
       while(imax-imin <3) ++imax;
              
       const int Ndim(3);
       const int Npts  = imax-imin+1;
       const double x0 = xvec_[imin+1];

       CLHEP::HepMatrix mat(Npts,Ndim);
       for (int i=0;i<Npts;++i) mat[i][0] = 1;
       for (int i=0;i<Npts;++i) mat[i][1] = xvec_[imin+i+1]-x0;
       for (int i=0;i<Npts;++i) mat[i][2] = (xvec_[imin+i+1]-x0)*(xvec_[imin+i+1]-x0);

       CLHEP::HepMatrix yv(Npts,1);
       for (int i=0;i<Npts;++i) yv[i][0] = yvec_[imin+i+1];

       CLHEP::HepMatrix matt = mat.T();
       CLHEP::HepMatrix xtx  = matt*mat;       
       CLHEP::HepMatrix pfit = xtx.inverse()*matt*yv;

       double level = 0.1*alpha/pulseCache_.factor();
       double t0 = 0.5*(sqrt(pfit[0][1]*pfit[0][1]-4*pfit[0][2]*(pfit[0][0]-level)) - pfit[0][1])/pfit[0][2] + x0;
       
       return t0;
       */
   }
   
    
    
    


    
    
    
   //-------------------------------------------------------------
   void FixedFastProcessor::buildXRange(std::vector<unsigned int>& peakLoc)
   {      
	std::set<unsigned int> tempX;
	for (unsigned int ipeak : peakLoc)
	{             
             unsigned int is = (ipeak > pulseLowBuffer_) ? ipeak-pulseLowBuffer_ : 0;
	     unsigned int ie = (ipeak+pulseHighBuffer_ < xvec_.size()) ? ipeak+pulseHighBuffer_ :  xvec_.size();	     	     
	     for (unsigned int ip=is; ip<ie; ++ip) tempX.insert(ip); 
	}

	xindices_.clear();
	for (auto i : tempX) xindices_.push_back(i);  
   }
     
   //------------------------------------------------------------
   double FixedFastProcessor::meanParabol(unsigned int i1, unsigned int i2, unsigned int i3)
   {
       if (i1==0 || i3 == xvec_.size()) return xvec_[i1];
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

   //---------------------------------------
   void FixedFastProcessor::plot(std::string pname)
   {
       double dx = xvec_[1]-xvec_[0];
       
       TH1F h("test","Amplitude vs time",xvec_.size(),xvec_.front()-0.5*dx,xvec_.back()+0.5*dx);
       h.GetXaxis()->SetTitle("Time (ns)");
       h.GetYaxis()->SetTitle("Amplitude");
       for (unsigned int i=0;i<xvec_.size();++i) h.SetBinContent(i+1,yvec_[i]);

       int istart = xindices_.back();
       int iend   = xindices_.front();
       
       TF1 f("f",fitfunction2,xvec_[istart],xvec_[iend],nparTot_);
       for (unsigned int i=0;i<res_.size();++i)f.SetParameter(i,res_[i]);

       TCanvas c1("c1","c1");
       h.Draw();
       f.Draw("same");
       std::cout<<"Save file as "<<pname<<std::endl;

       c1.SaveAs(pname.c_str());
       
       return;       
   }


}



//minimize (o-e)^2/e instead of (o-e)^2/o
//void myfcn(int& npar, double* , double &f, double *par, int)
//if (val>1e-5) f += (y-val)*(y-val)/val;

//calcAlpha(double testTime)
//if ( y> 1e-5 ) { ytot += y; x2tot += yvec_[i]*yvec_[i]/y;}       
//return (ytot > 1e-5) ? sqrt(x2tot/ytot) : 0;

//calcChi2(double testTime, double alpha)
//if (y > 1e-5 ) difference += (yvec_[i]-alpha*y)*(yvec_[i]-alpha*y)/alpha/y;     
