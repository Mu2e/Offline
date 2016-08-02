// This is an hybrid signal extraction algorithm with pre-calculated shape function

// If there is one peak, the normalization for a given time can be found analytically. This analytical 
// solution is used as imput to find the starting time with a quasi-Netwon method. Since the Netwon 
// method can oscillate around the minimum, an additional logarithmic scaling is used to make sure 
// the minimum is found  

// If there are more than one peak, we use a generic gradient descent methos, namely minuit.


#include "CaloReco/inc/FixedFastProcessor.hh"
#include "CaloReco/inc/CaloPulseCache.hh"
#include "art/Framework/Services/Optional/TFileDirectory.h" 
#include "art/Framework/Services/Optional/TFileService.h"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/ConditionsService.hh"


#include "TFile.h"
#include "TMinuit.h"
#include "TH2.h"
#include "TF1.h"
#include "TCanvas.h"
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/Matrix.h"

#include <algorithm>
#include <string> 
#include <iostream>
#include <vector>
#include <iterator>




//anonymous namespace containing the data structures and functions required by Minuit

namespace {

   int                    nparTot_, nparFcn_;
   unsigned int           istart_,iend_;
   std::vector<double>    xvec_,yvec_;
   mu2e::CaloPulseCache*  pulseCachePtr_;
   
   //-----------------------------------------------
   double logn(double x, double *par)
   {
        
        int idx = int((x-par[1]+pulseCachePtr_->deltaT())/pulseCachePtr_->step());
        if (idx < 0 || idx > (pulseCachePtr_->cacheSize()-2)) return 0;
     
        double dcache = pulseCachePtr_->cache(idx+1)-pulseCachePtr_->cache(idx);
	double step   = pulseCachePtr_->step();
	return par[0]*(dcache/step*(x-par[1]+pulseCachePtr_->deltaT()-idx*step)+pulseCachePtr_->cache(idx));
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
           if (val>1e-5) {f += (y-val)*(y-val)/val;}
       }
   }
      
}









namespace mu2e {

   //-----------------------------------------------------------------------------
   FixedFastProcessor::FixedFastProcessor(fhicl::ParameterSet const& PSet) :

      WaveformProcessor(PSet),
      windowPeak_         (PSet.get<int>   ("windowPeak")),
      minPeakAmplitude_   (PSet.get<double>("minPeakAmplitude")),
      psdThreshold_       (PSet.get<double>("psdThreshold")),
      pulseHighBuffer_    (PSet.get<int>   ("pulseHighBuffer")),
      shiftTime_          (PSet.get<double>("shiftTime")),
      printLevel_         (PSet.get<int>   ("fitPrintLevel",-1)),
      fitStrategy_        (PSet.get<int>   ("fitStrategy",1)),
      diagLevel_          (PSet.get<int>   ("diagLevel",0)),
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
	          
       art::ServiceHandle<art::TFileService> tfs;
       art::TFileDirectory tfdir = tfs->mkdir("FastFixedDiag");
       _hTime    = tfdir.make<TH1F>("hTime",    "time",                  100, 0., 2000);
       _hTimeErr = tfdir.make<TH1F>("hTimeErr", "time error",            100, 0.,   10);
       _hEner    = tfdir.make<TH1F>("hEner",    "Amplitude",             100, 0., 5000);
       _hEnerErr = tfdir.make<TH1F>("hEnerErr", "Amplitude error",       100, 0.,  100);
       _hChi2    = tfdir.make<TH1F>("hChi2",    "Chi2/ndf",              100, 0.,   20);
       _hNpeak   = tfdir.make<TH1F>("hNpeak",   "Number of peak fitted",  10, 0.,   10);
       _hchi2Amp = tfdir.make<TH2F>("hchi2Amp", "Amp vs chi2/ndf",        50, 0.,   20, 100, 0, 5000);
       _hDelta   = tfdir.make<TH1F>("hDelta",   "Delta t",               100, -50.,   50);
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
       iend_ = xvec_.size();

       if (xvec_.size() < 2) return;

       
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
       ndf_ = (iend_- istart_) - nparFcn_*nPeaks_;
       


       if (diagLevel_ > 2) _hChi2->Fill(chi2/float(ndf_)); 
       if (diagLevel_ > 2) for (auto amp : resAmp_) _hchi2Amp->Fill(chi2/float(ndf_),amp);          
       
       return;
   }


      
   //---------------------------
   void FixedFastProcessor::reset()
   {
       xvec_.clear();
       yvec_.clear();
       res_.clear();
       resAmp_.clear();
       resAmpErr_.clear();
       resTime_.clear();
       resTimeErr_.clear();
       
       istart_  = 0;
       iend_    = 0;
       nparTot_ = 0;
       nPeaks_  = 0;
       chi2_    = 999;
       ndf_     = 0;
   }
      
      
   //----------------------------------------------------------------------------------------------------------------------
   void FixedFastProcessor::findPeak(double* parInit)
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

             parInit[nparTot_++] = pulseCache_.factor()*(yvec_[ipeak] - currentAmplitudeX); 
             parInit[nparTot_++] = loc; 
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
             double currentAMplitudeAtx = fitfunction(xvec_[ipeak],parInit);

             parInit[nparTot_++] = pulseCache_.factor()*(yvec_[ipeak] - currentAMplitudeAtx); 
             parInit[nparTot_++] = xvec_[ipeak]; 
             
             if (ipeak < istart) istart = ipeak;
             if (ipeak > iend  ) iend = ipeak; 
        }
	
        


        istart_ = (istart > 3) ? istart-3 : 0;
        iend_   = (iend+pulseHighBuffer_ < xvec_.size()) ? iend+pulseHighBuffer_ :  xvec_.size(); 

        return;           
   }
   

   
   
   
   //----------------------------------------------------------------------------------------------------------------------
   void FixedFastProcessor::doFit(double* parInit, double *sfpar, double *errsfpar, double& chi2)
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
                
        
        unsigned int nPeak = nparTot_/nparFcn_;
	for (unsigned int ip=0;ip<nPeak;++ip)
        {
             double p0 = parInit[nparFcn_*ip];
             double p1 = parInit[nparFcn_*ip+1];
             minuit.mnparm(nparFcn_*ip+0, "par 0",  p0,  0.001,      0,    1e6, ierr);
             minuit.mnparm(nparFcn_*ip+1, "par 1",  p1,  0.001,  p1-10,  p1+10, ierr);
        }

        
	arglist[0] = 2000;
        arglist[1] = 0.1;
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
        }

        
	if (nPeak > 1 && refit)
	{  
	    minuit.mnexcm("MIGRAD", arglist ,2,ierr);  
	   if (!minuit.fCstatu.Contains("CONVERGED")) minuit.mnexcm("MIGRAD", arglist ,2,ierr);   
	}    


        minuit.mnstat(chi2,edm,errdef,nvpar,nparx,istat);        
        for (int i=0;i<nparTot_;++i) minuit.GetParameter(i,sfpar[i],errsfpar[i]);    
        
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
           double p1   = (calcChi2(tmin+1e-5)-calcChi2(tmin-1e-5))/2e-5;
           double p2   = (calcChi2(tmin+1e-5)-2*calcChi2(tmin)+calcChi2(tmin-1e-5))/1e-5/1e-5;
           
           if (fabs(p2) < 1e-8) p2 = (p2>0) ? 1e-8 : -1e-8;
           tmin = tmin - p1/p2;            
           
           if (fabs(p1) < 1e-3) break;          
       }

       if (nTry==10)
       {
           tmin = refineMin(nTry, tmin, 1);
           tmin = refineMin(nTry, tmin, 0.1);
           tmin = refineMin(nTry, tmin, 0.01);         
       }
       
//can add code here for refit strategy
       
       
       
       double alpha = calcAlpha(tmin);
       double v11   = 0.5*(calcChi2(tmin+1e-5,alpha)-2*calcChi2(tmin,alpha)+calcChi2(tmin-1e-5,alpha))/1e-5/1e-5;
       double v22   = 0.5*(calcChi2(tmin,alpha+1e-5)-2*calcChi2(tmin,alpha)+calcChi2(tmin,alpha-1e-5))/1e-5/1e-5;
       double v12   = 0.5*(calcChi2(tmin+1e-5,alpha+1e-5)-calcChi2(tmin,alpha+1e-5)-calcChi2(tmin+1e-5,alpha)+calcChi2(tmin,alpha))/1e-5/1e-5;
       double det   = v11*v22-v12*v12;

       sfpar[0]    = alpha;
       sfpar[1]    = tmin;
       errsfpar[0] = sqrt(v11/det);
       errsfpar[1] = sqrt(v22/det);

       chi2 = calcChi2(tmin,alpha);
       
       //sfpar[1] = refitTime(tmin,alpha);
      
       return;           
   }  

   
   //------------------------------------------------------------
   double FixedFastProcessor::refineMin(int& nTry, double tmin, double step)
   {
       double t0(tmin);
       double min0  = calcChi2(t0);
       double min1  = calcChi2(t0+0.0001);              
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
     
   //------------------------------------------------------------
   double FixedFastProcessor::meanParabol(int i1, int i2, int i3)
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


   //--------------------------------------------
   double FixedFastProcessor::calcAlpha(double testTime)
   {

       double ytot(0),x2tot(0);
       for (unsigned int i=istart_;i<iend_;++i)
       {
           int idx = int((xvec_[i]-testTime+pulseCache_.deltaT())/pulseCache_.step());         

           if (idx < 0) continue;
           if (idx > (pulseCache_.cacheSize()-2)) break;

           double y = (pulseCache_.cache()[idx+1]-pulseCache_.cache()[idx])/pulseCache_.step()*(xvec_[i]-testTime+pulseCache_.deltaT()-idx*pulseCache_.step())+pulseCache_.cache()[idx];           
           if ( y> 1e-5 ) { ytot += y; x2tot += yvec_[i]*yvec_[i]/y;}       
       }

       return (ytot > 1e-5) ? sqrt(x2tot/ytot) : 0;
   }

   //--------------------------------------------
   double FixedFastProcessor::calcChi2(double testTime, double alpha)
   {
       if (alpha<1e-5) alpha = calcAlpha(testTime);
       
       double difference(0);
       for (unsigned int i=istart_;i<iend_;++i)
       {
           int idx = int((xvec_[i]-testTime+pulseCache_.deltaT())/pulseCache_.step());
           if (idx < 0) continue;
           if (idx > (pulseCache_.cacheSize()-2)) break;

           double y = (pulseCache_.cache()[idx+1]-pulseCache_.cache()[idx])/pulseCache_.step()*(xvec_[i]-testTime+pulseCache_.deltaT()-idx*pulseCache_.step())+pulseCache_.cache()[idx];
           if (y > 0 ) difference += (yvec_[i]-alpha*y)*(yvec_[i]-alpha*y)/alpha/y;     
       }

       return difference;
   }
   


   
   //--------------------------------------------
   double FixedFastProcessor::refitTime(double tmin, double alpha)
   {

       int imax = (tmin-xvec_[0])/(xvec_[1]-xvec_[0]);
       int imin(imax);
       
       
       double yimax = yvec_[imax];
       while( yvec_[imax] > 0.8*yimax && imax > 0) -- imax;
       while( yvec_[imin] > minPeakAmplitude_) --imin;
       
       
       
       istart_ = imin;
       iend_   = imax;
       
       int dummy(0);
       double tmin2 = refineMin(dummy, tmin, 0.002); 
       return tmin2;
       
       
       
       /*
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
















   //---------------------------------------
   void FixedFastProcessor::plot(std::string pname)
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
       std::cout<<"Save file as "<<pname<<std::endl;

       c1.SaveAs(pname.c_str());
       
       return;       
   }


}
