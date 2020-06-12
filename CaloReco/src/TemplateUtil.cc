#include "CaloReco/inc/TemplateUtil.hh"
#include "Mu2eUtilities/inc/CaloPulseShape.hh"

#include "TMinuit.h"
#include "TF1.h"
#include "TH2.h"
#include "TCanvas.h"

#include <algorithm>
#include <vector>
#include <sstream>



//An anonymous namespace to use Minuit
namespace 
{
   unsigned              npTot_,npFcn_;
   std::vector<double>   xvec_,yvec_;
   mu2e::CaloPulseShape* pulseCachePtr_;
      
   double logn(double x, double *par) {return par[0]*pulseCachePtr_->evaluate(x-par[1]); }

   double fitfunction(double x, double *par)
   {   
       double result(0);
       for (unsigned i=0;i<npTot_;i+=npFcn_) result += logn(x,&par[i]);
       return result;
   }      
   double fitfunctionPlot(double* x, double *par) {return fitfunction(x[0],par);}
   
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
        
   
   TemplateUtil::TemplateUtil(double minPeakAmplitude, double digiSampling) : 
      pulseCache_(CaloPulseShape(digiSampling)),
      minPeakAmplitude_(minPeakAmplitude),
      fitStrategy_(1),
      diagLevel_(0),
      printLevel_(-1),
      param_(),
      paramErr_(),
      nParTot_(2),
      nParFcn_(2),
      chi2_(999.0)
   {
      pulseCachePtr_ = &pulseCache_;
      npTot_ = nParTot_;
      npFcn_ = nParFcn_;
   }       
   

   void   TemplateUtil::initialize  ()                                                                 {pulseCache_.buildShapes();}
   void   TemplateUtil::setXYVector (const std::vector<double>& xvec, const std::vector<double>& yvec) {xvec_ = xvec; yvec_ = yvec;}
   void   TemplateUtil::setPar      (const std::vector<double>& par)                                   {param_ = par; nParTot_= par.size(); npTot_ = par.size();}
   void   TemplateUtil::reset       ()                                                                 {param_.clear(); paramErr_.clear(); nParTot_=0; npTot_ = 0;}

   
   //-----------------------------------------------------------------------------------------------------
   void TemplateUtil::fit() 
   {
       if (param_.empty() || param_.size()>49) return;
       
       if (param_.size()==nParFcn_)  fitNewton();
       else                          fitMinuit();             
   }
   
   
   //-----------------------------------------------------------------------------------------------------
   void TemplateUtil::refitMin()
   {
       if (status_!=0) return;
       
       double tmin = refineMin(param_[1],0.1);
       if (status_==0) param_[1] = tmin;
   }

  
   //-----------------------------------------------------------------------------------------------------
   void TemplateUtil::fitMinuit()
   {      
       status_ = 1;
       if (xvec_.empty() || nParTot_ > 99) return;

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

       for (unsigned ip=0;ip<nParTot_;++ip)
       {
 	    std::stringstream sss; 
	    sss<<"par "<<ip;
            minuit.mnparm(ip, sss.str().c_str(),  param_[ip],  0.001,  0,  1e6, ierr);
       }

       // Perform first fit with initial model
       arglist[0] = 2000;
       arglist[1] = 0.1;
       minuit.mnexcm("MIGRAD", arglist ,2,ierr);  
       //if (!minuit.fCstatu.Contains("CONVERGED")) minuit.mnexcm("MIGRAD", arglist ,2,ierr);          

       
       // Remove small or "close" components and redo the fit with simplified model
       std::vector<double> tempPar(nParTot_,0),tempErr(nParTot_,1);
       for (unsigned i=0;i<nParTot_;++i) minuit.GetParameter(i,tempPar[i],tempErr[i]);    
       
       bool refit(false);
       for (unsigned ip=0; ip<nParTot_/nParFcn_; ++ip)
       {    
           if (selectComponent(tempPar,ip)) continue;           
	   minuit.mnparm(nParFcn_*ip,   "fixed par", 0, 0.01, -1e6, 1e6, ierr);
	   minuit.mnparm(nParFcn_*ip+1, "fixed par", 0, 0.01, 0, 1e6, ierr);
           minuit.FixParameter(nParFcn_*ip);
           minuit.FixParameter(nParFcn_*ip+1);
           refit = true;
       }

       if (refit) minuit.mnexcm("MIGRAD", arglist ,2, ierr);  

       
       // Save the final results - exclude low components
       minuit.mnstat(chi2_,edm,errdef,nvpar,nparx,istat);        
       
       param_.clear();
       paramErr_.clear();

       double val(0),err(0);
       for (unsigned i=0;i<nParTot_; i+=nParFcn_)
       {
           minuit.GetParameter(i,val,err);
           if (val<1e-3) continue;
           param_.push_back(val);
           paramErr_.push_back(err);
	     
           for (unsigned j=1;j<nParFcn_;++j)
           {
               minuit.GetParameter(i+j,val,err);
               param_.push_back(val);
               paramErr_.push_back(err);
	   } 
       }
       
       nParTot_ = param_.size();
       npTot_   = nParTot_;
       status_  = 0;
   }
   
   //----------------------------------------------------------------------------------
   bool TemplateUtil::selectComponent(const std::vector<double>& tempPar, unsigned ip)
   {
        // first check if component is too small
       if (tempPar[nParFcn_*ip] < minPeakAmplitude_) return false;

       //remove peaks close in time with smaller amplitude
       for (unsigned ip2=0; ip2<nParTot_/nParFcn_; ++ip2)
       {
          if (ip==ip2) continue;
          double dt = std::abs(tempPar[nParFcn_*ip2+1]-tempPar[nParFcn_*ip+1]);
          if (dt <20 && tempPar[nParFcn_*ip2]<tempPar[nParFcn_*ip]) return false;
       }
       
       return true;
   }


   //-----------------------------------------------------------------------------------------------------
   //Newtonian method to find the minimum - for a single peak only
   void TemplateUtil::fitNewton()
   {    
       status_=1;
       if (xvec_.empty() || nParTot_ != nParFcn_) return;

       int nTryMax(40);
       double EDM(1e-3), step(0.0001); 

       int nTry(0);
       double delta(0),tmin(param_[1]);
       while (nTry<=nTryMax)
       {
	   double c2 = calcChi2(tmin+step);
	   double c1 = calcChi2(tmin);
	   double p1 = (c2-c1)/step;
           	   
	   if (fabs(p1) < 1e-5) {status_=0; break;}

	   delta = c1/p1;
	   tmin -= delta;

	   ++nTry;          	    
           if (fabs(delta) < EDM ) {status_=0; break;}
           if (fabs(delta) > 100)  {status_=2; break;}  
       }

       //try binary search if we oscillate around the minimum
       if (nTry>nTryMax) tmin = refineMin(tmin, 0.1); 
       if (status_==2)   tmin = refineMin(param_[1], 0.1); 

       param_.clear();
       paramErr_.clear();
       if (status_==0) 
       {  
          double alpha = calcAlpha(tmin);
          double v11   = 0.5*(calcChi2(tmin+step,alpha)-2*calcChi2(tmin,alpha)+calcChi2(tmin-step,alpha))/step/step;
          double v22   = 0.5*(calcChi2(tmin,alpha+step)-2*calcChi2(tmin,alpha)+calcChi2(tmin,alpha-step))/step/step;
          double v12   = 0.5*(calcChi2(tmin+step,alpha+step)-calcChi2(tmin,alpha+step)-calcChi2(tmin+step,alpha)+calcChi2(tmin,alpha))/step/step;
          double det   = v11*v22-v12*v12;
	  double ealpha = (det > 1e-5) ? sqrt(v11/det) : 1000;
	  double etmin  = (det > 1e-5) ? sqrt(v22/det) : 1000;

          param_.push_back(alpha);
          param_.push_back(tmin);
          paramErr_.push_back(ealpha);
          paramErr_.push_back(etmin);
          chi2_ = calcChi2(tmin,alpha);
       }
   }  

   //------------------------------------------------------------
   double TemplateUtil::refineMin(double tmin, double stepInit)
   {        
       double t0(tmin), step(stepInit);
       unsigned nSteps(0);

       status_=1;
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
       if (diagLevel_ > 2) std::cout<<"TemplateUtil::refineMin  NSteps="<<nSteps<<" time "<<param_[1]<<" -> "<<t0<<std::endl;   

       if (nSteps>100) status_=2;
       else            status_=0; 

       return t0;
   }

   //--------------------------------------------
   double TemplateUtil::calcAlpha(double testTime)
   {       
       double ytot(0),x2tot(0);
       for (unsigned i=0;i<xvec_.size();++i)
       {
	   double y = pulseCachePtr_->evaluate(xvec_[i]-testTime);
	   if (yvec_[i]> 0) {ytot += y*yvec_[i]; x2tot += y*y;}       
       }
       return (x2tot > 1e-5) ? ytot/x2tot : 0.0;
   }

   //----------------------------------------------------------
   // NOTE: SHOULD WE USE THE FLOORED VALUE OR NOT?
   double TemplateUtil::calcChi2(double testTime, double alpha)
   {       
       //use pre-cached alpha if available
       if (alpha < 1e-5) alpha = calcAlpha(testTime);  

       double difference(0);
       for (unsigned i=0;i<xvec_.size();++i)
       {
           double y = alpha*pulseCachePtr_->evaluate(xvec_[i]-testTime);
           //float y = floor(alpha*pulseCachePtr_->evaluate(xvec_[i],testTime));
           if (yvec_[i] > 1e-5 ) difference += (yvec_[i]-y)*(yvec_[i]-y)/yvec_[i];     
       }
       return difference;
   }

   //----------------------------------------------------------------------
   double TemplateUtil::eval_fcn(double x)
   {       
       if (param_.size()<nParFcn_) return 0.0;
       return fitfunction(x,&param_[0]);       
   }
   //------------------------------------------------------------
   double TemplateUtil::eval_logn(double x, int ioffset)
   {
       if (param_.size() < ioffset+nParFcn_) return 0.0;
       return logn(x,&param_[ioffset]);       
   }

   //------------------------------------------------------------
   void TemplateUtil::calcTimeCorr(std::vector<double>& par)
   {	
	int Npt(0);
	double sx(0),sy(0),sxy(0),sx2(0);
	for (double toffset=0;toffset<4.99;toffset += 0.25)
	{
	    std::vector<double> t1 = pulseCachePtr_->digitizedPulse(toffset);
	    std::vector<double> wf(t1.size(),0);
	    for (unsigned i=0;i<t1.size();++i) wf[i] += (100*t1[i]);

	    unsigned imax = std::distance(wf.begin(),std::max_element(wf.begin(),wf.end()));
	    double tmin   = 1.5-(0.5*wf[imax-1]+1.5*wf[imax]+2.5*wf[imax+1])/(wf[imax-1]+wf[imax]+wf[imax+1]);

	    sx  += wf[imax+1]/wf[imax];
	    sx2 += wf[imax+1]*wf[imax+1]/wf[imax]/wf[imax];
	    sy  += tmin;
	    sxy += tmin*wf[imax+1]/wf[imax];    
	    ++Npt;
	}
        double m = (Npt*sxy-sx*sy)/(Npt*sx2-sx*sx);	
	par.push_back((sy-m*sx)/Npt);
	par.push_back(m);
   }




   //---------------------------------------------------------------------------
   void TemplateUtil::plotFit(const std::string& pname) const
   {
       if (xvec_.empty()) return;
       double dx = xvec_[1]-xvec_[0];

       TH1F h("test","Amplitude vs time",xvec_.size(),xvec_.front()-0.5*dx,xvec_.back()+0.5*dx);
       for (unsigned int i=0;i<xvec_.size();++i)h.SetBinContent(i+1,yvec_[i]);
       h.GetXaxis()->SetTitle("Time (ns)");
       h.GetYaxis()->SetTitle("Amplitude"); 
       h.SetStats(0);

       TF1 f("f",fitfunctionPlot,xvec_.front(),xvec_.back(),param_.size());
       for (unsigned i=0;i<param_.size();++i) f.SetParameter(i,param_[i]);

       TCanvas c1("TemplateUtil::Fit","TemplateUtil::Fit");
       h.Draw();
       if (!param_.empty()) f.Draw("same");
       std::cout<<"Save file as "<<pname<<std::endl;
       c1.SaveAs(pname.c_str());

       return;       
   }




}
