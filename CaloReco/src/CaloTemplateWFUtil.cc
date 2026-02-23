#include "Offline/CaloReco/inc/CaloTemplateWFUtil.hh"
#include "Offline/Mu2eUtilities/inc/CaloPulseShape.hh"

#include "TMinuit.h"
#include "TF1.h"
#include "TH2.h"
#include "TCanvas.h"

#include <algorithm>
#include <vector>
#include <sstream>

//
// The fit function has been modified to include the correlations in the waveform (signal bins are fully correlated)
// Instead of using the full covariance matrix, one can obtain similr results if the expected numnber of events
// in a bin is replaced by the background estimate, and the uncertainty on the efficiency is modified to account for
// the signal (see doc-db 36707 for a full explanation)


//An anonymous namespace to use Minuit
namespace
{
    unsigned              npTot_(0),npFcn_(0),npBkg_(0),x0_(0),x1_(0);
    std::vector<double>   xvec_{},yvec_{};
    mu2e::CaloPulseShape* pulseCachePtr_=(nullptr);

    double logn(double x, double *par) {return par[0]*pulseCachePtr_->evaluate(x-par[1]); }

    double fitfunction(double x, double *par)
    {
        double result(par[0]);
        for (unsigned i=npBkg_; i<npTot_; i+=npFcn_) result += logn(x,&par[i]);
        return result;
    }
    double fitfunctionPlot(double* x, double *par) {return fitfunction(x[0],par);}

    void myfcn(int& npar, double* , double &f, double *par, int)
    {
        f=0;
        for (unsigned i=x0_;i<x1_;++i)
        {
            double x = xvec_[i];
            double y = yvec_[i];
            double val = fitfunction(x, par);
            // modified fit function
            if (fabs(par[0]) > 1e-5) f += (y-val)*(y-val)/par[0];
        }
    }
}







namespace mu2e {


   CaloTemplateWFUtil::CaloTemplateWFUtil(const std::string& pulseFileName, const std::string& pulseHistName,
                                          double minPeakAmplitude, double digiSampling, double minDTPeaks, int printLevel) :
      pulseCache_(CaloPulseShape(pulseFileName, pulseHistName, digiSampling)),
      minPeakAmplitude_(minPeakAmplitude),
      minDTPeaks_(minDTPeaks),
      fitStrategy_(1),
      diagLevel_(0),
      printLevel_(printLevel),
      param_(),
      paramErr_(),
      nParTot_(3),
      nParFcn_(2),
      nParBkg_(1),
      chi2_(999.0)
   {
      pulseCachePtr_ = &pulseCache_;
      npTot_ = nParTot_;
      npFcn_ = nParFcn_;
      npBkg_ = nParBkg_;
   }


   //-----------------------------------------------------------------------------------------------------
   void   CaloTemplateWFUtil::initialize ()                                                                 {pulseCache_.buildShapes();}
   void   CaloTemplateWFUtil::reset      ()                                                                 {param_.clear(); paramErr_.clear(); nParTot_=0; npTot_ = 0;}
   void   CaloTemplateWFUtil::setXYVector(const std::vector<double>& xvec, const std::vector<double>& yvec) {xvec_ = xvec; yvec_ = yvec; x0_=0; x1_ = xvec_.size();}
   void   CaloTemplateWFUtil::setPar     (const std::vector<double>& par)                                   {param_ = par; nParTot_ = npTot_ = par.size();}

   //-----------------------------------------------------------------------------------------------------
   void CaloTemplateWFUtil::fit()
   {
       status_ = 0;
       if (param_.empty() || param_.size()>49 || xvec_.empty()) return;
       if (nParTot_ < nParBkg_  || (nParTot_-nParBkg_)%nParFcn_ !=0) return;

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
             std::string sss = "par " + std::to_string(ip);
            minuit.mnparm(ip, sss.c_str(),  param_[ip],  0.001,  0,  1e6, ierr);
       }

       // Perform first fit with initial model
       //
       arglist[0] = 2000;
       arglist[1] = 0.1;
       minuit.mnexcm("MIGRAD", arglist ,2,ierr);
       //if (!minuit.fCstatu.Contains("CONVERGED")) minuit.mnexcm("MIGRAD", arglist ,2,ierr);


       // Remove small or "duplicate" components and redo the fit with simplified model if there is more than one peak
       // A duplicate peak is defined as a peak shortly after a previous peak with a smaller amplitude
       //
       if (nParTot_ > nParFcn_+nParBkg_)
       {
           bool refit(false);
           std::vector<double> tempPar(nParTot_,0),tempErr(nParTot_,1);
           for (unsigned i=0;i<nParTot_;++i) minuit.GetParameter(i,tempPar[i],tempErr[i]);

           for (unsigned ip=nParBkg_; ip<nParTot_; ip += nParFcn_)
           {
               if (selectComponent(tempPar,tempErr,ip)) continue;
               minuit.mnparm(ip,   "fixed par", 0, 0.01, -1e6, 1e6, ierr);
               minuit.mnparm(ip+1, "fixed par", 0, 0.01, -1e6, 1e6, ierr);
               minuit.FixParameter(ip);
               minuit.FixParameter(ip+1);
               refit = true;
           }

           if (refit) minuit.mnexcm("MIGRAD", arglist ,2, ierr);
       }


       // Save the results - exclude low components
       //
       param_.clear();
       paramErr_.clear();

       unsigned i(0);
       while (i<nParTot_)
       {
           double val(0),err(0);
           minuit.GetParameter(i,val,err);

           //if the amplitude is too small, jump to the next peak
           if (val<1 && i >=nParBkg_ && (i-nParBkg_)%nParFcn_==0) {i+=nParFcn_;continue;}

           param_.push_back(val);
           paramErr_.push_back(err);
           ++i;
       }

       minuit.mnstat(chi2_,edm,errdef,nvpar,nparx,istat);

       //recalculate the chi2 removing the baseline to better reject the noise ?
       //chi2_=0;
       //for (unsigned i=x0_;i<x1_;++i)
       //{
       //    double val = fitfunction(xvec_[i], &param_[0]);
       //    if (yvec_[i]>1e-5) chi2_ += (yvec_[i]-val)*(yvec_[i]-val)/(yvec_[i]-param_[0]);
       //}

       nParTot_ = param_.size();
       npTot_   = nParTot_;
       status_  = istat;
   }



   //-----------------------------------------------------------------------------------------------------
   void CaloTemplateWFUtil::refitEdge()
   {
       status_ = 0;
       if (param_.size()<nParBkg_+nParFcn_ || xvec_.empty()) return;

       unsigned imax(0),ilow(0);
       while (xvec_[imax]<param_[2]) ++imax;
       for (unsigned i=imax;i>0;--i) if ((yvec_[i]-param_[0])/(yvec_[imax]-param_[0])>0.1) ilow = i;
       if (imax < ilow+4) return; //need at least 4 points to fit
       x0_ = 0;
       //x0_ = ilow;
       x1_ = imax;

       int ierr(0),nvpar(999), nparx(999), istat(999);
       double arglist[2]={0,0}, edm(999), errdef(999),chi(9999),val(0),err(0);

       TMinuit minuit(nParBkg_+nParFcn_);
       minuit.SetFCN(myfcn);

       arglist[0] = printLevel_;
       minuit.mnexcm("SET PRI", arglist ,1,ierr);
       arglist[0] = 1;
       minuit.mnexcm("SET NOW", arglist, 1, ierr);
       arglist[0] = fitStrategy_;
       minuit.mnexcm("SET STR", arglist ,1,ierr);

       for (unsigned ip=0;ip<nParBkg_+nParFcn_;++ip) minuit.mnparm(ip, "par",  param_[ip],  0.001, 0,  1e6, ierr);

       arglist[0] = 2000;
       arglist[1] = 0.1;
       minuit.mnexcm("MIGRAD", arglist ,2,ierr);

       minuit.GetParameter(nParBkg_+1,val,err);
       minuit.mnstat(chi,edm,errdef,nvpar,nparx,istat);

       param_[nParBkg_+1]    = val;
       paramErr_[nParBkg_+1] = err;
       status_               = istat;

       x0_     = 0;
       x1_     = xvec_.size();
   }

   //----------------------------------------------------------------------------------
   bool CaloTemplateWFUtil::selectComponent(const std::vector<double>& tempPar, const std::vector<double>& tempErr, unsigned ip)
   {
       // first check if component is too small, error too large or out of time
       if (tempPar[ip] < minPeakAmplitude_)                               return false;
       if (tempPar[ip+1] < xvec_.front() || tempPar[ip+1] > xvec_.back()) return false;
       if (tempErr[ip] >1e3)                                              return false;

       //remove peaks close in time with smaller amplitude
       for (unsigned ip2=nParBkg_; ip2<npTot_; ip2 += nParFcn_)
       {
           if (ip==ip2) continue;
           double dt = std::abs(tempPar[ip2+1]-tempPar[ip+1]);
           if (dt <minDTPeaks_ && tempPar[ip2]>tempPar[ip]) return false;
       }

       return true;
   }




   //----------------------------------------------------------------------
   double CaloTemplateWFUtil::eval_fcn(double x)
   {
       if (param_.size()<nParFcn_) return 0.0;
       return fitfunction(x,&param_[0]);
   }
   //------------------------------------------------------------
   double CaloTemplateWFUtil::eval_logn(double x, int ioffset)
   {
       if (param_.size() < ioffset+nParFcn_) return 0.0;
       return logn(x,&param_[ioffset]);
   }
   //------------------------------------------------------------
   double CaloTemplateWFUtil::maxAmplitude()
   {
       double maxAmplitudeFound(0.0);
       for (unsigned i=nParBkg_; i<param_.size(); i += nParFcn_) maxAmplitudeFound = std::max(maxAmplitudeFound,param_[i]);
       return maxAmplitudeFound;
   }

   //------------------------------------------------------------
   double CaloTemplateWFUtil::peakNorm(const std::vector<double>& xvalues, const std::vector<double>& yvalues, double x0, unsigned i0, unsigned i1)
   {
      double s1(0),s2(0);
      for (unsigned i=i0;i<=i1;++i)
      {
         double ff = pulseCachePtr_->evaluate(xvalues[i]-x0);

         s1 += ff*ff;
         s2 += yvalues[i]*ff;
      }
      if (std::abs(s1) < 1e-6) return 1e6;
      return s2/s1;
   }

   //------------------------------------------------------------
   double CaloTemplateWFUtil::sumSquare(const std::vector<double>& xvalues, const std::vector<double>& yvalues, double x0, unsigned i0, unsigned i1)
   {
      double A = peakNorm(xvalues, yvalues, x0, i0, i1);
      double chi2(0);
      for (unsigned i=i0;i<=i1;++i)
      {
         double cc = A*pulseCachePtr_->evaluate(xvalues[i]-x0)-yvalues[i];
         chi2 += cc*cc;
      }
      return chi2;
   }

   //------------------------------------------------------------
   double CaloTemplateWFUtil::peakToFunc(unsigned ip, double xmax, double ymax)
   {
      if (ip+1 > param_.size()) return 1e6;
      return ymax*pulseCache_.evaluate(param_[ip+1]-xmax)/param_[ip];
   }


   //---------------------------------------------------------------------------
   void CaloTemplateWFUtil::plotFit(const std::string& pname) const
   {
       if (xvec_.empty()) return;
       double dx = xvec_[1]-xvec_[0];

       TH1F h("test","Amplitude vs time",x1_-x0_,xvec_[x0_]-0.5*dx,xvec_[x1_-1]+0.5*dx);
       for (unsigned i=x0_;i<x1_;++i) h.SetBinContent(i+1-x0_,yvec_[i]);
       h.GetXaxis()->SetTitle("Time (ns)");
       h.GetYaxis()->SetTitle("Amplitude");
       h.SetStats(0);
       h.SetMinimum(0);

       TF1 f("f",fitfunctionPlot,xvec_[x0_],xvec_[x1_-1],param_.size());
       for (unsigned i=0;i<param_.size();++i) f.SetParameter(i,param_[i]);

       std::vector<TF1*> f2(nParTot_);
       int nPeaks(0);
       for (unsigned i=nParBkg_;i<nParTot_;i+=nParFcn_)
       {
          f2[nPeaks] = new TF1("f2",fitfunctionPlot,xvec_[x0_],xvec_[x1_-1],param_.size());
          for (unsigned j=0;j<nParTot_;++j) f2[nPeaks]->SetParameter(i,0);
          f2[nPeaks]->SetParameter(i,param_[i]);
          f2[nPeaks]->SetParameter(i+1,param_[i+1]);
          f2[nPeaks]->SetLineStyle(2);
          f2[nPeaks]->SetLineColor(92);
          ++nPeaks;
       }

       TCanvas c1("cc1","cc1");
         h.Draw("");
         if (!param_.empty()) f.Draw("same");
         for (int i=0; i<nPeaks;++i) f2[i]->Draw("same");
         std::cout<<"Save file as "<<pname<<std::endl;
       c1.SaveAs(pname.c_str());

       for (int i=0;i<nPeaks;++i) delete f2[i];

       return;
   }

}
