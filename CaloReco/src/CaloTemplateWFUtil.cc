#include "Offline/CaloReco/inc/CaloTemplateWFUtil.hh"
#include "Offline/Mu2eUtilities/inc/CaloPulseUtil.hh"

#include "TF1.h"
#include "TH2.h"
#include "TCanvas.h"

#include <algorithm>
#include <cmath>
#include <vector>
#include <sstream>

//
// The fit function has been modified to include the correlations in the waveform (signal bins are fully correlated)
// Instead of using the full covariance matrix, one can obtain similr results if the expected numnber of events
// in a bin is replaced by the background estimate, and the uncertainty on the efficiency is modified to account for
// the signal (see doc-db 36707 for a full explanation)


//An anonymous namespace for the fit function and its minimisation
namespace
{
    unsigned              npTot_(0),npFcn_(0),npBkg_(0),x0_(0),x1_(0);
    std::vector<double>   xvec_{},yvec_{};
    mu2e::CaloPulseUtil* pulseCachePtr_=(nullptr);

    constexpr double parMin(0.0), parMax(1e6);
    constexpr int    maxIterations(30);
    constexpr double initialDamping(100.0);

    double logn(double x, double *par) {return par[0]*pulseCachePtr_->evaluate(x-par[1]); }

    double fitfunction(double x, double *par)
    {
        double result(par[0]);
        for (unsigned i=npBkg_; i<npTot_; i+=npFcn_) result += logn(x,&par[i]);
        return result;
    }
    double fitfunctionPlot(double* x, double *par) {return fitfunction(x[0],par);}

    bool isTime(unsigned j) {return j >= npBkg_ && (j-npBkg_)%npFcn_ == 1;}

    // fit objective over [x0_,x1_): sum of (y-f)^2/baseline, with f built from the first npar parameters
    double objective(const std::vector<double>& par, unsigned npar)
    {
        double f(0);
        for (unsigned i=x0_;i<x1_;++i)
        {
            double val(par[0]);
            for (unsigned j=npBkg_; j+1<npar; j+=npFcn_) val += par[j]*pulseCachePtr_->evaluate(xvec_[i]-par[j+1]);
            if (std::abs(par[0]) > 1e-5) f += (yvec_[i]-val)*(yvec_[i]-val)/par[0];
        }
        return f;
    }

    // analytic gradient and Hessian of the objective (baseline above 1e-5); fixed parameters get identity rows
    void derivatives(const std::vector<double>& par, unsigned npar, const std::vector<bool>& free,
                     std::vector<double>& grad, std::vector<double>& hess)
    {
        grad.assign(npar,0.0);
        hess.assign(npar*npar,0.0);
        std::vector<double> jac(npar,0.0), slope(npar,0.0), rj(npar,0.0), cross(npar,0.0);
        double b(par[0]), sumSq(0);
        jac[0] = 1.0;
        for (unsigned i=x0_;i<x1_;++i)
        {
            double val(par[0]);
            for (unsigned j=npBkg_; j+1<npar; j+=npFcn_)
            {
                double p = pulseCachePtr_->evaluate(xvec_[i]-par[j+1],slope[j]);
                val     += par[j]*p;
                jac[j]   = p;
                jac[j+1] = -par[j]*slope[j];
            }
            double r = yvec_[i]-val;
            sumSq += r*r;
            for (unsigned a=0;a<npar;++a)
            {
                rj[a] += r*jac[a];
                for (unsigned c=0;c<=a;++c) hess[a*npar+c] += jac[a]*jac[c];
            }
            for (unsigned j=npBkg_; j+1<npar; j+=npFcn_) cross[j] += r*slope[j];
        }

        for (unsigned a=0;a<npar;++a)
        {
            grad[a] = -2.0*rj[a]/b - (a==0 ? sumSq/(b*b) : 0.0);
            for (unsigned c=0;c<=a;++c)
            {
                double h = 2.0*hess[a*npar+c]/b;
                if (c==0) h += 2.0*rj[a]/(b*b);
                if (a==0) h += 2.0*rj[c]/(b*b);
                if (a==0 && c==0) h += 2.0*sumSq/(b*b*b);
                if (isTime(a) && c+1==a) h += 2.0*cross[c]/b;
                hess[a*npar+c] = hess[c*npar+a] = h;
            }
        }

        for (unsigned a=0;a<npar;++a)
        {
            if (free[a]) continue;
            grad[a] = 0.0;
            for (unsigned c=0;c<npar;++c) hess[a*npar+c] = hess[c*npar+a] = 0.0;
            hess[a*npar+a] = 1.0;
        }
    }

    // solves m x = x (x holds the right-hand side) by Gaussian elimination with partial pivoting; false if singular
    bool solve(std::vector<double> m, std::vector<double>& x, unsigned n)
    {
        for (unsigned k=0;k<n;++k)
        {
            unsigned piv(k);
            for (unsigned i=k+1;i<n;++i) if (std::abs(m[i*n+k]) > std::abs(m[piv*n+k])) piv = i;
            if (!(std::abs(m[piv*n+k]) > 1e-300)) return false;
            if (piv != k)
            {
                for (unsigned c=0;c<n;++c) std::swap(m[k*n+c],m[piv*n+c]);
                std::swap(x[k],x[piv]);
            }
            for (unsigned i=k+1;i<n;++i)
            {
                double f = m[i*n+k]/m[k*n+k];
                for (unsigned c=k;c<n;++c) m[i*n+c] -= f*m[k*n+c];
                x[i] -= f*x[k];
            }
        }
        for (unsigned k=n;k-->0;)
        {
            for (unsigned c=k+1;c<n;++c) x[k] -= m[k*n+c]*x[c];
            x[k] /= m[k*n+k];
        }
        return true;
    }

    // Levenberg-Marquardt minimisation within [parMin,parMax]; each step stays in a trust region (baseline and
    // amplitudes by half their value, times by 5 ns) so the iterates do not collapse to a baseline-only solution
    bool minimize(std::vector<double>& par, unsigned npar, const std::vector<bool>& free,
                  std::vector<double>& err, double& fmin)
    {
        std::vector<double> grad, hess, step(npar,0.0), trial(par);
        double lambda(initialDamping);
        fmin = objective(par,npar);

        for (int it=0; it<maxIterations; ++it)
        {
            derivatives(par,npar,free,grad,hess);
            std::vector<double> damped(hess);
            for (unsigned j=0;j<npar;++j) damped[j*npar+j] += lambda*(std::abs(hess[j*npar+j])+1e-9);
            for (unsigned j=0;j<npar;++j) step[j] = -grad[j];
            if (!solve(damped,step,npar)) {lambda *= 10.0; continue;}

            for (unsigned j=0;j<npar;++j)
            {
                if (!free[j]) {trial[j] = par[j]; continue;}
                double cap = isTime(j) ? 5.0 : 0.5*std::abs(par[j]) + (j==0 ? 1.0 : 10.0);
                trial[j] = std::clamp(par[j]+std::clamp(step[j],-cap,cap), parMin, parMax);
            }

            double ftrial = objective(trial,npar);
            if (std::abs(trial[0]) > 1e-5 && ftrial <= fmin)
            {
                bool converged = (fmin-ftrial) <= 1e-12*fmin;
                par     = trial;
                fmin    = ftrial;
                lambda *= 0.3;
                if (converged) break;
            }
            else lambda *= 10.0;
        }

        err.assign(npar,0.0);
        if (!(std::abs(par[0]) > 1e-5)) return false;
        derivatives(par,npar,free,grad,hess);
        for (unsigned j=0;j<npar;++j)
        {
            if (!free[j]) continue;
            std::vector<double> col(npar,0.0);
            col[j] = 1.0;
            if (!solve(hess,col,npar)) return false;
            err[j] = std::sqrt(std::abs(2.0*col[j]));
        }
        return true;
    }
}







namespace mu2e {


   CaloTemplateWFUtil::CaloTemplateWFUtil(const CaloPulseUtil::Config& configPulseCache, double minPeakAmplitude,
                                          double minDTPeaks, int printLevel) :
      pulseCache_(CaloPulseUtil(configPulseCache)),
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
   void   CaloTemplateWFUtil::initialize ()                                                                 {pulseCache_.buildCache();}
   void   CaloTemplateWFUtil::reset      ()                                                                 {param_.clear(); paramErr_.clear(); nParTot_=0; npTot_ = 0;}
   void   CaloTemplateWFUtil::setXYVector(const std::vector<double>& xvec, const std::vector<double>& yvec) {xvec_ = xvec; yvec_ = yvec; x0_=0; x1_ = xvec_.size();}
   void   CaloTemplateWFUtil::setPar     (const std::vector<double>& par)                                   {param_ = par; nParTot_ = npTot_ = par.size();}

   //-----------------------------------------------------------------------------------------------------
   void CaloTemplateWFUtil::fit()
   {
       status_ = 0;
       if (param_.empty() || param_.size()>49 || xvec_.empty()) return;
       if (nParTot_ < nParBkg_  || (nParTot_-nParBkg_)%nParFcn_ !=0) return;

       std::vector<double> par(param_), err;
       for (auto& val : par) val = std::clamp(val, parMin, parMax);
       std::vector<bool> free(nParTot_,true);
       double fmin(0);

       // Perform first fit with initial model
       //
       bool ok = minimize(par, nParTot_, free, err, fmin);


       // Remove small or "duplicate" components and redo the fit with simplified model if there is more than one peak
       // A duplicate peak is defined as a peak shortly after a previous peak with a smaller amplitude
       //
       if (nParTot_ > nParFcn_+nParBkg_)
       {
           bool refit(false);
           std::vector<double> tempPar(par),tempErr(err);

           for (unsigned ip=nParBkg_; ip<nParTot_; ip += nParFcn_)
           {
               if (selectComponent(tempPar,tempErr,ip)) continue;
               par[ip]    = par[ip+1]  = 0.0;
               free[ip]   = free[ip+1] = false;
               refit = true;
           }

           if (refit) ok = minimize(par, nParTot_, free, err, fmin);
       }


       // Save the results - exclude low components
       //
       param_.clear();
       paramErr_.clear();

       unsigned i(0);
       while (i<nParTot_)
       {
           //if the amplitude is too small, jump to the next peak
           if (par[i]<1 && i >=nParBkg_ && (i-nParBkg_)%nParFcn_==0) {i+=nParFcn_;continue;}

           param_.push_back(par[i]);
           paramErr_.push_back(err[i]);
           ++i;
       }

       chi2_    = fmin;
       nParTot_ = param_.size();
       npTot_   = nParTot_;
       status_  = ok ? 3 : 0;
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

       // fit the baseline and the first peak only
       unsigned npar(nParBkg_+nParFcn_);
       std::vector<double> par(param_.begin(),param_.begin()+npar), err;
       for (auto& val : par) val = std::clamp(val, parMin, parMax);
       std::vector<bool> free(npar,true);
       double fmin(0);
       bool ok = minimize(par, npar, free, err, fmin);

       param_[nParBkg_+1]    = par[nParBkg_+1];
       paramErr_[nParBkg_+1] = err[nParBkg_+1];
       status_               = ok ? 3 : 0;

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
