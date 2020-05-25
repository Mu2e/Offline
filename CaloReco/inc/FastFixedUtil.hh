#ifndef FastFixedUtil_HH
#define FastFixedUtil_HH

#include "CaloReco/inc/CaloPulseCache.hh"
#include <vector>


namespace mu2e {


  class FastFixedUtil  {
     
     public:
     
        FastFixedUtil(unsigned nParFcn, double minPeakAmplitude, double minDiffTime, int printLevel, int fitStrategy, int diagLevel);
        ~FastFixedUtil() {};
	
        void   initialize();
        void   setXYVector(const std::vector<double>& xvec, const std::vector<double>& yvec);
        void   fitMinuit(const std::vector<double>& parInit);
        void   fitNewton(double tminInit);
        double refineMin(double t0init, double stepInit);
        
        double eval_logn(double x, int ioffset, const std::vector<double>& par);
        double eval_fcn(double x, const std::vector<double>& par);
        void   plot(std::string pname, const std::vector<double>& param);
 
        double                      chi2()        const {return chi2_;}
        const std::vector<double>&  sfpar()       const {return sfpar_;}
        const std::vector<double>&  esfpar()      const {return esfpar_;}


    private:
              
      double meanDistance(unsigned ip, const std::vector<double>& tempPar);
      double calcAlpha(double testTime);
      double calcChi2(double testTime, double alpha=0.0);


      CaloPulseCache      pulseCache_;
      int                 diagLevel_;
      int                 printLevel_;
      int                 fitStrategy_;
      double              minPeakAmplitude_;
      double              minDiffTime_;
      bool                isFitDone_;   
      std::vector<double> sfpar_;
      std::vector<double> esfpar_;
      double              chi2_;
  };

}
#endif
