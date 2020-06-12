#ifndef TemplateUtil_HH
#define TemplateUtil_HH

#include "Mu2eUtilities/inc/CaloPulseShape.hh"
#include <vector>
#include <string>


namespace mu2e {

  class TemplateUtil  {
     
     public:     
        TemplateUtil(double minPeakAmplitude, double digiSampling);
        
        void                        initialize    (); 
	void                        setXYVector   (const std::vector<double>& xvec, const std::vector<double>& yvec);
	void                        setPar        (const std::vector<double>& par);
	void                        reset         ();
        
        void                        fit           ();
        void                        refitMin      ();        
        void                        plotFit       (const std::string& pname) const;
        double                      eval_fcn      (double x); 
        double                      eval_logn     (double x, int ioffset);  
        void                        calcTimeCorr  (std::vector<double>& par);
 
	void                        setStrategy   (int val) {fitStrategy_ = val;}
	void                        setPrintLevel (int val) {printLevel_  = val;}
	void                        setFitStartegy(int val) {fitStrategy_ = val;}
	void                        setDiagLevel  (int val) {diagLevel_   = val;}
        
        unsigned                    status        ()                const {return status_;}
        double                      chi2          ()                const {return chi2_;}
        const std::vector<double>&  par           ()                const {return param_;}
        const std::vector<double>&  parErr        ()                const {return paramErr_;}
	unsigned                    nParFcn       ()                const {return nParFcn_;}
	unsigned                    nPeaks        ()                const {return param_.size()/nParFcn_;}
	double                      fromPeakToT0  (double timePeak) const {return pulseCache_.fromPeakToT0(timePeak);} 

     private:              
        void   fitMinuit      ();
        void   fitNewton      ();
	double calcAlpha      (double testTime);
	double calcChi2       (double testTime, double alpha=0.0);
        double refineMin      (double tinit,    double stepInit); 
        bool   selectComponent(const std::vector<double>& tempPar, unsigned ip);       

	CaloPulseShape      pulseCache_;
	double              minPeakAmplitude_;
	int                 fitStrategy_;
	int                 diagLevel_;
	int                 printLevel_;
	std::vector<double> param_;
	std::vector<double> paramErr_;
	unsigned            nParTot_;
	unsigned            nParFcn_;
	double              chi2_;
	unsigned            status_;

  };
  
}
#endif
