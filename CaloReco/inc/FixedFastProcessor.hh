#ifndef FixedFastProcessor_HH
#define FixedFastProcessor_HH


#include "CaloReco/inc/WaveformProcessor.hh"
#include "CaloReco/inc/CaloPulseCache.hh"
#include "fhiclcpp/ParameterSet.h"
#include <vector>
#include "TH2.h"


namespace mu2e {


  class FixedFastProcessor : public WaveformProcessor {
     

     public:
     
                    FixedFastProcessor(fhicl::ParameterSet const& param);
        virtual    ~FixedFastProcessor() {};
	

        virtual void   initialize();
        virtual void   reset();
        virtual void   extract(std::vector<double> &xInput, std::vector<double> &yInput);

        virtual int    nPeaks()                     const {return nPeaks_;}
        virtual double chi2()                       const {return chi2_;}
        virtual int    ndf()                        const {return ndf_;}
        virtual double amplitude(unsigned int i)    const {return resAmp_.at(i);}
        virtual double amplitudeErr(unsigned int i) const {return resAmpErr_.at(i);}
        virtual double time(unsigned int i)         const {return resTime_.at(i);}
        virtual double timeErr(unsigned int i)      const {return resTimeErr_.at(i);}  
        virtual bool   isPileUp(unsigned int i)     const {return nPeaks_ > 1;}


        virtual void   plot(std::string pname);
	


    private:
       
       int                 windowPeak_ ;
       double              minPeakAmplitude_;
       double              psdThreshold_;
       unsigned int        pulseLowBuffer_;
       unsigned int        pulseHighBuffer_;
       unsigned int        minDiffTime_;
       double              shiftTime_;
       int                 printLevel_;
       int                 fitStrategy_;
       int                 diagLevel_;
       
       CaloPulseCache      pulseCache_;                     
       int                 nPeaks_;
       double              chi2_;
       int                 ndf_;
       std::vector<double> res_;
       std::vector<double> resAmp_;
       std::vector<double> resAmpErr_;
       std::vector<double> resTime_;
       std::vector<double> resTimeErr_;
       

       TH1F* _hTime;
       TH1F* _hTimeErr;
       TH1F* _hEner;
       TH1F* _hEnerErr;
       TH1F* _hChi2;
       TH1F* _hNpeak;
       TH1F* _hRescale;       
       TH2F* _hchi2Amp;
       TH1F* _hDelta;

    
       void   findPeak(double* parInit);
       void   buildXRange(std::vector<unsigned int>& peakLoc);
       void   doFit(double* parInit, double *sfpar, double *errsfpar, double& chi2);
       void   doNewton(double* parInit, double *sfpar, double *errsfpar, double& chi2);
       double refineMin(int& nTry,double tmin, double step);

       double meanParabol(unsigned int i1, unsigned int i2, unsigned int i3);
       double calcAlpha(double testTime);
       double calcChi2(double testTime, double alpha = -1.0);
       double refitTime(double tmin, double alpha);

  };

}
#endif
