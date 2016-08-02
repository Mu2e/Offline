#ifndef LogNormalProcessor_HH
#define LogNormalProcessor_HH


#include "CaloReco/inc/WaveformProcessor.hh"
#include "fhiclcpp/ParameterSet.h"
#include <vector>
#include "TH1.h"


namespace mu2e {


  class LogNormalProcessor : public WaveformProcessor {
     

     public:
     
        enum shapeFixMode {All=0, None=1, Close=2};

                    LogNormalProcessor(fhicl::ParameterSet const& param);
        virtual    ~LogNormalProcessor() {};
	

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


        virtual void plot(std::string pname);
	


    private:
       
       int                 windowPeak_;
       double              minPeakAmplitude_;
       bool                fixShapeSig_;
       double              psdThreshold_;
       int                 pulseHighBuffer_;
       double              timeFraction_;
       double              shiftTime_;
       int                 printLevel_;
       int                 fitStrategy_;
       int                 diagLevel_;
       
       double              peakFactor_;
              
       int                 nPeaks_;
       double              chi2_;
       int                 ndf_;
       std::vector<double> res_;
       std::vector<double> resAmp_;
       std::vector<double> resAmpErr_;
       std::vector<double> resTime_;
       std::vector<double> resTimeErr_;

    
       void      findPeak(double* parInit);
       void      doFit(double* parInit, double *sfpar, double *errsfpar, double& chi2);
       double    findTime(int ipeak);
       double    meanParabol(int i1, int i2, int i3);

       TH1F* _hTime;
       TH1F* _hTimeErr;
       TH1F* _hEner;
       TH1F* _hEnerErr;
       TH1F* _hChi2;
       TH1F* _hNpeak;
       TH1F* _hRescale;       
       TH1F* _hDelta;

  };

}
#endif
