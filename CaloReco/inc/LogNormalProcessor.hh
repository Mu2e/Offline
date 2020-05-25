#ifndef LogNormalProcessor_HH
#define LogNormalProcessor_HH


#include "CaloReco/inc/WaveformProcessor.hh"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include <vector>
#include "TH1.h"


namespace mu2e {


  class LogNormalProcessor : public WaveformProcessor {
     

     public:
     
        enum shapeFixMode {All=0, None=1, Close=2};

        struct Config {
          using Name    = fhicl::Name;
          using Comment = fhicl::Comment;        
          fhicl::Atom<int>    windowPeak{       Name("windowPeak"),       Comment("Number of bins around central vlue to inspect")};
          fhicl::Atom<double> minPeakAmplitude{ Name("minPeakAmplitude"), Comment("Minimum peak amplitude")};
          fhicl::Atom<bool>   fixShapeSig{      Name("fixShapeSig"),      Comment("Fix pulse shape to template")};
          fhicl::Atom<double> psdThreshold{     Name("psdThreshold"),     Comment("Pulse Shape discrimination threshold")};          
          fhicl::Atom<int>    pulseHighBuffer{  Name("pulseHighBuffer"),  Comment("Buffer after last bin of waveform")};
          fhicl::Atom<double> timeFraction{     Name("timeFraction"),     Comment("CFD for timing measurement")};
          fhicl::Atom<double> shiftTime{        Name("shiftTime"),        Comment("Time bwtween beginning and maximum value of pusle")};
          fhicl::Atom<int>    fitPrintLevel{    Name("fitPrintLevel"),    Comment("minuit fit print level")};
          fhicl::Atom<int>    fitStrategy{      Name("fitStrategy"),      Comment("Minuit fit strategy")};
          fhicl::Atom<int>    diagLevel{        Name("diagLevel"),        Comment("Diagnosis level")};
        };


                    LogNormalProcessor(const Config& config);
        virtual    ~LogNormalProcessor() {};	       

        virtual void   initialize();
	virtual void   reset();
        virtual void   extract(const std::vector<double> &xInput, const std::vector<double> &yInput);

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
       int                 fitPrintLevel_;
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
