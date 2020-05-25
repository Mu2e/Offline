#ifndef FixedFastProcessor_HH
#define FixedFastProcessor_HH

#include "CaloReco/inc/WaveformProcessor.hh"
#include "CaloReco/inc/CaloPulseCache.hh"
#include "CaloReco/inc/FastFixedUtil.hh"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "TH2.h"

#include <vector>

namespace mu2e {


  class FixedFastProcessor : public WaveformProcessor 
  {
     public:
        struct Config
        {
           using Name = fhicl::Name;
           using Comment = fhicl::Comment;        
           fhicl::Atom<unsigned>     windowPeak      { Name("windowPeak"),       Comment("Number of bins around central vlue to inspect") };
           fhicl::Atom<double>       minPeakAmplitude{ Name("minPeakAmplitude"), Comment("Minimum peak amplitude") };
           fhicl::Atom<double>       psdThreshold    { Name("psdThreshold"),     Comment("Pulse Shape discrimination threshold") };
           fhicl::Atom<unsigned int> pulseLowBuffer  { Name("pulseLowBuffer"),   Comment("Buffer before first bin of waveform") };
           fhicl::Atom<unsigned int> pulseHighBuffer { Name("pulseHighBuffer"),  Comment("Buffer after last bin of waveform") };
           fhicl::Atom<unsigned int> minDiffTime     { Name("minDiffTime"),      Comment("Minimum time difference between two consecutive peaks") };
           fhicl::Atom<bool>         refitLeadingEdge{ Name("refitLeadingEdge"), Comment("Refit the leading edge to extract peak time") };
           fhicl::Atom<double>       shiftTime       { Name("shiftTime"),        Comment("Time bwtween beginning and maximum value of pusle") };
           fhicl::Atom<int>          fitPrintLevel   { Name("fitPrintLevel"),    Comment("minuit fit print level") };
           fhicl::Atom<int>          fitStrategy     { Name("fitStrategy"),      Comment("Minuit fit strategy") };
           fhicl::Atom<int>          diagLevel       { Name("diagLevel"),        Comment("Diagnosis level") };
        };

     
                    FixedFastProcessor(const Config& config);
        virtual    ~FixedFastProcessor() {};	

        virtual void   initialize();
        virtual void   reset();
        virtual void   extract(const std::vector<double>& xInput, const std::vector<double>& yInput);

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
       void   findPeak(const std::vector<double>& xvec, const std::vector<double>& yvec, std::vector<double>& parInit, std::vector<unsigned>& xindices);
       double meanParabol(double x1, double x2, double x3, double y1, double y2, double y3);
       void   buildXRange(const std::vector<unsigned>& peakLoc, const std::vector<double>& xvec, std::vector<unsigned>& xindices);
       void   refitTime(const std::vector<double>& xvec, const std::vector<double>& yvec, std::vector<double>& parFinal);
       
       unsigned            windowPeak_ ;
       double              minPeakAmplitude_;
       double              psdThreshold_;
       unsigned            pulseLowBuffer_;
       unsigned            pulseHighBuffer_;
       unsigned            minDiffTime_;
       bool                refitLeadingEdge_;
       double              shiftTime_;
       int                 printLevel_;
       int                 fitStrategy_;
       int                 diagLevel_;       
       CaloPulseCache      pulseCache_;
       unsigned            nParFcn_;
       FastFixedUtil       fmutil_;
                            
       int                 nPeaks_;
       double              chi2_;
       int                 ndf_;
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
       TH1F* _hpsd;
  };

}
#endif
