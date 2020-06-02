#ifndef TemplateProcessor_HH
#define TemplateProcessor_HH

// This is an signal extraction method based on the waveform template. 
//
// Each peak in the waveform is described by two parameters: amplitide and peak time
// For a single peak, the amplitude can be found analytically for a given start time, and a 
// quasi-Netwon method can be used to fit the waveform. 
// If there are more than one peak, we use a generic gradient descent method, namely minuit.
//
// There is an additional option to refit the leding edge of the first peak to improve 
// timing accuracy
//
// A peak inside a waveform is defined to be a pile-up if it is on top of another peak. 
// --> the first peak is never a pile-up
//

#include "CaloReco/inc/WaveformProcessor.hh"
#include "CaloReco/inc/TemplateUtil.hh"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "TH2.h"
#include <vector>



namespace mu2e {

  class TemplateProcessor : public WaveformProcessor 
  {
     public:
        struct Config
        {
            using Name    = fhicl::Name;
            using Comment = fhicl::Comment;        
            fhicl::Atom<unsigned>     windowPeak        { Name("windowPeak"),        Comment("Number of bins around central value to inspect") };
            fhicl::Atom<double>       minPeakAmplitude  { Name("minPeakAmplitude"),  Comment("Minimum peak amplitude") };
            fhicl::Atom<unsigned int> pulseLowBuffer    { Name("pulseLowBuffer"),    Comment("Buffer before first bin of waveform") };
            fhicl::Atom<unsigned int> pulseHighBuffer   { Name("pulseHighBuffer"),   Comment("Buffer after last bin of waveform") };
            fhicl::Atom<unsigned int> minDeltaPeakBin   { Name("minDeltaPeakBin"),   Comment("Minimum number of waveform bins between two consecutive peaks") };
            fhicl::Atom<bool>         refitLeadingEdge  { Name("refitLeadingEdge"),  Comment("Refit the leading edge to extract peak time") };
            fhicl::Atom<double>       timeCorr          { Name("timeCorr"),          Comment("Time between beginning and maximum value of pusle") };
            fhicl::Atom<double>       digiSampling      { Name("digiSampling"),      Comment("Digitization time sampling") }; 
            fhicl::Atom<int>          pulseIntegralSteps{ Name("pulseIntegralSteps"),Comment("Numer of time sub-division for CaloPulseChape") }; 
            fhicl::Atom<int>          fitPrintLevel     { Name("fitPrintLevel"),     Comment("minuit fit print level") };
            fhicl::Atom<int>          fitStrategy       { Name("fitStrategy"),       Comment("Minuit fit strategy") };
            fhicl::Atom<int>          diagLevel         { Name("diagLevel"),         Comment("Diagnosis level") };
        };

     
        TemplateProcessor(const Config& config);

        virtual void   initialize() override;
        virtual void   reset() override;
        virtual void   extract(const std::vector<double>& xInput, const std::vector<double>& yInput) override;
        virtual void   plot(std::string pname) const override;

        virtual int    nPeaks()                     const override {return fmutil_.nPeaks();}
        virtual double chi2()                       const override {return chi2_;}
        virtual int    ndf()                        const override {return ndf_;}
        virtual double amplitude(unsigned int i)    const override {return resAmp_.at(i);}
        virtual double amplitudeErr(unsigned int i) const override {return resAmpErr_.at(i);}
        virtual double time(unsigned int i)         const override {return resTime_.at(i);}
        virtual double timeErr(unsigned int i)      const override {return resTimeErr_.at(i);}  
        virtual bool   isPileUp(unsigned int i)     const override {return i > 1;}   // fmutil_.nPeaks() > 1 as alternative?


    private:
       void   findPeak        (const std::vector<double>& xvec, const std::vector<double>& yvec, std::vector<double>& parInit, std::vector<unsigned>& xindices);
       double meanParabol     (double x1, double x2, double x3, double y1, double y2, double y3);
       void   buildXRange     (const std::vector<unsigned>& peakLoc, const std::vector<double>& xvec, std::vector<unsigned>& xindices);
       void   refitLeadingEdge(const std::vector<double>& xvec, const std::vector<double>& yvec);
      
       unsigned            windowPeak_ ;
       double              minPeakAmplitude_;
       unsigned            pulseLowBuffer_;
       unsigned            pulseHighBuffer_;
       unsigned            minDeltaPeakBin_;
       bool                refitLeadingEdge_;
       double              timeCorr_;
       int                 diagLevel_;       
       TemplateUtil        fmutil_;                           
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
  };

}
#endif
