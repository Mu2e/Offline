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
            fhicl::Atom<unsigned> windowPeak        { Name("windowPeak"),       Comment("Number of bins around central value to inspect") };
            fhicl::Atom<double>   minPeakAmplitude  { Name("minPeakAmplitude"), Comment("Minimum peak amplitude") };
            fhicl::Atom<double>   minDTPeaks        { Name("minDTPeaks"),       Comment("Minimum time difference between consecutive peaks") };
            fhicl::Atom<unsigned> numNoiseBins      { Name("numNoiseBins"),     Comment("Number of bins to estimate noise") };
            fhicl::Atom<unsigned> minDeltaPeakBin   { Name("minDeltaPeakBin"),  Comment("Minimum number of waveform bins between two consecutive peaks") };
            fhicl::Atom<bool>     doSecondaryPeak   { Name("doSecondaryPeak"),  Comment("Extract secondary peaks") }; 
            fhicl::Atom<double>   psdThreshold      { Name("psdThreshold"),     Comment("Pulse shape discrimination threshold for secondary peaks") }; 
            fhicl::Atom<bool>     refitLeadingEdge  { Name("refitLeadingEdge"), Comment("Refit the leading edge to extract peak time") };
            fhicl::Atom<double>   timeCorr          { Name("timeCorr"),         Comment("Time between beginning and maximum value of pusle") };
            fhicl::Atom<double>   digiSampling      { Name("digiSampling"),     Comment("Digitization time sampling") }; 
            fhicl::Atom<int>      fitPrintLevel     { Name("fitPrintLevel"),    Comment("minuit fit print level") };
            fhicl::Atom<int>      fitStrategy       { Name("fitStrategy"),      Comment("Minuit fit strategy") };
            fhicl::Atom<int>      diagLevel         { Name("diagLevel"),        Comment("Diagnosis level") };
        };

     
        TemplateProcessor(const Config& config);

        virtual void   initialize() ;
        virtual void   reset() ;
        virtual void   extract(const std::vector<double>& xInput, const std::vector<double>& yInput) ;
        virtual void   plot   (const std::string& pname) const ;

        virtual int      nPeaks()                     const  {return resAmp_.size();}
        virtual double   chi2()                       const  {return chi2_;}
        virtual int      ndf()                        const  {return ndf_;}
        virtual double   amplitude(unsigned int i)    const  {return resAmp_.at(i);}
        virtual double   amplitudeErr(unsigned int i) const  {return resAmpErr_.at(i);}
        virtual double   time(unsigned int i)         const  {return resTime_.at(i);}
        virtual double   timeErr(unsigned int i)      const  {return resTimeErr_.at(i);}  
        virtual bool     isPileUp(unsigned int i)     const  {return i > 1;}   // resAmp_.size() > 1 as alternative?

    private:
       void   setPrimaryPeakPar  (const std::vector<double>& xvec, const std::vector<double>& yvec, std::vector<unsigned>& peakLocation);
       void   setSecondaryPeakPar(const std::vector<double>& xvec, const std::vector<double>& yvec, std::vector<unsigned>& peakLocation);
       double estimatePeakTime   (double x1, double x2, double x3, double y1, double y2, double y3);
       bool   checkPeakDist      (unsigned i, const std::vector<unsigned>& pealLocation);		              
       void   dump               (const std::string& name, const std::vector<double>& val) const; 

       unsigned            windowPeak_ ;
       double              minPeakAmplitude_;
       unsigned            numNoiseBins_;
       unsigned            minDeltaPeakBin_;
       bool                doSecondaryPeak_;
       double              psdThreshold_;
       bool                refitLeadingEdge_;
       double              timeCorr_;
       int                 diagLevel_;       
       TemplateUtil        fmutil_;                           
       double              chi2_;
       int                 ndf_;
       std::vector<double> shiftPar_;
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
