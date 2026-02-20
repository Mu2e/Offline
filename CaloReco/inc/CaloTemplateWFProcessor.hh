#ifndef CaloTemplateWFProcessor_HH
#define CaloTemplateWFProcessor_HH

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

#include "Offline/CaloReco/inc/CaloWaveformProcessor.hh"
#include "Offline/CaloReco/inc/CaloTemplateWFUtil.hh"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "TH2.h"
#include <vector>



namespace mu2e {

  class CaloTemplateWFProcessor : public CaloWaveformProcessor
  {
     public:
        struct Config
        {
            using Name    = fhicl::Name;
            using Comment = fhicl::Comment;
            fhicl::Atom<std::string> pulseFileName     { Name("pulseFileName"),    Comment("Calo pulse file name") };
            fhicl::Atom<std::string> pulseHistName     { Name("pulseHistName"),    Comment("Calo pulse hist name") };
            fhicl::Atom<unsigned>    windowPeak        { Name("windowPeak"),       Comment("Number of bins around central value to inspect") };
            fhicl::Atom<double>      minPeakAmplitude  { Name("minPeakAmplitude"), Comment("Minimum peak amplitude") };
            fhicl::Atom<double>      minDTPeaks        { Name("minDTPeaks"),       Comment("Minimum time difference between consecutive peaks") };
            fhicl::Atom<unsigned>    numNoiseBins      { Name("numNoiseBins"),     Comment("Number of bins to estimate noise") };
            fhicl::Atom<double>      psdThreshold      { Name("psdThreshold"),     Comment("Pulse shape discrimination threshold for secondary peaks") };
            fhicl::Atom<double>      chiThreshold      { Name("chiThreshold"),     Comment("Min chi2 for refit strategy") };
            fhicl::Atom<bool>        refitLeadingEdge  { Name("refitLeadingEdge"), Comment("Refit the leading edge to extract peak time") };
            fhicl::Atom<double>      digiSampling      { Name("digiSampling"),     Comment("Digitization time sampling") };
            fhicl::Atom<int>         fitPrintLevel     { Name("fitPrintLevel"),    Comment("minuit fit print level") };
            fhicl::Atom<int>         fitStrategy       { Name("fitStrategy"),      Comment("Minuit fit strategy") };
            fhicl::Atom<int>         diagLevel         { Name("diagLevel"),        Comment("Diagnosis level") };
        };


        CaloTemplateWFProcessor(const Config& config);

        virtual void     initialize  () override;
        virtual void     reset       () override;
        virtual void     extract     (const std::vector<double>& xInput, const std::vector<double>& yInput) override;
        virtual void     plot        (const std::string& pname) const override;

        virtual int      nPeaks      ()               const override {return resAmp_.size();}
        virtual double   chi2        ()               const override {return chi2_;}
        virtual int      ndf         ()               const override {return ndf_;}
        virtual double   amplitude   (unsigned int i) const override {return resAmp_.at(i);}
        virtual double   amplitudeErr(unsigned int i) const override {return resAmpErr_.at(i);}
        virtual double   time        (unsigned int i) const override {return resTime_.at(i);}
        virtual double   timeErr     (unsigned int i) const override {return resTimeErr_.at(i);}
        virtual bool     isPileUp    (unsigned int i) const override {return i > 1;}   // resAmp_.size() > 1 as alternative?


    private:
       void   initHistos         ();
       void   setPrimaryPeakPar1 (const std::vector<double>& xvec, const std::vector<double>& yvec);
       void   setPrimaryPeakPar2 (const std::vector<double>& xvec, const std::vector<double>& yvec);
       void   findRisingPeak     (int ipeak, std::vector<double>& parInit, const std::vector<double>& xvec, const std::vector<double>& yvec, std::vector<double>& ywork);
       void   setSecondaryPeakPar(const std::vector<double>& xvec, const std::vector<double>& yvec);
       double estimatePeakTime   (const std::vector<double>& xvec, const std::vector<double>& ywork, int ic);
       bool   checkPeakDist      (double x0);
       void   dump               (const std::string& name, const std::vector<double>& val) const;

       unsigned            windowPeak_ ;
       double              minPeakAmplitude_;
       unsigned            numNoiseBins_;
       double              minDTPeaks_;
       double              psdThreshold_;
       double              chiThreshold_;
       bool                refitLeadingEdge_;
       int                 diagLevel_;
       CaloTemplateWFUtil  fmutil_;
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
       TH2F* _hchi2Peak;
  };

}
#endif
