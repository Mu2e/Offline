#ifndef STMRawWFProcessor_HH
#define STMRawWFProcessor_HH


#include "Offline/STMReco/inc/STMWaveformProcessor.hh"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include <vector>


namespace mu2e {


  class STMRawWFProcessor : public STMWaveformProcessor {


     public:

        struct Config
        {
           using Name = fhicl::Name;
           using Comment = fhicl::Comment;
           fhicl::Atom<int>    windowPeak{       Name("windowPeak"),       Comment("Number of bins around central vlue to inspect")};
           fhicl::Atom<double> minPeakAmplitude{ Name("minPeakAmplitude"), Comment("Minimum peak amplitude")};
           fhicl::Atom<double> shiftTime{        Name("shiftTime"),        Comment("Time between beginning and maximum value of pusle")};
           fhicl::Atom<double> scaleFactor{      Name("scaleFactor"),      Comment("Factor to convert bin height to signal amplitude")};
           fhicl::Atom<int>    diagLevel{        Name("diagLevel"),        Comment("Diagnosis level")};
        };


        STMRawWFProcessor(const Config& config);

        virtual void   initialize() override;
        virtual void   reset() override;
        virtual void   extract(const std::vector<double> &xInput, const std::vector<double> &yInput) override;
        virtual void   plot   (const std::string& pname) const override;

        virtual int    nPeaks()                     const override {return nPeaks_;}
        virtual double chi2()                       const override {return chi2_;}
        virtual int    ndf()                        const override {return ndf_;}
        virtual double amplitude(unsigned int i)    const override {return resAmp_.at(i);}
        virtual double amplitudeErr(unsigned int i) const override {return resAmpErr_.at(i);}
        virtual double time(unsigned int i)         const override {return resTime_.at(i);}
        virtual double timeErr(unsigned int i)      const override {return resTimeErr_.at(i);}
        virtual bool   isPileUp(unsigned int i)     const override {return nPeaks_ > 1;}

    private:

       int                 windowPeak_;
       double              minPeakAmplitude_;
       double              shiftTime_;
       double              scaleFactor_;
       int                 diagLevel_;

       int                 nPeaks_;
       double              chi2_;
       int                 ndf_;
       std::vector<double> res_;
       std::vector<double> resAmp_;
       std::vector<double> resAmpErr_;
       std::vector<double> resTime_;
       std::vector<double> resTimeErr_;
       std::vector<double> xvec_;
       std::vector<double> yvec_;

  };

}
#endif
