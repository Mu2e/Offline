#ifndef CaloNoiseSimGenerator_HH
#define CaloNoiseSimGenerator_HH
//
// Generate long noise waveform to use for calorimeter digitization
//
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "Offline/SeedService/inc/SeedService.hh"

#include "Offline/CaloMC/inc/CaloWFExtractor.hh"
#include "Offline/Mu2eUtilities/inc/CaloPulseShape.hh"

#include "CLHEP/Random/RandPoissonQ.h"
#include "CLHEP/Random/RandGaussQ.h"
#include "CLHEP/Random/RandFlat.h"


namespace mu2e {

  class CaloNoiseSimGenerator
  {
     public:
        struct Config
        {
            using Name    = fhicl::Name;
            using Comment = fhicl::Comment;
            fhicl::Atom<std::string> pulseFileName  { Name("pulseFileName"),  Comment("Calo pulse file name") };
            fhicl::Atom<std::string> pulseHistName  { Name("pulseHistName"),  Comment("Calo pulse hist name") };
            fhicl::Atom<double>      elecNphotPerNs { Name("elecNphotPerNs"), Comment("Electronics noise number of PE / ns ") };
            fhicl::Atom<double>      rinNphotPerNs  { Name("rinNphotPerNs"),  Comment("RIN noise number of PE / ns ") };
            fhicl::Atom<double>      darkNphotPerNs { Name("darkNphotPerNs"), Comment("SiPM Dark noise number of PE / ns ") };
            fhicl::Atom<double>      digiSampling   { Name("digiSampling"),   Comment("Digitization time sampling") };
            fhicl::Atom<double>      pePerMeV       { Name("readoutPEPerMeV"),Comment("Number of pe / MeV for Readout") };
            fhicl::Atom<double>      ADCToMeV       { Name("ADCToMeV"),       Comment("ADC to MeV conversion factor") };
            fhicl::Atom<unsigned>    noiseWFSize    { Name("noiseWFSize"),    Comment("Noise WF size") };
            fhicl::Atom<unsigned>    nMaxFragment   { Name("nMaxFragment"),   Comment("maximum number of wf generated for extracting noise fragments ") };
            fhicl::Atom<int>         minPeakADC     { Name("minPeakADC"),     Comment("Minimum ADC hits of local peak to digitize") };
            fhicl::Atom<int>         diagLevel      { Name("diagLevel"),      Comment("Diag Level"),0 };
        };


        CaloNoiseSimGenerator(const Config& config, CLHEP::HepRandomEngine& engine, int iRO);

        void                         initialize(const CaloWFExtractor& wfExtractor);
        void                         refresh();

        void                         addSampleNoise(std::vector<double>& wfVector, unsigned istart, unsigned ilength);
        void                         addSaltAndPepper(std::vector<double>& wfVector);
        void                         plotNoise(const std::string& name);

        const std::vector<double>&   noise()    const {return waveform_;}
        double                       pedestal() const {return pedestal_;}


     private:
        using vvd = std::vector<std::vector<double>>;

        void                  generateWF(std::vector<double>& wfVector);
        void                  generateFragments(const CaloWFExtractor& wfExtractor);

        unsigned              iRO_;
        std::vector<double>   waveform_;
        int                   pedestal_;
        vvd                   digiNoise_;
        double                digiNoiseProb_;
        double                digiSampling_;
        double                noiseRinDark_;
        double                noiseElec_;
        double                minPeakADC_;
        double                pePerMeV_;
        double                MeVToADC_;
        CLHEP::RandPoissonQ   randPoisson_;
        CLHEP::RandGaussQ     randGauss_;
        CLHEP::RandFlat       randFlat_;
        unsigned              nMaxFragment_;
        CaloPulseShape        pulseShape_;
        int                   diagLevel_;
   };

}
#endif
