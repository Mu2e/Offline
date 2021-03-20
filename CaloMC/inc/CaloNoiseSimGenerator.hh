#ifndef CaloNoiseSimGenerator_HH
#define CaloNoiseSimGenerator_HH
//
// Generate long noise waveform to use for calorimeter digitization 
//
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "SeedService/inc/SeedService.hh"

#include "CaloMC/inc/CaloWFExtractor.hh"
#include "CaloMC/inc/CaloNoiseARFitter.hh"
#include "Mu2eUtilities/inc/CaloPulseShape.hh"

#include "CLHEP/Random/RandPoissonQ.h"
#include "CLHEP/Random/RandGaussQ.h"
#include "CLHEP/Random/RandFlat.h"
#include <vector>
#include <string>


namespace mu2e {

  class CaloNoiseSimGenerator 
  {
     public:
        struct Config
        {
            using Name    = fhicl::Name;
            using Comment = fhicl::Comment;        
            fhicl::Atom<double>   elecNphotPerNs { Name("elecNphotPerNs"), Comment("Electronics noise number of PE / ns ") }; 
            fhicl::Atom<double>   rinNphotPerNs  { Name("rinNphotPerNs"),  Comment("RIN noise number of PE / ns ") }; 
            fhicl::Atom<double>   darkNphotPerNs { Name("darkNphotPerNs"), Comment("SiPM Dark noise number of PE / ns ") }; 
            fhicl::Atom<double>   digiSampling   { Name("digiSampling"),   Comment("Digitization time sampling") }; 
            fhicl::Atom<unsigned> noiseWFSize    { Name("noiseWFSize"),    Comment("Noise WF size") };
            fhicl::Atom<bool>     enableAR       { Name("enableAR"),       Comment("Enable AR noise generation ") }; 
            fhicl::Atom<double>   nparAR         { Name("nParAR"),         Comment("Number parameters for AR fit ") }; 
            fhicl::Atom<unsigned> nMaxFragment   { Name("nMaxFragment"),   Comment("maximum number of wf generated for extracting noise fragments ") }; 
            fhicl::Atom<int>      diagLevel      { Name("diagLevel"),      Comment("Diag Level"),0 };
        };

     
        CaloNoiseSimGenerator(const Config& config, CLHEP::HepRandomEngine& engine, int iRO);

        void                         initialize(const CaloWFExtractor& wfExtractor);
        void                         refresh();

        void                         addFullNoise(std::vector<double>& wfVector, bool doAR);
        void                         addSampleNoise(std::vector<double>& wfVector, unsigned istart, unsigned ilength);
        void                         addSaltAndPepper(std::vector<double>& wfVector);
        void                         plotNoise(std::string name);

        const std::vector<double>&   noise()    const {return waveform_;}
        int                          pedestal() const {return pedestal_;}


     private:
        using vvd = std::vector<std::vector<double>>;
 
        void                  generateWF(std::vector<double>& wfVector);
        void                  generateFragments(const CaloWFExtractor& wfExtractor);
        void                  initAR();

        unsigned              iRO_;
        std::vector<double>   waveform_;
        int                   pedestal_;
        vvd                   digiNoise_;
        double                digiNoiseProb_;
        double                digiSampling_;
        double                noiseRinDark_;
        double                noiseElec_;
        CLHEP::RandPoissonQ   randPoisson_;
        CLHEP::RandGaussQ     randGauss_;
        CLHEP::RandFlat       randFlat_;
        unsigned              nMaxFragment_;
        bool                  enableAR_;
        unsigned              nparFitAR_;
        CaloNoiseARFitter     ARFitter_;
        CaloPulseShape        pulseShape_;
        int                   diagLevel_;
   };

}
#endif
