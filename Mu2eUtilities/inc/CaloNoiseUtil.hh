#ifndef Mu2eUtilities_CaloNoiseUtil_hh
#define Mu2eUtilities_CaloNoiseUtil_hh
//
// Cache and provide noise waveforms for readouts
//
// Call prepare cache before using it - peToADC needed if data are regenerated
// When regenerating, throw if histoID is already used for a different peToADC - each peToADC requires a __different__ noise histoID
// Return a view of the vector to add noise for efficiency reasons - invalid if cache is modified in the meantime
// TODO: fix tolerance on PeToADC ratio
//
#include "fhiclcpp/types/Atom.h"
#include "Offline/Mu2eUtilities/inc/CaloPulseUtil.hh"
#include "CLHEP/Random/RandPoissonQ.h"
#include "CLHEP/Random/RandGaussQ.h"
#include "CLHEP/Random/RandFlat.h"

#include <map>
#include <vector>
#include <span>


namespace mu2e {

  class CaloNoiseUtil
  {
     public:
        struct Config
        {
            using Name    = fhicl::Name;
            using Comment = fhicl::Comment;
            using CPG     = CaloPulseUtil::Config;
            fhicl::Table<CPG>        pulseCache     { Name("pulseCache"),     Comment("Pulse cache maker config") };
            fhicl::Atom<bool>        generate       { Name("generate"),       Comment("Regenerate waveform (true) or use histogram (false)") };
            fhicl::Atom<bool>        dumpGenerated  { Name("dumpGenerated"),  Comment("Dump generated waveform") };
            fhicl::Atom<std::string> histoFileName  { Name("histoFileName"),  Comment("Calo noise histo file name") };
            fhicl::Atom<std::string> histoPrefix    { Name("histoPrefix"),    Comment("Noise histogram prefix") };
            fhicl::Atom<double>      elecNphotPerNs { Name("elecNphotPerNs"), Comment("Electronics noise number of PE / ns ") };
            fhicl::Atom<double>      rinNphotPerNs  { Name("rinNphotPerNs"),  Comment("RIN noise number of PE / ns ") };
            fhicl::Atom<double>      darkNphotPerNs { Name("darkNphotPerNs"), Comment("SiPM Dark noise number of PE / ns ") };
            fhicl::Atom<double>      digiSampling   { Name("digiSampling"),   Comment("Digitization time sampling") };
        };


        CaloNoiseUtil(const Config& config, CLHEP::HepRandomEngine& engine);

        void              prepare(int histoID, double peToADC);
        std::span<double> noiseSegment(int histoID, size_t istart, size_t ilength);
        double            pedestal(int histoID) const;
        void              printCache() const;
        void              dumpNoise(const std::string& name, const std::vector<double>& wave);


     private:
        void fillCache(int histoBaseID);
        void generateCache(int histoID, double peToADC);

        bool                  generate_;
        std::string           fileName_;
        std::string           prefix_;
        double                digiSampling_;
        double                noiseRinDark_;
        double                noiseElec_;
        CLHEP::RandPoissonQ   randPoisson_;
        CLHEP::RandGaussQ     randGauss_;
        CLHEP::RandFlat       randFlat_;
        CaloPulseUtil         pulseCache_;
        bool                  dumpGenerated_;
        std::map<int,double>  pedestal_;
        std::map<int,double>  peToADC_;
        std::map<int,std::vector<double>> noiseMap_;

        static constexpr int base = 10000;
   };

}
#endif
