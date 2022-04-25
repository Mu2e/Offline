//
// Generate a single time sample randomly distributed
// between two times for cosmic reasampling
//
// Original author: Stefano Roberto Soleti, roberto@lbl.gov, 2021

#include <string>
#include <memory>

#include "fhiclcpp/types/OptionalTable.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "canvas/Persistency/Common/Ptr.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"

#include "Offline/MCDataProducts/inc/SimTimeOffset.hh"
#include "Offline/MCDataProducts/inc/GenParticle.hh"
#include "Offline/SeedService/inc/SeedService.hh"

#include "CLHEP/Random/RandFlat.h"

namespace mu2e {

  class CosmicTimeOffset : public art::EDProducer {

    public:

      struct Config {
        using Name=fhicl::Name;
        using Comment=fhicl::Comment;

        fhicl::OptionalAtom<std::string> cosmicModuleLabel{Name("cosmicModuleLabel"), Comment("Name of cosmic module label")};
        fhicl::Atom<int> verbosityLevel{ Name("verbosityLevel"), Comment("Levels 0, 1, >1"), 0 };
        fhicl::Atom<float> intervalStart{ Name("intervalStart"), Comment("Start time of the time offset window"), 250 };
        fhicl::Atom<float> intervalEnd{ Name("intervalEnd"), Comment("end time of the time offset window"), 1700 };
      };

      using Parameters = art::EDProducer::Table<Config>;
      explicit CosmicTimeOffset(const Parameters& conf);

      virtual void beginRun(art::Run&   r) override;
      virtual void produce (art::Event& e) override;

    private:
      art::RandomNumberGenerator::base_engine_t& engine_;

      std::string cosmicModuleLabel_;
      bool addTimeOffset_;
      int  verbosityLevel_;
      float intervalStart_;
      float intervalEnd_;

      CLHEP::RandFlat flatTime_;
  };

  //================================================================
  CosmicTimeOffset::CosmicTimeOffset(const Parameters& conf)
    : EDProducer{conf}
  , engine_(createEngine(art::ServiceHandle<SeedService>()->getSeed()))
    , addTimeOffset_(false)
    , verbosityLevel_(conf().verbosityLevel())
    , intervalStart_(conf().intervalStart())
    , intervalEnd_(conf().intervalEnd())
    , flatTime_(engine_, intervalStart_, intervalEnd_)
    {
      produces<SimTimeOffset>();

      addTimeOffset_ = conf().cosmicModuleLabel(cosmicModuleLabel_);

      if (addTimeOffset_) {
        consumes<GenParticleCollection>(cosmicModuleLabel_);
        produces<GenParticleCollection>();
      }
    }

  //================================================================
  void CosmicTimeOffset::beginRun(art::Run& run) {
    if ( verbosityLevel_ > 0 ) {
      std::ostringstream timeSpectrum;
      timeSpectrum << "Genarating random time offset between "
        << intervalStart_
        << " ns and "
        << intervalEnd_ << " ns\n";
      mf::LogInfo("Info") << "Cosmic time distribution\n" << timeSpectrum.str();
    }
  }

  //================================================================
  void CosmicTimeOffset::produce(art::Event& event) {
    // Generate and record offset
    const float timeOffset = flatTime_.fire();
    std::unique_ptr<SimTimeOffset> toff(new SimTimeOffset(timeOffset));
    if( verbosityLevel_ > 1) std::cout << "CosmicTimeOffset " << toff->timeOffset_ << std::endl;
    event.put(std::move(toff));

    if (addTimeOffset_ ) {
      // Add offset to all particles
      art::Handle<mu2e::GenParticleCollection> cosmicParticles;
      event.getByLabel(cosmicModuleLabel_, cosmicParticles);
      const GenParticleCollection &particles(*cosmicParticles);
      std::unique_ptr<GenParticleCollection> offsetParticles(new GenParticleCollection);

      for (GenParticleCollection::const_iterator i = particles.begin(); i != particles.end(); ++i) {
        GenParticle const &particle = *i;
        offsetParticles->push_back(GenParticle(particle.pdgId(), particle.generatorId(),
              particle.position(), particle.momentum(),
              particle.time() + timeOffset));
      }

      event.put(std::move(offsetParticles));
    }
  }

} // end namespace mu2e

DEFINE_ART_MODULE(mu2e::CosmicTimeOffset)
