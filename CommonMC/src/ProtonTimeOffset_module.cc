//  
// generate a single time sample from the proton bunch distribution
// for proton beam reasampling
//

#include <string>
#include <memory>

#include "fhiclcpp/types/OptionalTable.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "CLHEP/Random/RandGaussQ.h"

#include "canvas/Persistency/Common/Ptr.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"

#include "Offline/MCDataProducts/inc/SimTimeOffset.hh"
#include "Offline/SeedService/inc/SeedService.hh"
#include "Offline/Mu2eUtilities/inc/ProtonPulseRandPDF.hh"

namespace mu2e {

  class ProtonTimeOffset : public art::EDProducer {

  public:

    struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;

      fhicl::OptionalTable<ProtonPulseRandPDF::Config> randPDFparameters { Name("randPDFparameters") };
      fhicl::Atom<int> verbosityLevel{ Name("verbosityLevel"), Comment("Levels 0, 1, >1"), 0 };
    };

    using Parameters = art::EDProducer::Table<Config>;
    explicit ProtonTimeOffset(const Parameters& conf);

    virtual void beginRun(art::Run&   r) override;
    virtual void produce (art::Event& e) override;

  private:
    art::RandomNumberGenerator::base_engine_t& engine_;
    ProtonPulseRandPDF::Config protonPulseConf_;
    int  verbosityLevel_;
    std::unique_ptr<ProtonPulseRandPDF>  protonPulse_;
  };

  //================================================================
  ProtonTimeOffset::ProtonTimeOffset(const Parameters& conf)
    : EDProducer{conf}
    , engine_(createEngine(art::ServiceHandle<SeedService>()->getSeed()) )
    , verbosityLevel_(conf().verbosityLevel())
  {
    produces<SimTimeOffset>();
    conf().randPDFparameters(protonPulseConf_);
  }

  //================================================================
  void ProtonTimeOffset::beginRun(art::Run& run) {
    protonPulse_.reset( new ProtonPulseRandPDF( engine_, protonPulseConf_ ) );
    if ( verbosityLevel_ > 0 ) {
      std::ostringstream timeSpectrum;
      std::cout << " Size of proton pulse: " << protonPulse_->getTimes().size() << std::endl;
      for ( std::size_t i(0) ; i < protonPulse_->getTimes().size(); i++ ) {
        timeSpectrum << "   POT time: "
                     << protonPulse_->getTimes().at(i)
                     << "     "
                     << protonPulse_->getSpectrum().at(i) << "\n";
      }
      mf::LogInfo("Info") << "Longitudinal POT time distribution\n" << timeSpectrum.str();
    }
  }

  //================================================================
  void ProtonTimeOffset::produce(art::Event& event) {
// Generate and record offset
    std::unique_ptr<SimTimeOffset> toff(new SimTimeOffset(protonPulse_->fire()));
    if( verbosityLevel_ > 1) std::cout<<"ProtonTimeOffset "<< toff->timeOffset_ << std::endl;
    event.put(std::move(toff));
  } 

} // end namespace mu2e

DEFINE_ART_MODULE(mu2e::ProtonTimeOffset)
