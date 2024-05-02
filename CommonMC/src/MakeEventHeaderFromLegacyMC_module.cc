//
// Read legacy MC data products and create an mu2e::EventHeader data product, that is populated
// from information found in the EventWindowMarker and ProtonBunchTime data products.  This emulates
// what data will look like.
//
// The plan is to prepare for data by migrating trigger and reco algorithms to using this
// data product instead of EventWindowMarker and ProtonBunchTime.
//
// The mu2e::EventHEadevent header combines information from the heartbeat packet, defined in,
//    https://mu2e-docdb.fnal.gov/cgi-bin/sso/ShowDocument?docid=4914
// and the propossed CFO Event Window Data Record, defined on page 22 of
//   https://mu2e-docdb.fnal.gov/cgi-bin/sso/ShowDocument?docid=19095
//
//  Rob Kutschke, 2024
//

#include "Offline/DataProducts/inc/EventWindowMarker.hh"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/PhysicsParams.hh"
#include "Offline/MCDataProducts/inc/ProtonBunchTimeMC.hh"
#include "Offline/Mu2eUtilities/inc/EventHeaderFacade.hh"

#include "artdaq-core-mu2e/Data/EventHeader.hh"

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"

#include <iostream>

namespace mu2e {

  class MakeEventHeaderFromLegacyMC : public art::EDProducer {

  public:

    struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;

      fhicl::Atom<art::InputTag> protonBunchTimeMCTag{Name("protonBunchTimeMCTag"), Comment("Input tag for a ProtonBunchTimeMC data product.")};
      fhicl::Atom<art::InputTag> eventWindowMarkerTag{Name("eventWindowMarkerTag"), Comment("Input tag for an EventWindowMarker data product.")};
      fhicl::Atom<unsigned>                  maxPrint{Name("maxPrint"),             Comment("Maximum number of events to print."), 0};
    };
    typedef art::EDProducer::Table<Config> Parameters;

    explicit MakeEventHeaderFromLegacyMC(const Parameters& conf);

    void produce( art::Event& event) override;

  private:

    // A copy of the run time configuration.
    Config _conf;

    // Initialized from run-time configuration
    art::ProductToken<ProtonBunchTimeMC>   _protonBunchTimeMCToken;
    art::ProductToken<EventWindowMarker>   _eventWindowMarkerToken;
    unsigned                               _maxPrint;

    // Intialized from GlobalConstants service
    // Duration of one tick of the DAQ clock
    double _nominalDAQClockTick; // ns

    // Duration of one tick of the clock that measures RF0
    double _nominalRF0ClockTick; // ns

    // Counter to limit printout.
    unsigned _nEvents = 0;

  };

  MakeEventHeaderFromLegacyMC::MakeEventHeaderFromLegacyMC(const Parameters& conf)
    : art::EDProducer(conf),
      _conf(conf()),
      _protonBunchTimeMCToken{consumes<ProtonBunchTimeMC>(conf().protonBunchTimeMCTag())},
    _eventWindowMarkerToken{consumes<EventWindowMarker>(conf().eventWindowMarkerTag())},
    _maxPrint(conf().maxPrint()),
    _nominalDAQClockTick{GlobalConstantsHandle<PhysicsParams>()->getNominalDAQClockTick()},
    _nominalRF0ClockTick{GlobalConstantsHandle<PhysicsParams>()->getNominalRF0ClockTick()}
   {
     produces<EventHeader>();
   }

  void MakeEventHeaderFromLegacyMC::produce( art::Event& event){

    auto const& pbtmc = event.getProduct(_protonBunchTimeMCToken);
    auto const& ewm   = event.getProduct(_eventWindowMarkerToken);

    const uint16_t eventDuration{uint16_t(std::round(ewm.eventLength()/_nominalDAQClockTick))};

    // Fixme: this may need more work when new bits are assigned.
    const uint8_t flags{uint8_t(ewm.spillType())};

    // Fixme: do we want to mock this up?
    const EWT ewt{0};

    // Convert MC truth time offset to the unsigned binary format.
    // Fixme: need to understand and document the actual conventions for this parameter.
    //        currently reverse engineering what was done in MDC2020.
    const uint8_t rfmTDC_estimated = uint8_t( std::floor(-pbtmc.pbtime_ / _nominalRF0ClockTick) );

    // Fixme: Good enough for now.  Our legacy MC info does not have separate estimated and measured values.
    const uint8_t rfmTDC_measured{rfmTDC_estimated};

    // Fixme: Need to set real bits here.
    const uint32_t mode{0};

    auto header = std::make_unique<EventHeader>( ewt, mode, rfmTDC_estimated, flags, eventDuration, rfmTDC_measured);

    if ( _nEvents < _maxPrint ) {
      ++_nEvents;
      std::cout << "\nEvent:  " << event.id() << std::endl;
      std::cout << "  EWM:             Spill type: " << ewm.spillType() << "  Event length: " << ewm.eventLength() <<  std::endl;
      std::cout << "  ProtonBunchTimeMC:  pbtime_: " << pbtmc.pbtime_ << std::endl;
      std::cout << "  Header:                      " << *header << std::endl;
      EventHeaderFacade f(*header, *GlobalConstantsHandle<PhysicsParams>());
      std::cout << "  Readback: " << std::endl;
      std::cout << "     Event duration:   " << f.eventDuration() << std::endl;
      std::cout << "     pbtime measured:  " << f.rf0OffsetMeasured()
                << "     residual:         " << f.rf0OffsetMeasured()- pbtmc.pbtime_
                << std::endl;
      std::cout << "     pbtime estimated: " << f.rf0OffsetEstimated() << std::endl;
    }

    event.put(std::move(header));

  } // end analyze

} // end namespace mu2e

DEFINE_ART_MODULE(mu2e::MakeEventHeaderFromLegacyMC)
