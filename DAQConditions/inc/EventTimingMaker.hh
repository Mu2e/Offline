#ifndef DAQConditions_EventTimingMaker_hh
#define DAQConditions_EventTimingMaker_hh

//
// construct a EventTiming conditions entity
// from fcl or database
//

#include "Offline/DAQConditions/inc/EventTiming.hh"
#include "Offline/DAQConfig/inc/EventTimingConfig.hh"
#include "Offline/DbTables/inc/DAQTiming.hh"


namespace mu2e {

  class EventTimingMaker {
  public:
    EventTimingMaker(EventTimingConfig const& config):_config(config) {}
    EventTiming::ptr_t fromFcl();
    EventTiming::ptr_t fromDb(DAQTiming::cptr_t dt_p);

  private:

    // this object needs to be thread safe,
    // _config should only be initialized once
    const EventTimingConfig _config;

  };
}


#endif

