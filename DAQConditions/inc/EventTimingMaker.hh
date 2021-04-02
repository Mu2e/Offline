#ifndef DAQConditions_EventTimingMaker_hh
#define DAQConditions_EventTimingMaker_hh

//
// construct a EventTiming conditions entity
// from fcl or database
//

#include "DAQConditions/inc/EventTiming.hh"
#include "DAQConfig/inc/EventTimingConfig.hh"


namespace mu2e {

  class EventTimingMaker {
  public:
    EventTimingMaker(EventTimingConfig const& config):_config(config) {}
    EventTiming::ptr_t fromFcl();
  
  private:

    // this object needs to be thread safe, 
    // _config should only be initialized once
    const EventTimingConfig _config;

  };
}


#endif

