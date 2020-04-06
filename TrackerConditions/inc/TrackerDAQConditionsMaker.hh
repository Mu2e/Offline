#ifndef TrackerConditions_TrackerDAQConditionsMaker_hh
#define TrackerConditions_TrackerDAQConditionsMaker_hh

//
// construct a TrackerDAQConditions conditions entity
// from fcl or database
//

#include "TrackerConditions/inc/TrackerDAQConditions.hh"
#include "TrackerConfig/inc/TrackerDAQConditionsConfig.hh"
#include "DbTables/inc/TrkDRACtoStraw.hh"
#include "DbTables/inc/TrkROCtoPanel.hh"


namespace mu2e {

  class TrackerDAQConditionsMaker {
  public:
    TrackerDAQConditionsMaker(TrackerDAQConditionsConfig const& config):_config(config) {}
    TrackerDAQConditions::ptr_t fromFcl();
    TrackerDAQConditions::ptr_t fromDb(TrkDRACtoStraw::cptr_t tdts,
				   TrkROCtoPanel::cptr_t trtp );
  
  private:

    // this object needs to be thread safe, 
    // _config should only be initialized once
    const TrackerDAQConditionsConfig _config;

  };
}


#endif

