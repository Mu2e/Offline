#ifndef TrackerConditions_StrawDriftMaker_hh
#define TrackerConditions_StrawDriftMaker_hh

//
// construct a StrawDrift conditions entity
// from fcl or database
//


#include <iostream>
#include <vector>
#include <string>
#include "TrackerConditions/inc/StrawDrift.hh"
#include "TrackerConfig/inc/StrawDriftConfig.hh"


namespace mu2e {
  class StrawDriftMaker {
  public:
    StrawDriftMaker(StrawDriftConfig const& config):_config(config) {}
    StrawDrift::ptr_t fromFcl();
    StrawDrift::ptr_t fromDb( /* db tables will go here*/ );
  
  private:

    // this object needs to be thread safe, 
    // _config should only be initialized once
    const StrawDriftConfig _config;

  };
}
#endif

