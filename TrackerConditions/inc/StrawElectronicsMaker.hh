#ifndef TrackerConditions_StrawElectronicsMaker_hh
#define TrackerConditions_StrawElectronicsMaker_hh

//
// construct a StrawElectronics conditions entity
// from fcl or database
//

#include "TrackerConditions/inc/StrawElectronics.hh"
#include "TrackerConditions/inc/StrawElectronicsConfig.hh"


namespace mu2e {

  class StrawElectronicsMaker {
  public:
    StrawElectronicsMaker(StrawElectronicsConfig const& config):_config(config) {}
    StrawElectronics::ptr_t fromFcl();
    StrawElectronics::ptr_t fromDb( /* db tables will go here*/ );
  
  private:

    // this object needs to be thread safe, 
    // _config should only be initialized once
    const StrawElectronicsConfig _config;

  };
}


#endif

