#ifndef TrackerConditions_StrawPhysicsMaker_hh
#define TrackerConditions_StrawPhysicsMaker_hh

//
// construct a StrawPhysics conditions entity
// from fcl or database
//


#include <iostream>
#include <vector>
#include <string>
#include "TrackerConditions/inc/StrawPhysics.hh"
#include "TrackerConfig/inc/StrawPhysicsConfig.hh"


namespace mu2e {
  class StrawPhysicsMaker {
  public:
    StrawPhysicsMaker(StrawPhysicsConfig const& config):_config(config) {}
    StrawPhysics::ptr_t fromFcl(StrawDrift::cptr_t strawDrift);
    StrawPhysics::ptr_t fromDb( /* db tables will go here*/ );
  
  private:

    // this object needs to be thread safe, 
    // _config should only be initialized once
    const StrawPhysicsConfig _config;

  };
}
#endif

