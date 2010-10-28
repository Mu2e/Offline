//
// Some physical parameters
//
// 

// Mu2e include files
#include "ConditionsService/inc/PhysicsParams.hh"
#include "Mu2eUtilities/inc/SimpleConfig.hh"

namespace mu2e {

  PhysicsParams::PhysicsParams( SimpleConfig const& config ){
    
    // Throws if the entity is not given in the config file.
    decayTime  = config.getDouble("physicsParams.decayTime");
    
  }
  
}
