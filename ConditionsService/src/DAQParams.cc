//
// Parameters of the DAQ system.
//
//

// Mu2e include files
#include "ConditionsService/inc/DAQParams.hh"
#include "ConfigTools/inc/SimpleConfig.hh"

namespace mu2e {

  DAQParams::DAQParams( SimpleConfig const& config ){

    // These throw if the entity is not given in the config file.
    t0 = config.getDouble("DAQParams.t0");

  }

}
