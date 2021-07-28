//
// Parameters of the DAQ system.
//
//

// Mu2e include files
#include "Offline/ConditionsService/inc/DAQParams.hh"
#include "Offline/ConfigTools/inc/SimpleConfig.hh"

namespace mu2e {

  DAQParams::DAQParams( SimpleConfig const& config ){

    // These throw if the entity is not given in the config file.
    t0 = config.getDouble("DAQParams.t0");

  }

}
