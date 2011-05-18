//
// Parameters of the DAQ system.
//
// $Id: DAQParams.cc,v 1.2 2011/05/18 02:27:15 wb Exp $
// $Author: wb $
// $Date: 2011/05/18 02:27:15 $
//

// Mu2e include files
#include "ConditionsService/inc/DAQParams.hh"
#include "Mu2eUtilities/inc/SimpleConfig.hh"

namespace mu2e {

  DAQParams::DAQParams( SimpleConfig const& config ){

    // These throw if the entity is not given in the config file.
    t0 = config.getDouble("DAQParams.t0");

  }

}
