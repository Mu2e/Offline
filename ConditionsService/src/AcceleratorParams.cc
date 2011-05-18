//
// Some parameters of the accelerator complex.
//
// $Id: AcceleratorParams.cc,v 1.2 2011/05/18 02:27:15 wb Exp $
// $Author: wb $
// $Date: 2011/05/18 02:27:15 $
//

// Mu2e include files
#include "ConditionsService/inc/AcceleratorParams.hh"
#include "Mu2eUtilities/inc/SimpleConfig.hh"

namespace mu2e {

  AcceleratorParams::AcceleratorParams( SimpleConfig const& config ){

    // Throws if the entity is not given in the config file.
    deBuncherPeriod  = config.getDouble("acceleratorParams.deBuncherPeriod");

  }

}
