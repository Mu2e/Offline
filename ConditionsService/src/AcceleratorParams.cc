//
// Some parameters of the accelerator complex.
//
// $Id: AcceleratorParams.cc,v 1.3 2012/07/15 22:06:16 kutschke Exp $
// $Author: kutschke $
// $Date: 2012/07/15 22:06:16 $
//

// Mu2e include files
#include "ConditionsService/inc/AcceleratorParams.hh"
#include "ConfigTools/inc/SimpleConfig.hh"

namespace mu2e {

  AcceleratorParams::AcceleratorParams( SimpleConfig const& config ){

    // Throws if the entity is not given in the config file.
    deBuncherPeriod  = config.getDouble("acceleratorParams.deBuncherPeriod");

  }

}
