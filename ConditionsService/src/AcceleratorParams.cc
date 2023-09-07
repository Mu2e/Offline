//
// Some parameters of the accelerator complex.
//
//

// Mu2e include files
#include "Offline/ConditionsService/inc/AcceleratorParams.hh"
#include "Offline/ConfigTools/inc/SimpleConfig.hh"

#include <iostream>

namespace mu2e {

  AcceleratorParams::AcceleratorParams( SimpleConfig const& config ){

    // Throws if the entity is not given in the config file.
    deBuncherPeriod     = config.getDouble("acceleratorParams.deBuncherPeriod"   );

  }

}
