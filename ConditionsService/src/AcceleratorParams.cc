//
// Some parameters of the accelerator complex.
//
// $Id: AcceleratorParams.cc,v 1.4 2014/02/28 21:11:19 knoepfel Exp $
// $Author: knoepfel $
// $Date: 2014/02/28 21:11:19 $
//

// Mu2e include files
#include "ConditionsService/inc/AcceleratorParams.hh"
#include "ConfigTools/inc/SimpleConfig.hh"

#include <iostream>

namespace mu2e {

  AcceleratorParams::AcceleratorParams( SimpleConfig const& config ){

    // Throws if the entity is not given in the config file.
    deBuncherPeriod     = config.getDouble("acceleratorParams.deBuncherPeriod"   );
    intrinsicExtinction = config.getDouble("acceleratorParams.intrinsicExt"      );
    limitingHalfWidth   = config.getDouble("acceleratorParams.limitingHalfWidth" );

  }

}
