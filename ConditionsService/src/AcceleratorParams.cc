//
// Some parameters of the accelerator complex.
//
// $Id: AcceleratorParams.cc,v 1.5 2014/04/01 15:03:16 knoepfel Exp $
// $Author: knoepfel $
// $Date: 2014/04/01 15:03:16 $
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

    potPulse  = config.getString("acceleratorParams.potPulse" );  
    acDipole  = config.getString("acceleratorParams.acDipole" );
    outOfTime = config.getString("acceleratorParams.outOfTime");

  }

}
