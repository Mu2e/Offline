//
// Some parameters for the CRV.
//
//

#include "ConditionsService/inc/CrvParams.hh"
#include "ConfigTools/inc/SimpleConfig.hh"

#include <iostream>

namespace mu2e 
{
  CrvParams::CrvParams( SimpleConfig const& config )
  {
    digitizationPeriod   = config.getDouble("crv.digitizationPeriod");
    pedestal             = config.getDouble("crv.pedestal");
    calibrationFactor    = config.getDouble("crv.calibrationFactor");
    calibrationFactorPulseHeight = config.getDouble("crv.calibrationFactorPulseHeight");
  }
}
