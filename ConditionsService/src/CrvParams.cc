//
// Some parameters for the CRV.
//
// $Id: AcceleratorParams.cc,v 1.6 2014/04/14 18:12:55 knoepfel Exp $
// $Author: knoepfel $
// $Date: 2014/04/14 18:12:55 $
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
