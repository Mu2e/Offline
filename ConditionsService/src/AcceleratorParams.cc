//
// Some parameters of the accelerator complex.
//
// $Id: AcceleratorParams.cc,v 1.1 2009/11/13 23:07:51 kutschke Exp $
// $Author: kutschke $
// $Date: 2009/11/13 23:07:51 $
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
