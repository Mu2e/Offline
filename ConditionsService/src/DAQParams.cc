//
// Parameters of the DAQ system.
//
// $Id: DAQParams.cc,v 1.1 2009/11/13 23:07:51 kutschke Exp $
// $Author: kutschke $
// $Date: 2009/11/13 23:07:51 $
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
