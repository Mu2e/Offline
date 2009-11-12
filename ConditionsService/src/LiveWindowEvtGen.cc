//
// The live time window for Event Generators.
//
// $Id: LiveWindowEvtGen.cc,v 1.2 2009/11/12 01:35:23 kutschke Exp $
// $Author: kutschke $
// $Date: 2009/11/12 01:35:23 $
// 

// Mu2e include files
#include "ConditionsService/inc/LiveWindowEvtGen.hh"
#include "Mu2eUtilities/inc/SimpleConfig.hh"

namespace mu2e {

  LiveWindowEvtGen::LiveWindowEvtGen( SimpleConfig const& config ){
    
    // These throw if the entity is not given in the config file.
    t0   = config.getDouble("liveWindowEvtGen.t0");
    tend = config.getDouble("liveWindowEvtGen.tend");
    
  }
  
}
