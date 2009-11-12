//
// The live time window of the experiment.
//
// $Id: LiveWindowEvtGen.cc,v 1.1 2009/11/12 00:51:08 kutschke Exp $
// $Author: kutschke $
// $Date: 2009/11/12 00:51:08 $
// 

// Mu2e include files
#include "ConditionsService/inc/LiveWindowEvtGen.hh"
#include "Mu2eUtilities/inc/SimpleConfig.hh"

namespace mu2e {

  LiveWindowEvtGen::LiveWindowEvtGen( SimpleConfig const& config ){
    
    t0   = config.getDouble("liveWindowEvtGen.t0");
    tend = config.getDouble("liveWindowEvtGen.tend");
    
  }
  
}
