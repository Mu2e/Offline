#ifndef LiveWindowEvtGen_H
#define LiveWindowEvtGen_H
//
// The time window during which event generators will create events.  
// This is distinct from live window over which the digitizer simulation
// is asked to work and both are distinct from the settings of the
// real experiment for any particular run.
//
// $Id: LiveWindowEvtGen.hh,v 1.2 2009/11/12 01:35:23 kutschke Exp $
// $Author: kutschke $
// $Date: 2009/11/12 01:35:23 $
//
// Original author Rob Kutschke
//

// C++ includes.
#include <iostream>

// Mu2e includes.
#include "ConditionsService/inc/ConditionsEntity.hh"

namespace mu2e
{
  class SimpleConfig;

  struct LiveWindowEvtGen: public ConditionsEntity{

    // Start and end of the counting window, in ns, relative 
    // to the arrival time of the protons on the production target.
    double t0;     
    double tend;

    LiveWindowEvtGen ( SimpleConfig const& config );

    // Accept compiler generated:
    // copy c'tor, d'tor and assignment operator.

  private:

    // We want to discourage multi-phase construction.
    LiveWindowEvtGen ();

  };

  // Shift left (printing) operator.
  inline std::ostream& operator<<(std::ostream& ost,
                                  const LiveWindowEvtGen& lw ){
    ost << "( "
	<< lw.t0 << ", "
	<< lw.tend 
	<< " )";

    return ost;
  }


}

#endif
