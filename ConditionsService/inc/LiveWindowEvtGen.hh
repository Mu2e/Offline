#ifndef LiveWindowEvtGen_H
#define LiveWindowEvtGen_H
//
// The time window during which the event generator will
// create events.  This can be different from the live window
// over which the digitizers work.
//
// $Id: LiveWindowEvtGen.hh,v 1.1 2009/11/12 00:51:08 kutschke Exp $
// $Author: kutschke $
// $Date: 2009/11/12 00:51:08 $
//
// Original author Rob Kutschke
//

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

  private:

    // We want to discourage multi-phase construction.
    LiveWindowEvtGen ();

  };
}

#endif
