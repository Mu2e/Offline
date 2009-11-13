#ifndef DAQParams_H
#define DAQParams_H
//
// Parameters of the DAQ system.
//
// $Id: DAQParams.hh,v 1.1 2009/11/13 23:07:51 kutschke Exp $
// $Author: kutschke $
// $Date: 2009/11/13 23:07:51 $
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

  struct DAQParams: public ConditionsEntity{

    // Start of the daq live window, in ns, relative to the arrival time of 
    // the protons on the production target.
    double t0;

    DAQParams ( SimpleConfig const& config );

    // Accept compiler generated:
    // copy c'tor, d'tor and assignment operator.

  private:

    // We want to discourage multi-phase construction.
    DAQParams ();

  };

  // Shift left (printing) operator.
  inline std::ostream& operator<<(std::ostream& ost,
                                  const DAQParams& daqpar ){
    ost << "( "
	<< daqpar.t0 << ", "
	<< " )";

    return ost;
  }


}

#endif
