#ifndef ConditionsService_DAQParams_hh
#define ConditionsService_DAQParams_hh
//
// Parameters of the DAQ system.
//
// $Id: DAQParams.hh,v 1.6 2012/02/24 20:55:48 gandr Exp $
// $Author: gandr $
// $Date: 2012/02/24 20:55:48 $
//
// Original author Rob Kutschke
//

// C++ includes.
#include <iostream>

// Mu2e includes.
#include "Mu2eInterfaces/inc/ConditionsEntity.hh"

namespace mu2e
{
  class SimpleConfig;

  struct DAQParams: virtual public ConditionsEntity{

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

#endif /* ConditionsService_DAQParams_hh */
