#ifndef AcceleratorParams_H
#define AcceleratorParams_H
//
// Some parameters of the accelerator complex.
//
// $Id: AcceleratorParams.hh,v 1.1 2009/11/13 23:07:51 kutschke Exp $
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

  struct AcceleratorParams: public ConditionsEntity{

    // The nominal debuncher orbital period.
    double deBuncherPeriod;

    AcceleratorParams ( SimpleConfig const& config );

    // Accept compiler generated:
    // copy c'tor, d'tor and assignment operator.

  private:

    // We want to discourage multi-phase construction.
    AcceleratorParams ();

  };

  // Shift left (printing) operator.
  inline std::ostream& operator<<(std::ostream& ost,
                                  const AcceleratorParams& lw ){
    ost << "( "
	<< lw.deBuncherPeriod << ", "
	<< " )";

    return ost;
  }


}

#endif
