#ifndef ConditionsService_CrvParams_hh
#define ConditionsService_CrvParams_hh
//
// Some parameters of the CRV.
//
// $Id: AcceleratorParams.hh,v 1.8 2014/04/14 18:12:55 knoepfel Exp $
// $Author: knoepfel $
// $Date: 2014/04/14 18:12:55 $
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

  struct CrvParams: virtual public ConditionsEntity{

    // CRV digi parameters
    double digitizationPeriod;

    CrvParams ( SimpleConfig const& config );

    // Accept compiler generated:
    // copy c'tor, d'tor and assignment operator.

  private:

    // We want to discourage multi-phase construction.
    CrvParams ();

  };

  // Shift left (printing) operator.
  inline std::ostream& operator<<(std::ostream& ost,
                                  const CrvParams& lw ){
    ost << "( "
        << lw.digitizationPeriod 
        << " )";

    return ost;
  }

}

#endif /* ConditionsService_CrvParams_hh */
