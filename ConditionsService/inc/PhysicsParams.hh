#ifndef ConditionsService_PhysicsParams_hh
#define ConditionsService_PhysicsParams_hh
//
// Some physical parameters.
//
// Original author Gianni Onorato
//

// C++ includes.
#include <iostream>

// Mu2e includes.
#include "ConditionsService/inc/ConditionsEntity.hh"

namespace mu2e
{
  class SimpleConfig;

  struct PhysicsParams: public ConditionsEntity{

    // Nominal decay time for bound state in alluminum nucleus
    double decayTime;

    PhysicsParams( SimpleConfig const& config );

    // Accept compiler generated:
    // copy c'tor, d'tor and assignment operator.

  private:

    // We want to discourage multi-phase construction.
    PhysicsParams ();

  };

  // Shift left (printing) operator.
  inline std::ostream& operator<<(std::ostream& ost,
                                  const PhysicsParams& lw ){
    ost << "( "
        << lw.decayTime << ", "
        << " )";

    return ost;
  }

}

#endif /* ConditionsService_PhysicsParams_hh */
