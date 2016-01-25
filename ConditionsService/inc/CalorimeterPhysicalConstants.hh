#ifndef ConditionsService_CalorimeterPhysicalConstants_hh
#define ConditionsService_CalorimeterPhysicalConstants_hh
//
// Physical Constants for calorimeter calibrations.
//

// C++ includes.
#include <iostream>
#include <string>

// Mu2e includes.
#include "Mu2eInterfaces/inc/ConditionsEntity.hh"

namespace mu2e
{
  class SimpleConfig;


  struct CalorimeterPhysicalConstants: virtual public ConditionsEntity{


    CalorimeterPhysicalConstants ( SimpleConfig const& config );

    double radiationLength()    const{return _radiationLength;}
    double criticalEnergyPos()  const{return _criticalEnergyPos;}
    double criticalEnergyNeg()  const{return _criticalEnergyNeg;}
    double density()            const{return _density;}


  private:

    double _radiationLength;
    double _criticalEnergyPos;
    double _criticalEnergyNeg;
    double _density;


    // We want to discourage multi-phase construction.
    CalorimeterPhysicalConstants();

  };


}

#endif /* ConditionsService_CalorimeterPhysicalConstants_hh */
