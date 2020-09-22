#ifndef ConditionsService_CrvParams_hh
#define ConditionsService_CrvParams_hh
//
// Some parameters of the CRV.
//
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
    double digitizationPeriod;     //ns

    // CRV calibration parameters 
    // This is only a placeholder. 
    // In the future, these values need to be a function of time and SiPM number.
    double pedestal;               //ADC
    double calibrationFactor;      //ADC*ns/PE
    double calibrationFactorPulseHeight; //ADC/PE  (for calibrations using the pulse height)

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
        << lw.digitizationPeriod <<"ns, " 
        << lw.pedestal <<"ADC, "
        << lw.calibrationFactor <<"ADC*ns/PE, "
        << lw.calibrationFactorPulseHeight <<"ADC/PE"
        << " )";

    return ost;
  }

}

#endif /* ConditionsService_CrvParams_hh */
