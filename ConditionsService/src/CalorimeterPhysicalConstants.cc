//
// Parameters for calorimeter crystal internal properties.
//
//

// Mu2e include files
#include "ConditionsService/inc/CalorimeterPhysicalConstants.hh"
#include "ConfigTools/inc/SimpleConfig.hh"

#include "CLHEP/Units/SystemOfUnits.h"

#include <cmath>



namespace mu2e {

  CalorimeterPhysicalConstants::CalorimeterPhysicalConstants( SimpleConfig const& config )
  {
    //
    // initialize with absurd values
	_radiationLength = 0 ;	
	_criticalEnergyPos = 0;
	_criticalEnergyNeg = 0;
        _density = 0;

	_radiationLength   = config.getDouble("calorimeter.radiationLengthCsI")*CLHEP::cm ; //cm
	_criticalEnergyPos = config.getDouble("calorimeter.criticalEnergyPosCsI")*CLHEP::MeV; // e^-
	_criticalEnergyNeg = config.getDouble("calorimeter.criticalEnergyNegCsI")*CLHEP::MeV; // e^-
	_density           = config.getDouble("calorimeter.densityCsI")*CLHEP::g/CLHEP::cm3 ; // g/cc 

  }

}
