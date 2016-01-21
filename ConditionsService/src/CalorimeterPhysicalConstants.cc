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

    _materialType = config.getString("calorimeter.crystalMaterial","G4_CESIUM_IODIDE");
    if (_materialType == "G4_CESIUM_IODIDE")
      {
	_radiationLength = 1.860*CLHEP::cm ; //cm
	_criticalEnergyPos = 11.17*CLHEP::MeV; // e^-
	_criticalEnergyNeg = 10.80*CLHEP::MeV; // e^+
	_density = 4.51*CLHEP::g/CLHEP::cm3;// gm/cc
      }

  }

}
