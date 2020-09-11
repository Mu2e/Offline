//
// Parameters for calorimeter calibrations.
//
//

// Mu2e include files
#include "ConditionsService/inc/CalorimeterCalibrations.hh"
#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "ConfigTools/inc/SimpleConfig.hh"

#include <cmath>



namespace mu2e {

     CalorimeterCalibrations::CalorimeterCalibrations( SimpleConfig const& config )
     {
        ConfigFileLookupPolicy configFile;
	
	_pulseFileName = configFile(config.getString("calorimeter.pulseFileName"));
        _pulseHistName = config.getString("calorimeter.pulseHistName");

        // Here we should eventually interface to some database
        _LRUpar0 = config.getDouble("CrystalNonUniformity_0");
 
        // Here we should eventually interface to some database
        _BirkCorrHadron = config.getDouble("BirkCorrHadron");

        //RO photo-statistic number
        _peMeV = config.getDouble("ROphotostatistic");//p.e. / MeV

        //value of the sigma used to do the Gaussian smearing due to the electronic noise
        _ROnoise = config.getDouble("ReadOutElectronicNoise");//MeV

	//conversion factor between ADC counts and MeV for a specific RO
	_ADC2MeV  = config.getDouble("ADC2MeVConversionFactor");

	//conversion factor between ADC counts and MeV for Fast clustering
	_Peak2MeV  = config.getDouble("Peak2MeVConversionFactor");
      
      }
}

