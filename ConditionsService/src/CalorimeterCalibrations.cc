//
// Parameters for calorimeter calibrations.
//
// $Id: CalorimeterCalibrations.cc,v 1.1 2013/03/05 20:33:25 aluca Exp $
// $Author: aluca $
// $Date: 2013/03/05 20:33:25 $
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
      
      }
}

