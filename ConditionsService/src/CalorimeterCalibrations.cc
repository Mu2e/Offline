//
// Parameters for calorimeter calibrations.
//
// $Id: CalorimeterCalibrations.cc,v 1.1 2013/03/05 20:33:25 aluca Exp $
// $Author: aluca $
// $Date: 2013/03/05 20:33:25 $
//

// Mu2e include files
#include "ConditionsService/inc/CalorimeterCalibrations.hh"
#include "ConfigTools/inc/SimpleConfig.hh"
//#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
//#include "SeedService/inc/SeedService.hh"
//#include "art/Framework/Services/Registry/ServiceHandle.h"
//#include "art/Framework/Core/EngineCreator.h"
//#include "EventGenerator/inc/GeneratorBase.hh"

#include "CLHEP/Units/SystemOfUnits.h"
//#include "CLHEP/Random/RandomEngine.h"

#include <cmath>

// Mu2e includes.
//#include "GeometryService/inc/GeometryService.hh"
//#include "GeometryService/inc/GeomHandle.hh"


namespace mu2e {

CalorimeterCalibrations::CalorimeterCalibrations( SimpleConfig const& config )
{

        // Here we should eventually interface to some database
        _LRUpar0 = config.getDouble("CrystalNonUniformity_0", 0.45) / CLHEP::cm * 0.01;
        _LRUpar0Err = config.getDouble("CrystalNonUniformityErr_0", 0.0) / CLHEP::cm * 0.01;

        //non-linearity parameters
        _linpar0 = config.getDouble("CrystalNonLinearity_0", 10.35);
        _linpar0Err = config.getDouble("CrystalNonLinearity_0Err", 0.0);

        _linpar1 = config.getDouble("CrystalNonLinearity_1", 0.9919);
        _linpar1Err = config.getDouble("CrystalNonLinearity_1Err", 0.0);

        _linpar2 = config.getDouble("CrystalNonLinearity_2", 0.01812);
        _linpar2Err = config.getDouble("CrystalNonLinearity_2Err", 0.0);

        _linpar3 = config.getDouble("CrystalNonLinearity_3",  0.06654);
        _linpar3Err = config.getDouble("CrystalNonLinearity_3Err",  0.0);

        //RO photo-statistic number
        _ROpe = config.getDouble("ROphotostatistic", 1660.0);//p.e. / MeV
	_ROpeErr = config.getDouble("ROphotostatisticErr", 224.0);//p.e. / MeV
        _ROfano = config.getDouble("ROfanoFactor", 1.3);

        //value of the sigma used to do the Gaussian smearing due to the electronic noise
        _ROnoise = config.getDouble("ReadOutElectronicNoise", 0.030);//MeV


	//conversion factor between ADC counts and MeV for a specific RO
	_ADC2MeV  = config.getDouble("ADC2MeVConversionFactor", 1.);
      }

}
