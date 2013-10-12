// Mu2e includes
#include "ConfigTools/inc/SimpleConfig.hh"
#include "ConfigTools/inc/checkForStale.hh"
#include "DetectorSolenoidGeom/inc/DetectorSolenoid.hh"
#include "InternalNeutronAbsorberGeom/inc/InternalNeutronAbsorberMaker.hh"
#include "InternalNeutronAbsorberGeom/inc/InternalNeutronAbsorber.hh"

// C++ includes
#include <algorithm>
#include <iostream>
#include <vector>

// CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"

// Other includes
#include "cetlib/exception.h"

namespace mu2e {

  std::unique_ptr<InternalNeutronAbsorber> InternalNeutronAbsorberMaker::make(const SimpleConfig& c, const DetectorSolenoid& ds ) {

    checkForStale( "neutronabsorber", c );

    std::unique_ptr<InternalNeutronAbsorber> ina ( new InternalNeutronAbsorber() );
    
    ina->_rOut        = c.getDouble("intneutronabs.rOut");
    ina->_r4          = c.getDouble("intneutronabs.r4");
    ina->_rIn3        = c.getDouble("intneutronabs.rIn3");
    ina->_rIn2        = c.getDouble("intneutronabs.rIn2");
    ina->_rIn1        = c.getDouble("intneutronabs.rIn1");
    
    ina->_halfLength1 = c.getDouble("intneutronabs.halfLength1");
    ina->_halfLength2 = c.getDouble("intneutronabs.halfLength2");
    ina->_halfLength3 = c.getDouble("intneutronabs.halfLength3");
    ina->_halfLength4 = c.getDouble("intneutronabs.halfLength4");

    ina->_position    = CLHEP::Hep3Vector( ds.position().x(),0,c.getDouble("intneutronabs.z0"));

    ina->_mat1        = c.getString("intneutronabs.material1Name");
    ina->_mat2        = c.getString("intneutronabs.material2Name");
    ina->_mat3        = c.getString("intneutronabs.material3Name");
    ina->_mat4        = c.getString("intneutronabs.material4Name");
    
    return ina;

  } // make()

} // namespace mu2e
