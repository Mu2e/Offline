// Mu2e includes
#include "ConfigTools/inc/SimpleConfig.hh"
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

    // Check for deprecated variables
    std::vector<std::string> variables;
    c.getNames( variables );
    std::for_each( variables.begin(),
                   variables.end(),
                   [](std::string var) {
                     if ( var.find("neutronabsorber.") != std::string::npos) {
                       throw cet::exception("GEOM") <<
                         " Variable with name \"neutronabsorber.*\" now deprecated.\n" ;                       
                     }
                   } );
    
    std::unique_ptr<InternalNeutronAbsorber> ina ( new InternalNeutronAbsorber() );
    
    ina->_rOut        = c.getDouble("intneutronabs.rOut");
    ina->_rIn2        = c.getDouble("intneutronabs.rIn2");
    ina->_rIn1        = c.getDouble("intneutronabs.rIn1");
    
    ina->_halfLength1 = c.getDouble("intneutronabs.halfLength1");
    ina->_halfLength2 = c.getDouble("intneutronabs.halfLength2");

    ina->_position    = CLHEP::Hep3Vector( ds.position().x(),0,c.getDouble("intneutronabs.z0"));

    ina->_mat1        = c.getString("intneutronabs.material1Name");
    ina->_mat2        = c.getString("intneutronabs.material2Name");
    
    return ina;

  } // make()

} // namespace mu2e
