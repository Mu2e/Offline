// Mu2e includes
#include "ConfigTools/inc/SimpleConfig.hh"
#include "DetectorSolenoidGeom/inc/DetectorSolenoid.hh"
#include "ExternalNeutronAbsorberGeom/inc/ExternalNeutronAbsorberMaker.hh"
#include "ExternalNeutronAbsorberGeom/inc/ExternalNeutronAbsorber.hh"

// C++ includes
#include <vector>

// CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"

// Other includes
#include "cetlib/exception.h"

namespace mu2e {

  std::unique_ptr<ExternalNeutronAbsorber> ExternalNeutronAbsorberMaker::make(const SimpleConfig& c, const DetectorSolenoid& ds ) {
    
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

    std::unique_ptr<ExternalNeutronAbsorber> ena ( new ExternalNeutronAbsorber() );
    
    ena->_halfLengthZ   = c.getDouble("extneutronabs.halfLengthZ");
    ena->_halfLengthXY  = c.getDouble("extneutronabs.halfLengthXY");
    ena->_halfThickness = c.getDouble("extneutronabs.halfThickness");
    
    ena->_position      = CLHEP::Hep3Vector( ds.position().x(),0,c.getDouble("extneutronabs.z0"));

    ena->_material      = c.getString("extneutronabs.materialName");
    
    return ena;

  } // make()

} // namespace mu2e
