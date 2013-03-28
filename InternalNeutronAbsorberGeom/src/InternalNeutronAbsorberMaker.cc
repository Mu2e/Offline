// Mu2e includes
#include "ConfigTools/inc/SimpleConfig.hh"
#include "DetectorSolenoidGeom/inc/DetectorSolenoid.hh"
#include "InternalNeutronAbsorberGeom/inc/InternalNeutronAbsorberMaker.hh"
#include "InternalNeutronAbsorberGeom/inc/InternalNeutronAbsorber.hh"

// C++ includes
#include <vector>

// CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"

// Other includes
#include "cetlib/exception.h"

namespace mu2e {

  std::unique_ptr<InternalNeutronAbsorber> InternalNeutronAbsorberMaker::make(const SimpleConfig& c, const DetectorSolenoid& ds ) {
    
    std::unique_ptr<InternalNeutronAbsorber> ina ( new InternalNeutronAbsorber() );
    
    ina->_rOut                      = c.getDouble("intneutronabs.rOut");
    
    // Set abs1 parameters
    c.getVectorDouble("intneutronabs.abs1rIn", ina->_rInAbs1Vec );
    ina->_halfLengthAbs1            = c.getDouble("intneutronabs.abs1HalfLength");
    ina->_positionAbs1              = CLHEP::Hep3Vector( ds.position().x(),0,c.getDouble("intneutronabs.abs1Z0"));
    c.getVectorString("intneutronabs.Abs1materialName", ina->_matAbs1Vec);

    // Set abs2 parameters
    ina->_rInAbs2                   = c.getDouble("intneutronabs.abs2rIn");
    ina->_matAbs2                   = c.getString("intneutronabs.abs2materialName");
    
    return ina;

  } // make()

} // namespace mu2e
