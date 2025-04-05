// Mu2e includes
#include "Offline/ConfigTools/inc/SimpleConfig.hh"
#include "Offline/ConfigTools/inc/checkForStale.hh"
#include "Offline/DetectorSolenoidGeom/inc/DetectorSolenoid.hh"
#include "Offline/BeamlineGeom/inc/TSdA.hh"
#include "Offline/GeometryService/inc/TSdAMaker.hh"

// C++ includes
#include <algorithm>
#include <iostream>
#include <vector>

// CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"

// Other includes
#include "cetlib_except/exception.h"

namespace mu2e {

  std::unique_ptr<TSdA> TSdAMaker::make(const SimpleConfig& c, const DetectorSolenoid& ds ) {

    checkForStale( "intneutronabs", c );

    std::unique_ptr<TSdA> tsda ( new TSdA() );

    tsda->_r4          = c.getDouble("tsda.r4");

    tsda->_halfLength4 = c.getDouble("tsda.halfLength4");

    tsda->_position    = CLHEP::Hep3Vector( ds.position().x(),0,c.getDouble("tsda.z0"));

    tsda->_mat4        = c.getString("tsda.materialName");

    tsda->_version     = c.getInt("tsda.version",  1);
    tsda->_build       = c.getInt("tsda.build"  , -3);

    return tsda;

  } // make()

} // namespace mu2e
