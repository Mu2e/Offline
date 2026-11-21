// Mu2e includes
#include "Offline/ConfigTools/inc/SimpleConfig.hh"
#include "Offline/ConfigTools/inc/checkForStale.hh"
#include "Offline/DetectorSolenoidGeom/inc/DetectorSolenoid.hh"
#include "Offline/BeamlineGeom/inc/TSdA.hh"
#include "Offline/BeamlineGeom/inc/Beamline.hh"
#include "Offline/BeamlineGeom/inc/StraightSection.hh"
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

  std::unique_ptr<TSdA> TSdAMaker::make(const SimpleConfig& c, const DetectorSolenoid& ds, const Beamline& beamline ) {

    checkForStale( "intneutronabs", c );

    std::unique_ptr<TSdA> tsda ( new TSdA() );

    // inner radius: the TS5 inner cryostat radius if not provided
    tsda->_rIn         = c.getDouble("tsda.rin",0.0);
    if ( tsda->_rIn < 1.0e-06 ) {
      StraightSection const* ts5out = beamline.getTS().getTSCryo<StraightSection>( TransportSolenoid::TSRegion::TS5,
                                                                                    TransportSolenoid::TSRadialPart::IN );
      tsda->_rIn = ts5out->rIn();
    }

    tsda->_r4          = c.getDouble("tsda.r4");

    tsda->_halfLength4 = c.getDouble("tsda.halfLength4");

    tsda->_position    = CLHEP::Hep3Vector( ds.position().x(),0,c.getDouble("tsda.z0"));

    tsda->_mat4        = c.getString("tsda.materialName");

    tsda->_version     = c.getInt("tsda.version", 1);

    return tsda;

  } // make()

} // namespace mu2e
