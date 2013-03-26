#include "BeamlineGeom/inc/Beamline.hh"
#include "DetectorSolenoidGeom/inc/DetectorSolenoidMaker.hh"
#include "DetectorSolenoidGeom/inc/DetectorSolenoid.hh"

#include "cetlib/exception.h"

#include "CLHEP/Units/SystemOfUnits.h"

#include "ConfigTools/inc/SimpleConfig.hh"

#include <vector>

namespace mu2e {

  std::unique_ptr<DetectorSolenoid> DetectorSolenoidMaker::make(const SimpleConfig& c, const Beamline& bl ) {

    std::unique_ptr<DetectorSolenoid> ds ( new DetectorSolenoid() );
    
    ds->_rIn        = c.getDouble("toyDS.rIn");
    ds->_rOut       = c.getDouble("toyDS.rOut");
    ds->_halfLength = c.getDouble("toyDS.halfLength");

    ds->_ds1HalfLength = c.getDouble("toyDS1.halfLength");
    ds->_ds2HalfLength = c.getDouble("toyDS2.halfLength");
    ds->_frontHalfLength = c.getDouble("toyDS.frontHalfLength");

    // Position is computed on the fly, relative to the TS torus
    // radius, and the lengths of TS5 and the vacuum volumes
    // specified0; assumption is made that the front frace is flesh
    // with the edge of the DS
    double dsPosZ      = bl.getTS().torusRadius() +
      2.*bl.getTS().getTS5().getHalfLength() -
      2.*ds->halfLengthDs1()-
      2.*ds->frontHalfLength()+
      ds->halfLength();

    // Different test
    dsPosZ = c.getDouble("toyDS.z0");

    // for x component: +1(-1)*solenoidOffset for PS (DS)
    ds->_position           = CLHEP::Hep3Vector(-bl.solenoidOffset(), 0, dsPosZ );

    ds->_materialName       = c.getString("toyDS.materialName");
    ds->_insideMaterialName = c.getString("toyDS.insideMaterialName");

    return ds;

  } // make()

} // namespace mu2e
