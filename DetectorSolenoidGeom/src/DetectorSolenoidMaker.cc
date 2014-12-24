// Mu2e includes
#include "BeamlineGeom/inc/Beamline.hh"
#include "BeamlineGeom/inc/StraightSection.hh"
#include "ConfigTools/inc/SimpleConfig.hh"
#include "DetectorSolenoidGeom/inc/DetectorSolenoidMaker.hh"
#include "DetectorSolenoidGeom/inc/DetectorSolenoid.hh"

// Framework includes
#include "cetlib/exception.h"

// CLHEP includes
#include "CLHEP/Units/SystemOfUnits.h"

// C++ includes
#include <iostream>
#include <vector>

using namespace std;

namespace mu2e {

  std::unique_ptr<DetectorSolenoid> DetectorSolenoidMaker::make(const SimpleConfig& c, const Beamline& bl ) {

    std::unique_ptr<DetectorSolenoid> ds ( new DetectorSolenoid() );

    // DS cryostat
    ds->_materialName = c.getString("ds.materialName");
    ds->_insideMaterialName = c.getString("ds.insideMaterialName");
    ds->_rIn1         = c.getDouble("ds.rIn.in");
    ds->_rIn2         = c.getDouble("ds.rIn.out");
    ds->_rOut1        = c.getDouble("ds.rOut.in");
    ds->_rOut2        = c.getDouble("ds.rOut.out");
    ds->_halfLength   = c.getDouble("ds.halfLength");
    ds->_endWallHalfLength = c.getDouble("ds.endWallHalfLength");
    ds->_frontHalfLength   = c.getDouble("ds.frontHalfLength");

    // DS shield
    ds->_shield_materialName       = c.getString("dsShield.materialName");
    ds->_shield_insideMaterialName = c.getString("dsShield.insideMaterialName");
    ds->_shield_zOffset            = c.getDouble("dsShield.zOffset");
    ds->_shield_halfLength         = c.getDouble("dsShield.halfLength");
    ds->_shield_endWallHalfLength  = c.getDouble("dsShield.endWallHalfLength");
    ds->_shield_rIn1               = c.getDouble("dsShield.rIn.in");
    ds->_shield_rIn2               = c.getDouble("dsShield.rIn.out");
    ds->_shield_rOut1              = c.getDouble("dsShield.rOut.in");
    ds->_shield_rOut2              = c.getDouble("dsShield.rOut.out");

    // DS solenoid coils
    ds->_coil_materialName = c.getString("dsCoil.materialName"); 
    ds->_coil_rIn          = c.getDouble("dsCoil.rIn"); 
    c.getVectorDouble("dsCoil.rOut"     , ds->_coil_rOut     , ds->nCoils() ); 
    c.getVectorDouble("dsCoil.zLength"  , ds->_coil_zLength  , ds->nCoils() ); 
    c.getVectorDouble("dsCoil.zPosition", ds->_coil_zPosition, ds->nCoils() ); 

    // DS coil support system
    ds->_support_materialName = c.getString("dsSupport.materialName"); 
    ds->_support_rIn          = c.getDouble("dsSupport.rIn"); 
    ds->_support_rOut         = c.getDouble("dsSupport.rOut"); 
    ds->_support_halfLength   = c.getDouble("dsSupport.halfLength"); 

    // Vacuum volumes
    ds->_vacuumMaterialName = c.getString("ds.vacuumMaterialName");
    ds->_vacuumPressure     = c.getDouble("ds.vacuumPressure");
    ds->_ds1HalfLength      = c.getDouble("ds1.halfLength");
    ds->_ds2HalfLength      = c.getDouble("ds2.halfLength");

    StraightSection const * ts5 = bl.getTS().getTSCryo<StraightSection>(TransportSolenoid::TSRegion::TS5,TransportSolenoid::TSRadialPart::IN );

    ds->_locationDs23Split  =  bl.getTS().torusRadius() 
      + 2.*ts5->getHalfLength() 
      + 2.*ds->vac_halfLengthDs2();

    // Position is computed on the fly, relative to the TS torus
    // radius, and the lengths of TS5 and the vacuum volumes
    // specified0; assumption is made that the front frace is flush
    // with the edge of the DS
    double dsPosZ      = bl.getTS().torusRadius() 
      + 2.*ts5->getHalfLength() 
      - 2.*ds->vac_halfLengthDs1()
      - 2.*ds->frontHalfLength()
      + ds->halfLength();

    // for x component: +1(-1)*solenoidOffset for PS (DS)
    ds->_position = CLHEP::Hep3Vector(-bl.solenoidOffset(), 0, dsPosZ );

    return ds;

  } // make()

} // namespace mu2e
