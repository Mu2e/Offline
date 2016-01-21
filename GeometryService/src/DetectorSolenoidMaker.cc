// Mu2e includes
#include "BeamlineGeom/inc/Beamline.hh"
#include "BeamlineGeom/inc/StraightSection.hh"
#include "ConfigTools/inc/SimpleConfig.hh"
#include "DetectorSolenoidGeom/inc/DetectorSolenoid.hh"
#include "GeometryService/inc/DetectorSolenoidMaker.hh"

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

    // Rings (David Norvil Brown, May 2015)
    ds->_rInRingSide = c.getDouble("ds.rInRingSide");
    ds->_rOutRingSide = c.getDouble("ds.rOutRingSide");
    ds->_thickRingSide = c.getDouble("ds.thickRingSide");
    ds->_rInRing = c.getDouble("ds.rInRing");
    ds->_rOutRing = c.getDouble("ds.rOutRing");
    ds->_lengthRing = c.getDouble("ds.lengthRing");
    ds->_RingMaterial = c.getString("ds.RingMaterialType");
    c.getVectorDouble("ds.xRing", ds->_xRing, 2);
    c.getVectorDouble("ds.yRing", ds->_yRing, 2);
    c.getVectorDouble("ds.zRing", ds->_zRing, 2);

    // Rails
    int nPtRail = c.getInt("ds.nPtRail");
    c.getVectorDouble("ds.outlineU",ds->_uOutlineRail, nPtRail);
    c.getVectorDouble("ds.outlineV",ds->_vOutlineRail, nPtRail);
    ds->_lengthRail2 = c.getDouble("ds.lengthRail2");
    ds->_lengthRail3 = c.getDouble("ds.lengthRail3");
    ds->_RailMaterial = c.getString("ds.RailMaterialType");
    ds->_n2RailCenter = c.getHep3Vector("ds.n2RailCenter");
    ds->_s2RailCenter = c.getHep3Vector("ds.s2RailCenter");
    ds->_n3RailCenter = c.getHep3Vector("ds.n3RailCenter");
    ds->_s3RailCenter = c.getHep3Vector("ds.s3RailCenter");
    // Bearing blocks riding rails
    int nPtBBlock = c.getInt("ds.nPtBBlock");
    c.getVectorDouble("ds.outlineBBlockU",ds->_uOutlineBBlock, nPtBBlock);
    c.getVectorDouble("ds.outlineBBlockV",ds->_vOutlineBBlock, nPtBBlock);
    ds->_lengthBBlock2 = c.getDouble("ds.lengthBBlock2");
    ds->_lengthBBlock3 = c.getDouble("ds.lengthBBlock3");
    ds->_BBlockMaterial = c.getString("ds.BBlockMaterialType");
    int nBlocks = c.getInt("ds.nBBlocks",0);
    std::vector<double> xBBCs;
    xBBCs.reserve(nBlocks);
    std::vector<double> yBBCs;
    yBBCs.reserve(nBlocks);
    std::vector<double> zBBCs;
    xBBCs.reserve(nBlocks);
    c.getVectorDouble("ds.xCentersBBlock",xBBCs,nBlocks);
    c.getVectorDouble("ds.yCentersBBlock",yBBCs,nBlocks);
    c.getVectorDouble("ds.zCentersBBlock",zBBCs,nBlocks);
    for ( int iBB = 0; iBB < nBlocks; iBB++ ) {
      CLHEP::Hep3Vector center(xBBCs[iBB],yBBCs[iBB],zBBCs[iBB]);
      if ( zBBCs[iBB]*CLHEP::mm > 8340.0 ) {
	ds->_BBlockCenters3.push_back(center);
      } else {
	ds->_BBlockCenters2.push_back(center);
      }
    }

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
