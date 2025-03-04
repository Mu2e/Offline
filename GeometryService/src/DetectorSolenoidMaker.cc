// Mu2e includes
#include "Offline/BeamlineGeom/inc/Beamline.hh"
#include "Offline/BeamlineGeom/inc/StraightSection.hh"
#include "Offline/ConfigTools/inc/SimpleConfig.hh"
#include "Offline/DetectorSolenoidGeom/inc/DetectorSolenoid.hh"
#include "Offline/GeometryService/inc/DetectorSolenoidMaker.hh"

// Framework includes
#include "cetlib_except/exception.h"

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

    // experimental inner lining for shielding
    ds->_hasInnerLining = c.getBool("ds.hasInnerLining",false);
    ds->_innerLiningThickness = c.getDouble("ds.innerLiningThickness",0.0)*CLHEP::mm;
    ds->_innerLiningMaterial = c.getString("ds.innerLiningMaterial","None");

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
    ds->_coil_rIn          = c.getDouble("dsCoil.rIn");
    c.getVectorDouble("dsCoil.rOut"     , ds->_coil_rOut     , ds->nCoils() );
    c.getVectorDouble("dsCoil.zLength"  , ds->_coil_zLength  , ds->nCoils() );
    c.getVectorDouble("dsCoil.zPosition", ds->_coil_zPosition, ds->nCoils() );
    ds->_coilVersion        = c.getInt("dsCoil.version",1);
    if ( ds->_coilVersion == 1 ) {
      ds->_coil_materialName = c.getString("dsCoil.materialName");
    } else {
      c.getVectorString("dsCoil.materialNameVector", ds->_coil_mats , ds->nCoils() );
      // DS coil spacers
      ds->_spacer_materialName = c.getString("dsSpacer.materialName");
      ds->_spacer_rIn          = c.getDouble("dsSpacer.rIn");
      c.getVectorDouble("dsSpacer.rOut"     , ds->_spacer_rOut     , ds->nSpacers() );
      c.getVectorDouble("dsSpacer.zLength"  , ds->_spacer_zLength  , ds->nSpacers() );
      c.getVectorDouble("dsSpacer.zPosition", ds->_spacer_zPosition, ds->nSpacers() );
    }

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
    ds->_widthCoupler = c.getDouble("ds.widthCoupler",0.0);
    ds->_heightCoupler = c.getDouble("ds.heightCoupler",0.0);
    ds->_yCenterCoupler = c.getDouble("ds.yCenterCoupler",0.0);
    ds->_couplerScheme = c.getInt("ds.couplerScheme",-1);

    // MBS spherical support structure
    ds->_hasMBSS = c.getBool("ds.hasMBSSupport",false);
    if ( ds->_hasMBSS ) {
      ds->_lengthMBSS = c.getDouble("ds.MBSSupport.length");
      int nPtMBSS = c.getInt("ds.MBSSupport.nVertices");
      c.getVectorDouble("ds.MBSSupport.outlineU",ds->_uOutlineMBSS, nPtMBSS );
      c.getVectorDouble("ds.MBSSupport.outlineV",ds->_vOutlineMBSS, nPtMBSS );
      ds->_locationMBSS = c.getHep3Vector("ds.MBSSupport.location");
      ds->_materialMBSS = c.getString("ds.MBSSupport.material");
    }

    // Cable Runs
    ds->_hasCableRunCal = c.getBool("ds.hasCableRunCal",false);
    ds->_hasCableRunTrk = c.getBool("ds.hasCableRunTrk",false);
    if ( ds->_hasCableRunCal ) {
      ds->_cableRunVersion = c.getInt("ds.CableRun.version",1);
      if ( ds->_cableRunVersion > 1 ) {
        ds->_upRInCableRunCal  = c.getDouble("ds.CableRunCal.UpRin");
        ds->_upROutCableRunCal = c.getDouble("ds.CableRunCal.UpRout");
        ds->_upHL2CableRunCal  = c.getDouble("ds.CableRunCal.UpHL2");
        ds->_upZC2CableRunCal  = c.getDouble("ds.CableRunCal.UpZC2");
      }
      if ( ds->_cableRunVersion > 2 ) {
        //cable core parameters
        ds->_rCableRunCalCoreFract    = c.getDouble("ds.CableRunCalCore.RadiusFraction");
        ds->_rdCableRunCalCoreFract   = c.getDouble("ds.CableRunCalCore.dRadiusFraction");
        ds->_dPhiCableRunCalCoreFract = c.getDouble("ds.CableRunCalCore.dPhiFraction");
        ds->_materialCableRunCalCore  = c.getString("ds.CableRunCalCore.material");

        //IFB cabling
        ds->_calR1CableRunIFB    = c.getDouble("ds.CableRunIFB.CalR1");
        ds->_calR2CableRunIFB    = c.getDouble("ds.CableRunIFB.CalR2");
        ds->_calPhi0CableRunIFB  = c.getDouble("ds.CableRunIFB.CalPhi0");
        ds->_calDPhiCableRunIFB  = c.getDouble("ds.CableRunIFB.CalDPhi");
        ds->_calREndCableRunIFB  = c.getDouble("ds.CableRunIFB.CalREnd");
        ds->_calEndWCableRunIFB  = c.getDouble("ds.CableRunIFB.CalEndW");
        ds->_calPhiECableRunIFB  = c.getDouble("ds.CableRunIFB.CalPhiE");
        //IFB patch panel
        ds->_calPR1CableRunIFB   = c.getDouble("ds.CableRunIFB.CalPR1");
        ds->_calPR2CableRunIFB   = c.getDouble("ds.CableRunIFB.CalPR2");
        ds->_calPPhi0CableRunIFB = c.getDouble("ds.CableRunIFB.CalPPhi0");
        ds->_calPDPhiCableRunIFB = c.getDouble("ds.CableRunIFB.CalPDPhi");
        ds->_calPZInCableRunIFB  = c.getDouble("ds.CableRunIFB.CalPZIn");
        ds->_calPZHLCableRunIFB  = c.getDouble("ds.CableRunIFB.CalPZHL");
        ds->_calPZOutCableRunIFB = c.getDouble("ds.CableRunIFB.CalPZOut");
        ds->_calPMatCableRunIFB  = c.getString("ds.CableRunIFB.CalPMat");
        //cabling at bottom of IFB cabling
        ds->_calBCXCableRunIFB   = c.getDouble("ds.CableRunIFB.CalBCX");
        ds->_calBLCableRunIFB    = c.getDouble("ds.CableRunIFB.CalBL");

        ds->_trkR1CableRunIFB    = c.getDouble("ds.CableRunIFB.TrkR1");
        ds->_trkR2CableRunIFB    = c.getDouble("ds.CableRunIFB.TrkR2");
        ds->_trkPhi0CableRunIFB  = c.getDouble("ds.CableRunIFB.TrkPhi0");
        ds->_trkDPhiCableRunIFB  = c.getDouble("ds.CableRunIFB.TrkDPhi");
        ds->_trkREndCableRunIFB  = c.getDouble("ds.CableRunIFB.TrkREnd");
        ds->_trkEndWCableRunIFB  = c.getDouble("ds.CableRunIFB.TrkEndW");
        ds->_trkPhiECableRunIFB  = c.getDouble("ds.CableRunIFB.TrkPhiE");
        //IFB patch panel
        ds->_trkPR1CableRunIFB   = c.getDouble("ds.CableRunIFB.TrkPR1");
        ds->_trkPR2CableRunIFB   = c.getDouble("ds.CableRunIFB.TrkPR2");
        ds->_trkPPhi0CableRunIFB = c.getDouble("ds.CableRunIFB.TrkPPhi0");
        ds->_trkPDPhiCableRunIFB = c.getDouble("ds.CableRunIFB.TrkPDPhi");
        ds->_trkPZInCableRunIFB  = c.getDouble("ds.CableRunIFB.TrkPZIn");
        ds->_trkPZHLCableRunIFB  = c.getDouble("ds.CableRunIFB.TrkPZHL");
        ds->_trkPZOutCableRunIFB = c.getDouble("ds.CableRunIFB.TrkPZOut");
        ds->_trkPMatCableRunIFB  = c.getString("ds.CableRunIFB.TrkPMat");
        //cabling at bottom of IFB cabling
        ds->_trkBCXCableRunIFB   = c.getDouble("ds.CableRunIFB.TrkBCX");
        ds->_trkBLCableRunIFB    = c.getDouble("ds.CableRunIFB.TrkBL");

        ds->_zHLCableRunIFB      = c.getDouble("ds.CableRunIFB.ZHL");
        ds->_materialCalCableRunIFB = c.getString("ds.CalCableRunIFB.Material");
        ds->_materialTrkCableRunIFB = c.getString("ds.TrkCableRunIFB.Material");
        ds->_zCCableRunIFB       = c.getDouble("ds.CableRunIFB.ZC");
      }
      ds->_lengthCableRunCal = c.getDouble("ds.CableRunCal.length");
      ds->_rInCableRunCal    = c.getDouble("ds.CableRunCal.Rin");
      ds->_rOutCableRunCal   = c.getDouble("ds.CableRunCal.Rout");
      ds->_dPhiCableRunCal   = c.getDouble("ds.CableRunCal.dPhi");
      ds->_zCCableRunCal     = c.getDouble("ds.CableRunCal.zC"  );
      ds->_phi0CableRunCal   = c.getDouble("ds.CableRunCal.phi0");
      ds->_materialCableRunCal = c.getString("ds.CableRunCal.material");

    }
    if ( ds->_hasCableRunTrk ) {
      ds->_lengthCableRunTrk = c.getDouble("ds.CableRunTrk.length");
      ds->_rInCableRunTrk    = c.getDouble("ds.CableRunTrk.Rin");
      ds->_rOutCableRunTrk   = c.getDouble("ds.CableRunTrk.Rout");
      ds->_dPhiCableRunTrk   = c.getDouble("ds.CableRunTrk.dPhi");
      ds->_zCCableRunTrk     = c.getDouble("ds.CableRunTrk.zC"  );
      ds->_phi0CableRunTrk   = c.getDouble("ds.CableRunTrk.phi0");
      ds->_materialCableRunTrk = c.getString("ds.CableRunTrk.material");

      if ( ds->_cableRunVersion > 2 ) {
        ds->_rCableRunTrkCoreFract    = c.getDouble("ds.CableRunTrkCore.RadiusFraction");
        ds->_rdCableRunTrkCoreFract   = c.getDouble("ds.CableRunTrkCore.dRadiusFraction");
        ds->_dPhiCableRunTrkCoreFract = c.getDouble("ds.CableRunTrkCore.dPhiFraction");
        ds->_materialCableRunTrkCore  = c.getString("ds.CableRunTrkCore.material");
      }

    }

    // Service pipes
    bool hasServicePipes = c.getBool("ds.hasServicePipes",false);
    ds->_hasServicePipes = hasServicePipes;
    if ( hasServicePipes ) {
      ds->_servicePipeRIn = c.getDouble("ds.servicePipeRIn");
      ds->_servicePipeROut = c.getDouble("ds.servicePipeROut");
      ds->_servicePipeHL  = c.getDouble("ds.servicePipeHL");
      ds->_servicePipeMat = c.getString("ds.servicePipeMat");
      ds->_servicePipeFillMat = c.getString("ds.servicePipeFillMat");
      ds->_servicePipeZC = c.getDouble("ds.servicePipeZC");
      ds->_servicePipeYC = c.getDouble("ds.servicePipeYC");
      c.getVectorDouble("ds.servicePipeXCs",ds->_servicePipeXCs);
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
