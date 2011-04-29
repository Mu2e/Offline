//
// Free function to create Muon Beam Stop and some elements of the Cryostat in G4
//
// $Id: constructMBS.cc,v 1.3 2011/04/29 17:43:23 genser Exp $
// $Author: genser $
// $Date: 2011/04/29 17:43:23 $
//
// Original author KLG 
//
// Notes: 
//
// The initial implementaion is described in Mu2e Document 1519 

// Note the interdependence of the position of the CryoSeal on
// the position of the neutron absorber to avoid overlaps, it should
// probably be more related to the SolenoidCoil once the dimensions of
// those components are reconciled

// clhep includes
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Units/PhysicalConstants.h"

// Mu2e includes.

#include "Mu2eG4/inc/constructMBS.hh"
#include "BeamlineGeom/inc/Beamline.hh"
#include "G4Helper/inc/VolumeInfo.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "CosmicRayShieldGeom/inc/CRSSteelShield.hh"
#include "G4Helper/inc/G4Helper.hh"
#include "Mu2eG4/inc/findMaterialOrThrow.hh"
#include "Mu2eG4/inc/nestBox.hh"
#include "Mu2eG4/inc/nestTubs.hh"
#include "Mu2eG4/inc/finishNesting.hh"

// G4 includes
#include "G4Material.hh"
#include "G4Color.hh"
#include "G4VSolid.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4SubtractionSolid.hh"
#include "G4VPhysicalVolume.hh"

using namespace std;

namespace mu2e {

  void constructMBS(SimpleConfig const * const _config){

    // BSTS is the main support tube; the other z positions are dependent upon it

    double const BSTSInnerRadius   = _config->getDouble("mbs.BSTSInnerRadius");
    double const BSTSOuterRadius   = _config->getDouble("mbs.BSTSOuterRadius");
    double const BSTSHLength       = _config->getDouble("mbs.BSTSHLength");
    string const BSTSMaterialName  = _config->getString("mbs.BSTSMaterialName");
    double const BSTSZ             = _config->getDouble("mbs.BSTSZ");
    double const SPBSInnerRadius   = _config->getDouble("mbs.SPBSInnerRadius");
    double const SPBSOuterRadius   = _config->getDouble("mbs.SPBSOuterRadius");
    double const SPBSHLength       = _config->getDouble("mbs.SPBSHLength");
    string const SPBSMaterialName  = _config->getString("mbs.SPBSMaterialName");
    double const SPBSZ             = BSTSZ + BSTSHLength - SPBSHLength;
    double const BSTCInnerRadius   = _config->getDouble("mbs.BSTCInnerRadius");
    double const BSTCOuterRadius   = _config->getDouble("mbs.BSTCOuterRadius");
    double const BSTCHLength       = _config->getDouble("mbs.BSTCHLength");
    string const BSTCMaterialName  = _config->getString("mbs.BSTCMaterialName");
    double const BSTCZ             = BSTSZ - BSTSHLength + BSTCHLength;
    double const BSBSInnerRadius   = _config->getDouble("mbs.BSBSInnerRadius");
    double const BSBSOuterRadius   = _config->getDouble("mbs.BSBSOuterRadius");
    double const BSBSHLength       = BSTSHLength - BSTCHLength;
    string const BSBSMaterialName  = _config->getString("mbs.BSBSMaterialName");
    double const BSBSZ             = BSTSZ + BSTSHLength - BSBSHLength;
    double const CLV2InnerRadius   = _config->getDouble("mbs.CLV2InnerRadius");
    double const CLV2OuterRadius   = _config->getDouble("mbs.CLV2OuterRadius");
    double const CLV2HLength       = _config->getDouble("mbs.CLV2HLength");
    string const CLV2MaterialName  = _config->getString("mbs.CLV2MaterialName");
    double const CLV2Z             = BSTSZ + BSTSHLength - CLV2HLength;

    bool const MBSisVisible        = _config->getBool("mbs.visible",true);
    bool const MBSisSolid          = _config->getBool("mbs.solid", false);
    int  const verbosityLevel           = _config->getInt("mbs.verbosityLevel", 0);

    bool const forceAuxEdgeVisible = _config->getBool("g4.forceAuxEdgeVisible",false);
    bool const doSurfaceCheck      = _config->getBool("g4.doSurfaceCheck",false);
    bool const placePV             = true;

    // Access to the G4HelperService.
    G4Helper* _helper = &(*(edm::Service<G4Helper>()));
    
    GeomHandle<Beamline> beamg;
    double solenoidOffset = -beamg->solenoidOffset(); 

 
    cout << __func__ << " verbosityLevel                          : " << verbosityLevel  << endl;

    if ( verbosityLevel > 0) {
      cout << __func__ << " solenoidOffset                   : " << solenoidOffset  << endl;
    }

     // mother volumes
    VolumeInfo const & detSolDownstreamVacInfo = _helper->locateVolInfo("ToyDS3Vacuum");
    VolumeInfo const & hallInfo =                _helper->locateVolInfo("HallAir");

    if ( verbosityLevel > 0) {
      double pzhl = static_cast<G4Tubs*>(detSolDownstreamVacInfo.solid->GetConstituentSolid(0))->GetZHalfLength();
      double pZOffset = detSolDownstreamVacInfo.centerInMu2e()[CLHEP::Hep3Vector::Z];
      cout << __func__ << " ToyDS3Vacuum   Offset in Mu2e    : " << 
        detSolDownstreamVacInfo.centerInMu2e() << endl;
      cout << __func__ << " ToyDS3Vacuum Z extent in Mu2e    : " << 
        pZOffset - pzhl << ", " << pZOffset + pzhl << endl;
    }

    // BSTS see WBS 5.8 for the naming conventions

    CLHEP::Hep3Vector BSTSOffsetInMu2e  = CLHEP::Hep3Vector(solenoidOffset,0.,BSTSZ);

    if ( verbosityLevel > 0) {
      cout << __func__ << " BSTSOffsetInMu2e                 : " << BSTSOffsetInMu2e << endl;
    }

    // now local offset in mother volume
    CLHEP::Hep3Vector BSTSOffset =  BSTSOffsetInMu2e - detSolDownstreamVacInfo.centerInMu2e();

    TubsParams BSTSParams ( BSTSInnerRadius,
                            BSTSOuterRadius,
                            BSTSHLength);

    VolumeInfo BSTSInfo  = nestTubs("BSTS",
                                    BSTSParams,
                                    findMaterialOrThrow(BSTSMaterialName),
                                    0,
                                    BSTSOffset,
                                    detSolDownstreamVacInfo,
                                    0, 
                                    MBSisVisible,
                                    G4Colour::Gray(),
                                    MBSisSolid,
                                    forceAuxEdgeVisible,
                                    placePV,
                                    doSurfaceCheck
                                    );

    if ( verbosityLevel > 0) {
      double zhl         = static_cast<G4Tubs*>(BSTSInfo.solid)->GetZHalfLength();
      double BSTSOffsetInMu2eZ = BSTSOffsetInMu2e[CLHEP::Hep3Vector::Z];
      cout << __func__ << " BSTSOffsetZ           in Mu2e    : " << 
        BSTSOffsetInMu2eZ << endl;
      cout << __func__ << " BSTS         Z extent in Mu2e    : " << 
        BSTSOffsetInMu2eZ - zhl << ", " << BSTSOffsetInMu2eZ + zhl << endl;
    }

    // SPBS

    CLHEP::Hep3Vector SPBSOffsetInMu2e  = CLHEP::Hep3Vector(solenoidOffset,0.,SPBSZ);

    if ( verbosityLevel > 0) {
      cout << __func__ << " SPBSOffsetInMu2e                 : " << SPBSOffsetInMu2e << endl;
    }

    // now local offset in mother volume
    CLHEP::Hep3Vector SPBSOffset =  SPBSOffsetInMu2e - detSolDownstreamVacInfo.centerInMu2e();

    TubsParams SPBSParams ( SPBSInnerRadius,
                            SPBSOuterRadius,
                            SPBSHLength);

    VolumeInfo SPBSInfo  = nestTubs("SPBS",
                                    SPBSParams,
                                    findMaterialOrThrow(SPBSMaterialName),
                                    0,
                                    SPBSOffset,
                                    detSolDownstreamVacInfo,
                                    0, 
                                    MBSisVisible,
                                    G4Colour::Blue(),
                                    MBSisSolid,
                                    forceAuxEdgeVisible,
                                    placePV,
                                    doSurfaceCheck
                                    );

    if ( verbosityLevel > 0) {
      double zhl         = static_cast<G4Tubs*>(SPBSInfo.solid)->GetZHalfLength();
      double SPBSOffsetInMu2eZ = SPBSOffsetInMu2e[CLHEP::Hep3Vector::Z];
      cout << __func__ << " SPBS         Z extent in Mu2e    : " << 
        SPBSOffsetInMu2eZ - zhl << ", " << SPBSOffsetInMu2eZ + zhl << endl;
    }

    // BSTC

    CLHEP::Hep3Vector BSTCOffsetInMu2e  = CLHEP::Hep3Vector(solenoidOffset,0.,BSTCZ);

    if ( verbosityLevel > 0) {
      cout << __func__ << " BSTCOffsetInMu2e                 : " << BSTCOffsetInMu2e << endl;
    }

    // now local offset in mother volume
    CLHEP::Hep3Vector BSTCOffset =  BSTCOffsetInMu2e - detSolDownstreamVacInfo.centerInMu2e();

    TubsParams BSTCParams ( BSTCInnerRadius,
                            BSTCOuterRadius,
                            BSTCHLength);

    G4Colour  orange  (.75, .55, .0);
    VolumeInfo BSTCInfo  = nestTubs("BSTC",
                                    BSTCParams,
                                    findMaterialOrThrow(BSTCMaterialName),
                                    0,
                                    BSTCOffset,
                                    detSolDownstreamVacInfo,
                                    0, 
                                    MBSisVisible,
                                    orange,
                                    MBSisSolid,
                                    forceAuxEdgeVisible,
                                    placePV,
                                    doSurfaceCheck
                                    );


    if ( verbosityLevel > 0) {
      double zhl         = static_cast<G4Tubs*>(BSTCInfo.solid)->GetZHalfLength();
      double BSTCOffsetInMu2eZ = BSTCOffsetInMu2e[CLHEP::Hep3Vector::Z];
      cout << __func__ << " BSTC         Z extent in Mu2e    : " << 
        BSTCOffsetInMu2eZ - zhl << ", " << BSTCOffsetInMu2eZ + zhl << endl;
    }

    // BSBS

    CLHEP::Hep3Vector BSBSOffsetInMu2e  = CLHEP::Hep3Vector(solenoidOffset,0.,BSBSZ);

    if ( verbosityLevel > 0) {
      cout << __func__ << " BSBSOffsetInMu2e                 : " << BSBSOffsetInMu2e << endl;
    }

    // now local offset in mother volume
    CLHEP::Hep3Vector BSBSOffset =  BSBSOffsetInMu2e - detSolDownstreamVacInfo.centerInMu2e();

    TubsParams BSBSParams ( BSBSInnerRadius,
                            BSBSOuterRadius,
                            BSBSHLength);

    VolumeInfo BSBSInfo  = nestTubs("BSBS",
                                    BSBSParams,
                                    findMaterialOrThrow(BSBSMaterialName),
                                    0,
                                    BSBSOffset,
                                    detSolDownstreamVacInfo,
                                    0, 
                                    MBSisVisible,
                                    G4Colour::Yellow(),
                                    MBSisSolid,
                                    forceAuxEdgeVisible,
                                    placePV,
                                    doSurfaceCheck
                                    );


    if ( verbosityLevel > 0) {
      double zhl         = static_cast<G4Tubs*>(BSBSInfo.solid)->GetZHalfLength();
      double BSBSOffsetInMu2eZ = BSBSOffsetInMu2e[CLHEP::Hep3Vector::Z];
      cout << __func__ << " BSBS         Z extent in Mu2e    : " << 
        BSBSOffsetInMu2eZ - zhl << ", " << BSBSOffsetInMu2eZ + zhl << endl;
    }

    // CLV2

    CLHEP::Hep3Vector CLV2OffsetInMu2e  = CLHEP::Hep3Vector(solenoidOffset,0.,CLV2Z);

    if ( verbosityLevel > 0) {
      cout << __func__ << " CLV2OffsetInMu2e                 : " << CLV2OffsetInMu2e << endl;
    }

    // now local offset in mother volume
    CLHEP::Hep3Vector CLV2Offset =  CLV2OffsetInMu2e - detSolDownstreamVacInfo.centerInMu2e();

    TubsParams CLV2Params ( CLV2InnerRadius,
                            CLV2OuterRadius,
                            CLV2HLength);

    VolumeInfo CLV2Info  = nestTubs("CLV2",
                                    CLV2Params,
                                    findMaterialOrThrow(CLV2MaterialName),
                                    0,
                                    CLV2Offset,
                                    detSolDownstreamVacInfo,
                                    0, 
                                    MBSisVisible,
                                    orange,
                                    MBSisSolid,
                                    forceAuxEdgeVisible,
                                    placePV,
                                    doSurfaceCheck
                                    );

    if ( verbosityLevel > 0) {
      double zhl         = static_cast<G4Tubs*>(CLV2Info.solid)->GetZHalfLength();
      double CLV2OffsetInMu2eZ = CLV2OffsetInMu2e[CLHEP::Hep3Vector::Z];
      cout << __func__ << " CLV2         Z extent in Mu2e    : " << 
        CLV2OffsetInMu2eZ - zhl << ", " << CLV2OffsetInMu2eZ + zhl << endl;
    }

    // constructing endplug closing the cryostat and enclosing the MBS

    // the hollow disk aka CryoSeal

    double CryoSealInnerRadius = SPBSOuterRadius; // may need to be BSTSOuterRadius TBD
    double CryoSealOuterRadius = _config->getDouble("toyDS.rIn");
    double CryoSealHLength     = _config->getDouble("mbs.CryoSealHLength");

    TubsParams CryoSealParams ( CryoSealInnerRadius,
                                CryoSealOuterRadius,
                                CryoSealHLength);

    double CryoSealZ = _config->getDouble("neutronabsorber.internalZ01") +
      _config->getDouble("neutronabsorber.internalHalfLengthZ01") + 
      2.*_config->getDouble("neutronabsorber.internalHalfLengthZ02") + CryoSealHLength;

    CLHEP::Hep3Vector CryoSealOffsetInMu2e = CLHEP::Hep3Vector(solenoidOffset,0.,CryoSealZ);

    CLHEP::Hep3Vector CryoSealOffset = CryoSealOffsetInMu2e - hallInfo.centerInMu2e();

    string const CryoSealMaterialName  = _config->getString("mbs.CryoSealMaterialName");

    VolumeInfo CryoSealInfo  = nestTubs("DSCryoSeal",
                                        CryoSealParams,
                                        findMaterialOrThrow(CryoSealMaterialName),
                                        0,
                                        CryoSealOffset,
                                        hallInfo,
                                        0, 
                                        MBSisVisible,
                                        G4Colour::Green(),
                                        MBSisSolid,
                                        forceAuxEdgeVisible,
                                        placePV,
                                        doSurfaceCheck
                                        );


    if ( verbosityLevel > 0) {
      double zhl         = static_cast<G4Tubs*>(CryoSealInfo.solid)->GetZHalfLength();
      double CryoSealOffsetInMu2eZ = CryoSealOffsetInMu2e[CLHEP::Hep3Vector::Z];
      cout << __func__ << " CryoSeal     Z extent in Mu2e    : " << 
        CryoSealOffsetInMu2eZ - zhl << ", " << CryoSealOffsetInMu2eZ + zhl << endl;
    }

    // now the endplug itself, the hollow part first

    double EndPlug1InnerRadius = SPBSOuterRadius;
    double EndPlug1OuterRadius = CryoSealOuterRadius;

    // from the end of the BSTS to the Iron   x 0.5

    double EndPlug1DSZ = BSTSZ + BSTSHLength;

    if ( verbosityLevel > 0) {

      cout << __func__ << " EndPlug1DSZ : " << EndPlug1DSZ << endl;

    }

    GeomHandle<CosmicRayShield> CosmicRayShieldGeomHandle;

    CRSSteelShield const & dssshield = 
      CosmicRayShieldGeomHandle->getCRSSteelShield("CRSSteelDownstreamShield");

    double EndPlug1USZ = dssshield.getHalfLengths()[2] + dssshield.getGlobalOffset()[CLHEP::Hep3Vector::Z];

    if ( verbosityLevel > 0) {

      cout << __func__ << " EndPlug1USZ : " << EndPlug1USZ << endl;

    }

    double EndPlug1HLength   = (EndPlug1DSZ - EndPlug1USZ)*0.5;

    TubsParams EndPlug1Params ( EndPlug1InnerRadius,
                               EndPlug1OuterRadius,
                               EndPlug1HLength );

    double EndPlug1Z = (EndPlug1DSZ + EndPlug1USZ)*0.5;

    if ( verbosityLevel > 0) {

      cout << __func__ << " EndPlug1Z  : " << EndPlug1Z << endl;

    }

    CLHEP::Hep3Vector EndPlug1OffsetInMu2e = CLHEP::Hep3Vector(solenoidOffset,0.,EndPlug1Z);

    CLHEP::Hep3Vector EndPlug1Offset = EndPlug1OffsetInMu2e - hallInfo.centerInMu2e();

    string const EndPlug1MaterialName  = _config->getString("mbs.EndPlugMaterialName");

    VolumeInfo EndPlug1Info  = nestTubs("DSEndPlug1",
                                        EndPlug1Params,
                                        findMaterialOrThrow(EndPlug1MaterialName),
                                        0,
                                        EndPlug1Offset,
                                        hallInfo,
                                        0, 
                                        MBSisVisible,
                                        G4Colour::Gray(),
                                        MBSisSolid,
                                        forceAuxEdgeVisible,
                                        placePV,
                                        doSurfaceCheck
                                        );


    if ( verbosityLevel > 0) {
      double zhl         = static_cast<G4Tubs*>(EndPlug1Info.solid)->GetZHalfLength();
      double EndPlug1OffsetInMu2eZ = EndPlug1OffsetInMu2e[CLHEP::Hep3Vector::Z];
      cout << __func__ << " EndPlug1     Z extent in Mu2e    : " << 
        EndPlug1OffsetInMu2eZ - zhl << ", " << EndPlug1OffsetInMu2eZ + zhl << endl;
    }

    // the end plug end disk 

    double EndPlugDiskInnerRadius = 0.0;
    double EndPlugDiskOuterRadius = CryoSealOuterRadius;

    double EndPlugDiskHLength     = _config->getDouble("mbs.EndPlugDiskHLength");

    TubsParams EndPlugDiskParams ( EndPlugDiskInnerRadius,
                                   EndPlugDiskOuterRadius,
                                   EndPlugDiskHLength );

    double EndPlugDiskZ = EndPlug1Z + EndPlug1HLength + EndPlugDiskHLength;

    if ( verbosityLevel > 0) {

      cout << __func__ << " EndPlugDiskZ  : " << EndPlugDiskZ << endl;

    }

    CLHEP::Hep3Vector EndPlugDiskOffsetInMu2e = CLHEP::Hep3Vector(solenoidOffset,0.,EndPlugDiskZ);

    CLHEP::Hep3Vector EndPlugDiskOffset = EndPlugDiskOffsetInMu2e - hallInfo.centerInMu2e();

    string const EndPlugDiskMaterialName  = _config->getString("mbs.EndPlugMaterialName");

    VolumeInfo EndPlugDiskInfo  = nestTubs("DSEndPlugDisk",
                                           EndPlugDiskParams,
                                           findMaterialOrThrow(EndPlugDiskMaterialName),
                                           0,
                                           EndPlugDiskOffset,
                                           hallInfo,
                                           0, 
                                           MBSisVisible,
                                           G4Colour::Gray(),
                                           MBSisSolid,
                                           forceAuxEdgeVisible,
                                           placePV,
                                           doSurfaceCheck
                                           );


    if ( verbosityLevel > 0) {
      double zhl         = static_cast<G4Tubs*>(EndPlugDiskInfo.solid)->GetZHalfLength();
      double EndPlugDiskOffsetInMu2eZ = EndPlugDiskOffsetInMu2e[CLHEP::Hep3Vector::Z];
      cout << __func__ << " EndPlugDisk  Z extent in Mu2e    : " << 
        EndPlugDiskOffsetInMu2eZ - zhl << ", " << EndPlugDiskOffsetInMu2eZ + zhl << endl;
    }


  } // end of constructMBS;

}
