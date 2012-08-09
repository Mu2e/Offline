//
// Free function to create Muon Beam Stop and some elements of the Cryostat in G4
//
// $Id: constructMBS.cc,v 1.11 2012/08/09 22:22:25 kutschke Exp $
// $Author: kutschke $
// $Date: 2012/08/09 22:22:25 $
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
#include "MBSGeom/inc/MBS.hh"
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

    MBS const & mbsgh = *(GeomHandle<MBS>());

    Tube const & pBSTSParams = *mbsgh.getBSTSPtr();
    Tube const & pSPBSParams = *mbsgh.getSPBSPtr();
    Tube const & pBSTCParams = *mbsgh.getBSTCPtr();
    Tube const & pBSBSParams = *mbsgh.getBSBSPtr();
    Tube const & pCLV2Params = *mbsgh.getCLV2Ptr();

    bool const MBSisVisible        = _config->getBool("mbs.visible",true);
    bool const MBSisSolid          = _config->getBool("mbs.solid", false);
    int  const verbosityLevel      = _config->getInt("mbs.verbosityLevel", 0);

    bool const forceAuxEdgeVisible = _config->getBool("g4.forceAuxEdgeVisible",false);
    bool const doSurfaceCheck      = _config->getBool("g4.doSurfaceCheck",false);
    bool const placePV             = true;

    // Access to the G4HelperService.
    G4Helper* _helper = &(*(art::ServiceHandle<G4Helper>()));

    GeomHandle<Beamline> beamg;
    double solenoidOffset = -beamg->solenoidOffset();

    if ( verbosityLevel > 0) {
      cout << __func__ << " verbosityLevel                   : " << verbosityLevel  << endl;
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

    CLHEP::Hep3Vector BSTSOffsetInMu2e = pBSTSParams.originInMu2e();

    if ( verbosityLevel > 0) {
      cout << __func__ << " BSTSOffsetInMu2e                 : " << BSTSOffsetInMu2e << endl;
    }

    // now local offset in mother volume
    CLHEP::Hep3Vector BSTSOffset =  BSTSOffsetInMu2e - detSolDownstreamVacInfo.centerInMu2e();

    // Use BSTS as mother volume for MBS

    /*
    TubsParams MBSMotherParams ( 0,
				  BSTSOuterRadius,
				  BSTSHLength);
    G4Material* vacuumMaterial     = findMaterialOrThrow(_config->getString("toyDS.insideMaterialName"));
    */

    /*
    VolumeInfo MBSMotherInfo  = nestTubs("MBSMother",
					 MBSMotherParams,
					 vacuumMaterial,
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
    */

    VolumeInfo BSTSInfo  = nestTubs("BSTS",
                                    pBSTSParams.getTubsParams(),
                                    findMaterialOrThrow(pBSTSParams.materialName()),
                                    0,
				    BSTSOffset,
				    //G4ThreeVector(0,0,0),
				    detSolDownstreamVacInfo,
				    //MBSMotherInfo,
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
    // This one is placed directly into DS3Vacuum; 
    // Note that its radius is lager than theone of BSTS
    
    CLHEP::Hep3Vector SPBSOffsetInMu2e = pSPBSParams.originInMu2e();

    if ( verbosityLevel > 0) {
      cout << __func__ << " SPBSOffsetInMu2e                 : " << SPBSOffsetInMu2e << endl;
    }

    // now local offset in mother volume
    CLHEP::Hep3Vector SPBSOffset = SPBSOffsetInMu2e - detSolDownstreamVacInfo.centerInMu2e();

    VolumeInfo SPBSInfo  = nestTubs("SPBS",
                                    pSPBSParams.getTubsParams(),
                                    findMaterialOrThrow(pSPBSParams.materialName()),
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

    CLHEP::Hep3Vector BSTCOffsetInMu2e = pBSTCParams.originInMu2e();
    // now local offset in mother volume
    CLHEP::Hep3Vector BSTCOffset =  BSTCOffsetInMu2e - detSolDownstreamVacInfo.centerInMu2e(); // - MBSMotherInfo.centerInMu2e(); 

    if ( verbosityLevel > 0) {
      cout << __func__ << " BSTCOffsetInMu2e                 : " << BSTCOffsetInMu2e << endl;
      cout << __func__ << " BSTCOffsetInMBS                  : " << BSTCOffset << endl;
    }

    G4Colour  orange  (.75, .55, .0);
    VolumeInfo BSTCInfo  = nestTubs("BSTC",
                                    pBSTCParams.getTubsParams(),
                                    findMaterialOrThrow(pBSTCParams.materialName()),
                                    0,
                                    BSTCOffset,
                                    //MBSMotherInfo,
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

    CLHEP::Hep3Vector BSBSOffsetInMu2e = pBSBSParams.originInMu2e();

    // now local offset in mother volume
    CLHEP::Hep3Vector BSBSOffset =  BSBSOffsetInMu2e - detSolDownstreamVacInfo.centerInMu2e(); //-MBSMotherInfo.centerInMu2e();

    if ( verbosityLevel > 0) {
      cout << __func__ << " BSBSOffsetInMu2e                 : " << BSBSOffsetInMu2e << endl;
      cout << __func__ << " BSBSOffsetInMBS                  : " << BSBSOffset << endl;
    }

    VolumeInfo BSBSInfo  = nestTubs("BSBS",
                                    pBSBSParams.getTubsParams(),
                                    findMaterialOrThrow(pBSBSParams.materialName()),
                                    0,
                                    BSBSOffset,
                                    //MBSMotherInfo,
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

    CLHEP::Hep3Vector CLV2OffsetInMu2e = pCLV2Params.originInMu2e();

    // now local offset in mother volume
    CLHEP::Hep3Vector CLV2Offset = CLV2OffsetInMu2e - detSolDownstreamVacInfo.centerInMu2e();//- MBSMotherInfo.centerInMu2e(); 

    if ( verbosityLevel > 0) {
      cout << __func__ << " CLV2OffsetInMu2e                 : " << CLV2OffsetInMu2e << endl;
      cout << __func__ << " CLV2OffsetInMBS                  : " << CLV2Offset << endl;
    }

    VolumeInfo CLV2Info  = nestTubs("CLV2",
                                    pCLV2Params.getTubsParams(),
                                    findMaterialOrThrow(pCLV2Params.materialName()),
                                    0,
                                    CLV2Offset,
                                    detSolDownstreamVacInfo,
                                    //MBSMotherInfo,
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

    bool const hasCryoSeal = _config->getBool("mbs.hasCryoSeal", true);

    if (hasCryoSeal) {

      double const CryoSealInnerRadius = _config->getDouble("mbs.CryoSealInnerRadius");
      double const CryoSealOuterRadius = _config->getDouble("mbs.CryoSealOuterRadius");
      double const CryoSealHLength     = _config->getDouble("mbs.CryoSealHLength");

      TubsParams CryoSealParams ( CryoSealInnerRadius,
                                  CryoSealOuterRadius,
                                  CryoSealHLength);

      double CryoSealZ           = _config->getDouble("mbs.CryoSealZ");

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

    }

    bool const hasEndPlug = _config->getBool("mbs.hasEndPlug", true);

    // now the endplug itself, the hollow part first

    if (hasEndPlug) {

      double const EndPlugTubeInnerRadius = _config->getDouble("mbs.EndPlugTubeInnerRadius");
      double const EndPlugTubeOuterRadius = _config->getDouble("mbs.EndPlugTubeOuterRadius");
      double const EndPlugTubeHLength     = _config->getDouble("mbs.EndPlugTubeHLength");

      TubsParams EndPlugTubeParams ( EndPlugTubeInnerRadius,
                                     EndPlugTubeOuterRadius,
                                     EndPlugTubeHLength );

      double const EndPlugTubeZ           = _config->getDouble("mbs.EndPlugTubeZ");

      if ( verbosityLevel > 0) {

        double EndPlugTubeDSZ = EndPlugTubeZ + EndPlugTubeHLength;
        double EndPlugTubeUSZ = EndPlugTubeZ - EndPlugTubeHLength;
        cout << __func__ << " EndPlugTubeDSZ : " << EndPlugTubeDSZ << endl;
        cout << __func__ << " EndPlugTubeUSZ : " << EndPlugTubeUSZ << endl;
        cout << __func__ << " EndPlugTubeZ   : " << EndPlugTubeZ << endl;

      }

      CLHEP::Hep3Vector EndPlugTubeOffsetInMu2e = CLHEP::Hep3Vector(solenoidOffset,0.,EndPlugTubeZ);

      CLHEP::Hep3Vector EndPlugTubeOffset = EndPlugTubeOffsetInMu2e - hallInfo.centerInMu2e();

      string const EndPlugTubeMaterialName  = _config->getString("mbs.EndPlugMaterialName");

      VolumeInfo EndPlugTubeInfo  = nestTubs("DSEndPlugTube",
                                             EndPlugTubeParams,
                                             findMaterialOrThrow(EndPlugTubeMaterialName),
                                             0,
                                             EndPlugTubeOffset,
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
        double zhl         = static_cast<G4Tubs*>(EndPlugTubeInfo.solid)->GetZHalfLength();
        double EndPlugTubeOffsetInMu2eZ = EndPlugTubeOffsetInMu2e[CLHEP::Hep3Vector::Z];
        cout << __func__ << " EndPlugTube     Z extent in Mu2e    : " <<
          EndPlugTubeOffsetInMu2eZ - zhl << ", " << EndPlugTubeOffsetInMu2eZ + zhl << endl;
      }

      // the end plug end disk

      double const EndPlugDiskInnerRadius  = _config->getDouble("mbs.EndPlugDiskInnerRadius");
      double const EndPlugDiskOuterRadius  = _config->getDouble("mbs.EndPlugDiskOuterRadius");
      double const EndPlugDiskHLength      = _config->getDouble("mbs.EndPlugDiskHLength");

      TubsParams EndPlugDiskParams ( EndPlugDiskInnerRadius,
                                     EndPlugDiskOuterRadius,
                                     EndPlugDiskHLength );

      double EndPlugDiskZ = EndPlugTubeZ + EndPlugTubeHLength + EndPlugDiskHLength;

      if ( verbosityLevel > 0) {

        cout << __func__ << " EndPlugDiskZ  : " << EndPlugDiskZ << endl;

      }

      CLHEP::Hep3Vector EndPlugDiskOffsetInMu2e = CLHEP::Hep3Vector(solenoidOffset,0.,EndPlugDiskZ);

      CLHEP::Hep3Vector EndPlugDiskOffset  = EndPlugDiskOffsetInMu2e - hallInfo.centerInMu2e();

      string const EndPlugDiskMaterialName = _config->getString("mbs.EndPlugMaterialName");

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

    }

  } // end of constructMBS;

}
