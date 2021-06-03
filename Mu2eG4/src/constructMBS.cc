//
// Free function to create Muon Beam Stop and some elements of the Cryostat in G4
//
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

// Update by D. No. Brown, Louisville, 30/10/2015.  Updated version is 
// "Version 2", based on v7 of Doc 1351.  In "Version 3", Anthony Palladino
// widened the downstream beam exit for the STM design.
// Version 4 adds pump-out holes in the bottom of the MBS to Version 2.
// Version 5 does the same thing to Version 3.  
// At the time of the latter modification, version 4 is default.  Once 
// the design of the CRV and STM are finalized, version 5 will be default.

// art includes
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"

// clhep includes
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Units/PhysicalConstants.h"

// Mu2e includes.

#include "Mu2eG4/inc/constructMBS.hh"
#include "DetectorSolenoidGeom/inc/DetectorSolenoid.hh"
#include "Mu2eG4Helper/inc/VolumeInfo.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/G4GeometryOptions.hh"
#include "CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "MBSGeom/inc/MBS.hh"
#include "Mu2eG4Helper/inc/Mu2eG4Helper.hh"
#include "Mu2eG4/inc/findMaterialOrThrow.hh"
#include "Mu2eG4/inc/nestBox.hh"
#include "Mu2eG4/inc/nestTubs.hh"
#include "Mu2eG4/inc/nestPolycone.hh"
#include "Mu2eG4/inc/finishNesting.hh"

// G4 includes
#include "Geant4/G4Material.hh"
#include "Geant4/G4Box.hh"
#include "Geant4/G4Color.hh"
#include "Geant4/G4VSolid.hh"
#include "Geant4/G4Tubs.hh"
#include "Geant4/G4Polycone.hh"
#include "Geant4/G4Cons.hh"
#include "Geant4/G4SubtractionSolid.hh"
#include "Geant4/G4VPhysicalVolume.hh"

#include <iostream>
using namespace std;

namespace mu2e {

  void constructMBS(SimpleConfig const & _config){

    // BSTS is the main support tube; the other z positions are dependent upon it

    MBS const & mbsgh = *(GeomHandle<MBS>());
    int MBSversion = mbsgh.getVersion();

    Polycone const & pMBSMParams = *mbsgh.getMBSMPtr();
    Polycone const & pBSTSParams = *mbsgh.getBSTSPtr();
    Tube const & pSPBSCParams    = *mbsgh.getSPBSCPtr();
    Polycone const & pBSTCParams = *mbsgh.getBSTCPtr();
    Polycone const & pBSBSParams = *mbsgh.getBSBSPtr();
    Polycone const & pCLV2Params = *mbsgh.getCLV2Ptr();
    Tube const & pSPBSSup1Params = *mbsgh.getSPBSSup1Ptr();
    Tube const & pSPBSSup2Params = *mbsgh.getSPBSSup2Ptr();
    Tube const & pSPBSLParams    = *mbsgh.getSPBSLPtr();
    Tube const & pSPBSRParams    = *mbsgh.getSPBSRPtr();
    Tube const & pCLV2ABSParams  = *mbsgh.getCLV2ABSPtr();
    Tube const & pCalShieldRingParams  = *mbsgh.getCalRingShieldPtr();


    const auto geomOptions = art::ServiceHandle<GeometryService>()->geomOptions();
    geomOptions->loadEntry( _config, "mbs", "mbs" );

    const bool MBSisVisible        = geomOptions->isVisible("mbs"); 
    const bool MBSisSolid          = geomOptions->isSolid("mbs"); 
    const bool forceAuxEdgeVisible = geomOptions->forceAuxEdgeVisible("mbs"); 
    const bool doSurfaceCheck      = geomOptions->doSurfaceCheck("mbs"); 
    const bool placePV             = geomOptions->placePV("mbs"); 
    int  const verbosityLevel      = _config.getInt("mbs.verbosityLevel", 0);



    // Access to the Mu2eG4HelperService.
    Mu2eG4Helper* _helper = &(*(art::ServiceHandle<Mu2eG4Helper>()));

    // Fetch DS geometry handle
    GeomHandle<DetectorSolenoid> ds;

    if ( verbosityLevel > 0) {
      cout << __func__ << " verbosityLevel                   : " << verbosityLevel  << endl;
      cout << __func__ << " solenoidOffset                   : " << ds->position().x()  << endl;
    }

     // mother volumes
    std::string theDS3("DS3Vacuum");
    if ( _config.getBool("inGaragePosition",false) ) theDS3 = "garageFakeDS3Vacuum";
    VolumeInfo const & detSolDownstreamVacInfo = _helper->locateVolInfo(theDS3);

    if ( verbosityLevel > 0) {
      cout << __func__ << " z-extent of DS3Vacuum portion in DS in Mu2e  : " <<
        ds->vac_zLocDs23Split() << ", " << ds->cryoZMax() << endl;
    }


    // MBSM

    CLHEP::Hep3Vector MBSMOffsetInMu2e = pMBSMParams.originInMu2e();
    // now local offset in mother volume
    CLHEP::Hep3Vector MBSMOffset =  MBSMOffsetInMu2e - detSolDownstreamVacInfo.centerInMu2e(); // - MBSMotherInfo.centerInMu2e();

    if ( verbosityLevel > 0) {
      cout << __func__ << " MBSMotherOffsetInMu2e                 : " << MBSMOffsetInMu2e << endl;
      cout << __func__ << " MBSMotherOffsetInMBS                  : " << MBSMOffset << endl;
    }


    VolumeInfo MBSMotherInfo  = nestPolycone("MBSMother",
                                    pMBSMParams.getPolyconsParams(),
                                    findMaterialOrThrow(pMBSMParams.materialName()),
                                    0,
                                    MBSMOffset,
                                    detSolDownstreamVacInfo,
                                    0,
                                    MBSisVisible,
                                    0,
                                    MBSisSolid,
                                    forceAuxEdgeVisible,
                                    placePV,
                                    doSurfaceCheck
                                    );


    if ( verbosityLevel > 0) {
      G4Polycone *tmpMBSMInfo  =  static_cast<G4Polycone*>(MBSMotherInfo.solid);
      double zhl = tmpMBSMInfo->GetCorner(tmpMBSMInfo->GetNumRZCorner()-1).z-tmpMBSMInfo->GetCorner(0).z;
      zhl*=0.5;
      double MBSMOffsetInMu2eZ = MBSMOffsetInMu2e[CLHEP::Hep3Vector::Z];
      cout << __func__ << " MBSM         Z extent in Mu2e    : " <<
        MBSMOffsetInMu2eZ - zhl << ", " << MBSMOffsetInMu2eZ + zhl << endl;
    }

    // BSTS see WBS 5.8 for the naming conventions
    // Stainless steel pipe - in Version 1, a single tube rendered as 
    // a polycone.  In Version 2, a polycone with 3 segments of varying
    // thickness.

    CLHEP::Hep3Vector BSTSOffsetInMu2e = pBSTSParams.originInMu2e();

    if ( verbosityLevel > 0) {
      cout << __func__ << " BSTSOffsetInMu2e                 : " << BSTSOffsetInMu2e << endl;
    }

    // now local offset in mother volume
    CLHEP::Hep3Vector BSTSOffset =  BSTSOffsetInMu2e - MBSMOffsetInMu2e;


    // Now build steel.  For versions 1-3, no holes.  For 4 and up, holes.
    if ( MBSversion < 4 ) {
      // Use BSTS as mother volume for MBS (DNB:  ??)
      VolumeInfo BSTSInfo  = nestPolycone("BSTS",
			       pBSTSParams.getPolyconsParams(),
			       findMaterialOrThrow(pBSTSParams.materialName()),
			       0,
			       BSTSOffset,
			       //G4ThreeVector(0,0,0),
			       //detSolDownstreamVacInfo,
			       MBSMotherInfo,
			       0,
			       MBSisVisible,
			       G4Colour::Gray(),
			       MBSisSolid,
			       forceAuxEdgeVisible,
			       placePV,
			       doSurfaceCheck
			       );
    } else {
      // Now we'll build the version with pump-out holes
      // Get information for pump-out holes
      double const BSTSHoleXDim     = mbsgh.getHoleXDimInSteel();
      double const BSTSHoleYDim     = mbsgh.getHoleYDimInSteel();
      double const BSTSHoleZDim     = mbsgh.getHoleZDimInSteel();
      std::vector<CLHEP::Hep3Vector> BSTSHoleCenters 
	= mbsgh.getHoleCentersInSteel();

      // Make the steel
      VolumeInfo BSTSInfo( "BSTS", BSTSOffset, MBSMotherInfo.centerInWorld );

      G4Polycone * steel = new G4Polycone ( "BSTS",
		   pBSTSParams.getPolyconsParams().phi0(),
		   pBSTSParams.getPolyconsParams().phiTotal(),
		   pBSTSParams.getPolyconsParams().numZPlanes(),
		   &pBSTSParams.getPolyconsParams().zPlanes()[0],
		   &pBSTSParams.getPolyconsParams().rInner()[0],
		   &pBSTSParams.getPolyconsParams().rOuter()[0] );
		   
      
      G4SubtractionSolid* aSolid = 0;
      
      // Now loop over the holes and make them.
      for ( unsigned int jHole = 0; jHole < BSTSHoleCenters.size(); jHole++ ) {
	std::ostringstream hname;
	hname << "BSTSHole" << jHole+1;
	G4Box* aHoleBox = new G4Box( hname.str(),
				     BSTSHoleXDim/2.0*CLHEP::mm,
				     BSTSHoleYDim/2.0*CLHEP::mm,
				     BSTSHoleZDim/2.0*CLHEP::mm );

	if ( 0 == aSolid ) {
	  aSolid = new G4SubtractionSolid( "BSTS",
					   steel, aHoleBox,
					   0,
					   BSTSHoleCenters[jHole] );
	} else {
	  G4SubtractionSolid * bSolid 
	    = new G4SubtractionSolid ( "BSTS",
				       aSolid, aHoleBox,
				       0,
				       BSTSHoleCenters[jHole] );
	  aSolid = bSolid;
	}
      }

      BSTSInfo.solid = aSolid;
      finishNesting( BSTSInfo,
		     findMaterialOrThrow(pBSTSParams.materialName()),
		     0,
		     BSTSInfo.centerInParent,
		     MBSMotherInfo.logical,
		     0,
		     MBSisVisible,
		     G4Colour::Gray(),
		     MBSisSolid,
		     forceAuxEdgeVisible,
		     placePV,
		     doSurfaceCheck
		     );
      
    }


    if ( MBSversion == 1 ) {

      // SPBSSup1
      // This one is placed directly into DS3Vacuum;
      // Note that its radius is lager than that of BSTS

      CLHEP::Hep3Vector SPBSSup1OffsetInMu2e = pSPBSSup1Params.originInMu2e();

      if ( verbosityLevel > 0) {
	cout << __func__ << " SPBSSup1OffsetInMu2e                 : " << SPBSSup1OffsetInMu2e << endl;
      }

      // now local offset in mother volume
      CLHEP::Hep3Vector SPBSSup1Offset = SPBSSup1OffsetInMu2e - MBSMOffsetInMu2e;

      VolumeInfo SPBSSup1Info  = nestTubs("SPBSSup1",
					  pSPBSSup1Params.getTubsParams(),
					  findMaterialOrThrow(pSPBSSup1Params.materialName()),
					  0,
					  SPBSSup1Offset,
					  //detSolDownstreamVacInfo,
					  MBSMotherInfo,
					  0,
					  MBSisVisible,
					  G4Colour::Blue(),
					  MBSisSolid,
					  forceAuxEdgeVisible,
					  placePV,
					  doSurfaceCheck
					  );

      if ( verbosityLevel > 0) {
	double zhl         = static_cast<G4Tubs*>(SPBSSup1Info.solid)->GetZHalfLength();
	double SPBSSup1OffsetInMu2eZ = SPBSSup1OffsetInMu2e[CLHEP::Hep3Vector::Z];
	cout << __func__ << " SPBSSup1         Z extent in Mu2e    : " <<
	  SPBSSup1OffsetInMu2eZ - zhl << ", " << SPBSSup1OffsetInMu2eZ + zhl << endl;
      }

      // SPBSSup2
      // This one is placed directly into DS3Vacuum;
      // Note that its radius is lager than that of BSTS

      CLHEP::Hep3Vector SPBSSup2OffsetInMu2e = pSPBSSup2Params.originInMu2e();

      if ( verbosityLevel > 0) {
	cout << __func__ << " SPBSSup2OffsetInMu2e                 : " << SPBSSup2OffsetInMu2e << endl;
      }

      // now local offset in mother volume
      CLHEP::Hep3Vector SPBSSup2Offset = SPBSSup2OffsetInMu2e - MBSMOffsetInMu2e;

      VolumeInfo SPBSSup2Info  = nestTubs("SPBSSup2",
					  pSPBSSup2Params.getTubsParams(),
					  findMaterialOrThrow(pSPBSSup2Params.materialName()),
					  0,
					  SPBSSup2Offset,
					  //detSolDownstreamVacInfo,
					  MBSMotherInfo,
					  0,
					  MBSisVisible,
					  G4Colour::Blue(),
					  MBSisSolid,
					  forceAuxEdgeVisible,
					  placePV,
					  doSurfaceCheck
					  );

      if ( verbosityLevel > 0) {
	double zhl         = static_cast<G4Tubs*>(SPBSSup2Info.solid)->GetZHalfLength();
	double SPBSSup2OffsetInMu2eZ = SPBSSup2OffsetInMu2e[CLHEP::Hep3Vector::Z];
	cout << __func__ << " SPBSSup2         Z extent in Mu2e    : " <<
	  SPBSSup2OffsetInMu2eZ - zhl << ", " << SPBSSup2OffsetInMu2eZ + zhl << endl;
      }


      // SPBSL
      // This one is placed directly into DS3Vacuum;
      // Note that its radius is lager than that of BSTS

      CLHEP::Hep3Vector SPBSLOffsetInMu2e = pSPBSLParams.originInMu2e();

      if ( verbosityLevel > 0) {
	cout << __func__ << " SPBSLOffsetInMu2e                 : " << SPBSLOffsetInMu2e << endl;
      }

      // now local offset in mother volume
      CLHEP::Hep3Vector SPBSLOffset = SPBSLOffsetInMu2e - MBSMOffsetInMu2e;

      VolumeInfo SPBSLInfo  = nestTubs("SPBSL",
				       pSPBSLParams.getTubsParams(),
				       findMaterialOrThrow(pSPBSLParams.materialName()),
				       0,
				       SPBSLOffset,
				       //detSolDownstreamVacInfo,
				       MBSMotherInfo,
				       0,
				       MBSisVisible,
				       G4Colour::Blue(),
				       MBSisSolid,
				       forceAuxEdgeVisible,
				       placePV,
				       doSurfaceCheck
				       );

      if ( verbosityLevel > 0) {
	double zhl         = static_cast<G4Tubs*>(SPBSLInfo.solid)->GetZHalfLength();
	double SPBSLOffsetInMu2eZ = SPBSLOffsetInMu2e[CLHEP::Hep3Vector::Z];
	cout << __func__ << " SPBSL         Z extent in Mu2e    : " <<
	  SPBSLOffsetInMu2eZ - zhl << ", " << SPBSLOffsetInMu2eZ + zhl << endl;
      }


      // SPBSR
      // This one is placed directly into DS3Vacuum;
      // Note that its radius is larger than that of BSTS

      CLHEP::Hep3Vector SPBSROffsetInMu2e = pSPBSRParams.originInMu2e();

      if ( verbosityLevel > 0) {
	cout << __func__ << " SPBSROffsetInMu2e                 : " << SPBSROffsetInMu2e << endl;
      }

      // now local offset in mother volume
      CLHEP::Hep3Vector SPBSROffset = SPBSROffsetInMu2e - MBSMOffsetInMu2e;

      VolumeInfo SPBSRInfo  = nestTubs("SPBSR",
				       pSPBSRParams.getTubsParams(),
				       findMaterialOrThrow(pSPBSRParams.materialName()),
				       0,
				       SPBSROffset,
				       //detSolDownstreamVacInfo,
				       MBSMotherInfo,
				       0,
				       MBSisVisible,
				       G4Colour::Blue(),
				       MBSisSolid,
				       forceAuxEdgeVisible,
				       placePV,
                                    doSurfaceCheck
				       );

      if ( verbosityLevel > 0) {
	double zhl         = static_cast<G4Tubs*>(SPBSRInfo.solid)->GetZHalfLength();
	double SPBSROffsetInMu2eZ = SPBSROffsetInMu2e[CLHEP::Hep3Vector::Z];
	cout << __func__ << " SPBSR         Z extent in Mu2e    : " <<
	  SPBSROffsetInMu2eZ - zhl << ", " << SPBSROffsetInMu2eZ + zhl << endl;
      }


    } // end of version 1 specific items


    // SPBSC
    // This one is placed directly into DS3Vacuum;
    // Note that its radius is larger than the one of BSTS

    CLHEP::Hep3Vector SPBSCOffsetInMu2e = pSPBSCParams.originInMu2e();

    if ( verbosityLevel > 0) {
      cout << __func__ << " SPBSCOffsetInMu2e                 : " << SPBSCOffsetInMu2e << endl;
    }

    // now local offset in mother volume
    CLHEP::Hep3Vector SPBSCOffset = SPBSCOffsetInMu2e - MBSMOffsetInMu2e;

    VolumeInfo SPBSCInfo  = nestTubs("SPBSC",
                                    pSPBSCParams.getTubsParams(),
                                    findMaterialOrThrow(pSPBSCParams.materialName()),
                                    0,
                                    SPBSCOffset,
                                    //detSolDownstreamVacInfo,
                                    MBSMotherInfo,
                                    0,
                                    MBSisVisible,
                                    G4Colour::Blue(),
                                    MBSisSolid,
                                    forceAuxEdgeVisible,
                                    placePV,
                                    doSurfaceCheck
                                    );

    if ( verbosityLevel > 0) {
      double zhl         = static_cast<G4Tubs*>(SPBSCInfo.solid)->GetZHalfLength();
      double SPBSCOffsetInMu2eZ = SPBSCOffsetInMu2e[CLHEP::Hep3Vector::Z];
      cout << __func__ << " SPBSC         Z extent in Mu2e    : " <<
        SPBSCOffsetInMu2eZ - zhl << ", " << SPBSCOffsetInMu2eZ + zhl << endl;
    }

    // =========================
    // BSTC - This is the upstream inner HDPE liner
    // =========================

    CLHEP::Hep3Vector BSTCOffsetInMu2e = pBSTCParams.originInMu2e();
    // now local offset in mother volume
    CLHEP::Hep3Vector BSTCOffset =  BSTCOffsetInMu2e - MBSMOffsetInMu2e; // - MBSMotherInfo.centerInMu2e();
    //BSTCOffset.setZ(0);

    if ( verbosityLevel > 0) {
      cout << __func__ << " BSTCOffsetInMu2e                 : " << BSTCOffsetInMu2e << endl;
      cout << __func__ << " BSTCOffsetInMBS                  : " << BSTCOffset << endl;
    }

    G4Colour  orange  (.75, .55, .0);
    if ( MBSversion < 4 ) {
      VolumeInfo BSTCInfo  = nestPolycone("BSTC",
					  pBSTCParams.getPolyconsParams(),
					  findMaterialOrThrow(pBSTCParams.materialName()),
					  0,
					  BSTCOffset,
					  //detSolDownstreamVacInfo,
					  MBSMotherInfo,
					  0,
					  MBSisVisible,
					  orange,
					  MBSisSolid,
					  forceAuxEdgeVisible,
					  placePV,
					  doSurfaceCheck
					  );
      
      
      if ( verbosityLevel > 0) {
	G4Polycone *tmpBSTCInfo  =  static_cast<G4Polycone*>(BSTCInfo.solid);
	double zhl = tmpBSTCInfo->GetCorner(tmpBSTCInfo->GetNumRZCorner()-1).z-tmpBSTCInfo->GetCorner(0).z;
	zhl*=0.5;
	double BSTCOffsetInMu2eZ = BSTCOffsetInMu2e[CLHEP::Hep3Vector::Z];
	cout << __func__ << " BSTC         Z extent in Mu2e    : " <<
	  BSTCOffsetInMu2eZ - zhl << ", " << BSTCOffsetInMu2eZ + zhl << endl;
      }

    } else {
      // Now we'll build the version with pump-out holes
      // Get information for pump-out holes
      double const BSTCHoleXDim   = mbsgh.getHoleXDimInUpPoly();
      double const BSTCHoleYDim   = mbsgh.getHoleYDimInUpPoly();
      double const BSTCHoleZDim   = mbsgh.getHoleZDimInUpPoly();
      std::vector<CLHEP::Hep3Vector> BSTCHoleCenters 
      	= mbsgh.getHoleCentersInUpstreamPoly();

      // Make the steel
      VolumeInfo BSTCInfo( "BSTC", BSTCOffset, MBSMotherInfo.centerInWorld );

      G4Polycone * hdpe = new G4Polycone ( "BSTC",
		   pBSTCParams.getPolyconsParams().phi0(),
		   pBSTCParams.getPolyconsParams().phiTotal(),
		   pBSTCParams.getPolyconsParams().numZPlanes(),
		   &pBSTCParams.getPolyconsParams().zPlanes()[0],
		   &pBSTCParams.getPolyconsParams().rInner()[0],
		   &pBSTCParams.getPolyconsParams().rOuter()[0] );
		   
      
      G4SubtractionSolid* aSolid = 0;
      
      // Now loop over the holes and make them.
      for ( unsigned int jHole = 0; jHole < BSTCHoleCenters.size(); jHole++ ) {
	std::ostringstream hname;
	hname << "BSTCHole" << jHole+1;
	G4Box* aHoleBox = new G4Box( hname.str(),
				     BSTCHoleXDim/2.0*CLHEP::mm,
				     BSTCHoleYDim/2.0*CLHEP::mm,
				     BSTCHoleZDim/2.0*CLHEP::mm );

	if ( 0 == aSolid ) {
	  CLHEP::Hep3Vector realCenter = BSTCHoleCenters[jHole] - BSTCOffset;
	  aSolid = new G4SubtractionSolid( "BSTC",
	       			   hdpe, aHoleBox,
       				   0,
	       			   realCenter );
	} else {
	  CLHEP::Hep3Vector realCenter = BSTCHoleCenters[jHole] - BSTCOffset;
	  G4SubtractionSolid * bSolid 
	    = new G4SubtractionSolid ( "BSTC",
				       aSolid, aHoleBox,
				       0,
				       realCenter );
	  aSolid = bSolid;
	}
      }

      BSTCInfo.solid = aSolid;
      finishNesting( BSTCInfo,
		     findMaterialOrThrow(pBSTCParams.materialName()),
		     0,
		     BSTCInfo.centerInParent,
		     MBSMotherInfo.logical,
		     0,
		     MBSisVisible,
		     G4Colour::Gray(),
		     MBSisSolid,
		     forceAuxEdgeVisible,
		     placePV,
		     doSurfaceCheck
		     );
	       
    }


    // =========================
    // BSBS - This is the downstream inner HDPE liner
    // =========================

    CLHEP::Hep3Vector BSBSOffsetInMu2e = pBSBSParams.originInMu2e();

    // now local offset in mother volume
    CLHEP::Hep3Vector BSBSOffset =  BSBSOffsetInMu2e - MBSMOffsetInMu2e; //-MBSMotherInfo.centerInMu2e();


    if ( verbosityLevel > 0) {
      cout << __func__ << " BSBSOffsetInMu2e                 : " << BSBSOffsetInMu2e << endl;
      cout << __func__ << " BSBSOffsetInMBS                  : " << BSBSOffset << endl;
    }
    if ( MBSversion < 4 ) {
      VolumeInfo BSBSInfo  = nestPolycone("BSBS",
					  pBSBSParams.getPolyconsParams(),
					  findMaterialOrThrow(pBSBSParams.materialName()),
					  0,
					  BSBSOffset,
					  //detSolDownstreamVacInfo,
					  MBSMotherInfo,
					  0,
					  MBSisVisible,
					  G4Colour::Yellow(),
					  MBSisSolid,
					  forceAuxEdgeVisible,
					  placePV,
					  doSurfaceCheck
					  );
      
      if ( verbosityLevel > 0) {
	G4Polycone *tmpBSBSInfo  =  static_cast<G4Polycone*>(BSBSInfo.solid);
	double zhl = tmpBSBSInfo->GetCorner(tmpBSBSInfo->GetNumRZCorner()-1).z-tmpBSBSInfo->GetCorner(0).z;
	zhl*=0.5;
	double BSBSOffsetInMu2eZ = BSBSOffsetInMu2e[CLHEP::Hep3Vector::Z];
	cout << __func__ << " BSBS         Z extent in Mu2e    : " <<
	  BSBSOffsetInMu2eZ - zhl << ", " << BSBSOffsetInMu2eZ + zhl << endl;
      }
    } else { 
      // Now we'll build the version with pump-out holes
      // Get information for pump-out holes
      double const BSBSHoleXDim = mbsgh.getHoleXDimInDownPoly();
      double const BSBSHoleYDim = mbsgh.getHoleYDimInDownPoly();
      double const BSBSHoleZDim = mbsgh.getHoleZDimInDownPoly();
       std::vector<CLHEP::Hep3Vector> BSBSHoleCenters 
	 = mbsgh.getHoleCentersInDownstreamPoly();

      // Make the upstream HDPE
      VolumeInfo BSBSInfo( "BSBS", BSBSOffset, MBSMotherInfo.centerInWorld );

      G4Polycone * hdpe = new G4Polycone ( "BSBS",
		   pBSBSParams.getPolyconsParams().phi0(),
		   pBSBSParams.getPolyconsParams().phiTotal(),
		   pBSBSParams.getPolyconsParams().numZPlanes(),
		   &pBSBSParams.getPolyconsParams().zPlanes()[0],
		   &pBSBSParams.getPolyconsParams().rInner()[0],
		   &pBSBSParams.getPolyconsParams().rOuter()[0] );
		   
      
      G4SubtractionSolid* aSolid = 0;
      
      // Now loop over the holes and make them.
      for ( unsigned int jHole = 0; jHole < BSBSHoleCenters.size(); jHole++ ) {
	CLHEP::Hep3Vector realCenter = BSBSHoleCenters[jHole] - BSBSOffset;
	std::ostringstream hname;
	hname << "BSBSHole" << jHole+1;
	G4Box* aHoleBox = new G4Box( hname.str(),
				     BSBSHoleXDim/2.0*CLHEP::mm,
				     BSBSHoleYDim/2.0*CLHEP::mm,
				     BSBSHoleZDim/2.0*CLHEP::mm );

	if ( 0 == aSolid ) {
	  aSolid = new G4SubtractionSolid( "BSBS",
					   hdpe, aHoleBox,
					   0,
					   realCenter );
	} else {
	  G4SubtractionSolid * bSolid 
	    = new G4SubtractionSolid ( "BSBS",
				       aSolid, aHoleBox,
				       0,
				       realCenter );
	  aSolid = bSolid;
	}
      }

      BSBSInfo.solid = aSolid;
      finishNesting( BSBSInfo,
		     findMaterialOrThrow(pBSBSParams.materialName()),
		     0,
		     BSBSInfo.centerInParent,
		     MBSMotherInfo.logical,
		     0,
		     MBSisVisible,
		     G4Colour::Gray(),
		     MBSisSolid,
		     forceAuxEdgeVisible,
		     placePV,
		     doSurfaceCheck
		     );
	       
    }


    // ==============================
    // CLV2 - This is the end plug
    // ==============================

    CLHEP::Hep3Vector CLV2OffsetInMu2e = pCLV2Params.originInMu2e();

    // now local offset in mother volume
    CLHEP::Hep3Vector CLV2Offset = CLV2OffsetInMu2e - MBSMOffsetInMu2e;//- MBSMotherInfo.centerInMu2e();
    //CLV2Offset.setZ(0.0);

    if ( verbosityLevel > 0) {
      cout << __func__ << " CLV2OffsetInMu2e                 : " << CLV2OffsetInMu2e << endl;
      cout << __func__ << " CLV2OffsetInMBS                  : " << CLV2Offset << endl;
    }

    VolumeInfo CLV2Info  = nestPolycone("CLV2",
                                    pCLV2Params.getPolyconsParams(),
                                    findMaterialOrThrow(pCLV2Params.materialName()),
                                    0,
                                    CLV2Offset,
                                    //detSolDownstreamVacInfo,
                                    MBSMotherInfo,
                                    0,
                                    MBSisVisible,
                                    orange,
                                    MBSisSolid,
                                    forceAuxEdgeVisible,
                                    placePV,
                                    doSurfaceCheck
                                    );

    if ( verbosityLevel > 0) {
      G4Polycone *tmpCLV2Info  =  static_cast<G4Polycone*>(CLV2Info.solid);
      double zhl = tmpCLV2Info->GetCorner(tmpCLV2Info->GetNumRZCorner()-1).z-tmpCLV2Info->GetCorner(0).z;
      zhl*=0.5;
      double CLV2OffsetInMu2eZ = CLV2OffsetInMu2e[CLHEP::Hep3Vector::Z];
      cout << __func__ << " CLV2 number of RZ corners        : " << tmpCLV2Info->GetNumRZCorner() << endl;
      cout << __func__ << " CLV2         Z extent in Mu2e    : " <<
        CLV2OffsetInMu2eZ - zhl << ", " << CLV2OffsetInMu2eZ + zhl << endl;
    }
    
    // CLV2 Absorber :  Variable thickness plug in MBS axial hole
    if (_config.getBool("mbs.CLV2.absorber.build",false)){
      
      CLHEP::Hep3Vector CLV2AbsOffsetInMu2e = pCLV2ABSParams.originInMu2e();
      // now local offset in mother volume
      CLHEP::Hep3Vector CLV2AbsOffset = CLV2AbsOffsetInMu2e - MBSMOffsetInMu2e;
        
      VolumeInfo CLV2AbsorberInfo  = nestTubs("CLV2Absorber",
					      pCLV2ABSParams.getTubsParams(),
					      findMaterialOrThrow(pCLV2ABSParams.materialName()),
					      0,
					      CLV2AbsOffset,
					      //detSolDownstreamVacInfo,
					      MBSMotherInfo,
					      0,
					      MBSisVisible,
					      orange,
					      MBSisSolid,
					      forceAuxEdgeVisible,
					      placePV,
					      doSurfaceCheck
					      );
        
      if ( verbosityLevel > 0) {
	cout << __func__ << " CLV2AbsOffsetInMu2e                 : " << CLV2AbsOffsetInMu2e << endl;
	cout << __func__ << " CLV2AbsOffsetInMBS                  : " << CLV2AbsOffset << endl;       
	double zhl         = static_cast<G4Tubs*>(CLV2AbsorberInfo.solid)->GetZHalfLength();
	double CLV2AbsOffsetInMu2eZ = CLV2AbsOffsetInMu2e[CLHEP::Hep3Vector::Z];
	cout << __func__ << " CLV2Absorber         Z extent in Mu2e    : " <<
          CLV2AbsOffsetInMu2eZ - zhl << ", " << CLV2AbsOffsetInMu2eZ + zhl << endl;
      }        
    }

    //Adding a shield at the front of the MBS to protect the calorimeter
    if(MBSversion == 6) {
      CLHEP::Hep3Vector RingOffset = pCalShieldRingParams.originInMu2e() - detSolDownstreamVacInfo.centerInMu2e(); // MBSMOffsetInMu2e;
       
      VolumeInfo MBSCalShieldRing  = nestTubs("CalShieldRing",
					      pCalShieldRingParams.getTubsParams(),
					      findMaterialOrThrow(pCalShieldRingParams.materialName()),
					      0,
					      RingOffset,
					      detSolDownstreamVacInfo, // MBSMotherInfo,
					      0,
					      MBSisVisible,
					      G4Colour::Blue(),
					      MBSisSolid,
					      forceAuxEdgeVisible,
					      placePV,
					      doSurfaceCheck
					      );
    }

  } // end of constructMBS;

}
