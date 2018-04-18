//
// Free function to create DS. (Detector Solenoid)
//
// $Id: constructDS.cc,v 1.22 2014/09/19 19:14:59 knoepfel Exp $
// $Author: knoepfel $
// $Date: 2014/09/19 19:14:59 $
//
// Original author KLG based on Mu2eWorld constructDS
//

// Mu2e includes.
#include "BeamlineGeom/inc/Beamline.hh"
#include "BeamlineGeom/inc/StraightSection.hh"
#include "DetectorSolenoidGeom/inc/DetectorSolenoid.hh"
#include "DetectorSolenoidGeom/inc/DetectorSolenoidShielding.hh"
#include "G4Helper/inc/G4Helper.hh"
#include "G4Helper/inc/VolumeInfo.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeomPrimitives/inc/PolyconsParams.hh"
#include "MBSGeom/inc/MBS.hh"
#include "Mu2eG4/inc/findMaterialOrThrow.hh"
#include "Mu2eG4/inc/constructDS.hh"
#include "Mu2eG4/inc/nestBox.hh"
#include "Mu2eG4/inc/nestTubs.hh"
#include "Mu2eG4/inc/nestPolycone.hh"
#include "Mu2eG4/inc/nestExtrudedSolid.hh"
#include "Mu2eG4/inc/finishNesting.hh"
#include "Mu2eG4/inc/MaterialFinder.hh"
#include "Mu2eG4/inc/SensitiveDetectorName.hh"
#include "Mu2eG4/inc/SensitiveDetectorHelper.hh"

// G4 includes
#include "G4ThreeVector.hh"
#include "G4Material.hh"
#include "G4Color.hh"
#include "G4ExtrudedSolid.hh"
#include "G4Polycone.hh"
#include "G4Tubs.hh"
#include "G4SubtractionSolid.hh"
#include "G4VSensitiveDetector.hh"
#include "G4SDManager.hh"
// CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/TwoVector.h"

using namespace std;

namespace mu2e {

  void constructDS( const VolumeInfo& parent,
                    const SimpleConfig& _config,
                    const SensitiveDetectorHelper& sdHelper
                    ){
    MaterialFinder materialFinder(_config);

    // Load flags
    int const verbosityLevel = _config.getInt("ds.verbosityLevel",0);
    bool const inGaragePosition = _config.getBool("inGaragePosition",false);

    G4GeometryOptions* geomOptions = art::ServiceHandle<GeometryService>()->geomOptions();
    geomOptions->loadEntry( _config, "DS"         , "ds"          );
    geomOptions->loadEntry( _config, "DSCoil"     , "dsCoil"      );
    geomOptions->loadEntry( _config, "DSRing"     , "dsRing"      );
    geomOptions->loadEntry( _config, "DSSpacer"   , "dsSpacer"    );
    geomOptions->loadEntry( _config, "DSSupport"  , "dsSupport"   );
    geomOptions->loadEntry( _config, "DSThShield" , "dsThShield"  );
    geomOptions->loadEntry( _config, "DSVacuum"   , "dsVacuum"    );
    geomOptions->loadEntry( _config, "DSShielding", "dsShielding" );

    // Fetch parent (hall) position
    G4ThreeVector _hallOriginInMu2e = parent.centerInMu2e();

    // Fetch DS geom. object
   GeomHandle<DetectorSolenoid> ds;
    CLHEP::Hep3Vector const & dsP ( ds->position() );

    // DS cryostat
    G4Material*   dsCryoMaterial = findMaterialOrThrow( ds->material() );

    // - inner cryo shell
    G4ThreeVector dsInnerCryoPosition( dsP.x(), dsP.y(), dsP.z());
    TubsParams    dsInnerCryoParams  ( ds->rIn1(), ds->rIn2(), ds->halfLength() );
    nestTubs( "DSInnerCryoShell",
              dsInnerCryoParams,
              dsCryoMaterial,
              0,
              dsInnerCryoPosition-_hallOriginInMu2e,
              parent,
              0,
              G4Color::Magenta(),
	      "DS"
              );
    // ***
    // Lining for inner cryo shell - for test of shielding
    // added 11 June 2017
    // ***
    if ( ds->hasInnerLining() ) { // only if specifically enabled
      nestTubs( "DSInnerCryoLining",
		TubsParams( ds->rIn1()-ds->innerLiningThickness(),ds->rIn1(),ds->halfLength()-2.*ds->endWallHalfLength()),
		findMaterialOrThrow(ds->innerLiningMaterial()),
		0,
		dsInnerCryoPosition-_hallOriginInMu2e,
		parent,
		0,
		G4Color::Magenta(),
		"DS"
		);
    }

    // - outer cryo shell
    G4ThreeVector dsOuterCryoPosition( dsP.x(), dsP.y(), dsP.z());
    TubsParams    dsOuterCryoParams  ( ds->rOut1(), ds->rOut2(), ds->halfLength() );
    nestTubs( "DSOuterCryoShell",
              dsOuterCryoParams,
              dsCryoMaterial,
              0,
              dsOuterCryoPosition-_hallOriginInMu2e,
              parent,
              0,
              G4Color::Magenta(),
	      "DS"
              );

    // DNB (Lou) May 2017.  Making a vacuum volume inside the cryostat
    G4Tubs * dsCryoVacTub = new G4Tubs( "DSCryoVacuumTube",
					ds->rIn2(), ds->rOut1(),
					ds->halfLength() - 2.*ds->endWallHalfLength(),
					0.0, 360.0*CLHEP::degree);

    CLHEP::Hep3Vector dsCryoVacLocationInMu2e( dsP.x(), dsP.y(), dsP.z() );
    VolumeInfo dsCryoVacMother( "DSCryoVacuumRegion",
				dsCryoVacLocationInMu2e - _hallOriginInMu2e,
				parent.centerInWorld );
    dsCryoVacMother.solid = dsCryoVacTub;
				
    finishNesting ( dsCryoVacMother,
		    findMaterialOrThrow("DSVacuum"),
		    0, dsCryoVacMother.centerInParent,
		    parent.logical, 0, G4Colour::White(),
		    "DS" );
		    

    // - end walls
    TubsParams    dsEndWallParams    ( ds->rIn2(), ds->rOut1(), ds->endWallHalfLength() );
    G4ThreeVector dsUpEndWallPosition( dsP.x(), dsP.y(), 
                                       dsP.z() - ds->halfLength() + ds->endWallHalfLength());
    nestTubs( "DSUpEndWallShell",
              dsEndWallParams,
              dsCryoMaterial,
              0,
              dsUpEndWallPosition-_hallOriginInMu2e,
              parent,
              0,
              G4Color::Magenta(),
	      "DS"
              );

    G4ThreeVector dsDownEndWallPosition( dsP.x(), dsP.y(), 
                                         dsP.z() + ds->halfLength() - ds->endWallHalfLength());
    nestTubs( "DSDownEndWallShell",
              dsEndWallParams,
              dsCryoMaterial,
              0,
              dsDownEndWallPosition-_hallOriginInMu2e,
              parent,
              0,
              G4Color::Magenta(),
	      "DS"
              );

    // - upstream face
    GeomHandle<Beamline> beamg;
    const StraightSection * ts5out = beamg->getTS().getTSCryo<StraightSection>( TransportSolenoid::TSRegion::TS5,
                                                                                TransportSolenoid::TSRadialPart::OUT );

    double dsFrontZ0 = dsP.z() - ds->halfLength() + ds->frontHalfLength();
    TubsParams    dsFrontParams  ( ts5out->rOut(), dsInnerCryoParams.innerRadius(), ds->frontHalfLength() );
    G4ThreeVector dsFrontPosition( dsP.x(), dsP.y(), dsFrontZ0);

    nestTubs( "DSFront",
              dsFrontParams,
              dsCryoMaterial,
              0,
              dsFrontPosition-_hallOriginInMu2e,
              parent,
              0,
              G4Color::Blue(),
	      "DS"
              );

    // DS thermal shield
    G4Material*   dsThShieldMaterial = findMaterialOrThrow( ds->shield_material() );

    // - inner shield shell
    G4ThreeVector dsInnerThShieldPosition( dsP.x(), dsP.y(), dsP.z());
    TubsParams    dsInnerThShieldParams  ( ds->shield_rIn1(), ds->shield_rIn2(), ds->shield_halfLength() );
    nestTubs( "DSInnerThShieldShell",
              dsInnerThShieldParams,
              dsThShieldMaterial,
              0,
              dsInnerThShieldPosition-dsCryoVacLocationInMu2e,
              dsCryoVacMother,
              0,
              G4Color::Cyan(),
	      "DSThShield"
              );

    // - outer shield shell
    G4ThreeVector dsOuterThShieldPosition( dsP.x(), dsP.y(), dsP.z());
    TubsParams    dsOuterThShieldParams  ( ds->shield_rOut1(), ds->shield_rOut2(), ds->shield_halfLength() );
    nestTubs( "DSOuterThShieldShell",
              dsOuterThShieldParams,
              dsThShieldMaterial,
              0,
              dsOuterThShieldPosition-dsCryoVacLocationInMu2e,
              dsCryoVacMother,
              0,
              G4Color::Cyan(),
	      "DSThShield"
              );

    // SKIPPING END WALLS OF THERMAL SHIELD
    // DUE TO INCONSISTENCY IN DRAWINGS!

    // DS coils
    G4Material* dsCoilMaterial;
    for ( int i(0); i < ds->nCoils() ; i++ ) {
      TubsParams coilParams( ds->coil_rIn(), 
                             ds->coil_rOut().at(i),
                             ds->coil_zLength().at(i)*0.5 );
      if ( ds->coilVersion() == 1 ) {
	dsCoilMaterial = findMaterialOrThrow( ds->coil_material() );
      } else {
	dsCoilMaterial = findMaterialOrThrow( ds->coil_materials().at(i) );
      }

      G4ThreeVector coilPosition( dsP.x(), dsP.y(), 
                                  ds->coil_zPosition().at(i) + coilParams.zHalfLength()); 

      ostringstream coilname;
      coilname << "DSCoil_" << i+1;

      nestTubs( coilname.str(),
                coilParams,
                dsCoilMaterial,
                0,
                coilPosition-dsCryoVacLocationInMu2e,
                dsCryoVacMother,
                0,
                G4Color::Green(),
		"DSCoil"
                );
    }

    // DS coil spacers
    if ( ds->coilVersion() > 1 ) {  // spacers added to version 2 and up
      G4Material* dsSpacerMaterial = findMaterialOrThrow( ds->spacer_material() );
      for ( int i(0); i < ds->nSpacers() ; i++ ) {
	TubsParams spacerParams( ds->spacer_rIn(), 
				 ds->spacer_rOut().at(i),
				 ds->spacer_zLength().at(i)*0.5 );

	G4ThreeVector spacerPosition( dsP.x(), dsP.y(), 
				      ds->spacer_zPosition().at(i) 
				      + spacerParams.zHalfLength()); 

	ostringstream spacername;
	spacername << "DSSpacer_" << i+1;
	
	nestTubs( spacername.str(),
		  spacerParams,
		  dsSpacerMaterial,
		  0,
		  spacerPosition-dsCryoVacLocationInMu2e,
		  dsCryoVacMother,
		  0,
		  G4Color::Green(),
		  "DSSpacer"
		  );
      }
    }


    // DS coils support system
    G4Material*   dsSupportMaterial = findMaterialOrThrow( ds->support_material() );
    G4ThreeVector dsSupportPosition( dsP.x(), dsP.y(), dsP.z());
    TubsParams    dsSupportParams  ( ds->support_rIn(), ds->support_rOut(), ds->support_halfLength() );
    nestTubs( "DSCoilSupport",
              dsSupportParams,
              dsSupportMaterial,
              0,
              dsSupportPosition-dsCryoVacLocationInMu2e,
              dsCryoVacMother,
              0,
              G4Color::Blue(),
	      "DSSupport"
              );



    // Build DS Rings.  Added by David Norvil Brown, May 2015
    double rirs = ds->rInRingSide();
    double rors = ds->rOutRingSide();
    double trs = ds->thickRingSide();
    double rir = ds->rInRing();
    double ror = ds->rOutRing();
    double lr = ds->lengthRing();
    G4Material* ringMaterial = findMaterialOrThrow(ds->RingMaterial());
    std::vector<double> xr = ds->xRing();
    std::vector<double> yr = ds->yRing();
    std::vector<double> zr = ds->zRing();
    
    for ( unsigned int iRing = 0; iRing < xr.size(); iRing++ ) {

      // Let's build a mother volume first
      std::ostringstream ringMotherName;
      ringMotherName << "DSRingMother" << iRing;
      CLHEP::HepRotation* ringRotat = new CLHEP::HepRotation(CLHEP::HepRotation::IDENTITY);

      double motherx = xr[iRing];
      double mothery = yr[iRing];
      double motherz = zr[iRing];

      VolumeInfo motherVol = nestTubs( ringMotherName.str(),
				       TubsParams( rirs, rors, trs + lr/2.0 ),
				       findMaterialOrThrow("G4_AIR"),
				       ringRotat, 
				       CLHEP::Hep3Vector(motherx,mothery,motherz) - _hallOriginInMu2e,
				       parent, 0, G4Color::Blue(),
				       "DSRing" );
 
     std::ostringstream leftName;
      leftName << "DSleftSideRing" << iRing;


      nestTubs( leftName.str(),
		TubsParams( rirs, rors, trs/2.0 ),
		ringMaterial,
                ringRotat,
		CLHEP::Hep3Vector(0.0,0.0,-lr/2.0-trs/2.0),
		motherVol,
		0,
		G4Color::Blue(),
		"DSRing"
		);

      std::ostringstream centerName;
      centerName << "DScenterRing" << iRing;

      nestTubs( centerName.str(),
		TubsParams( rir, ror, lr/2.0 ),
		ringMaterial,
                ringRotat,
		CLHEP::Hep3Vector(0.0,0.0,0.0),
		motherVol,
		0,
		G4Color::Blue(),
		"DSRing"
		);

      std::ostringstream rightName;
      rightName << "DSrightSideRing" << iRing;

      nestTubs( rightName.str(),
		TubsParams( rirs, rors, trs/2.0 ),
		ringMaterial,
                ringRotat,
		CLHEP::Hep3Vector(0.0,0.0,lr/2.0+trs/2.0),
		motherVol,
		0,
		G4Color::Blue(),
		"DSRing"
		);

    } // finished inserting Rings


    // DS vacuum volumes
    G4Material* vacuumMaterial = findMaterialOrThrow( ds->vacuumMaterial() );
    TubsParams ds1VacParams    ( ts5out->rOut(), ds->rIn1(), ds->vac_halfLengthDs1()   );
    TubsParams ds2VacParams    ( 0.            , ds->rIn1(), ds->vac_halfLengthDs2()   );

    // Compute/set positions of vacuum volumes in Mu2e coordinates.
    // - DS position is fixed by TS torus radius, and half lengths of 
    //   front face, DS1, and TS5
    double ds1Z0     = dsFrontZ0 + ds->frontHalfLength() + ds->vac_halfLengthDs1();
    double ds2Z0     = ds->vac_zLocDs23Split() - ds->vac_halfLengthDs2();
    //    double ds2HalfLength     = _config.getDouble("ds2.halfLength");
    
    if ( verbosityLevel > 0 ) {
      cout << __func__ << " DS2 vacuum extent: " 
           << " [ " << ds2Z0 - ds->vac_halfLengthDs2() << " , " 
           << ds->vac_zLocDs23Split() << " ] " << endl;
    }

    G4ThreeVector ds1Position( dsP.x(), dsP.y(), ds1Z0 );
    G4ThreeVector ds2Position( dsP.x(), dsP.y(), ds2Z0 );

    nestTubs( "DS1Vacuum",
              ds1VacParams,
              vacuumMaterial,
              0,
              ds1Position-_hallOriginInMu2e,
              parent,
              0,
              G4Colour::Green(),
	      "DSVacuum"
              );

    VolumeInfo ds2VacInfo = 
      nestTubs( "DS2Vacuum",
		ds2VacParams,
		vacuumMaterial,
		0,
		ds2Position-_hallOriginInMu2e,
		parent,
		0,
		G4Colour::Yellow(),
		"DSVacuum"
		);

    // Polycone geometry allows for MBS to extend beyond solenoid
    // physical boundaries
    
    GeomHandle<DetectorSolenoidShielding> dss;

    // Define polycone parameters
    vector<double> tmp_zPlanesDs3 = { 
      ds->vac_zLocDs23Split(), 
      dss->getVPSPmain()->zBegin(),
      dss->getVPSPCryoSeal()->zBegin(),
      dss->getVPSPCryoSeal()->zEnd(),
      dss->getVPSPCryoSeal()->zEnd(),
      dss->getVPSPendSeal()->zBegin(), 
      dss->getVPSPendSeal()->zBegin(),
      dss->getVPSPendFlange()->zEnd(),
      dss->getVPSPendFlange()->zEnd(),
      dss->getIFBmain()->zEnd(),
      dss->getIFBendPlug()->zBegin(),
      dss->getIFBendPlug()->zEnd() };

    vector<double> tmp_rOuterDs3  = {
      ds->rIn1(),
      dss->getVPSPmain()->innerRadius(),
      dss->getVPSPCryoSeal()->outerRadius(),
      dss->getVPSPCryoSeal()->outerRadius(),
      dss->getVPSPmain()->outerRadius(),
      dss->getVPSPmain()->outerRadius(),
      dss->getVPSPendSeal()->outerRadius(),
      dss->getVPSPendFlange()->outerRadius(),
      dss->getVPSPendFlange()->innerRadius(),
      dss->getIFBmain()->outerRadius(),
      dss->getIFBendPlug()->outerRadius(),
      dss->getIFBendPlug()->outerRadius()
    };

    assert( tmp_zPlanesDs3.size() == tmp_rOuterDs3.size() );

    vector<double> tmp_rInnerDs3 ( tmp_rOuterDs3.size(), 0. );

    PolyconsParams    ds3PolyParams    ( tmp_zPlanesDs3, tmp_rInnerDs3, tmp_rOuterDs3 );
    CLHEP::Hep3Vector ds3positionInMu2e( dsP.x(), dsP.y(), 0.);
    
    VolumeInfo dsShieldParent = nestPolycone( "DS3Vacuum",
					      ds3PolyParams,
					      vacuumMaterial,
					      0,
					      ds3positionInMu2e - parent.centerInMu2e(),
					      parent,
					      0,
					      G4Colour::Yellow(),
					      "DSVacuum"
					      );

    if ( inGaragePosition ) {
      double zOffGarage = _config.getDouble("garage.zOffset",14000.0);
      CLHEP::Hep3Vector relPosFake(0.,0., zOffGarage);
      G4Material*  airMaterial = findMaterialOrThrow( _config.getString("hall.insideMaterialName","G4_AIR") );

      VolumeInfo dsShieldParent = nestPolycone( "garageFakeDS3Vacuum",
						ds3PolyParams,
						airMaterial,
						0,
						ds3positionInMu2e - parent.centerInMu2e() + relPosFake,
						parent,
						0,
						G4Colour::Yellow(),
						"DSVacuum"
						);
    }
      

    // Construct shielding downstream of DS
    for ( const auto & shield : dss->getTubes() ) {

      nestTubs( shield->getName(),
                shield->getTubsParams(),
                findMaterialOrThrow(shield->materialName()),
                0,
                shield->originInMu2e()-dsShieldParent.centerInMu2e(),
                dsShieldParent,
                0,
                G4Colour::Blue(),
		"DSShielding"
                );

    }

    // ************ End Shielding, begin Rails ************
    // D. No. Brown
    // Add rails in DS2Vacuum

    //    std::vector<CLHEP::Hep2Vector> railOutline;
    std::vector<double> uRailOutline = ds->uOutlineRail();
    std::vector<double> vRailOutline = ds->vOutlineRail();

    const bool forceAuxEdgeVisible = _config.getBool("g4.forceAuxEdgeVisible",false);
    const bool doSurfaceCheck      = _config.getBool("g4.doSurfaceCheck",false)|| _config.getBool("ds.doSurfaceCheck",false);
    const bool placePV             = true;
    CLHEP::HepRotation* nRailRotat = new CLHEP::HepRotation(CLHEP::HepRotation::IDENTITY);
    CLHEP::HepRotation* sRailRotat = new CLHEP::HepRotation(CLHEP::HepRotation::IDENTITY);
    sRailRotat->rotateY(180.0*CLHEP::degree);

    VolumeInfo RailN2 = nestExtrudedSolid
                     ( "NorthRailDS2", ds->lengthRail2()/2.0*CLHEP::mm,
		       uRailOutline, vRailOutline, 
		       findMaterialOrThrow(ds->RailMaterial()),
		       nRailRotat, ds->n2RailCenter(),
		       ds2VacInfo.logical, 0, _config.getBool("ds.visible"),
		       G4Colour::Blue(), _config.getBool("ds.solid"),
		       forceAuxEdgeVisible, placePV, doSurfaceCheck );

    VolumeInfo RailS2 = nestExtrudedSolid
                     ( "SouthRailDS2", ds->lengthRail2()/2.0*CLHEP::mm,
		       uRailOutline, vRailOutline, 
		       findMaterialOrThrow(ds->RailMaterial()),
		       sRailRotat, ds->s2RailCenter(),
		       ds2VacInfo.logical, 0, _config.getBool("ds.visible"),
		       G4Colour::Blue(), _config.getBool("ds.solid"),
		       forceAuxEdgeVisible, placePV, doSurfaceCheck );

    // And now in DS3Vacuum

     VolumeInfo RailN3 = nestExtrudedSolid
                      ( "NorthRailDS3", ds->lengthRail3()/2.0*CLHEP::mm,
 		       uRailOutline, vRailOutline, 
 		       findMaterialOrThrow(ds->RailMaterial()),
 		       nRailRotat, ds->n3RailCenter(),
 		       dsShieldParent, 0, _config.getBool("ds.visible"),
 		       G4Colour::Blue(), _config.getBool("ds.solid"),
 		       forceAuxEdgeVisible, placePV, doSurfaceCheck );

     VolumeInfo RailS3 = nestExtrudedSolid
                      ( "SouthRailDS3", ds->lengthRail3()/2.0*CLHEP::mm,
 		       uRailOutline, vRailOutline, 
 		       findMaterialOrThrow(ds->RailMaterial()),
 		       sRailRotat, ds->s3RailCenter(),
 		       dsShieldParent, 0, _config.getBool("ds.visible"),
 		       G4Colour::Blue(), _config.getBool("ds.solid"),
 		       forceAuxEdgeVisible, placePV, doSurfaceCheck );
    
     // Now put bearing blocks on rails
     // D. No. Brown, Jan 2016
     // First in DS2Vacuum region
     std::vector<double> uBBlockOutline = ds->uOutlineBBlock();
     std::vector<double> vBBlockOutline = ds->vOutlineBBlock();
     CLHEP::HepRotation* BBRotat = new CLHEP::HepRotation( CLHEP::HepRotation::IDENTITY);
     double lenBB2 = ds->lengthBBlock2()/2.0*CLHEP::mm;
     double lenBB3 = ds->lengthBBlock3()/2.0*CLHEP::mm;
     std::vector<CLHEP::Hep3Vector> BBCenters2 = ds->BBlockCenters2();
     std::vector<CLHEP::Hep3Vector> BBCenters3 = ds->BBlockCenters3();
     int cScheme = ds->couplerScheme();

     // First in DS2Vacuum region
     int nB2 = BBCenters2.size();
     CLHEP::Hep3Vector DS2Offset(0,0,ds2Z0);
     int coupleCounter = 0;
     double widCoupler = ds->widthCoupler();
     double hCoupler = ds->heightCoupler();
     double yCoupler = ds->yCenterCoupler();
     for ( int iB2 = 0; iB2 < nB2; iB2++ ) {
       std::stringstream sstm;
       sstm << "BearingBlock_DS2_" << iB2+1;
       VolumeInfo BBlock2 = nestExtrudedSolid
	 ( sstm.str().c_str(), lenBB2,
	   uBBlockOutline, vBBlockOutline, 
	   findMaterialOrThrow(ds->BBlockMaterial()),
	   BBRotat, BBCenters2[iB2] - DS2Offset,
	   ds2VacInfo.logical, 
	   0, _config.getBool("ds.visible"),
	   G4Colour::Blue(), _config.getBool("ds.solid"),
	   forceAuxEdgeVisible, placePV, doSurfaceCheck );
       // Now add Couplers
       if (iB2 < nB2 - 2 && (cScheme == 0 || (cScheme == 1 && BBCenters2[iB2].x() > 0.0) || (cScheme == 2 && BBCenters2[iB2].x() < 0) ) ) {
	 coupleCounter++;
	 std::stringstream couplerName;
	 couplerName << "Coupler_DS2_" << coupleCounter;
	 double lenCoupler = BBCenters2[iB2+2].z() - BBCenters2[iB2].z() 
	   - 2.*lenBB2 - 0.2; // The 0.2 is to avoid accidental overlaps. 
	 // The 2.* is because lenBB2 is actually halfLength
	 CLHEP::Hep3Vector cenCoupler( (BBCenters2[iB2] + BBCenters2[iB2+2] )*0.5 + CLHEP::Hep3Vector(0.0,yCoupler,0.0));
	 std::vector<double> halfDims = { widCoupler/2.0, 
					  hCoupler/2.0, 
					  lenCoupler/2.0 };
	 nestBox( couplerName.str(),
		  halfDims,
		  findMaterialOrThrow(ds->BBlockMaterial()),
		  BBRotat,
		  cenCoupler - DS2Offset,
		  ds2VacInfo.logical,
		  0, _config.getBool("ds.visible"),
		  G4Colour::Blue(), _config.getBool("ds.solid"),
		  forceAuxEdgeVisible, placePV, doSurfaceCheck );
		  
       } // end of if for adding coupler if not last bearing block
     }

     // Now in DS3Vacuum region
     int nB3 = BBCenters3.size();
     coupleCounter = 0;
     for ( int iB3 = 0; iB3 < nB3; iB3++ ) {
       std::stringstream sstm;
       sstm << "BearingBlock_DS3_" << iB3+1;

       VolumeInfo BBlock3 = nestExtrudedSolid
	 ( sstm.str().c_str(), lenBB3,
	   uBBlockOutline, vBBlockOutline, 
	   findMaterialOrThrow(ds->BBlockMaterial()),
	   BBRotat, BBCenters3[iB3],
	   dsShieldParent, 0, _config.getBool("ds.visible"),
	   G4Colour::Blue(), _config.getBool("ds.solid"),
	   forceAuxEdgeVisible, placePV, doSurfaceCheck );
       // Now add Couplers
       if (iB3 < nB3 - 2 && (cScheme == 0 || (cScheme == 1 && BBCenters3[iB3].x() > 0.0) || (cScheme == 2 && BBCenters3[iB3].x() < 0) ) ) {
	 coupleCounter++;
	 std::stringstream couplerName;
	 couplerName << "Coupler_DS3_" << coupleCounter;
	 double lenCoupler = BBCenters3[iB3+2].z() - BBCenters3[iB3].z() 
	   - 2.*lenBB3 - 0.2; // The 0.2 is to avoid accidental overlaps.
	 // The 2.* is because lenBB3 is actually halfLength
	 CLHEP::Hep3Vector cenCoupler( (BBCenters3[iB3] + BBCenters3[iB3+2] )*0.5 + CLHEP::Hep3Vector(0.0,yCoupler,0.0));
	 std::vector<double> halfDims = { widCoupler/2.0, 
					  hCoupler/2.0, 
					  lenCoupler/2.0 };
	 nestBox( couplerName.str(),
		  halfDims,
		  findMaterialOrThrow(ds->BBlockMaterial()),
		  BBRotat,
		  cenCoupler,
		  dsShieldParent,
		  0, _config.getBool("ds.visible"),
		  G4Colour::Blue(), _config.getBool("ds.solid"),
		  forceAuxEdgeVisible, placePV, doSurfaceCheck );
		  
       } // end of if for adding coupler if not last bearing block

 
     }


     //************* End of Rails and Bearing Blocks ************
     //************* Begin MBS Support block *****************
     if ( ds->hasMBSS() ) {
       std::vector<double> uOutlineMBSS = ds->uOutlineMBSS();
       std::vector<double> vOutlineMBSS = ds->vOutlineMBSS();
       VolumeInfo MBSS = nestExtrudedSolid
	 ( "MBSSphericalSupport", ds->lengthMBSS()/2.0*CLHEP::mm,
	   uOutlineMBSS, vOutlineMBSS, 
	   findMaterialOrThrow(ds->MBSSmaterial()),
	   0, ds->MBSSlocation(),
	   dsShieldParent, 0, _config.getBool("ds.visible"),
	   G4Colour::Blue(), _config.getBool("ds.solid"),
	   forceAuxEdgeVisible, placePV, doSurfaceCheck );
     }  // end of if ( ds->hasMBSS() )

     // End of MBS spherical shielding, begin cable runs for Cal and Tracker
     // Each is modeled as a thin wedge of a ring

     bool cableRunSensitive = _config.getBool("ds.CableRun.sensitive",false);

     if ( ds->hasCableRunCal() ) {


       TubsParams  calCableRunParams  ( ds->rInCableRunCal(), 
					ds->rOutCableRunCal(), 
					ds->lengthCableRunCal(),
					ds->phi0CableRunCal()*CLHEP::degree,
					ds->dPhiCableRunCal()*CLHEP::degree);

       CLHEP::Hep3Vector calCableRunLoc( 0.0, 0.0, ds->zCCableRunCal() );

       VolumeInfo ccrTemp = nestTubs( "CalCableRun",
				      calCableRunParams,
				      findMaterialOrThrow(ds->calCableRunMaterial()),
				      0,
				      calCableRunLoc,
				      dsShieldParent,
				      0,
				      G4Color::Magenta(),
				      "DS"
				      );

       if ( cableRunSensitive && sdHelper.enabled(StepInstanceName::DSCableRun) ) {
	 G4VSensitiveDetector* cableRunSD = G4SDManager::GetSDMpointer()->
	   FindSensitiveDetector(SensitiveDetectorName::DSCableRun());
	 if(cableRunSD) ccrTemp.logical->SetSensitiveDetector(cableRunSD);
       }

       if ( ds->cableRunVersion() > 1 ) {

	 // Now the part between the Calorimeter Disks
	 TubsParams  upCalCableRunParm1( ds->upRInCableRunCal(), 
					 ds->upROutCableRunCal(), 
					 ds->upHL1CableRunCal(),
					 ds->phi0CableRunCal()*CLHEP::degree,
					 ds->dPhiCableRunCal()*CLHEP::degree);

	 CLHEP::Hep3Vector upCalCableRunLoc1( 0.0, 0.0,ds->upZC1CableRunCal());

	 VolumeInfo ccrTempUG1 = nestTubs( "CalCableRunUpGap1",
					   upCalCableRunParm1,
					   findMaterialOrThrow(ds->calCableRunMaterial()),
					   0,
					   upCalCableRunLoc1,
					   dsShieldParent,
					   0,
					   G4Color::Magenta(),
					   "DS"
					   );

       if ( cableRunSensitive && sdHelper.enabled(StepInstanceName::DSCableRun) ) {
	 G4VSensitiveDetector* cableRunSD = G4SDManager::GetSDMpointer()->
	   FindSensitiveDetector(SensitiveDetectorName::DSCableRun());
	 if(cableRunSD) ccrTempUG1.logical->SetSensitiveDetector(cableRunSD);
       }


	 TubsParams  upCalCableRunParm2( ds->upRInCableRunCal(), 
					 ds->upROutCableRunCal(), 
					 ds->upHL2CableRunCal(),
					 ds->phi0CableRunCal()*CLHEP::degree,
					 ds->dPhiCableRunCal()*CLHEP::degree);

	 CLHEP::Hep3Vector upCalCableRunLoc2( 0.0, 0.0,ds->upZC2CableRunCal());
	 
	 VolumeInfo ccrTmpUG2 = nestTubs( "CalCableRunUpGap2",
					  upCalCableRunParm2,
					  findMaterialOrThrow(ds->calCableRunMaterial()),
					  0,
					  upCalCableRunLoc2,
					  dsShieldParent,
					  0,
					  G4Color::Magenta(),
					  "DS"
					  );

       if ( cableRunSensitive && sdHelper.enabled(StepInstanceName::DSCableRun) ) {
	 G4VSensitiveDetector* cableRunSD = G4SDManager::GetSDMpointer()->
	   FindSensitiveDetector(SensitiveDetectorName::DSCableRun());
	 if(cableRunSD) ccrTmpUG2.logical->SetSensitiveDetector(cableRunSD);
       }

	 // And last but not least the connector between the top of the Cal
	 // and the top of the MBS
	 // Implement this as a Polycone
	 std::vector<double> zs = { ds->upZC2CableRunCal() + ds->upHL2CableRunCal() + 4.0, ds->zCCableRunCal() - ds->lengthCableRunCal() };
	 std::vector<double> rins = { ds->upRInCableRunCal(), ds->rInCableRunCal()};
	 std::vector<double> routs= { ds->upROutCableRunCal(), ds->rOutCableRunCal()};
	 PolyconsParams myPars( zs, rins, routs,
				ds->phi0CableRunCal()*CLHEP::degree,
				ds->dPhiCableRunCal()*CLHEP::degree );
	 
	 VolumeInfo ccrTmpF = nestPolycone ( "calCableRunFall",
					     myPars,
					     findMaterialOrThrow(ds->calCableRunMaterial()),
					     0,
					     G4ThreeVector(0,0,0),
					     dsShieldParent,
					     0,
					     G4Colour::Magenta(),
					     "DS" );

       if ( cableRunSensitive && sdHelper.enabled(StepInstanceName::DSCableRun) ) {
	 G4VSensitiveDetector* cableRunSD = G4SDManager::GetSDMpointer()->
	   FindSensitiveDetector(SensitiveDetectorName::DSCableRun());
	 if(cableRunSD) ccrTmpF.logical->SetSensitiveDetector(cableRunSD);
       }


       } // end of if ( CableRunVersion > 1 )
     } // end of if ( ds->hasCableRunCal() )

     if ( ds->hasCableRunTrk() ) {

       TubsParams  trkCableRun1Params ( ds->rInCableRunTrk(), 
					ds->rOutCableRunTrk(), 
					ds->lengthCableRunTrk(),
					ds->phi0CableRunTrk()*CLHEP::degree,
					ds->dPhiCableRunTrk()*CLHEP::degree);

       CLHEP::Hep3Vector trkCableRunLoc( 0.0, 0.0, ds->zCCableRunTrk() );

       VolumeInfo tcrTmp1 = nestTubs( "TrkCableRun1",
				      trkCableRun1Params,
				      findMaterialOrThrow(ds->trkCableRunMaterial()),
				      0,
				      trkCableRunLoc,
				      dsShieldParent,
				      0,
				      G4Color::Magenta(),
				      "DS"
				      );

       if ( cableRunSensitive && sdHelper.enabled(StepInstanceName::DSCableRun) ) {
	 G4VSensitiveDetector* cableRunSD = G4SDManager::GetSDMpointer()->
	   FindSensitiveDetector(SensitiveDetectorName::DSCableRun());
	 if(cableRunSD) tcrTmp1.logical->SetSensitiveDetector(cableRunSD);
       }

       // Now the second one
       TubsParams  trkCableRun2Params ( ds->rInCableRunTrk(), 
					ds->rOutCableRunTrk(), 
					ds->lengthCableRunTrk(),
					(180.0 - ds->phi0CableRunTrk()
					 - ds->dPhiCableRunTrk())
					*CLHEP::degree,
					ds->dPhiCableRunTrk()*CLHEP::degree);

       VolumeInfo tcrTmp2=nestTubs( "TrkCableRun2",
				    trkCableRun2Params,
				    findMaterialOrThrow(ds->trkCableRunMaterial()),
				    0,
				    trkCableRunLoc,
				    dsShieldParent,
				    0,
				    G4Color::Magenta(),
				    "DS"
				    );

       if ( cableRunSensitive && sdHelper.enabled(StepInstanceName::DSCableRun) ) {
	 G4VSensitiveDetector* cableRunSD = G4SDManager::GetSDMpointer()->
	   FindSensitiveDetector(SensitiveDetectorName::DSCableRun());
	 if(cableRunSD) tcrTmp2.logical->SetSensitiveDetector(cableRunSD);
       }

       if ( ds->cableRunVersion() > 1 ) {
	 // Now the part between the Calorimeter Disks
	 TubsParams  upTrkCableRunParm1( ds->rInCableRunTrk(), 
					 ds->rOutCableRunTrk(), 
					 ds->upHL1CableRunCal(),
					 ds->phi0CableRunTrk()*CLHEP::degree,
					 ds->dPhiCableRunTrk()*CLHEP::degree);

	 CLHEP::Hep3Vector upTrkCableRunLoc1( 0.0, 0.0,ds->upZC1CableRunCal());

	 VolumeInfo tcrTmpG1=nestTubs( "TrkCableRunGap1",
				       upTrkCableRunParm1,
				       findMaterialOrThrow(ds->trkCableRunMaterial()),
				       0,
				       upTrkCableRunLoc1,
				       dsShieldParent,
				       0,
				       G4Color::Magenta(),
				       "DS"
				       );

       if ( cableRunSensitive && sdHelper.enabled(StepInstanceName::DSCableRun) ) {
	 G4VSensitiveDetector* cableRunSD = G4SDManager::GetSDMpointer()->
	   FindSensitiveDetector(SensitiveDetectorName::DSCableRun());
	 if(cableRunSD) tcrTmpG1.logical->SetSensitiveDetector(cableRunSD);
       }

	 TubsParams  upTrkCableRunParm1a( ds->rInCableRunTrk(), 
					  ds->rOutCableRunTrk(), 
					  ds->upHL1CableRunCal(),
					  (180.0 - ds->phi0CableRunTrk()
					   - ds->dPhiCableRunTrk())
					  *CLHEP::degree,
					  ds->dPhiCableRunTrk()*CLHEP::degree);

	 VolumeInfo tcrTmpG1a=nestTubs( "TrkCableRunGap1a",
					upTrkCableRunParm1a,
					findMaterialOrThrow(ds->trkCableRunMaterial()),
					0,
					upTrkCableRunLoc1,
					dsShieldParent,
					0,
					G4Color::Magenta(),
					"DS"
					);

       if ( cableRunSensitive && sdHelper.enabled(StepInstanceName::DSCableRun) ) {
	 G4VSensitiveDetector* cableRunSD = G4SDManager::GetSDMpointer()->
	   FindSensitiveDetector(SensitiveDetectorName::DSCableRun());
	 if(cableRunSD) tcrTmpG1a.logical->SetSensitiveDetector(cableRunSD);
       }


	 TubsParams  upTrkCableRunParm2( ds->rInCableRunTrk(), 
					 ds->rOutCableRunTrk(), 
					 ds->upHL2CableRunCal(),
					 ds->phi0CableRunTrk()*CLHEP::degree,
					 ds->dPhiCableRunTrk()*CLHEP::degree);

	 CLHEP::Hep3Vector upTrkCableRunLoc2( 0.0, 0.0,ds->upZC2CableRunCal());
	 
	 VolumeInfo tcrTmpG2=nestTubs( "TrkCableRunGap2",
				       upTrkCableRunParm2,
				       findMaterialOrThrow(ds->trkCableRunMaterial()),
				       0,
				       upTrkCableRunLoc2,
				       dsShieldParent,
				       0,
				       G4Color::Magenta(),
				       "DS"
				       );

	 if ( cableRunSensitive && sdHelper.enabled(StepInstanceName::DSCableRun) ) {
	   G4VSensitiveDetector* cableRunSD = G4SDManager::GetSDMpointer()->
	     FindSensitiveDetector(SensitiveDetectorName::DSCableRun());
	   if(cableRunSD) tcrTmpG2.logical->SetSensitiveDetector(cableRunSD);
	 }


	 TubsParams  upTrkCableRunParm2a( ds->rInCableRunTrk(), 
					  ds->rOutCableRunTrk(), 
					  ds->upHL2CableRunCal(),
					  (180.0 - ds->phi0CableRunTrk()
					   - ds->dPhiCableRunTrk())
					  *CLHEP::degree,
					  ds->dPhiCableRunTrk()*CLHEP::degree);

	 VolumeInfo tcrTmpG2a= nestTubs( "TrkCableRunGap2a",
					 upTrkCableRunParm2a,
					 findMaterialOrThrow(ds->trkCableRunMaterial()),
					 0,
					 upTrkCableRunLoc2,
					 dsShieldParent,
					 0,
					 G4Color::Magenta(),
					 "DS"
					 );

	 if ( cableRunSensitive && sdHelper.enabled(StepInstanceName::DSCableRun) ) {
	   G4VSensitiveDetector* cableRunSD = G4SDManager::GetSDMpointer()->
	     FindSensitiveDetector(SensitiveDetectorName::DSCableRun());
	   if(cableRunSD) tcrTmpG2a.logical->SetSensitiveDetector(cableRunSD);
	 }


       } // end of adding gap runs for trk cable runs
     } // end of if ( ds->hasCableRunTrk() )

     if ( ds->hasServicePipes() ) {
       TubsParams pipeMomParms( 0.0, 
				ds->servicePipeROut(),
				ds->servicePipeHalfLength());//default phi,dphi
       TubsParams pipeSelfParms(ds->servicePipeRIn(),
				ds->servicePipeROut(),
				ds->servicePipeHalfLength() );
       std::vector<double>theXs = ds->servicePipeXCs();
       double theY = ds->servicePipeYC();
       double theZ = ds->servicePipeZC();

       for ( unsigned int iP = 0; iP < theXs.size(); iP++ ) {
	 CLHEP::Hep3Vector pipeLoc(theXs[iP],theY, theZ);
	 ostringstream pmName;
	 pmName << "DSservicePipeMother" << iP+1;
	 VolumeInfo DSPipeMother = nestTubs ( pmName.str(),
					      pipeMomParms,
					      findMaterialOrThrow(ds->servicePipeFillMat()),
					      0,
					      pipeLoc,
					      dsShieldParent,
					      0,
					      G4Color::Magenta(),
					      "DS"
					      );

	 // Now make it a real pipe
	 ostringstream pName;
	 pName << "DSservicePipe" << iP+1;
	 nestTubs ( pName.str(),
		    pipeSelfParms,
		    findMaterialOrThrow(ds->servicePipeMaterial()),
		    0,
		    CLHEP::Hep3Vector(0,0,0),
		    DSPipeMother,
		    0,
		    G4Color::Magenta(),
		    "DS"
		    );

       } // End of loop over service pipes
     } // end of if ( ds->hasServicePipes() )

  } // end of Mu2eWorld::constructDS;

}
