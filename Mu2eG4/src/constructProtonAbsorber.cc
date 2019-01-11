//
// Free function to create Proton Absorber
//
// $Id: constructProtonAbsorber.cc,v 1.28 2014/02/28 21:11:19 knoepfel Exp $
// $Author: knoepfel $
// $Date: 2014/02/28 21:11:19 $
//
// Original author KLG based on Mu2eWorld constructProtonAbs
//
// Notes:
// Construct the  Proton Absorber

// C++ includes
#include <iostream>
#include <string>
#include <sstream>
#include <vector>

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"

// Mu2e includes.
#include "Mu2eG4/inc/constructProtonAbsorber.hh"
#include "G4Helper/inc/VolumeInfo.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/VirtualDetector.hh"
#include "StoppingTargetGeom/inc/StoppingTarget.hh"
#include "G4Helper/inc/G4Helper.hh"
#include "GeomPrimitives/inc/TubsParams.hh"
#include "GeomPrimitives/inc/PolyhedraParams.hh"
#include "Mu2eG4/inc/findMaterialOrThrow.hh"
#include "Mu2eG4/inc/MaterialFinder.hh"
#include "Mu2eG4/inc/nestBox.hh"
#include "Mu2eG4/inc/nestCons.hh"
#include "Mu2eG4/inc/nestTubs.hh"
#include "Mu2eG4/inc/nestPolyhedra.hh"
#include "Mu2eG4/inc/finishNesting.hh"
#include "Mu2eG4/inc/HelicalProtonAbsorber.hh"
#include "MECOStyleProtonAbsorberGeom/inc/MECOStyleProtonAbsorber.hh"
#include "DetectorSolenoidGeom/inc/DetectorSolenoid.hh"
#include "Mu2eG4/inc/checkForOverlaps.hh"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Units/SystemOfUnits.h"

// G4 includes
#include "G4Material.hh"
#include "G4Color.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Tubs.hh"
#include "G4BooleanSolid.hh"
#include "G4VPhysicalVolume.hh"
#include "G4SubtractionSolid.hh"
#include "G4SDManager.hh"

using namespace std;

namespace mu2e {

  void constructProtonAbsorber( const SimpleConfig& _config ){

    if( !_config.getBool("hasProtonAbsorber", true) ) return;
    
    int  const verbosityLevel           = _config.getInt("protonabsorber.verbosityLevel", 0);
      
    // Access to the G4HelperService.
    G4Helper* _helper = &(*(art::ServiceHandle<G4Helper>()));
    
    VolumeInfo const & parent1Info  = _helper->locateVolInfo("DS2Vacuum");
    VolumeInfo const & parent2Info  = _helper->locateVolInfo("DS3Vacuum");
    
    // Fetch DS geometry
    GeomHandle<DetectorSolenoid> ds;


    const auto geomOptions = art::ServiceHandle<GeometryService>()->geomOptions();
    geomOptions->loadEntry( _config, "protonabsorber", "protonabsorber");

    const bool pabsIsVisible       = geomOptions->isVisible("protonabsorber"); 
    const bool pabsIsSolid         = geomOptions->isSolid("protonabsorber"); 
    const bool forceAuxEdgeVisible = geomOptions->forceAuxEdgeVisible("protonabsorber"); 
    const bool doSurfaceCheck      = geomOptions->doSurfaceCheck("protonabsorber"); 
    const bool placePV             = geomOptions->placePV("protonabsorber"); 


    //**** Helical proton absorber ****//
    if ( _config.getBool("protonabsorber.isHelical", false) ) {
      if ( verbosityLevel > 0) cout << __func__ << " : Proton Absorber is Helical type" << endl;
      MaterialFinder materialFinder(_config);
      G4Material* pabsMaterial = materialFinder.get("protonabsorber.materialName");
      double ds2HalfLen = ds->vac_halfLengthDs2();
      double helHalfPabsLength = _config.getDouble("protonabsorber.halfLength");
      double helPabsLength = 2.0*helHalfPabsLength;
      double lengthScaleFact = 130.0/274.0; //I don't know why but in the implementation of the absorber the length is divided into 130 steps of Z and that there is a cycle over 274 number of Z steps. ????? //FIXME
      double helPabsThickness = _config.getDouble("protonabsorber.thickness");
      //double numOfTurns = _config.getDouble("protonabsorber.NumOfTurns");
      double innerRadii[3], outerRadii[3], innerPhis[3], outerPhis[3];
      innerRadii[0] = _config.getDouble("protonabsorber.InnerRadius_0");
      innerRadii[1] = _config.getDouble("protonabsorber.InnerRadius_1");
      innerRadii[2] = _config.getDouble("protonabsorber.InnerRadius_2");
      
      outerRadii[0] = _config.getDouble("protonabsorber.OuterRadius_0");
      outerRadii[1] = _config.getDouble("protonabsorber.OuterRadius_1");
      outerRadii[2] = _config.getDouble("protonabsorber.OuterRadius_2");
      
      innerPhis[0] = _config.getDouble("protonabsorber.phiInner_0");
      innerPhis[1] = _config.getDouble("protonabsorber.phiInner_1");
      innerPhis[2] = _config.getDouble("protonabsorber.phiInner_2");
      
      outerPhis[0] = _config.getDouble("protonabsorber.phiOuter_0");
      outerPhis[1] = _config.getDouble("protonabsorber.phiOuter_1");
      outerPhis[2] = _config.getDouble("protonabsorber.phiOuter_2");
      
      //bool addSD   = _config.getBool("protonabsorber.saveStepPnts",false);
      
      HelicalProtonAbsorber* hpabs = 
        new HelicalProtonAbsorber( ds2HalfLen-helPabsLength,
                                   helPabsLength*lengthScaleFact,/*ScaleFact is a trick FIXME */
                                   innerRadii, outerRadii, innerPhis, outerPhis,
                                   helPabsThickness, /*numOfTurns,*/
                                   _config.getInt("protonabsorber.NumOfVanes"), 
                                   pabsMaterial, parent1Info.logical );
      
      if ( _config.getBool("g4.doSurfaceCheck",false) || _config.getBool("protonAbsorber.doSurfaceCheck",false) ) {
        checkForOverlaps( hpabs->GetPhys(), _config, verbosityLevel>0);
      }
      
      if ( _config.getBool("protonabsorber.visible",true) ) {
        AntiLeakRegistry & reg = _helper->antiLeakRegistry();
        hpabs->SetVisibility( _config.getBool("protonabsorber.solid",true),
                              forceAuxEdgeVisible,
                              G4Color::Green(), reg);
      }
      
      delete hpabs;
      
    } 

    //**** MECO-Style proton absorber ****//
    else if (! _config.getBool("protonabsorber.isShorterCone", false) ) {
      
      if ( verbosityLevel > 0) cout << __func__ << " : Proton Absorber is MECO-style Conical type" << endl;
      
      // smaller and larger outer radii
      double pabs1rOut0   = _config.getDouble("protonabsorber.OutRadius0");
      double pabs2rOut1   = _config.getDouble("protonabsorber.OutRadius1");
      double pabsZHalfLen = _config.getDouble("protonabsorber.halfLength");
      double thick        = _config.getDouble("protonabsorber.thickness");
      
      // adding virtual detector before and after target
      double vdHL = 0.;
      art::ServiceHandle<GeometryService> geom;
      if( geom->hasElement<VirtualDetector>() ) {
        GeomHandle<VirtualDetector> vdg;
        if( vdg->nDet()>0 ) vdHL = vdg->getHalfLength();
      }
      
      // subtract virtual detector thickness from the larger outer
      // radius of the proton absorber
      
      pabs2rOut1 -= 2.*vdHL;
      
      double pabs1rIn0  = pabs1rOut0 - thick;
      double pabs2rIn1  = pabs2rOut1 - thick;
      
      MaterialFinder materialFinder(_config);
      G4Material* pabsMaterial = materialFinder.get("protonabsorber.materialName");

      GeomHandle<StoppingTarget> target;

      // The proton absorber starts at the target end.
      // we add space for the virtual detector here
      double pabsStartInMu2eZ = target->centerInMu2e().z() + 0.5*target->cylinderLength() + 2.*vdHL;;

      // Need to split it at the DS2/DS3 boundary
      double pabs1EndInMu2eZ = ds->vac_zLocDs23Split();

      double pabs1len = pabs1EndInMu2eZ - pabsStartInMu2eZ;

      G4ThreeVector  pabs1Offset(0.0, 0.0, (pabs1EndInMu2eZ+pabsStartInMu2eZ)/2. - parent1Info.centerInMu2e().z());

      if ( verbosityLevel > 0) {
        cout << __func__ <<
          " pabs1len                          : " << pabs1len << endl;
      }

      // interpolating the outer radius of the DS2 part

      double pabs1rOut1 = ((pabs2rOut1 - pabs1rOut0)*(pabs1len/(2.0*pabsZHalfLen))) + pabs1rOut0;
      double pabs1rIn1  = pabs1rOut1 - thick;

      double pabs2len     = (2.0*pabsZHalfLen) - pabs1len;
      double pabs2ZOffset = (pabs2len*0.5) +  ds->vac_zLocDs23Split();

      // protonabs2 should touch protonabs1 and both of them should touch the ds2/ds3 boundary
      // it looks like the boolean solid center is in the center of ConstituentSolid(0)
      // note that the subtraction solid may be shorter, depending how the subtraction is done

      if ( verbosityLevel > 0) {
        double theZ  = parent1Info.centerInMu2e()[CLHEP::Hep3Vector::Z];
        cout << __func__ << " " << parent1Info.name << " Z offset in Mu2e    : " <<
          theZ << endl;
        cout << __func__ << " " << parent1Info.name << " Z extent in Mu2e    : " <<
          theZ - ds->vac_halfLengthDs2() << ", " << theZ + ds->vac_halfLengthDs2() << endl;
      }

      if ( verbosityLevel > 0) {
        cout << __func__ << " " << parent2Info.name << " Z extent in Mu2e    : " <<
          ds->vac_zLocDs23Split() << ", " << ds->cryoZMax() << endl;
      }

      G4ThreeVector pabs2Offset(0.0, 0.0, pabs2ZOffset);

      // proton absorber in DS2
      double pabs1Param[7] = { pabs1rIn0, pabs1rOut0, pabs1rIn1, pabs1rOut1, pabs1len*0.5,
                               0.0, 360.0*CLHEP::degree };

      if ( verbosityLevel > 0 ) {
        mf::LogInfo log("GEOM");
        log << "Constructing Proton Absorber -- \n";
        log << "Proton Abs Offset in DS2:  " << pabs1Offset <<"\n";
        log << "rIn,  rOut (-z): "<< pabs1rIn0 <<"  "<< pabs1rOut0<<"  ";
        log << "rIn,  rOut (+z): "<< pabs1rIn1 <<"  "<< pabs1rOut1<<"  ";
        log << "halflength: "<< pabs1len*0.5 <<"\n";
        log << "Proton Abs Offset in DS3:  " << pabs2Offset <<"\n";
        log << "rIn,  rOut (-z): "<< pabs1rIn1 <<"  "<< pabs1rOut1<<"  ";
        log << "rIn,  rOut (+z): "<< pabs2rIn1 <<"  "<< pabs2rOut1<<"  ";
        log << "halflength: "<< pabs2len*0.5 <<"\n";
      }

      VolumeInfo protonabs1Info = nestCons( "protonabs1",
                                            pabs1Param,
                                            pabsMaterial,
                                            0,
                                            pabs1Offset,
                                            parent1Info,
                                            0,
                                            pabsIsVisible,
                                            G4Color::White(),
                                            pabsIsSolid,
                                            forceAuxEdgeVisible,
                                            placePV,
                                            doSurfaceCheck
                                            );
        

      if ( verbosityLevel > 0) {
        double pzhl   = static_cast<G4Cons*>(protonabs1Info.solid)->GetZHalfLength();
        double pabs1Z = protonabs1Info.centerInMu2e()[CLHEP::Hep3Vector::Z];
        cout << __func__ << " protonabs1 Z offset in Mu2e    : " <<
          pabs1Z << endl;
        cout << __func__ << " protonabs1 Z extent in Mu2e    : " <<
          pabs1Z - pzhl  << ", " <<  pabs1Z + pzhl  << endl;
      }


      // proton absorber in DS3
      double pabs2Param[7] = { pabs1rIn1, pabs1rOut1, pabs2rIn1, pabs2rOut1, pabs2len*0.5,
                               0.0, 360.0*CLHEP::degree };

      VolumeInfo protonabs2Info = nestCons( "protonabs2",
                                            pabs2Param,
                                            pabsMaterial,
                                            0,
                                            pabs2Offset,
                                            parent2Info,
                                            0,
                                            pabsIsVisible,
                                            G4Color::Yellow(),
                                            pabsIsSolid,
                                            forceAuxEdgeVisible,
                                            placePV,
                                            doSurfaceCheck
                                            );
        
      if ( verbosityLevel > 0) {
        double pzhl   = static_cast<G4Cons*>(protonabs2Info.solid)->GetZHalfLength();
        double pabs2Z = protonabs2Info.centerInMu2e()[CLHEP::Hep3Vector::Z];
        cout << __func__ << " " << protonabs2Info.name << " Z offset in Mu2e    : " <<
          pabs2Z << endl;
        cout << __func__ << " " << protonabs2Info.name << " Z extent in Mu2e    : " <<
          pabs2Z - pzhl  << ", " <<  pabs2Z + pzhl  << endl;

        // we also check how the offsets are handled

        cout << __func__ << " " << protonabs2Info.name << " local input offset in G4                  : " <<
          pabs2Offset << endl;
        cout << __func__ << " " << protonabs2Info.name << " local GetTranslation()       offset in G4 : " <<
          protonabs2Info.physical->GetTranslation() << endl; // const &

      }

    } 

    //**** Modifiable MECO-Style proton absorber ****//
    else {
      if ( verbosityLevel > 0) cout << __func__ << " : Proton Absorber is modified Conical type" << endl;

      GeomHandle<MECOStyleProtonAbsorber> pabs;

      double ds2zcenter   = parent1Info.centerInMu2e().z();
      double ds3zcenter   = parent2Info.centerInMu2e().z();

      if ( verbosityLevel > 0) {
        cout << __func__ <<
          " ds2zcenter               : " << ds2zcenter << endl;
        cout << __func__ <<
          " ds3zcenter               : " << ds3zcenter << endl;
      }
      //G4Material* pabsMaterial =  findMaterialOrThrow(pabs->fillMaterial());
      MaterialFinder materialFinder(_config);
      G4Material* pabsMaterial = materialFinder.get("protonabsorber.materialName");

      if ( verbosityLevel > 0) {
        double theZ  = parent1Info.centerInMu2e()[CLHEP::Hep3Vector::Z];
        cout << __func__ << " " << parent1Info.name << " Z offset in Mu2e    : " <<
          theZ << endl;
        cout << __func__ << " " << parent1Info.name << " Z extent in Mu2e    : " <<
          theZ - ds->vac_halfLengthDs2() << ", " << theZ + ds->vac_halfLengthDs2() << endl;
      }

      double pabs1ZOffset = 0, pabs2ZOffset = 0, pabs3ZOffset = 0, pabs4ZOffset = 0;
      if (pabs->isAvailable(0)) pabs1ZOffset = CLHEP::mm * pabs->part(0).center().z() - ds2zcenter;
      if (pabs->isAvailable(1)) pabs2ZOffset = CLHEP::mm * pabs->part(1).center().z() - ds3zcenter;
      if (pabs->isAvailable(2)) pabs3ZOffset = CLHEP::mm * pabs->part(2).center().z() - ds2zcenter;
      if (pabs->isAvailable(3)) pabs4ZOffset = CLHEP::mm * pabs->part(3).center().z() - ds3zcenter;

      G4ThreeVector pabs1Offset(0.0, 0.0, pabs1ZOffset);
      G4ThreeVector pabs2Offset(0.0, 0.0, pabs2ZOffset);
      G4ThreeVector pabs3Offset(0.0, 0.0, pabs3ZOffset);
      G4ThreeVector pabs4Offset(0.0, 0.0, pabs4ZOffset);

      // proton absorber in DS2
      if (pabs->isAvailable(0)) {
        double pabs1rIn0  = pabs->part(0).innerRadiusAtStart();
        double pabs1rOut0 = pabs->part(0).outerRadiusAtStart();
        double pabs1rIn1  = pabs->part(0).innerRadiusAtEnd();
        double pabs1rOut1 = pabs->part(0).outerRadiusAtEnd();
        double pabs1len   = pabs->part(0).halfLength() *2.;
  
        double pabs1Param[7] = { pabs1rIn0, pabs1rOut0, pabs1rIn1, pabs1rOut1, pabs1len*0.5,
                                 0.0, 360.0*CLHEP::degree };
  
        VolumeInfo protonabs1Info = nestCons( "protonabs1",
                                              pabs1Param,
                                              pabsMaterial,
                                              0,
                                              pabs1Offset,
                                              parent1Info,
                                              0,
                                              pabsIsVisible,
                                              G4Color::White(),
                                              pabsIsSolid,
                                              forceAuxEdgeVisible,
                                              placePV,
                                              doSurfaceCheck
                                              );
          
 
        if ( verbosityLevel > 0) {
          double pzhl   = static_cast<G4Cons*>(protonabs1Info.solid)->GetZHalfLength();
          double pabs1Z = protonabs1Info.centerInMu2e()[CLHEP::Hep3Vector::Z];
          cout << __func__ << " protonabs1 Z offset in Mu2e    : " <<
            pabs1Z << endl;
          cout << __func__ << " protonabs1 Z extent in Mu2e    : " <<
            pabs1Z - pzhl  << ", " <<  pabs1Z + pzhl  << endl;
        }
      }
      else {
        if( verbosityLevel > 0 ) cout << __func__ << " protonabs1 disabled" << endl;
      }


      // proton absorber in DS3
      if (pabs->isAvailable(1)) {
        double pabs2rIn0  = pabs->part(1).innerRadiusAtStart();
        double pabs2rOut0 = pabs->part(1).outerRadiusAtStart();
        double pabs2rIn1  = pabs->part(1).innerRadiusAtEnd();
        double pabs2rOut1 = pabs->part(1).outerRadiusAtEnd();
        double pabs2len   = pabs->part(1).halfLength() * 2.;
  
        double pabs2Param[7] = { pabs2rIn0, pabs2rOut0, pabs2rIn1, pabs2rOut1, pabs2len*0.5,
                                 0.0, 360.0*CLHEP::degree };
  
        VolumeInfo protonabs2Info = nestCons( "protonabs2",
                                              pabs2Param,
                                              pabsMaterial,
                                              0,
                                              pabs2Offset,
                                              parent2Info,
                                              0,
                                              pabsIsVisible,
                                              G4Color::Yellow(),
                                              pabsIsSolid,
                                              forceAuxEdgeVisible,
                                              placePV,
                                              doSurfaceCheck
                                              );
          
          
        if ( verbosityLevel > 0) {
          double pzhl   = static_cast<G4Cons*>(protonabs2Info.solid)->GetZHalfLength();
          double pabs2Z = protonabs2Info.centerInMu2e()[CLHEP::Hep3Vector::Z];
          cout << __func__ << " " << protonabs2Info.name << " Z offset in Mu2e    : " <<
            pabs2Z << endl;
          cout << __func__ << " " << protonabs2Info.name << " Z extent in Mu2e    : " <<
            pabs2Z - pzhl  << ", " <<  pabs2Z + pzhl  << endl;
          // we also check how the offsets are handled
          cout << __func__ << " " << protonabs2Info.name << " local input offset in G4                  : " <<
            pabs2Offset << endl;
          cout << __func__ << " " << protonabs2Info.name << " local GetTranslation()       offset in G4 : " <<
            protonabs2Info.physical->GetTranslation() << endl; // const &
        }
      }
      else {
        if ( verbosityLevel > 0 ) cout << __func__ << " protonabs2 disabled" << endl;
      }

      /*
      mf::LogInfo log("GEOM");
      log << "Constructing Proton Absorber -- \n";
      if (pabs->isAvailable(0)) {
        log << "Proton Abs Offset in DS2:  " << pabs1Offset <<"\n";
        log << "rIn,  rOut (-z): "<< pabs1rIn0 <<"  "<< pabs1rOut0<<"  ";
        log << "rIn,  rOut (+z): "<< pabs1rIn1 <<"  "<< pabs1rOut1<<"  ";
        log << "halflength: "<< pabs1len*0.5 <<"\n";
      }
      if (pabs->isAvailable(1)) {
        log << "Proton Abs Offset in DS3:  " << pabs2Offset <<"\n";
        log << "rIn,  rOut (-z): "<< pabs2rIn0 <<"  "<< pabs2rOut0<<"  ";
        log << "rIn,  rOut (+z): "<< pabs2rIn1 <<"  "<< pabs2rOut1<<"  ";
        log << "halflength: "<< pabs2len*0.5 <<"\n";
      }*/
      //TODO
  
      // ****************** OPA! *******************
      // outer proton absorber
      // ****************** OPA! *******************

      if (pabs->isAvailable(2)) {
        double pabs3rIn0  = pabs->part(2).innerRadiusAtStart();
        double pabs3rOut0 = pabs->part(2).outerRadiusAtStart();
        double pabs3rIn1  = pabs->part(2).innerRadiusAtEnd();
        double pabs3rOut1 = pabs->part(2).outerRadiusAtEnd();
        double pabs3len   = pabs->part(2).halfLength() * 2.;
	int    pabs3nS    = pabs->part(2).nSides();
	double pabs3SlotWidth = pabs->slotWidth();
	double pabs3SlotLength = pabs->slotLength();
	double pabs3SlotOffset = pabs->slotOffset();

        G4Material* pabs3Material = materialFinder.get("protonabsorber.outerPAMaterialName");
  
	if ( 0 == pabs3nS ) {
	  // This is the "classic" implementation as a conical frustrum
	  double pabs3Param[7] = { pabs3rIn0, pabs3rOut0, pabs3rIn1, pabs3rOut1, pabs3len*0.5,
                                 0.0, 360.0*CLHEP::degree };
  
	  VolumeInfo protonabs3Info = nestCons( "protonabs3",
						pabs3Param,
						pabs3Material,
						0,
						pabs3Offset,
						parent1Info,
						0,
						pabsIsVisible,
						G4Color::Yellow(),
						pabsIsSolid,
						forceAuxEdgeVisible,
						placePV,
						doSurfaceCheck
						);
	} else {
	  // This is an implementation as a polyhedron - like a barrel of
	  // slats.  We also put slots in it for stopping target support
	  // wires to go through
	  double zs[2] = { -pabs3len/2.0, pabs3len/2.0 };
	  double rins[2] = { pabs3rIn0, pabs3rIn1 };
	  double rous[2] = { pabs3rOut0, pabs3rOut1 };

	  VolumeInfo protonabs3Info("protonabs3",
				    pabs3Offset,
				    parent1Info.centerInWorld);

	  double turn = 180.0*CLHEP::degree/((double)pabs3nS);

	  G4Polyhedra* anOPApolyhedron = new G4Polyhedra("protonabs3PH",
							 turn,
							 360.0*CLHEP::degree,
							 pabs3nS,
							 2,
							 zs,
							 rins,
							 rous);


	  G4Box* slot = new G4Box ( "slotBox",
				     pabs3SlotWidth/2.0,
				     400*CLHEP::mm, pabs3SlotLength/2.0 );

	  CLHEP::HepRotation* slot1Rotat = new CLHEP::HepRotation(CLHEP::HepRotation::IDENTITY);
	  slot1Rotat->rotateZ(120.0*CLHEP::degree);

	  CLHEP::Hep3Vector slot1spot(300.0*CLHEP::mm*sin(120.0*CLHEP::degree),
				      300.0*CLHEP::mm*cos(120.0*CLHEP::degree),
				      pabs3SlotOffset*CLHEP::mm);

	  G4SubtractionSolid* aSolid = new G4SubtractionSolid( "FirstStep",
							       anOPApolyhedron,
							       slot,
							       slot1Rotat,
							       slot1spot );

	  CLHEP::HepRotation* slot2Rotat = new CLHEP::HepRotation(CLHEP::HepRotation::IDENTITY);

	  CLHEP::Hep3Vector slot2spot(0., 300.0*CLHEP::mm, pabs3SlotOffset*CLHEP::mm );

	  G4SubtractionSolid* bSolid = new G4SubtractionSolid( "Step2",
							       aSolid,
							       slot,
							       slot2Rotat,
							       slot2spot );

	  CLHEP::HepRotation* slot3Rotat = new CLHEP::HepRotation(CLHEP::HepRotation::IDENTITY);
	  slot3Rotat->rotateZ(240.0*CLHEP::degree);

	  CLHEP::Hep3Vector slot3spot(300.0*CLHEP::mm*sin(240.0*CLHEP::degree),
				      300.0*CLHEP::mm*cos(240.0*CLHEP::degree),
				      pabs3SlotOffset*CLHEP::mm);

	  G4SubtractionSolid* cSolid = new G4SubtractionSolid( protonabs3Info.name,
							       bSolid,
							       slot,
							       slot3Rotat,
							       slot3spot );

	  protonabs3Info.solid = cSolid;

	  finishNesting( protonabs3Info,
			 pabs3Material,
			 0,
			 pabs3Offset,
			 parent1Info.logical,
			 0,
			 pabsIsVisible,
			 G4Color::Yellow(),
			 pabsIsSolid,
			 forceAuxEdgeVisible,
			 placePV,
			 doSurfaceCheck);
	}

   //      mf::LogInfo log("GEOM");
   //      log << "Constructing Proton Absorber -- \n";
   //      log << "Proton Abs Offset in DS2:  " << pabs1Offset <<"\n";
   //      log << "rIn,  rOut (-z): "<< pabs1rIn0 <<"  "<< pabs1rOut0<<"  ";
   //      log << "rIn,  rOut (+z): "<< pabs1rIn1 <<"  "<< pabs1rOut1<<"  ";
   //      log << "halflength: "<< pabs1len*0.5 <<"\n";
   //      log << "Proton Abs Offset in DS3:  " << pabs3Offset <<"\n";
   //      log << "rIn,  rOut (-z): "<< pabs3rIn0 <<"  "<< pabs3rOut0<<"  ";
   //      log << "rIn,  rOut (+z): "<< pabs3rIn1 <<"  "<< pabs3rOut1<<"  ";
   //      log << "halflength: "<< pabs3len*0.5 <<"\n";

   //      if ( verbosityLevel > 0) {
   //        double pzhl   = static_cast<G4Cons*>(protonabs2Info.solid)->GetZHalfLength();
   //        double pabs3Z = protonabs2Info.centerInMu2e()[CLHEP::Hep3Vector::Z];
   //        cout << __func__ << " " << protonabs2Info.name << " Z offset in Mu2e    : " <<
   //          pabs3Z << endl;
   //        cout << __func__ << " " << protonabs2Info.name << " Z extent in Mu2e    : " <<
   //          pabs3Z - pzhl  << ", " <<  pabs3Z + pzhl  << endl;
  
   //        // we also check how the offsets are handled
  
   //        cout << __func__ << " " << protonabs2Info.name << " local input offset in G4                  : " <<
   //          pabs3Offset << endl;
   //        cout << __func__ << " " << protonabs2Info.name << " local GetTranslation()       offset in G4 : " <<
   //          protonabs2Info.physical->GetTranslation() << endl; // const &
   //      }
   // //TODO turn on logging
  
      }
      else {
        if ( verbosityLevel > 0 ) cout << __func__ << " outer protonabs1 disabled" << endl;
      }

      
      if (pabs->isAvailable(3)) {
        double pabs4rIn0  = pabs->part(3).innerRadiusAtStart();
        double pabs4rOut0 = pabs->part(3).outerRadiusAtStart();
        double pabs4rIn1  = pabs->part(3).innerRadiusAtEnd();
        double pabs4rOut1 = pabs->part(3).outerRadiusAtEnd();
        double pabs4len   = pabs->part(3).halfLength() * 2.;
	int    pabs4nS    = pabs->part(3).nSides();
        G4Material* pabs4Material = materialFinder.get("protonabsorber.outerPAMaterialName");
	VolumeInfo protonabs4Info;

	if ( 0 == pabs4nS ) {
	  double pabs4Param[7] = { pabs4rIn0, pabs4rOut0, pabs4rIn1, pabs4rOut1, pabs4len*0.5,
				   0.0, 360.0*CLHEP::degree };
  
	  protonabs4Info = nestCons( "protonabs4",
				     pabs4Param,
				     pabs4Material,
				     0,
				     pabs4Offset,
				     parent2Info,
				     0,
				     pabsIsVisible,
				     G4Color::Yellow(),
				     pabsIsSolid,
				     forceAuxEdgeVisible,
				     placePV,
				     doSurfaceCheck
				     );
	} else {
	  std::vector<double> zs = { -pabs4len/2.0, pabs4len/2.0 };
	  std::vector<double> rins = { pabs4rIn0, pabs4rIn1 };
	  std::vector<double> rous = { pabs4rOut0, pabs4rOut1 };
	  double turn = 180.0*CLHEP::degree/((double)pabs4nS);

	  PolyhedraParams pabs4Param( pabs4nS,
				      zs,
				      rins,
				      rous,
				      turn );
	  
	  protonabs4Info = nestPolyhedra( "protonabs4",
					  pabs4Param,
					  pabs4Material,
					  0,
					  pabs4Offset,
					  parent2Info,
					  0,
					  pabsIsVisible,
					  G4Color::Yellow(),
					  pabsIsSolid,
					  forceAuxEdgeVisible,
					  placePV,
					  doSurfaceCheck
					  );
	}


        if ( verbosityLevel > 0) {
          double pzhl   = static_cast<G4Cons*>(protonabs4Info.solid)->GetZHalfLength();
          double pabs4Z = protonabs4Info.centerInMu2e()[CLHEP::Hep3Vector::Z];
          cout << __func__ << " " << protonabs4Info.name << " Z offset in Mu2e    : " <<
            pabs4Z << endl;
          cout << __func__ << " " << protonabs4Info.name << " Z extent in Mu2e    : " <<
            pabs4Z - pzhl  << ", " <<  pabs4Z + pzhl  << endl;
	  
          // we also check how the offsets are handled	  
          cout << __func__ << " " << protonabs4Info.name << " local input offset in G4                  : " <<
            pabs4Offset << endl;
          cout << __func__ << " " << protonabs4Info.name << " local GetTranslation()       offset in G4 : " <<
            protonabs4Info.physical->GetTranslation() << endl; // const &
        }
  
      }
      else {
        if ( verbosityLevel > 0 ) cout << __func__ << " outer protonabs2 disabled" << endl;
      }

      //***************************
      // Now build the support rings for the OPA
      //***************************

      if ( pabs->oPAnSupports() > 0 ) {
	G4Material* oPAsupportMaterial = 
	  findMaterialOrThrow( pabs->oPAsupportMaterial() );
	double pabs1EndInMu2eZ = ds->vac_zLocDs23Split();
       
	for (int iSup = 0; iSup < pabs->oPAnSupports(); iSup++ ) {
	  double zm = pabs->oPAsupportZMidpoint().at(iSup);
	  double rin= pabs->oPAsupportInnerRadii().at(iSup);
	  double rou= pabs->oPAsupportOuterRadii().at(iSup);
	  double hl = pabs->oPAsupportHalflength().at(iSup);
	  bool hasExtra = (pabs->oPAsupportExtra().at(iSup) > 0.0 );
	  double xRad = 0.0;
	  double dPhiX = 0.0;
	  if ( hasExtra ) {
	    xRad = pabs->oPAsupportXRad().at(iSup);
	    dPhiX = pabs->oPAsupportDPhiX().at(iSup);
	  }

	  CLHEP::Hep3Vector location(-3904.0,0.0,zm);
	  ostringstream myName;
	  myName << "OPAsupport_" << iSup+1;
	  ostringstream myName2;
	  myName2 << "OPAsupportExtra_" << iSup+1;

	  if ( iSup == 1 ) {
	    // Make notches in this one for the Stopping target support slats
	    VolumeInfo ring2Info(myName.str(),
				 location - parent1Info.centerInMu2e(),
				 parent1Info.centerInWorld);

	    G4Tubs* aRingTub = new G4Tubs("TheRing",
					  rin, rou, hl,
					  0.0, 360.0*CLHEP::degree);

	    G4Box* notch = new G4Box("notch", 35.0*CLHEP::mm,
				      25.0*CLHEP::mm, hl*1.1 );

	    CLHEP::HepRotation* notch1Rotat = new CLHEP::HepRotation(CLHEP::HepRotation::IDENTITY);
	    
	    CLHEP::Hep3Vector notch1Locat(0,rin+22.0*CLHEP::mm,0);

	    G4SubtractionSolid* aSolid = new G4SubtractionSolid("firstTry",
								aRingTub,
								notch,
								notch1Rotat,
								notch1Locat);
	      
	    CLHEP::HepRotation* notch2Rotat = new CLHEP::HepRotation(CLHEP::HepRotation::IDENTITY);
	    notch2Rotat->rotateZ(120.0*CLHEP::degree);

	    CLHEP::Hep3Vector notch2Locat((rin+22.0)*CLHEP::mm*sin(120.0*CLHEP::degree),
					  (rin+22.0)*CLHEP::mm*cos(120.0*CLHEP::degree),
					  0.0 );

	    G4SubtractionSolid* bSolid = new G4SubtractionSolid("secondTry",
							       aSolid,
							       notch,
							       notch2Rotat,
							       notch2Locat);

	    CLHEP::HepRotation* notch3Rotat = new CLHEP::HepRotation(CLHEP::HepRotation::IDENTITY);
	    notch3Rotat->rotateZ(240.0*CLHEP::degree);

	    CLHEP::Hep3Vector notch3Locat((rin+22.0)*CLHEP::mm*sin(240.0*CLHEP::degree),
					  (rin+22.0)*CLHEP::mm*cos(240.0*CLHEP::degree),
					  0.0 );

	    G4SubtractionSolid* cSolid = new G4SubtractionSolid(ring2Info.name,
							       bSolid,
							       notch,
							       notch3Rotat,
							       notch3Locat);

	    ring2Info.solid = cSolid;

	    finishNesting( ring2Info,
			   oPAsupportMaterial,
			   0,
			   location - parent1Info.centerInMu2e(),
			   parent1Info.logical,
			   0,
			   pabsIsVisible,
			   G4Color::Blue(),
			   pabsIsSolid,
			   forceAuxEdgeVisible,
			   placePV,
			   doSurfaceCheck);
	    if ( hasExtra ) {
	      nestTubs( myName2.str(),
			TubsParams( rou, rou+xRad, hl, (270-dPhiX/2.)*CLHEP::deg,
				    dPhiX*CLHEP::deg),
			oPAsupportMaterial,
			0,
			location - parent1Info.centerInMu2e(),
			parent1Info,
			0,
			pabsIsVisible,
			G4Color::Blue(),
			pabsIsSolid,
			forceAuxEdgeVisible,
			placePV,
			doSurfaceCheck );	      
	    } 
	  } else {
	    if ( zm < pabs1EndInMu2eZ ) {
	      nestTubs( myName.str(),
			TubsParams( rin, rou, hl ),
			oPAsupportMaterial,
			0,
			location - parent1Info.centerInMu2e(),
			parent1Info,
			0,
			pabsIsVisible,
			G4Color::Blue(),
			pabsIsSolid,
			forceAuxEdgeVisible,
			placePV,
			doSurfaceCheck );
	      if ( hasExtra ) {
		nestTubs( myName2.str(),
			  TubsParams( rou, rou+xRad, hl, (270-dPhiX/2.)*CLHEP::deg,
				      dPhiX*CLHEP::deg),
			  oPAsupportMaterial,
			  0,
			  location - parent1Info.centerInMu2e(),
			  parent1Info,
			  0,
			  pabsIsVisible,
			  G4Color::Blue(),
			  pabsIsSolid,
			  forceAuxEdgeVisible,
			  placePV,
			  doSurfaceCheck );	      
	      } // end of adding extra for ds2Vac area support ring
	      if ( iSup == 2 ) {
		// Add slats to support ST
		double sHigh = pabs->slatHeight();
		double sWide = pabs->slatWidth();
		double sLong = pabs->slatLength();  // full length, not half...
		if ( sHigh > 1e-02 && sWide > 1e-02 && sLong > 1e-02 ) {
		  double boxPars[3] = {sWide/2.0,sHigh/2.0,sLong/2.0};
		  CLHEP::Hep3Vector slat1Loc(-3904.0,rin+sHigh/2.0,zm-hl-sLong/2.0);
		  nestBox ( "STsupportSlat1",
			    boxPars,
			    oPAsupportMaterial,
			    0,
			    slat1Loc - parent1Info.centerInMu2e(),
			    parent1Info,
			    0,
			    pabsIsVisible,
			    G4Colour::Blue(),
			    pabsIsSolid,
			    forceAuxEdgeVisible,
			    placePV,
			    doSurfaceCheck);

		  CLHEP::HepRotation* slat2Rotat = new CLHEP::HepRotation(CLHEP::HepRotation::IDENTITY);
		  slat2Rotat->rotateZ(120.0*CLHEP::degree);

		  CLHEP::Hep3Vector slat2Loc(-3904.0 + (rin+sHigh/2.0)*sin(120.0*CLHEP::degree),
					     (rin+sHigh/2.0)*cos(120.0*CLHEP::degree),
					     zm - hl - sLong/2.0 );
		  nestBox ( "STsupportSlat2",
			    boxPars,
			    oPAsupportMaterial,
			    slat2Rotat,
			    slat2Loc - parent1Info.centerInMu2e(),
			    parent1Info,
			    0,
			    pabsIsVisible,
			    G4Colour::Blue(),
			    pabsIsSolid,
			    forceAuxEdgeVisible,
			    placePV,
			    doSurfaceCheck);

		  CLHEP::HepRotation* slat3Rotat = new CLHEP::HepRotation(CLHEP::HepRotation::IDENTITY);
		  slat3Rotat->rotateZ(240.0*CLHEP::degree);

		  CLHEP::Hep3Vector slat3Loc(-3904.0 + (rin+sHigh/2.0)*sin(240.0*CLHEP::degree),
					     (rin+sHigh/2.0)*cos(240.0*CLHEP::degree),
					     zm - hl - sLong/2.0 );
		  nestBox ( "STsupportSlat3",
			    boxPars,
			    oPAsupportMaterial,
			    slat3Rotat,
			    slat3Loc - parent1Info.centerInMu2e(),
			    parent1Info,
			    0,
			    pabsIsVisible,
			    G4Colour::Blue(),
			    pabsIsSolid,
			    forceAuxEdgeVisible,
			    placePV,
			    doSurfaceCheck);

		} // end if dimensions > 0 

	      } // end if iSub == 2  (adding slats for ST support)
	    } else {
	      nestTubs( myName.str(),
			TubsParams( rin, rou, hl ),
			oPAsupportMaterial,
			0,
			location - parent2Info.centerInMu2e(),
			parent2Info,
			0,
			pabsIsVisible,
			G4Color::Blue(),
			pabsIsSolid,
			forceAuxEdgeVisible,
			placePV,
			doSurfaceCheck );
	    } // end of if for placing in DS2Vac or DS3Vac
	  } // end of if else for support 2 (iSup == 1)
	}// end for loop over OPA supports

      } // end if nSupports > 0 for OPA


      //***************************
      // Now build the Degrader if requested
      //***************************

      if ( pabs->degraderBuild() ) {
	CLHEP::HepRotation* degraderRot = new CLHEP::HepRotation(CLHEP::HepRotation::IDENTITY);
	degraderRot->rotateZ(pabs->degraderRotation()*CLHEP::degree);
	
	// Make Frame
	std::vector<double> frameDims = pabs->degraderFrameDims();
	std::vector<double> filterDims = pabs->degraderFilterDims();
	std::vector<double> counterDims = pabs->degraderCounterwtDims();
	std::vector<double> rodDims = pabs->degraderRodDims();
	CLHEP::Hep3Vector locationInMu2e (-3904.0-frameDims.at(3), 0.0, 
					  pabs->degraderZ0() 
					  + counterDims.at(2) );

	// Make mother volume for degrader
	std::string motherName("Degrader");
	VolumeInfo degraderMother ( motherName,
				    locationInMu2e - parent1Info.centerInMu2e(), parent1Info.centerInWorld);

	// Make box for degrader mother volume.
	G4Box* motherBox = new G4Box ( "degraderOutline", 
	 			       frameDims.at(1)+frameDims.at(3),
	 			       frameDims.at(1),
	 			       counterDims.at(2) );
	degraderMother.solid = motherBox;

	// Now put degraderMother in DS2Vacuum
	finishNesting ( degraderMother,
			findMaterialOrThrow("DSVacuum"),
			degraderRot, degraderMother.centerInParent,
			parent1Info.logical, 0, false, G4Colour::Red(),
			false,
			forceAuxEdgeVisible,
			placePV,
			doSurfaceCheck );

	// Start cobbling pieces together by putting frame in mother
	CLHEP::Hep3Vector trans1(frameDims.at(3),0,frameDims.at(2) -
	 			 counterDims.at(2));
	nestTubs("degraderFrame",
	 	 TubsParams(frameDims.at(0),frameDims.at(1),frameDims.at(2)),
		 findMaterialOrThrow(pabs->degraderFrameMaterial()),
	 	 0, trans1, degraderMother,
	 	 0, pabsIsVisible, G4Color::Red(),
	 	 pabsIsSolid,
	 	 forceAuxEdgeVisible,
		 placePV,
	 	 doSurfaceCheck );
	// G4Tubs * frame = new G4Tubs( "degraderFrameTubs",
	// 			     frameDims.at(0),
	// 			     frameDims.at(1),
	// 			     frameDims.at(2) );

	// Now put filter in mother
	CLHEP::Hep3Vector trans1b(frameDims.at(3),0,2.0 * frameDims.at(2) +
	 			  filterDims.at(2) - counterDims.at(2));
	nestTubs("degraderFilter",
	 	 TubsParams(filterDims.at(0),filterDims.at(1),filterDims.at(2)),
	 	 findMaterialOrThrow(pabs->degraderFilterMaterial()),
	 	 0, trans1b, degraderMother,
	 	 0, pabsIsVisible, G4Color::Red(),
	 	 pabsIsSolid,
		 forceAuxEdgeVisible,
		 placePV,
		 doSurfaceCheck );
	// G4Tubs * filter = new G4Tubs( "degraderFilterTubs",
	// 			      filterDims.at(0),
	// 			      filterDims.at(1),
	// 			      filterDims.at(2) );


	// Now put counterweight in mother
	CLHEP::Hep3Vector trans2( -counterDims.at(3), 0.0, 0.0);
	nestTubs("degraderCounterweight",
		 TubsParams(counterDims.at(0),counterDims.at(1),counterDims.at(2)),
		 findMaterialOrThrow(pabs->degraderCountwtMaterial()),
		 0, trans2, degraderMother,
		 0, pabsIsVisible, G4Color::Red(),
		 pabsIsSolid,
		 forceAuxEdgeVisible,
		 placePV,
		 doSurfaceCheck );

	// Creat rod rep
	double lenRod = frameDims.at(3) + counterDims.at(3) - frameDims.at(1) 
	  - counterDims.at(1) - 0.1;
	std::vector<double> lwhs = {lenRod/2.0,rodDims.at(0)/2.0,rodDims.at(1)/2.0};
	// Now put rod in mother volume
	CLHEP::Hep3Vector trans3( (frameDims.at(3)-counterDims.at(3)
				   + counterDims.at(1) - frameDims.at(1))/2.0, 
				  0.0, rodDims.at(1)/2.0 - counterDims.at(2));
	nestBox ("degraderRod",
		 lwhs,
		 findMaterialOrThrow(pabs->degraderRodMaterial()),
		 0, trans3, degraderMother,
		 0, pabsIsVisible, G4Color::Red(),
		 pabsIsSolid,
		 forceAuxEdgeVisible,
		 placePV,
		 doSurfaceCheck );



      }// end if degraderBuild

      
      if ( pabs->buildSupports() ) {

        AntiLeakRegistry& reg = art::ServiceHandle<G4Helper>()->antiLeakRegistry();

        const InnerProtonAbsSupport* ipaSup = pabs->getIPAsupport();
	const int ipa_version = _config.getInt("protonabsorber.version", 1); // also need to know the version number in this file
	const double wireRotation = _config.getDouble("protonabsorber.ipa.wireRotationToVertical", 45); // will be used for v2 only

	// Build the end rings
	// this will only happen for v2 of the IPA because the deafult number of end rings is 0
        for ( std::size_t iEndRing(0); iEndRing < ipaSup->nEndRings(); iEndRing++) {

            Tube endRing = ipaSup->getEndRing( iEndRing );

            ostringstream endRingName ; endRingName << "IPAsupport_endring" << iEndRing;

            nestTubs( endRingName.str() ,
                      endRing.getTubsParams(),
                      findMaterialOrThrow( endRing.materialName() ),
                      0,
                      endRing.originInMu2e()-parent1Info.centerInMu2e(),
                      parent1Info,
                      0,
                      pabsIsVisible,
                      G4Color::Red(),
                      pabsIsSolid,
                      forceAuxEdgeVisible,
                      placePV,
                      doSurfaceCheck );

	} // end ring loop


        int iS(0);
        for ( std::size_t iSet(0); iSet < ipaSup->nSets(); iSet++) {
	  ++iS; // increment here so we have set 1 and set 2, one for each end of the IPA
          const double pabs1rOut0 = pabs->part(0).outerRadiusAtStart();
          const double pabs1rOut1 = pabs->part(0).outerRadiusAtEnd();
          const double pabs1len   = pabs->part(0).halfLength() *2.;

          const double zstartOfIPA = pabs->part(0).center().z()-pabs->part(0).halfLength();

          int iW(0);
          for ( std::size_t iWire(0); iWire < ipaSup->nWiresPerSet() ; iWire++ ) {

            Tube supportWire = ipaSup->getWire( iSet, iWire );

            ostringstream wirename ; wirename << "IPAsupport_set" << iS << "_wire" << ++iW ;

            const double rStartOfWire = pabs1rOut0+(supportWire.originInMu2e().z()-zstartOfIPA)/pabs1len*(pabs1rOut1-pabs1rOut0);

            CLHEP::Hep3Vector additionalOffset ( (supportWire.halfLength()+0.005+rStartOfWire) * std::cos(iWire * 360.*CLHEP::deg / ipaSup->nWiresPerSet() ),
                                                 (supportWire.halfLength()+0.005+rStartOfWire) * std::sin(iWire * 360.*CLHEP::deg / ipaSup->nWiresPerSet() ), 
                                                 0 );

            // Now get appropriate rotation angles
            G4RotationMatrix* supportRot = reg.add(G4RotationMatrix());

	    if (ipa_version == 1) {
	      supportRot->rotateY(-M_PI/2.);
	      supportRot->rotateZ(-iW*360.*CLHEP::deg / ipaSup->nWiresPerSet() );
	    }
	    else if (ipa_version == 2) {
	      CLHEP::Hep3Vector rotationAxis(0, 1, 0); // start off with rotating around the y-axis
	      rotationAxis.rotateZ((iW-1)*360.*CLHEP::deg / ipaSup->nWiresPerSet()); // each wire wants to be rotated arounf a slightly different axis
	      supportRot->rotate(wireRotation*CLHEP::deg, rotationAxis);
	    
	      CLHEP::Hep3Vector extraZOffset(0, 0, supportWire.halfLength() * std::cos(wireRotation*CLHEP::deg)); // because of the rotation from the vertical
	      CLHEP::Hep3Vector extraROffset(std::cos(iWire * 360.*CLHEP::deg / ipaSup->nWiresPerSet()) * supportWire.halfLength()*(std::cos((90-wireRotation*CLHEP::deg)))*std::tan(wireRotation*CLHEP::deg), 
					   std::sin(iWire * 360.*CLHEP::deg / ipaSup->nWiresPerSet()) * supportWire.halfLength()*(std::cos((90-wireRotation*CLHEP::deg)))*std::tan(wireRotation*CLHEP::deg), 
					   0);
	      additionalOffset -= 0.90*extraROffset; 

	      if (supportWire.originInMu2e().z() > pabs->part(0).center().z()) { // if the wires are further away from the target...
		// ...we need to rotate them again
		supportRot->rotate(90.*CLHEP::deg, rotationAxis);
		
		// and move them
		additionalOffset += extraZOffset;
	      }
	      else {
		additionalOffset -= extraZOffset;
	      }
	    }

            nestTubs( wirename.str() ,
                      supportWire.getTubsParams(),
                      findMaterialOrThrow( supportWire.materialName() ),
                      supportRot,
                      supportWire.originInMu2e()+additionalOffset-parent1Info.centerInMu2e(),
                      parent1Info,
                      0,
                      pabsIsVisible,
                      G4Color::Red(),
                      pabsIsSolid,
                      forceAuxEdgeVisible,
                      placePV,
                      doSurfaceCheck );

          } // wire loop
        } // set loop

      } // build ipa supports
    }
  }
}

