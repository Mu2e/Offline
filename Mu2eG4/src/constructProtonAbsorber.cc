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
#include "Mu2eG4/inc/findMaterialOrThrow.hh"
#include "Mu2eG4/inc/MaterialFinder.hh"
#include "Mu2eG4/inc/nestCons.hh"
#include "Mu2eG4/inc/nestTubs.hh"
#include "Mu2eG4/inc/SensitiveDetectorName.hh"
#include "Mu2eG4/inc/HelicalProtonAbsorber.hh"
#include "MECOStyleProtonAbsorberGeom/inc/MECOStyleProtonAbsorber.hh"
#include "DetectorSolenoidGeom/inc/DetectorSolenoid.hh"
#include "Mu2eG4/inc/checkForOverlaps.hh"

// G4 includes
#include "G4Material.hh"
#include "G4Color.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4BooleanSolid.hh"
#include "G4VPhysicalVolume.hh"
#include "G4SDManager.hh"

using namespace std;

namespace mu2e {

  void constructProtonAbsorber( SimpleConfig const & _config
                                ){

    if( !_config.getBool("hasProtonAbsorber", true) ) return;
    
    int  const verbosityLevel           = _config.getInt("protonabsorber.verbosityLevel", 0);
    
    G4VSensitiveDetector *paSD = G4SDManager::GetSDMpointer()->FindSensitiveDetector (SensitiveDetectorName::ProtonAbsorber());
    
    // Access to the G4HelperService.
    G4Helper* _helper = &(*(art::ServiceHandle<G4Helper>()));
    
    VolumeInfo const & parent1Info  = _helper->locateVolInfo("DS2Vacuum");
    VolumeInfo const & parent2Info  = _helper->locateVolInfo("DS3Vacuum");
    
    // Fetch DS geometry
    GeomHandle<DetectorSolenoid> ds;

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
      
      bool addSD   = _config.getBool("protonabsorber.saveStepPnts",false);
      
      HelicalProtonAbsorber* hpabs = 
        new HelicalProtonAbsorber( ds2HalfLen-helPabsLength,
                                   helPabsLength*lengthScaleFact,/*ScaleFact is a trick FIXME */
                                   innerRadii, outerRadii, innerPhis, outerPhis,
                                   helPabsThickness, /*numOfTurns,*/
                                   _config.getInt("protonabsorber.NumOfVanes"), 
                                   pabsMaterial, parent1Info.logical, (addSD) ? paSD : 0x0 );
      
      if ( _config.getBool("g4.doSurfaceCheck",false) ) {
        checkForOverlaps( hpabs->GetPhys(), _config, verbosityLevel>0);
      }
      
      if ( _config.getBool("protonabsorber.visible",true) ) {
        AntiLeakRegistry & reg = _helper->antiLeakRegistry();
        hpabs->SetVisibility( _config.getBool("protonabsorber.solid",true),
                              _config.getBool("g4.forceAuxEdgeVisible",false),
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
      bool pabsIsVisible = _config.getBool("protonabsorber.visible",true);
      bool pabsIsSolid   = _config.getBool("protonabsorber.solid",true);

      bool forceAuxEdgeVisible = _config.getBool("g4.forceAuxEdgeVisible",false);
      bool doSurfaceCheck      = _config.getBool("g4.doSurfaceCheck",false);
      bool const placePV       = true;

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
      if(paSD) protonabs1Info.logical->SetSensitiveDetector (paSD);

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
      if(paSD) protonabs2Info.logical->SetSensitiveDetector (paSD);
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

      bool pabsIsVisible = _config.getBool("protonabsorber.visible",true);
      bool pabsIsSolid   = _config.getBool("protonabsorber.solid",true);
  
      bool forceAuxEdgeVisible = _config.getBool("g4.forceAuxEdgeVisible",false);
      bool doSurfaceCheck      = _config.getBool("g4.doSurfaceCheck",false);
      bool const placePV       = true;
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
        if(paSD) protonabs1Info.logical->SetSensitiveDetector (paSD);
  
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
        if(paSD) protonabs2Info.logical->SetSensitiveDetector (paSD);
  
  
  
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
  

      // outer proton absorber
      if (pabs->isAvailable(2)) {
        double pabs3rIn0  = pabs->part(2).innerRadiusAtStart();
        double pabs3rOut0 = pabs->part(2).outerRadiusAtStart();
        double pabs3rIn1  = pabs->part(2).innerRadiusAtEnd();
        double pabs3rOut1 = pabs->part(2).outerRadiusAtEnd();
        double pabs3len   = pabs->part(2).halfLength() * 2.;
        G4Material* pabs3Material = materialFinder.get("protonabsorber.outerPAMaterialName");
  
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
        if(paSD) protonabs3Info.logical->SetSensitiveDetector (paSD);
  
/*        mf::LogInfo log("GEOM");
        log << "Constructing Proton Absorber -- \n";
        log << "Proton Abs Offset in DS2:  " << pabs1Offset <<"\n";
        log << "rIn,  rOut (-z): "<< pabs1rIn0 <<"  "<< pabs1rOut0<<"  ";
        log << "rIn,  rOut (+z): "<< pabs1rIn1 <<"  "<< pabs1rOut1<<"  ";
        log << "halflength: "<< pabs1len*0.5 <<"\n";
        log << "Proton Abs Offset in DS3:  " << pabs3Offset <<"\n";
        log << "rIn,  rOut (-z): "<< pabs3rIn0 <<"  "<< pabs3rOut0<<"  ";
        log << "rIn,  rOut (+z): "<< pabs3rIn1 <<"  "<< pabs3rOut1<<"  ";
        log << "halflength: "<< pabs3len*0.5 <<"\n";

        if ( verbosityLevel > 0) {
          double pzhl   = static_cast<G4Cons*>(protonabs2Info.solid)->GetZHalfLength();
          double pabs3Z = protonabs2Info.centerInMu2e()[CLHEP::Hep3Vector::Z];
          cout << __func__ << " " << protonabs2Info.name << " Z offset in Mu2e    : " <<
            pabs3Z << endl;
          cout << __func__ << " " << protonabs2Info.name << " Z extent in Mu2e    : " <<
            pabs3Z - pzhl  << ", " <<  pabs3Z + pzhl  << endl;
  
          // we also check how the offsets are handled
  
          cout << __func__ << " " << protonabs2Info.name << " local input offset in G4                  : " <<
            pabs3Offset << endl;
          cout << __func__ << " " << protonabs2Info.name << " local GetTranslation()       offset in G4 : " <<
            protonabs2Info.physical->GetTranslation() << endl; // const &
        }
  */ //TODO turn on logging
  
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
        G4Material* pabs4Material = materialFinder.get("protonabsorber.outerPAMaterialName");
  
        double pabs4Param[7] = { pabs4rIn0, pabs4rOut0, pabs4rIn1, pabs4rOut1, pabs4len*0.5,
                                 0.0, 360.0*CLHEP::degree };
  
        VolumeInfo protonabs4Info = nestCons( "protonabs4",
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
        if(paSD) protonabs4Info.logical->SetSensitiveDetector (paSD);
  
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
	    
	      CLHEP::Hep3Vector extraZOffset(0, 0, supportWire.halfLength() * std::cos(wireRotation*CLHEP::deg) + 1.0); // because of the rotation from the vertical (1.0 to avoid contact with the end-rings)
	      CLHEP::Hep3Vector extraROffset(std::cos(iWire * 360.*CLHEP::deg / ipaSup->nWiresPerSet()) * supportWire.halfLength()*(std::cos((90-wireRotation*CLHEP::deg)))*std::tan(wireRotation*CLHEP::deg), 
					   std::sin(iWire * 360.*CLHEP::deg / ipaSup->nWiresPerSet()) * supportWire.halfLength()*(std::cos((90-wireRotation*CLHEP::deg)))*std::tan(wireRotation*CLHEP::deg), 
					   0);
	      additionalOffset -= 0.97*extraROffset; // try to avoid moving back and overlapping with the end ring (0.97 to avoid contact with the endrings)

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

