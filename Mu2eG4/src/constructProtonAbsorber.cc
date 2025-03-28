//
// Free function to create Proton Absorber
//
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
#include <cmath>

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"

// Mu2e includes.
#include "Offline/Mu2eG4/inc/constructProtonAbsorber.hh"
#include "Offline/Mu2eG4Helper/inc/VolumeInfo.hh"
#include "Offline/GeometryService/inc/GeometryService.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/VirtualDetector.hh"
#include "Offline/StoppingTargetGeom/inc/StoppingTarget.hh"
#include "Offline/Mu2eG4Helper/inc/Mu2eG4Helper.hh"
#include "Offline/GeomPrimitives/inc/TubsParams.hh"
#include "Offline/GeomPrimitives/inc/PolyhedraParams.hh"
#include "Offline/Mu2eG4/inc/findMaterialOrThrow.hh"
#include "Offline/Mu2eG4/inc/MaterialFinder.hh"
#include "Offline/Mu2eG4/inc/nestBox.hh"
#include "Offline/Mu2eG4/inc/nestCons.hh"
#include "Offline/Mu2eG4/inc/nestTubs.hh"
#include "Offline/Mu2eG4/inc/nestPolyhedra.hh"
#include "Offline/Mu2eG4/inc/finishNesting.hh"
#include "Offline/Mu2eG4/inc/HelicalProtonAbsorber.hh"
#include "Offline/BeamlineGeom/inc/ProtonAbsorber.hh"
#include "Offline/DetectorSolenoidGeom/inc/DetectorSolenoid.hh"
#include "Offline/Mu2eG4/inc/checkForOverlaps.hh"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Units/SystemOfUnits.h"

// G4 includes
#include "Geant4/G4Material.hh"
#include "Geant4/G4Color.hh"
#include "Geant4/G4Box.hh"
#include "Geant4/G4Cons.hh"
#include "Geant4/G4Tubs.hh"
#include "Geant4/G4BooleanSolid.hh"
#include "Geant4/G4UnionSolid.hh"
#include "Geant4/G4VPhysicalVolume.hh"
#include "Geant4/G4SubtractionSolid.hh"
#include "Geant4/G4SDManager.hh"

using namespace std;

namespace mu2e {

  void constructProtonAbsorber( const SimpleConfig& _config ){

    if( !_config.getBool("hasProtonAbsorber", true) ) return;

    int  const verbosityLevel           = _config.getInt("protonabsorber.verbosityLevel", 0);

    // Access to the Mu2eG4HelperService.
    Mu2eG4Helper* _helper = &(*(art::ServiceHandle<Mu2eG4Helper>()));
    AntiLeakRegistry & reg = _helper->antiLeakRegistry();

    const bool inGaragePosition = _config.getBool("inGaragePosition",false); //offset detector train elements for extracted position
    const bool OPA_IPA_ST_Extracted = (inGaragePosition) ? _config.getBool("garage.extractOPA_IPA_ST") : false;
    const double zOffGarage = (inGaragePosition && OPA_IPA_ST_Extracted) ? _config.getDouble("garage.zOffset") : 0.;
    const CLHEP::Hep3Vector relPosFake(0.,0., zOffGarage);

    std::string theDS2("DS2Vacuum");
    std::string theDS3("DS3Vacuum");
    if (inGaragePosition && OPA_IPA_ST_Extracted) {
      theDS2 = "garageFakeDS2Vacuum";
      theDS3 = "garageFakeDS3Vacuum";
    }
    VolumeInfo const & parent1Info  = _helper->locateVolInfo(theDS2);
    VolumeInfo const & parent2Info  = _helper->locateVolInfo(theDS3);

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
      double pabsStartInMu2eZ = target->centerInMu2e().z() + 0.5*target->cylinderLength() + 2.*vdHL + zOffGarage;

      // Need to split it at the DS2/DS3 boundary
      double pabs1EndInMu2eZ = ds->vac_zLocDs23Split() + zOffGarage;

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
      double pabs2ZOffset = (pabs2len*0.5) +  ds->vac_zLocDs23Split() + zOffGarage;

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

      GeomHandle<ProtonAbsorber> pabs;

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
      if (pabs->isAvailable(0)) pabs1ZOffset = CLHEP::mm * pabs->part(0).center().z() - ds2zcenter + relPosFake.z();
      if (pabs->isAvailable(1)) pabs2ZOffset = CLHEP::mm * pabs->part(1).center().z() - ds3zcenter + relPosFake.z();
      if (pabs->isAvailable(2)) pabs3ZOffset = CLHEP::mm * pabs->part(2).center().z() - ds2zcenter + relPosFake.z();
      if (pabs->isAvailable(3)) pabs4ZOffset = CLHEP::mm * pabs->part(3).center().z() - ds3zcenter + relPosFake.z();

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

        if(verbosityLevel > 0) {
          cout << __func__ << " protonabs3 center in Mu2e: " << pabs3Offset - parent1Info.centerInMu2e() << endl
               << __func__ << " protonabs3 {r1_0, r2_0}: {" << pabs3rIn0 << ", "
               << pabs3rOut0 << "}" << endl
               << __func__ << " protonabs3 {r1_1, r2_1}: {" << pabs3rIn1 << ", "
               << pabs3rOut1 << "}" << endl
               << __func__ << " protonabs3 length: " << pabs3len << endl;
        }
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
          if(verbosityLevel > 0) {
            cout << __func__ << "protonabs3 implemented with " << pabs3nS << " slats" << endl;
          }
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

          CLHEP::HepRotation* slot1Rotat = reg.add(CLHEP::HepRotation(CLHEP::HepRotation::IDENTITY));
          slot1Rotat->rotateZ(120.0*CLHEP::degree);

          CLHEP::Hep3Vector slot1spot(300.0*CLHEP::mm*sin(120.0*CLHEP::degree),
                                      300.0*CLHEP::mm*cos(120.0*CLHEP::degree),
                                      pabs3SlotOffset*CLHEP::mm);

          G4SubtractionSolid* aSolid = new G4SubtractionSolid( "FirstStep",
                                                               anOPApolyhedron,
                                                               slot,
                                                               slot1Rotat,
                                                               slot1spot );

          CLHEP::HepRotation* slot2Rotat (nullptr);

          CLHEP::Hep3Vector slot2spot(0., 300.0*CLHEP::mm, pabs3SlotOffset*CLHEP::mm );

          G4SubtractionSolid* bSolid = new G4SubtractionSolid( "Step2",
                                                               aSolid,
                                                               slot,
                                                               slot2Rotat,
                                                               slot2spot );

          CLHEP::HepRotation* slot3Rotat = reg.add(CLHEP::HepRotation(CLHEP::HepRotation::IDENTITY));
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
        if ( verbosityLevel > 0 ) cout << __func__ << " outer protonabs3 disabled" << endl;
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
        if ( verbosityLevel > 0 ) cout << __func__ << " outer protonabs4 disabled" << endl;
      }

      //***************************
      // Now build the support rings for the OPA
      //***************************


      if ( pabs->oPAnSupports() > 0 ) {
        double pabs1EndInMu2eZ = ds->vac_zLocDs23Split();
        int opaVersion = pabs->outerPAVersion();
        // Get stopping target parameters for stopping target support slat parameters
        GeomHandle<StoppingTarget> target;

        for (int iSup = 0; iSup < pabs->oPAnSupports(); iSup++ ) {
          G4Material* oPAsupportMaterial = findMaterialOrThrow( pabs->oPAsupportMaterials().at(iSup) );
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

          CLHEP::Hep3Vector location(parent1Info.centerInMu2e()); //inherit x and y position from DS
          location.setZ(zm);

          ostringstream myName;
          myName << "OPAsupport_" << iSup+1;
          ostringstream myName2;
          myName2 << "OPAsupportExtra_" << iSup+1;

          if(verbosityLevel > 0) {
            cout << __func__ << ": Constructing " << myName.str() << " with center in Mu2e " << location << " and {r1, r2} = {"
                 << rin << ", " << rou << "}" << endl;
          }
          if ( iSup == 1 ) {
            // Make notches in this one for the Stopping target support slats
            VolumeInfo ring2Info(myName.str(),
                                 location - parent1Info.centerInMu2e() + relPosFake,
                                 parent1Info.centerInWorld);

            G4Tubs* aRingTub = new G4Tubs("TheRing",
                                          rin, rou, hl,
                                          0.0, 360.0*CLHEP::degree);
            const double notchWidth  = _config.getDouble("protonabsorber.oPASupportNotchWidth",
                                                         2.*35.);
            const double notchHeight = _config.getDouble("protonabsorber.oPASupportNotchHeight",
                                                         2.*25.);

            //extend the box so that it notches through at either end, not just the center
            const double notchExtension = (opaVersion > 2) ?
              (rin - sqrt(rin*rin-notchWidth*notchWidth)) : 0.;

            G4Box* notch = new G4Box("notch", notchWidth/2.*CLHEP::mm,
                                     (notchHeight+notchExtension)/2.*CLHEP::mm, hl*1.1 );

            CLHEP::HepRotation* notch1Rotat (nullptr);
            double notchOffset = notchHeight/2.;
            notchOffset -= notchExtension/2.; //so the extension is only on the bottom
            if(opaVersion < 3) notchOffset = 22.; //adding for backwards compatibility
            CLHEP::Hep3Vector notch1Locat(0,(rin+notchOffset)*CLHEP::mm,0);

            G4SubtractionSolid* aSolid = new G4SubtractionSolid("firstTry",
                                                                aRingTub,
                                                                notch,
                                                                notch1Rotat,
                                                                notch1Locat);

            CLHEP::HepRotation* notch2Rotat = reg.add(CLHEP::HepRotation(CLHEP::HepRotation::IDENTITY));
            notch2Rotat->rotateZ(120.0*CLHEP::degree);

            CLHEP::Hep3Vector notch2Locat((rin+notchOffset)*CLHEP::mm*sin(120.0*CLHEP::degree),
                                          (rin+notchOffset)*CLHEP::mm*cos(120.0*CLHEP::degree),
                                          0.0 );

            G4SubtractionSolid* bSolid = new G4SubtractionSolid("secondTry",
                                                                aSolid,
                                                                notch,
                                                                notch2Rotat,
                                                                notch2Locat);

            CLHEP::HepRotation* notch3Rotat = reg.add(CLHEP::HepRotation(CLHEP::HepRotation::IDENTITY));
            notch3Rotat->rotateZ(240.0*CLHEP::degree);

            CLHEP::Hep3Vector notch3Locat((rin+notchOffset)*CLHEP::mm*sin(240.0*CLHEP::degree),
                                          (rin+notchOffset)*CLHEP::mm*cos(240.0*CLHEP::degree),
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
                           location - parent1Info.centerInMu2e() + relPosFake,
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
                        location - parent1Info.centerInMu2e() + relPosFake,
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
              auto vi = nestTubs( myName.str(),
                                  TubsParams( rin, rou, hl ),
                                  oPAsupportMaterial,
                                  0,
                                  location - parent1Info.centerInMu2e() + relPosFake,
                                  parent1Info,
                                  0,
                                  pabsIsVisible,
                                  G4Color::Blue(),
                                  pabsIsSolid,
                                  forceAuxEdgeVisible,
                                  placePV,
                                  doSurfaceCheck );
              if(verbosityLevel > 0) cout << " Support final center = " << vi.centerInMu2e() << endl;
              if ( hasExtra ) {
                nestTubs( myName2.str(),
                          TubsParams( rou, rou+xRad, hl, (270-dPhiX/2.)*CLHEP::deg,
                                      dPhiX*CLHEP::deg),
                          oPAsupportMaterial,
                          0,
                          location - parent1Info.centerInMu2e() + relPosFake,
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
                if(verbosityLevel > 0)
                  cout << "Constructing stopping target supports at the OPA:" << endl;
                // Add slats to support ST
                for(int iSlat = 0; iSlat < pabs->nOPASupportSlats(); ++iSlat) {
                  if(verbosityLevel > 0)
                    cout << " iSlat = " << iSlat << endl;
                  const double slatAngle = pabs->slatAngles().at(iSlat)*CLHEP::degree; //taken in as degrees

                  CLHEP::HepRotation* slatRotat = reg.add(CLHEP::HepRotation(CLHEP::HepRotation::IDENTITY));
                  slatRotat->rotateZ(slatAngle);
                  const int slatIndex = pabs->slatTypes().at(iSlat); //get the type of slat
                  const double sHigh = pabs->slatHeights().at(slatIndex);
                  const double sWide = pabs->slatWidths().at(slatIndex);
                  const double sLong = pabs->slatLengths().at(slatIndex);
                  const double boxPars[3] = {sWide/2.0,sHigh/2.0,sLong/2.0};
                  double rCenter = rin - sHigh/2.; //radial distance from center axis to box center
                  //remove overlap
                  if(opaVersion > 2)
                    rCenter -= rin - sqrt(rin*rin - sWide/2.*sWide/2.) + 1.;
                  else
                    rCenter = rin + sHigh/2.0;
                  CLHEP::Hep3Vector slatLoc(rCenter*sin(slatAngle),
                                            rCenter*cos(slatAngle),
                                            0.);

                  if(opaVersion > 2) {
                    location.setZ(target->centerInMu2e().z()); //newer versions set center at target center
                    slatLoc += location; //add z offset
                    G4Material* slatMaterial = findMaterialOrThrow( pabs->slatMaterials().at(slatIndex) );
                    //create mother volume box
                    VolumeInfo slatMotherInfo = nestBox ( "STsupportSlatMother"+std::to_string(iSlat+1),
                                                          boxPars,
                                                          parent1Info.logical->GetMaterial(),
                                                          slatRotat,
                                                          slatLoc - parent1Info.centerInMu2e() + relPosFake,
                                                          parent1Info,
                                                          0,
                                                          pabsIsVisible,
                                                          G4Colour::Blue(),
                                                          pabsIsSolid,
                                                          forceAuxEdgeVisible,
                                                          placePV,
                                                          doSurfaceCheck);
                    const double sideThickness = pabs->slatSideThicknesses().at(slatIndex);
                    const double topThickness = pabs->slatTopThicknesses().at(slatIndex);
                    if(slatIndex == 0) {
                      //add bottom and sides of the slat
                      const double botPars[3] = {sWide/2-sideThickness-0.001, topThickness/2., sLong/2.};
                      nestBox ( "STsupportSlatBottom"+std::to_string(iSlat+1),
                                botPars,
                                slatMaterial,
                                0,
                                CLHEP::Hep3Vector(0., topThickness/2.0 - sHigh/2.0, 0.),
                                slatMotherInfo,
                                0,
                                pabsIsVisible,
                                G4Colour::Blue(),
                                pabsIsSolid,
                                forceAuxEdgeVisible,
                                placePV,
                                doSurfaceCheck);
                      const double sidePars[3] = {sideThickness/2., sHigh/2., sLong/2.};
                      nestBox ( "STsupportSlatLeft"+std::to_string(iSlat+1),
                                sidePars,
                                slatMaterial,
                                0,
                                CLHEP::Hep3Vector(sideThickness/2.0 - sWide/2.0, 0., 0.),
                                slatMotherInfo,
                                0,
                                pabsIsVisible,
                                G4Colour::Blue(),
                                pabsIsSolid,
                                forceAuxEdgeVisible,
                                placePV,
                                doSurfaceCheck);
                      nestBox ( "STsupportSlatRight"+std::to_string(iSlat+1),
                                sidePars,
                                slatMaterial,
                                0,
                                CLHEP::Hep3Vector(-sideThickness/2.0 + sWide/2.0, 0., 0.),
                                slatMotherInfo,
                                0,
                                pabsIsVisible,
                                G4Colour::Blue(),
                                pabsIsSolid,
                                forceAuxEdgeVisible,
                                placePV,
                                doSurfaceCheck);
                      const double fillH      = pabs->slatFillParameter1().at(slatIndex);
                      const double fillW      = pabs->slatFillParameter2().at(slatIndex);
                      const double fillL      = pabs->slatFillParameter3().at(slatIndex);
                      const double fillOffset = pabs->slatFillParameter4().at(slatIndex);
                      const double fillPars[] = {fillW/2., fillH/2., fillL/2.};
                      G4Material* slatFillMaterial = findMaterialOrThrow( pabs->slatFillMaterials().at(slatIndex) );
                      CLHEP::Hep3Vector fillLoc(fillOffset, -sHigh/2.+fillH/2.+topThickness + 0.001, 0.);
                      nestBox ( "STsupportSlatFill"+std::to_string(iSlat+1),
                                fillPars,
                                slatFillMaterial,
                                0,
                                fillLoc,
                                slatMotherInfo,
                                0,
                                pabsIsVisible,
                                G4Colour::Blue(),
                                pabsIsSolid,
                                forceAuxEdgeVisible,
                                placePV,
                                doSurfaceCheck);

                    } else if (slatIndex == 1) {
                      const double sidePars[3] = {sideThickness/2., sHigh/4., sLong/2.};
                      nestBox ( "STsupportSlatLeft"+std::to_string(iSlat+1),
                                sidePars,
                                slatMaterial,
                                0,
                                CLHEP::Hep3Vector(+sideThickness/2.0 - sWide/2.0, 0., 0.),
                                slatMotherInfo,
                                0,
                                pabsIsVisible,
                                G4Colour::Blue(),
                                pabsIsSolid,
                                forceAuxEdgeVisible,
                                placePV,
                                doSurfaceCheck);
                      const double fillR       = pabs->slatFillParameter1().at(slatIndex);
                      const double fillL       = pabs->slatFillParameter2().at(slatIndex);
                      const double fillOffset1 = pabs->slatFillParameter3().at(slatIndex);
                      const double fillOffset2 = pabs->slatFillParameter4().at(slatIndex);
                      TubsParams fillPars(0., fillR, fillL/2.);
                      G4Material* slatFillMaterial = findMaterialOrThrow( pabs->slatFillMaterials().at(slatIndex) );
                      CLHEP::Hep3Vector fillLoc1(-sWide/2.+fillR + sideThickness+ 0.001, fillOffset1, 0.);
                      nestTubs ( "STsupportSlatFill1_"+std::to_string(iSlat+1),
                                 fillPars,
                                 slatFillMaterial,
                                 0,
                                 fillLoc1,
                                 slatMotherInfo,
                                 0,
                                 pabsIsVisible,
                                 G4Colour::Blue(),
                                 pabsIsSolid,
                                 forceAuxEdgeVisible,
                                 placePV,
                                 doSurfaceCheck);
                      CLHEP::Hep3Vector fillLoc2(sWide/2.-fillR-0.001, fillOffset2, 0.);
                      nestTubs ( "STsupportSlatFill2_"+std::to_string(iSlat+1),
                                 fillPars,
                                 slatFillMaterial,
                                 0,
                                 fillLoc2,
                                 slatMotherInfo,
                                 0,
                                 pabsIsVisible,
                                 G4Colour::Blue(),
                                 pabsIsSolid,
                                 forceAuxEdgeVisible,
                                 placePV,
                                 doSurfaceCheck);
                      const double weightBoxH   = _config.getDouble("protonabsorber.oPASTWeightBox.height");
                      const double weightBoxW   = _config.getDouble("protonabsorber.oPASTWeightBox.width");
                      const double weightBoxL   = _config.getDouble("protonabsorber.oPASTWeightBox.length");
                      const double weightBoxR   = _config.getDouble("protonabsorber.oPASTWeightBox.radius");
                      const double weightBoxPhi = _config.getDouble("protonabsorber.oPASTWeightBox.angle")*CLHEP::degree;
                      const double weightW      = _config.getDouble("protonabsorber.oPASTWeight.width");
                      const double weightH      = _config.getDouble("protonabsorber.oPASTWeight.height");
                      G4Material* weightMaterial = findMaterialOrThrow(_config.getString("protonabsorber.oPASTWeight.material"));
                      CLHEP::Hep3Vector weightBoxLoc(weightBoxR*cos(weightBoxPhi), weightBoxR*sin(weightBoxPhi), 0.);
                      weightBoxLoc += location;
                      //intersects OPA support 1, 5, and 2, so have to break up the box into z regions
                      vector<int> suppIndices = {1, 5, 2};
                      vector<double> suppZStart;
                      vector<double> suppZEnd;
                      if(verbosityLevel > 0)
                        cout << " Weight box z range = {" << weightBoxLoc.z() - weightBoxL/2. << ","
                             << weightBoxLoc.z() + weightBoxL/2. << "}" <<  endl;
                      unsigned count = 0;
                      for(int index : suppIndices) {
                        suppZStart.push_back(pabs->oPAsupportZMidpoint().at(index) - pabs->oPAsupportHalflength().at(index));
                        suppZEnd  .push_back(pabs->oPAsupportZMidpoint().at(index) + pabs->oPAsupportHalflength().at(index));
                        if(verbosityLevel > 0)
                          cout << "  supp index = " << index << " z start = " << suppZStart[count]
                               << " z end = " << suppZEnd[count] << endl;
                        ++count;
                      }
                      //add a fake support at the end
                      suppZStart.push_back(weightBoxLoc.z()+weightBoxL/2.);
                      suppZEnd.push_back(weightBoxLoc.z()+weightBoxL/2.);
                      double currentZ = weightBoxLoc.z() - weightBoxL/2.;
                      for(unsigned step = 0; step < suppZStart.size(); ++step) {
                        double currL = suppZStart[step] - currentZ;
                        CLHEP::Hep3Vector weightMotherLoc(weightBoxLoc);
                        weightMotherLoc.setZ(suppZStart[step]-currL/2.);
                        const double weightMotherParams[] = {weightBoxW/2., weightBoxH/2., currL/2.-0.0005};
                        if(verbosityLevel > 0)
                          cout << " step: " << step << " Current Weight box z range = {" << weightMotherLoc.z() - currL/2. << ","
                               << weightMotherLoc.z() + currL/2. << "} total length = " << currL << endl
                               << " Weight box mother params: {" << weightMotherParams[0] << ", "
                               << weightMotherParams[1] << ", " << weightMotherParams[2] << "}\n";
                        VolumeInfo weightMotherInfo = nestBox ( "STsupportWeightMother"+std::to_string(step),
                                                                weightMotherParams,
                                                                parent1Info.logical->GetMaterial(),
                                                                0,
                                                                weightMotherLoc - parent1Info.centerInMu2e() + relPosFake,
                                                                parent1Info,
                                                                0,
                                                                pabsIsVisible,
                                                                G4Colour::Blue(),
                                                                pabsIsSolid,
                                                                forceAuxEdgeVisible,
                                                                placePV,
                                                                doSurfaceCheck);
                        double zCenter = 0.;
                        if(step == 0 || step == suppZStart.size() - 1) {//add front panel
                          const double frontPars[3] = {weightBoxW/2., weightBoxH/2., sideThickness/2.-0.001};
                          int side = (2*(step==0) - 1);
                          nestBox ( "STsupportSlatWeightBox"+std::string((step == 0) ? "Front" : "Back"),
                                    frontPars,
                                    slatMaterial,
                                    0,
                                    CLHEP::Hep3Vector(0., 0., side*(sideThickness/2.0 - currL/2.0)),
                                    weightMotherInfo,
                                    0,
                                    pabsIsVisible,
                                    G4Colour::Blue(),
                                    pabsIsSolid,
                                    forceAuxEdgeVisible,
                                    placePV,
                                    doSurfaceCheck);
                          zCenter += side*sideThickness/2.;
                          currL   -= sideThickness;
                          if(verbosityLevel > 0)
                            cout << " Extra panel params: {" << frontPars[0] << ", "
                                 << frontPars[1] << ", " << frontPars[2] << "} loc = "
                                 << "{0, 0, " << side*(sideThickness/2.0 - currL/2.0) << "}\n";
                        }
                        //add bottom and sides of the slat
                        const double botPars[3] = {weightBoxW/2-sideThickness, topThickness/2., currL/2.-0.001};
                        CLHEP::Hep3Vector botLoc(0., topThickness/2.0 - weightBoxH/2.0, zCenter);
                        nestBox ( "STsupportWeightBoxBottom"+std::to_string(step)+"_"+std::to_string(iSlat+1),
                                  botPars,
                                  slatMaterial,
                                  0,
                                  botLoc,
                                  weightMotherInfo,
                                  0,
                                  pabsIsVisible,
                                  G4Colour::Blue(),
                                  pabsIsSolid,
                                  forceAuxEdgeVisible,
                                  placePV,
                                  doSurfaceCheck);
                        if(verbosityLevel > 0)
                          cout << "Bottom panel params: {" << botPars[0] << ", "
                               << botPars[1] << ", " << botPars[2] << "} loc = "
                               << botLoc << endl;
                        const double weightPars[3] = {weightW/2, weightH/2., currL/2.-0.001};
                        nestBox ( "STsupportWeight"+std::to_string(step)+"_"+std::to_string(iSlat+1),
                                  weightPars,
                                  weightMaterial,
                                  0,
                                  CLHEP::Hep3Vector(0., topThickness - weightBoxH/2.0+weightH/2., zCenter),
                                  weightMotherInfo,
                                  0,
                                  pabsIsVisible,
                                  G4Colour::Blue(),
                                  pabsIsSolid,
                                  forceAuxEdgeVisible,
                                  placePV,
                                  doSurfaceCheck);
                        const double sidePars[3] = {sideThickness/2., weightBoxH/2., currL/2.-0.001};
                        nestBox ( "STsupportWeightBoxLeft"+std::to_string(step),
                                  sidePars,
                                  slatMaterial,
                                  0,
                                  CLHEP::Hep3Vector(sideThickness/2.0 - weightBoxW/2.0, 0., zCenter),
                                  weightMotherInfo,
                                  0,
                                  pabsIsVisible,
                                  G4Colour::Blue(),
                                  pabsIsSolid,
                                  forceAuxEdgeVisible,
                                  placePV,
                                  doSurfaceCheck);
                        nestBox ( "STsupportWeightBoxRight"+std::to_string(step),
                                  sidePars,
                                  slatMaterial,
                                  0,
                                  CLHEP::Hep3Vector(-sideThickness/2.0 + weightBoxW/2.0, 0., zCenter),
                                  weightMotherInfo,
                                  0,
                                  pabsIsVisible,
                                  G4Colour::Blue(),
                                  pabsIsSolid,
                                  forceAuxEdgeVisible,
                                  placePV,
                                  doSurfaceCheck);
                        currentZ = suppZEnd[step];

                      }
                    } else {
                      mf::LogError("GEOM") << "Undefined stopping target slat type = " << slatIndex
                                           << " in " << __PRETTY_FUNCTION__;
                    }
                  } else { //end if(opaVersion > 2)
                    // full lengths, not half lengths
                    slatLoc.setZ(slatLoc.z() - (hl+sLong/2.0) );
                    slatLoc += location; //add z offset
                    ostringstream slatName;
                    slatName << "STsupportSlat" << iSlat+1;

                    nestBox ( slatName.str(),
                              boxPars,
                              oPAsupportMaterial,
                              slatRotat,
                              slatLoc - parent1Info.centerInMu2e() + relPosFake,
                              parent1Info,
                              0,
                              pabsIsVisible,
                              G4Colour::Blue(),
                              pabsIsSolid,
                              forceAuxEdgeVisible,
                              placePV,
                              doSurfaceCheck);
                  }
                }
              } // end if iSub == 2  (adding slats for ST support)
            } else {
              auto vi = nestTubs( myName.str(),
                                  TubsParams( rin, rou, hl ),
                                  oPAsupportMaterial,
                                  0,
                                  location - parent2Info.centerInMu2e() + relPosFake,
                                  parent2Info,
                                  0,
                                  pabsIsVisible,
                                  G4Color::Blue(),
                                  pabsIsSolid,
                                  forceAuxEdgeVisible,
                                  placePV,
                                  doSurfaceCheck );
              if(verbosityLevel > 0) cout << " Support final center = " << vi.centerInMu2e() << endl;
            } // end of if for placing in DS2Vac or DS3Vac
          } // end of if else for support 2 (iSup == 1)
        } // end for loop over OPA supports
//-----------------------------------------------------------------------------
// in opa version 3 and above, add cross supports between support rings
// P.M.: need Z coordinate for the degrader converter positioning
//-----------------------------------------------------------------------------
        CLHEP::Hep3Vector motherBoxLoc(0,0,0);

        if(opaVersion > 2) {
          G4Material* crossSupMat = findMaterialOrThrow( pabs->crossSupportMaterial() );
          for(int iCross = 0; iCross < pabs->nCrossSupports(); ++iCross) {
            const double barThickness = pabs->crossSupportThicknesses().at(iCross);
            const double barWidth     = pabs->crossSupportWidth      ().at(iCross);
            const int oneIndex        = pabs->crossSupportOneIndex   ().at(iCross);
            const int twoIndex        = pabs->crossSupportTwoIndex   ().at(iCross);
            const double phi          = pabs->crossSupportPhis       ().at(iCross);
            const double height       = pabs->crossSupportHeights    ().at(iCross);
            const double radius       = pabs->crossSupportRadii      ().at(iCross);
            //get the width using the opa support ring positions
            const double z1  = pabs->oPAsupportZMidpoint().at(oneIndex);
            const double z2  = pabs->oPAsupportZMidpoint().at(twoIndex);
            const double hl1 = pabs->oPAsupportHalflength().at(oneIndex);
            const double hl2 = pabs->oPAsupportHalflength().at(twoIndex);

            const double smallGap = 1.e-3; //1 um buffers
            const double width = abs(z1-z2)-hl1-hl2-2.*smallGap; //offset by thickness of each so distance between edges
            const double motherBoxParams[] = {height/2., barWidth/2., width/2.};

            const double zMother = (z1>z2) ? (z1-hl1+z2+hl2)/2. : (z1+hl1+z2-hl2)/2.;
            const double xMother = radius*sin(phi*CLHEP::degree) + parent1Info.centerInMu2e().x();
            const double yMother = radius*cos(phi*CLHEP::degree) + parent1Info.centerInMu2e().y();

            motherBoxLoc.set(xMother,yMother,zMother);

            CLHEP::HepRotation* motherBoxRot = reg.add(CLHEP::HepRotation(CLHEP::HepRotation::IDENTITY));
            motherBoxRot->rotateZ(phi*CLHEP::degree);
            ostringstream motherBoxName;
            motherBoxName << "OPA_CrossSupportMother_" << iCross+1;
            VolumeInfo motherBoxInfo = nestBox( motherBoxName.str(),
                                                motherBoxParams,
                                                parent1Info.logical->GetMaterial(),
                                                motherBoxRot,
                                                motherBoxLoc-parent1Info.centerInMu2e() + relPosFake,
                                                parent1Info,
                                                0,
                                                pabsIsVisible,
                                                G4Colour::Blue(),
                                                pabsIsSolid,
                                                forceAuxEdgeVisible,
                                                placePV,
                                                doSurfaceCheck);

            /* build the individual bars in the cross support */

            //top and bottom bars are identical, cross entire structure
            const double topBarParams[] = {barThickness/2., barWidth/2., width/2.};
            CLHEP::Hep3Vector topBarLoc(height/2.-barThickness/2., 0., 0.);
            ostringstream topBarName;
            topBarName << "OPA_CrossSupport_" << iCross+1 << "_top";
            nestBox( topBarName.str(),
                     topBarParams,
                     crossSupMat,
                     0,
                     topBarLoc,
                     motherBoxInfo,
                     0,
                     pabsIsVisible,
                     G4Colour::Blue(),
                     pabsIsSolid,
                     forceAuxEdgeVisible,
                     placePV,
                     doSurfaceCheck);
            CLHEP::Hep3Vector botBarLoc(-height/2.+barThickness/2., 0., 0.);
            ostringstream botBarName;
            botBarName << "OPA_CrossSupport_" << iCross+1 << "_bottom";
            nestBox( botBarName.str(),
                     topBarParams, //same parameters as top
                     crossSupMat,
                     0,
                     botBarLoc,
                     motherBoxInfo,
                     0,
                     pabsIsVisible,
                     G4Colour::Blue(),
                     pabsIsSolid,
                     forceAuxEdgeVisible,
                     placePV,
                     doSurfaceCheck);
            //left and right bars are identical, cross until top/bottom bars
            const double leftBarParams[] = {height/2.-barThickness-smallGap, barWidth/2., barThickness/2.};
            CLHEP::Hep3Vector leftBarLoc(0., 0., -width/2.+barThickness/2.);
            ostringstream leftBarName;
            leftBarName << "OPA_CrossSupport_" << iCross+1 << "_left";
            nestBox( leftBarName.str(),
                     leftBarParams,
                     crossSupMat,
                     0,
                     leftBarLoc,
                     motherBoxInfo,
                     0,
                     pabsIsVisible,
                     G4Colour::Blue(),
                     pabsIsSolid,
                     forceAuxEdgeVisible,
                     placePV,
                     doSurfaceCheck);
            CLHEP::Hep3Vector rightBarLoc(0., 0., width/2.-barThickness/2.);
            ostringstream rightBarName;
            rightBarName << "OPA_CrossSupport_" << iCross+1 << "_right";
            nestBox( rightBarName.str(),
                     leftBarParams,
                     crossSupMat,
                     0,
                     rightBarLoc,
                     motherBoxInfo,
                     0,
                     pabsIsVisible,
                     G4Colour::Blue(),
                     pabsIsSolid,
                     forceAuxEdgeVisible,
                     placePV,
                     doSurfaceCheck);
            //first cross bar
            const double crossBarAngle = atan((height-barThickness*2.)/(width-barThickness*2.));
            double crossBarLength = sqrt((height-barThickness*2.)*(height-barThickness*2.) +
                                         (width-barThickness*2.)*(width-barThickness*2.));
            //avoid overlap
            crossBarLength -= (height > width) ? height/width*barThickness : width/height*barThickness;
            crossBarLength -= 1.; //remove mm to ensure no overlaps of corners
            CLHEP::Hep3Vector crossBar1Loc(0., 0., 0.);
            CLHEP::Hep3Vector crossBar2Loc(0., 0., 0.); //wrt bar 1
            CLHEP::HepRotation* crossBar1Rot = reg.add(CLHEP::HepRotation(CLHEP::HepRotation::IDENTITY));
            crossBar1Rot->rotateY(crossBarAngle);
            CLHEP::HepRotation* crossBar2Rot = reg.add(CLHEP::HepRotation(CLHEP::HepRotation::IDENTITY));
            crossBar2Rot->rotateY(-2.*crossBarAngle); //double to account for rotation wrt bar 1
            ostringstream crossBarName;
            crossBarName << "OPA_CrossSupport_" << iCross+1 << "_CrossBar";
            ostringstream crossBar1Name;
            crossBar1Name << "OPA_CrossSupport_" << iCross+1 << "_CrossBar_1";
            ostringstream crossBar2Name;
            crossBar2Name << "OPA_CrossSupport_" << iCross+1 << "_CrossBar_2";
            G4Box* crossBar1Box= new G4Box(crossBar1Name.str(), barThickness/2., barWidth/2., crossBarLength/2.);
            G4Box* crossBar2Box= new G4Box(crossBar2Name.str(), barThickness/2., barWidth/2., crossBarLength/2.);
            VolumeInfo crossBarInfo(crossBarName.str(), crossBar1Loc, parent1Info.centerInMu2e() + relPosFake);
            crossBarInfo.solid = new G4UnionSolid(crossBar1Name.str()
                                                  , crossBar1Box, crossBar2Box, crossBar2Rot , crossBar2Loc);
            finishNesting( crossBarInfo,
                           crossSupMat,
                           crossBar1Rot,
                           crossBar1Loc,
                           motherBoxInfo.logical,
                           0,
                           pabsIsVisible,
                           G4Colour::Blue(),
                           pabsIsSolid,
                           forceAuxEdgeVisible,
                           placePV,
                           doSurfaceCheck);



          } //end loop of cross supports
        } //end opa version check for cross support
      } // end if nSupports > 0 for OPA

      //***************************
      // Now build the Pion Degrader if requested
      //***************************

      if ( pabs->degraderBuild() ) {
        int degraderVersion = pabs->degraderVersion();

        CLHEP::HepRotation* degraderRot = reg.add(CLHEP::HepRotation(CLHEP::HepRotation::IDENTITY));
        degraderRot->rotateZ(pabs->degraderRotation()*CLHEP::degree);

        // Make Frame, 'frame' is a mounting plate connecting the shaft and the filter
        std::vector<double> frameDims  = pabs->degraderFrameDims();
//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
        std::vector<double> filterDims  = pabs->degraderFilterDims();
        std::vector<double> counterDims = pabs->degraderCounterwtDims();
        std::vector<double> rodDims     = pabs->degraderRodDims();
        std::vector<double> pivotPos    = pabs->degraderPivotPos();
        std::vector<double> armDims     = pabs->degraderSupportArmDims();
        std::vector<double> plateDims   = pabs->degraderSupportPlateDims();

        std::vector<double> filter2Dims   = pabs->degraderFilter2Dims  ();
        std::vector<double> converterDims = pabs->degraderConverterDims();
        double              converterDz   = pabs->degraderConverterDz  ();
//-----------------------------------------------------------------------------
// Calculate the location in x and y for the mother box location
// This is the mother box containing frame, filter, and rod.
// Counterweight has to be put directly into DS2
//-----------------------------------------------------------------------------
        double lengthPDMoBox, dptoc;
        if (degraderVersion == 3) {
//-----------------------------------------------------------------------------
// P.M. : in v3, the converter is paced into the DegraderMother volume
// in v1,v2 there is no converter at all
// in v4, the converter is placed into the DS2
//-----------------------------------------------------------------------------
          lengthPDMoBox = converterDims.at(1) + frameDims.at(3) + counterDims.at(3);
          dptoc         = converterDims.at(1) + frameDims.at(3) - lengthPDMoBox/2.0;
        }
        else {
          lengthPDMoBox = frameDims.at(1) + frameDims.at(3) + counterDims.at(3);
          dptoc         = frameDims.at(1) + frameDims.at(3) - lengthPDMoBox/2.0;
        }

        double xLocMoBox     = pivotPos.at(0) - dptoc*cos(pabs->degraderRotation()*CLHEP::degree);
        double yLocMoBox     = pivotPos.at(1) + dptoc*sin(pabs->degraderRotation()*CLHEP::degree);

        CLHEP::Hep3Vector locationInMu2e;

        if (degraderVersion ==  3) {
          locationInMu2e.set(xLocMoBox,yLocMoBox,pabs->degraderZ0()+2.0*filterDims.at(2)+frameDims.at(2));
        }
        else {
          locationInMu2e.set(xLocMoBox,yLocMoBox,pabs->degraderZ0()+frameDims[2]);
        }

        if (verbosityLevel > 0) {
          printf(">> degrader location in Mu2e: X,Y,Z= %10.3f %10.3f %10.3f\n",
                 locationInMu2e.x(),locationInMu2e.y(),locationInMu2e.z());
        }
//-----------------------------------------------------------------------------
// degrader mother volume
//-----------------------------------------------------------------------------
        std::string motherName("Degrader");
        VolumeInfo  degraderMother (motherName,
                                    locationInMu2e - parent1Info.centerInMu2e() + relPosFake, parent1Info.centerInWorld);
//-----------------------------------------------------------------------------
// Make box for degrader mother volume.
//-----------------------------------------------------------------------------
        if (degraderVersion == 3) {
          degraderMother.solid = new G4Box ( "degraderOutline",
                                             lengthPDMoBox/2.0,
                                             converterDims.at(1),  // is converter larger than the frame
                                             frameDims[2] + filterDims[2] + filter2Dims[2]+converterDims[2]);
        }
        else {
          degraderMother.solid = new G4Box ( "degraderOutline",
                                             lengthPDMoBox/2.0,
                                             frameDims.at(1),
                                             2.0*filterDims.at(2) + frameDims.at(2) + 1);
        }

        // Now put degraderMother in DS2Vacuum
        finishNesting ( degraderMother,
                        findMaterialOrThrow("DSVacuum"),
                        degraderRot,
                        degraderMother.centerInParent,
                        parent1Info.logical,
                        0, false, G4Colour::Red(),
                        false,
                        forceAuxEdgeVisible,
                        placePV,
                        doSurfaceCheck );
//-----------------------------------------------------------------------------
// Start cobbling pieces together by putting frame in mother
//-----------------------------------------------------------------------------
        double mother_dz2 = ((G4Box*) degraderMother.solid)->GetZHalfLength();

        if (degraderVersion == 3) {
//-----------------------------------------------------------------------------
// v3: two stopping disks - one of polyethylene and one of lead,
//                          make the sum equivalent to 4 mm of Ti
//     for a ~1 mm gold converter at R=25 cm, position of the degrader disk may be offset upstream in Z
//     - cant move the converter downstream
//     - the disk shaft will be longer by ~ 1cm, but should still fit in...
// converter goes first - for validation, it can be shifted by dz back
//-----------------------------------------------------------------------------
          //          CLHEP::Hep3Vector trans3b(frameDims.at(1) - lengthPDMoBox/2.0,
          CLHEP::Hep3Vector trans3b(converterDims.at(1) - lengthPDMoBox/2.0,
                                    0,
                                    mother_dz2-2*frameDims[2]-converterDims[2]-converterDz);
          nestTubs("degraderConverter",
                   TubsParams(converterDims.at(0),converterDims.at(1),converterDims.at(2)),
                   findMaterialOrThrow(pabs->degraderConverterMaterial()),
                   0, trans3b, degraderMother,
                   0, pabsIsVisible, G4Color::Red(),
                   pabsIsSolid,
                   forceAuxEdgeVisible,
                   placePV,
                   doSurfaceCheck );
//-----------------------------------------------------------------------------
// then goes filter2 (Pb)
//-----------------------------------------------------------------------------
          //          CLHEP::Hep3Vector trans2b(frameDims.at(1) - lengthPDMoBox/2.0,
          CLHEP::Hep3Vector trans2b(converterDims.at(1) - lengthPDMoBox/2.0,
                                    0,
                                    //                                    -(frameDims.at(2) + filter2Dims.at(2)));
                                    mother_dz2-2*frameDims[2]-2*converterDims[2]+filter2Dims[2]);
          nestTubs("degraderFilter2",
                   TubsParams(filter2Dims.at(0),filter2Dims.at(1),filter2Dims.at(2)),
                   findMaterialOrThrow(pabs->degraderFilter2Material()),
                   0, trans2b, degraderMother,
                   0, pabsIsVisible, G4Color::Red(),
                   pabsIsSolid,
                   forceAuxEdgeVisible,
                   placePV,
                   doSurfaceCheck );
          //          CLHEP::Hep3Vector trans1b(frameDims.at(1)-lengthPDMoBox/2.0,
          CLHEP::Hep3Vector trans1b(converterDims.at(1)-lengthPDMoBox/2.0,
                                    0,
                                    //                                    -(frameDims.at(2) + filterDims.at(2)));
                                    mother_dz2-2*frameDims[2]-2*converterDims[2]-2*filter2Dims[2]-filterDims[2]);
          nestTubs("degraderFilter",
                   TubsParams(filterDims.at(0),filterDims.at(1),filterDims.at(2)),
                   findMaterialOrThrow(pabs->degraderFilterMaterial()),
                   0, trans1b, degraderMother,
                   0, pabsIsVisible, G4Color::Red(),
                   pabsIsSolid,
                   forceAuxEdgeVisible,
                   placePV,
                   doSurfaceCheck );
        }
        else {
//-----------------------------------------------------------------------------
// v1,v2, v4
// P.M.: for v4, I don't have the right frame positioning, for now, do not build it
//-----------------------------------------------------------------------------
          if (degraderVersion < 3) {
            CLHEP::Hep3Vector trans1(frameDims.at(1) - lengthPDMoBox/2.0, 0, 0);

            nestTubs("degraderFrame",
                     TubsParams(frameDims.at(0),frameDims.at(1),frameDims.at(2),
                                -frameDims.at(4)/2 * CLHEP::degree,
                              frameDims.at(4)*CLHEP::degree ),
                     findMaterialOrThrow(pabs->degraderFrameMaterial()),
                     0, trans1, degraderMother,
                     0, pabsIsVisible, G4Color::Red(),
                     pabsIsSolid,
                     forceAuxEdgeVisible,
                     placePV,
                     doSurfaceCheck );
          }
//-----------------------------------------------------------------------------
// Now put filter in mother
//-----------------------------------------------------------------------------
          CLHEP::Hep3Vector trans1b(frameDims.at(1) - lengthPDMoBox/2.0,
                                    0,
                                    mother_dz2-2*frameDims.at(2)-filterDims.at(2));
          nestTubs("degraderFilter",
                   TubsParams(filterDims.at(0),filterDims.at(1),filterDims.at(2)),
                   findMaterialOrThrow(pabs->degraderFilterMaterial()),
                   0, trans1b, degraderMother,
                   0, pabsIsVisible, G4Color::Red(),
                   pabsIsSolid,
                   forceAuxEdgeVisible,
                   placePV,
                   doSurfaceCheck );
//-----------------------------------------------------------------------------
// in v4, add converter, directly into the DS2 volume
//-----------------------------------------------------------------------------
          if (degraderVersion == 4) {
            //          CLHEP::Hep3Vector trans3b(frameDims.at(1) - lengthPDMoBox/2.0,

            CLHEP::Hep3Vector trans(pabs3Offset.x(),
                                    pabs3Offset.y(),
                                    pabs3Offset.z()-pabs->part(2).halfLength()+converterDims[2]);

            nestTubs("degraderConverter",
                     TubsParams(converterDims.at(0),converterDims.at(1),converterDims.at(2)),
                     findMaterialOrThrow(pabs->degraderConverterMaterial()),
                     0,
                     trans,
                     parent1Info,
                     0, pabsIsVisible, G4Color::Red(),
                     pabsIsSolid,
                     forceAuxEdgeVisible,
                     placePV,
                     doSurfaceCheck );
          }
        }
//-----------------------------------------------------------------------------
// Create rod rep and put it into degraderMother, if goes right underneath the frame
// in v3, temorary, don't put in the rod - that doesn't matter
//-----------------------------------------------------------------------------
        if (degraderVersion < 3) {
          double lenRod = frameDims.at(3) + counterDims.at(3) - frameDims.at(1)- 0.1;
          std::vector<double> lwhs = {lenRod/2.0,rodDims.at(0)/2.0,rodDims.at(1)/2.0};

          CLHEP::Hep3Vector trans3(frameDims.at(1), 0.0, mother_dz2-rodDims.at(1)/2);

          nestBox ("degraderRod",
                   lwhs,
                   findMaterialOrThrow(pabs->degraderRodMaterial()),
                   0, trans3, degraderMother,
                   0, pabsIsVisible, G4Color::Red(),
                   pabsIsSolid,
                   forceAuxEdgeVisible,
                   placePV,
                   doSurfaceCheck );
        }
//-----------------------------------------------------------------------------
// Now put Counterweight directly into DS2Vacuum.
// Has to be separate from the rod, frame, and filter because its
// size would cause overlaps otherwise.
//-----------------------------------------------------------------------------
        double xLocMoBox2 = pivotPos.at(0) +
          (counterDims.at(3)+counterDims.at(2)) * cos(pabs->degraderRotation()
                                                      *CLHEP::degree);
        double yLocMoBox2 = pivotPos.at(1) -
          (counterDims.at(3) + counterDims.at(2))*sin(pabs->degraderRotation()
                                                      *CLHEP::degree);

        CLHEP::Hep3Vector location2InMu2e (xLocMoBox2, yLocMoBox2,
                                           pabs->degraderZ0()
                                           + 2.0*filterDims.at(2)
                                           + frameDims.at(2) );

        // Now put counterweight in DS2Vacuum
        CLHEP::HepRotation* cwtRot = reg.add(CLHEP::HepRotation(CLHEP::HepRotation::IDENTITY));
        cwtRot->rotateY(90.0*CLHEP::degree);
        cwtRot->rotateX(pabs->degraderRotation()*CLHEP::degree);

        nestTubs("degraderCounterweight",
                 TubsParams(counterDims.at(0),counterDims.at(1),counterDims.at(2)),
                 findMaterialOrThrow(pabs->degraderCountwtMaterial()),
                 cwtRot, location2InMu2e-parent1Info.centerInMu2e() + relPosFake,
                 parent1Info,
                 0, pabsIsVisible, G4Color::Red(),
                 pabsIsSolid,
                 forceAuxEdgeVisible,
                 placePV,
                 doSurfaceCheck );

        // Now put in bracket rods (four of them)
        if ( armDims.at(0) > 0 && armDims.at(3) > 0) {
          std::vector<double> lwhs2 = {armDims.at(0),armDims.at(1),armDims.at(2)};
          // Get locations
          CLHEP::Hep3Vector transa0(pivotPos.at(0) + armDims.at(3),
                                    pivotPos.at(1) + armDims.at(4),
                                    pabs->degraderZ0() + armDims.at(5));

          nestBox ("degraderSupportArm0",
                   lwhs2,
                   findMaterialOrThrow(pabs->degraderSupportMaterial()),
                   0, transa0 -parent1Info.centerInMu2e() + relPosFake, parent1Info,
                   0, pabsIsVisible, G4Color::Red(),
                   pabsIsSolid,
                   forceAuxEdgeVisible,
                   placePV,
                   doSurfaceCheck );

          CLHEP::Hep3Vector transa1(pivotPos.at(0) + armDims.at(3),
                                    pivotPos.at(1) - armDims.at(4),
                                    pabs->degraderZ0() + armDims.at(5));

          nestBox ("degraderSupportArm1",
                   lwhs2,
                   findMaterialOrThrow(pabs->degraderSupportMaterial()),
                   0, transa1 -parent1Info.centerInMu2e() + relPosFake, parent1Info,
                   0, pabsIsVisible, G4Color::Red(),
                   pabsIsSolid,
                   forceAuxEdgeVisible,
                   placePV,
                   doSurfaceCheck );

          CLHEP::Hep3Vector transa2(pivotPos.at(0) - armDims.at(3),
                                    pivotPos.at(1) - armDims.at(4),
                                    pabs->degraderZ0() + armDims.at(5));

          nestBox ("degraderSupportArm2",
                   lwhs2,
                   findMaterialOrThrow(pabs->degraderSupportMaterial()),
                   0, transa2 -parent1Info.centerInMu2e() + relPosFake, parent1Info,
                   0, pabsIsVisible, G4Color::Red(),
                   pabsIsSolid,
                   forceAuxEdgeVisible,
                   placePV,
                   doSurfaceCheck );

          CLHEP::Hep3Vector transa3(pivotPos.at(0) - armDims.at(3),
                                    pivotPos.at(1) + armDims.at(4),
                                    pabs->degraderZ0() + armDims.at(5));

          nestBox ("degraderSupportArm3",
                   lwhs2,
                   findMaterialOrThrow(pabs->degraderSupportMaterial()),
                   0, transa3 -parent1Info.centerInMu2e() + relPosFake, parent1Info,
                   0, pabsIsVisible, G4Color::Red(),
                   pabsIsSolid,
                   forceAuxEdgeVisible,
                   placePV,
                   doSurfaceCheck );

          std::vector<double> lwhs3 = {plateDims.at(0), plateDims.at(1),
                                       plateDims.at(2)};

          CLHEP::Hep3Vector transp1(pivotPos.at(0),
                                    pivotPos.at(1),
                                    pabs->degraderZ0() + armDims.at(5)
                                    - armDims.at(2)-plateDims.at(2) -0.1);

          nestBox ("degraderSupportPlate",
                   lwhs3,
                   findMaterialOrThrow(pabs->degraderSupportMaterial()),
                   0, transp1 - parent1Info.centerInMu2e() + relPosFake, parent1Info,
                   0, pabsIsVisible, G4Color::Red(),
                   pabsIsSolid,
                   forceAuxEdgeVisible,
                   placePV,
                   doSurfaceCheck );

        } //end of if on non-zero support bracket arm dims

      }// end of if degraderBuild


      if ( pabs->buildSupports() ) {

        AntiLeakRegistry& reg = art::ServiceHandle<Mu2eG4Helper>()->antiLeakRegistry();

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
                    endRing.originInMu2e()-parent1Info.centerInMu2e() + relPosFake,
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
          double wire_angle_offset = ipaSup->wireAngleOffset()*CLHEP::degree; //control rotating wires around IPA together
          for ( std::size_t iWire(0); iWire < ipaSup->nWiresPerSet() ; iWire++ ) {

            Tube supportWire = ipaSup->getWire( iSet, iWire );

            ostringstream wirename ; wirename << "IPAsupport_set" << iS << "_wire" << ++iW ;

            const double rStartOfWire = pabs1rOut0+(supportWire.originInMu2e().z()-zstartOfIPA)/pabs1len*(pabs1rOut1-pabs1rOut0);
            const double wirePhi = iWire * 360.*CLHEP::deg / ipaSup->nWiresPerSet() + wire_angle_offset ;
            CLHEP::Hep3Vector additionalOffset ( (supportWire.halfLength()+0.005+rStartOfWire) *
                                                 std::cos(wirePhi),
                                                 (supportWire.halfLength()+0.005+rStartOfWire) *
                                                 std::sin(wirePhi ),
                                                 0 );

            // Now get appropriate rotation angles
            G4RotationMatrix* supportRot = reg.add(G4RotationMatrix());

            if (ipa_version == 1) {
              supportRot->rotateY(-M_PI/2.);
              supportRot->rotateZ(-wirePhi );
            }
            else if (ipa_version >= 2) {
              CLHEP::Hep3Vector rotationAxis(0, 1, 0); // start off with rotating around the y-axis
              rotationAxis.rotateZ(wirePhi); // each wire wants to be rotated arounf a slightly different axis
              supportRot->rotate(wireRotation*CLHEP::deg, rotationAxis);

              CLHEP::Hep3Vector extraZOffset(0, 0, supportWire.halfLength() * std::cos(wireRotation*CLHEP::deg)); // because of the rotation from the vertical
              CLHEP::Hep3Vector extraROffset(std::cos(wirePhi)
                                             * supportWire.halfLength()*(std::cos((90-wireRotation*CLHEP::deg)))*std::tan(wireRotation*CLHEP::deg),
                                             std::sin(wirePhi)
                                             * supportWire.halfLength()*(std::cos((90-wireRotation*CLHEP::deg)))*std::tan(wireRotation*CLHEP::deg),
                                             0);
              if (ipa_version > 2) {
                //correct z offset bug that cos should be sin, since angle from verticle, theta = 0 --> 0 offset (angles are 45 degrees so not noticed)
                extraZOffset.setZ(supportWire.halfLength() * std::sin(wireRotation*CLHEP::deg));
                //correct using 90 instead of 90*CLHEP::deg
                extraROffset.setX(std::cos(wirePhi)
                                  * supportWire.halfLength()*(std::cos((90.-wireRotation)*CLHEP::deg))*std::tan(wireRotation*CLHEP::deg));
                extraROffset.setY(std::sin(wirePhi)
                                  * supportWire.halfLength()*(std::cos((90.-wireRotation)*CLHEP::deg))*std::tan(wireRotation*CLHEP::deg));
              }
              if (ipa_version == 2)
                additionalOffset -= 0.90*extraROffset;
              else
                additionalOffset -= 0.99*extraROffset;

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
                      supportWire.originInMu2e()+additionalOffset-parent1Info.centerInMu2e() + relPosFake,
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
