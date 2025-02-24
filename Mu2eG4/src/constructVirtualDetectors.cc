//
// Free function to create the virtual detectors
//
//
// Original author KLG based on Mu2eWorld constructVirtualDetectors

// C++ includes
#include <iostream>
#include <string>

// art includes
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"

// Mu2e includes.
#include "Offline/Mu2eG4/inc/constructVirtualDetectors.hh"

#include "Offline/ConfigTools/inc/SimpleConfig.hh"
#include "Offline/ConfigTools/inc/checkForStale.hh"

#include "Offline/BeamlineGeom/inc/Beamline.hh"
#include "Offline/CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "Offline/CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "Offline/DetectorSolenoidGeom/inc/DetectorSolenoid.hh"
#include "Offline/Mu2eG4Helper/inc/Mu2eG4Helper.hh"
#include "Offline/Mu2eG4Helper/inc/VolumeInfo.hh"
#include "Offline/Mu2eG4Helper/inc/AntiLeakRegistry.hh"
#include "Offline/GeomPrimitives/inc/Tube.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/GeometryService.hh"
#include "Offline/GeometryService/inc/G4GeometryOptions.hh"
#include "Offline/GeometryService/inc/VirtualDetector.hh"
#include "Offline/DataProducts/inc/VirtualDetectorId.hh"
#include "Offline/MECOStyleProtonAbsorberGeom/inc/MECOStyleProtonAbsorber.hh"
#include "Offline/Mu2eG4/inc/checkForOverlaps.hh"
#include "Offline/Mu2eG4/inc/findMaterialOrThrow.hh"
#include "Offline/Mu2eG4/inc/finishNesting.hh"
#include "Offline/Mu2eG4/inc/nestTubs.hh"
#include "Offline/Mu2eG4/inc/nestBox.hh"
#include "Offline/ProductionSolenoidGeom/inc/PSVacuum.hh"
#include "Offline/ProductionSolenoidGeom/inc/PSEnclosure.hh"
#include "Offline/ProtonBeamDumpGeom/inc/ProtonBeamDump.hh"
#include "Offline/ProductionTargetGeom/inc/ProductionTarget.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"

// G4 includes
#include "Geant4/G4Material.hh"
#include "Geant4/G4SDManager.hh"
#include "Geant4/G4Color.hh"
#include "Geant4/G4Tubs.hh"
#include "Geant4/G4Cons.hh"
#include "Geant4/G4SubtractionSolid.hh"
#include "Geant4/G4IntersectionSolid.hh"
#include "Geant4/G4Box.hh"

using namespace std;

namespace mu2e {

  // Construct the virtual detectors

  void constructVirtualDetectors( const SimpleConfig& _config ){

    // Place virtual detectors
    const auto geomOptions = art::ServiceHandle<GeometryService>()->geomOptions();
    geomOptions->loadEntry( _config, "vd", "vd");

    const bool vdIsVisible          = geomOptions->isVisible("vd");
    const bool vdIsSolid            = geomOptions->isSolid("vd");
    const bool forceAuxEdgeVisible  = geomOptions->forceAuxEdgeVisible("vd");
    const bool doSurfaceCheck       = geomOptions->doSurfaceCheck("vd");
    const bool placePV              = geomOptions->placePV("vd");
    int static const verbosityLevel = _config.getInt("vd.verbosityLevel",0);

    AntiLeakRegistry& reg = art::ServiceHandle<Mu2eG4Helper>()->antiLeakRegistry();

    GeomHandle<VirtualDetector> vdg;
    if( vdg->nDet()<=0 ) return;

    GeomHandle<Beamline> beamg;

    GeomHandle<DetectorSolenoid> ds;
    TransportSolenoid const&  ts = beamg->getTS();
    G4Material* downstreamVacuumMaterial = findMaterialOrThrow( ds->vacuumMaterial() );
    G4Material* upstreamVacuumMaterial   = findMaterialOrThrow(  ts.upstreamVacuumMaterial() );

    double rCol = ts.getColl51().rOut();
    double vdHalfLength = CLHEP::mm * vdg->getHalfLength();

    TubsParams vdParams(0,rCol,vdHalfLength);

    // Virtual Detectors Coll1_In, COll1_Out are placed inside TS1

    Mu2eG4Helper* _helper = &(*(art::ServiceHandle<Mu2eG4Helper>()));

    if(verbosityLevel>0) {
      VirtualDetectorId::printAll();
    }

    // FIXME: one should factorize some the code below; the main
    // things which change: parent and offset
    for( int vdId=VirtualDetectorId::Coll1_In;
         vdId<=VirtualDetectorId::Coll1_Out;
         ++vdId) if( vdg->exist(vdId) ) {
        VolumeInfo const & parent = _helper->locateVolInfo("TS1Vacuum");
        if ( verbosityLevel > 0) {
          cout << __func__ << " constructing " << VirtualDetector::volumeName(vdId)
               << " at " << vdg->getGlobal(vdId) << endl;
          cout << __func__ << "    VD parameters: " << vdParams << endl;
          cout << __func__ << "    VD rel. posit: " << vdg->getLocal(vdId) << endl;
        }

        VolumeInfo vd = nestTubs( VirtualDetector::volumeName(vdId),
                                  vdParams, upstreamVacuumMaterial, 0,
                                  vdg->getLocal(vdId),
                                  parent,
                                  vdId, vdIsVisible, G4Color::Red(), vdIsSolid,
                                  forceAuxEdgeVisible,
                                  placePV,
                                  false);

        doSurfaceCheck && checkForOverlaps(vd.physical, _config, verbosityLevel>0);

      }

    // ************************** DNB (Lou) Jan 2016 **********
    // Virtual Detector TS2_Bend is placed inside TS2

    int myvdId = VirtualDetectorId::TS2_Bend;
    if ( vdg->exist(myvdId) ) {
      VolumeInfo const & parent = _helper->locateVolInfo("TS2Vacuum");
      if ( verbosityLevel > 0 ) {
        cout << __func__ << " constructing TS2_Bend " <<  " at "
             << vdg->getGlobal(myvdId) << endl;
        cout << __func__ << "    VD parameters: " << vdParams << endl;
        cout << __func__ << "    VD rel. posit: " << vdg->getLocal(myvdId) << endl;

      }
      CLHEP::HepRotation* tsBendrot =
        new CLHEP::HepRotation(CLHEP::HepRotation::IDENTITY);
      tsBendrot->rotateX(90.0*CLHEP::deg);
      tsBendrot->rotateY(-45.0*CLHEP::deg);
      VolumeInfo myvd = nestTubs( VirtualDetector::volumeName(myvdId),
                                  vdParams, upstreamVacuumMaterial, tsBendrot,
                                  vdg->getLocal(myvdId),
                                  parent, myvdId, vdIsVisible, G4Color::Red(),
                                  vdIsSolid, forceAuxEdgeVisible,
                                  placePV, false);

      doSurfaceCheck && checkForOverlaps(myvd.physical,
                                         _config, verbosityLevel>0);

    }

    // Virtual Detector TS4_Bend is placed inside TS4
    myvdId = VirtualDetectorId::TS4_Bend;
    if ( vdg->exist(myvdId) ) {
      VolumeInfo const & parent = _helper->locateVolInfo("TS4Vacuum");
      if ( verbosityLevel > 0 ) {
        cout << __func__ << " constructing TS4_Bend " <<  " at "
             << vdg->getGlobal(myvdId) << endl;
        cout << __func__ << "    VD parameters: " << vdParams << endl;
        cout << __func__ << "    VD rel. posit: " << vdg->getLocal(myvdId) << endl;

      }
      CLHEP::HepRotation* tsBendrot =
        new CLHEP::HepRotation(CLHEP::HepRotation::IDENTITY);
      tsBendrot->rotateX(90.0*CLHEP::deg);
      tsBendrot->rotateY(-45.0*CLHEP::deg);
      VolumeInfo myvd = nestTubs( VirtualDetector::volumeName(myvdId),
                                  vdParams, downstreamVacuumMaterial, tsBendrot,
                                  vdg->getLocal(myvdId),
                                  parent, myvdId, vdIsVisible, G4Color::Red(),
                                  vdIsSolid, forceAuxEdgeVisible,
                                  placePV, false);

      doSurfaceCheck && checkForOverlaps(myvd.physical,
                                         _config, verbosityLevel>0);

    }

    //***************************
    // Virtual Detectors Coll31_In, Coll31_Out, Coll32_In, Coll32_Out are placed inside TS3

    for( int vdId=VirtualDetectorId::Coll31_In;
         vdId<=VirtualDetectorId::Coll32_Out;
         ++vdId) if( vdg->exist(vdId) ) {
        VolumeInfo const & parent = _helper->locateVolInfo("TS3Vacuum");
        if ( verbosityLevel > 0) {
          cout << __func__ << " constructing " << VirtualDetector::volumeName(vdId)
               << " at " << vdg->getGlobal(vdId) << endl;
        }
        VolumeInfo vd = nestTubs( VirtualDetector::volumeName(vdId),
                                  vdParams, upstreamVacuumMaterial, 0,
                                  vdg->getLocal(vdId),
                                  parent,
                                  vdId, vdIsVisible, G4Color::Red(), vdIsSolid,
                                  forceAuxEdgeVisible,
                                  placePV,
                                  false);

        doSurfaceCheck && checkForOverlaps(vd.physical, _config, verbosityLevel>0);

      }

    // Virtual Detectors Coll5_In, Coll5_Out are placed inside TS5

    for( int vdId=VirtualDetectorId::Coll5_In;
         vdId<=VirtualDetectorId::Coll5_Out;
         ++vdId) if( vdg->exist(vdId) ) {
        VolumeInfo const & parent = _helper->locateVolInfo("TS5Vacuum");
        if ( verbosityLevel > 0) {
          cout << __func__ << " constructing " << VirtualDetector::volumeName(vdId)
               << " at " << vdg->getGlobal(vdId) <<  " parent: " << parent.centerInMu2e() << endl;
        }
        VolumeInfo vd = nestTubs( VirtualDetector::volumeName(vdId),
                                  vdParams, downstreamVacuumMaterial, 0,
                                  vdg->getLocal(vdId),
                                  parent,
                                  vdId, vdIsVisible, G4Color::Red(), vdIsSolid,
                                  forceAuxEdgeVisible,
                                  placePV,
                                  false);

        doSurfaceCheck && checkForOverlaps(vd.physical, _config, verbosityLevel>0);

      }

    // Virtual Detectors Coll5_OutSurf surrounds the outer cylindrical surface of collimator in TS5

    int vdId = VirtualDetectorId::Coll5_OutSurf;
    if( vdg->exist(vdId) ) {

      if ( verbosityLevel > 0) {
        cout << __func__ << " constructing " << VirtualDetector::volumeName(vdId)  << endl;
      }

      // the detector is on the outer surface of the coll5
      // it is thin cylinder, NOT a thin disk

      VolumeInfo const & parent = _helper->locateVolInfo("TS5Vacuum");

      double coll5OuterRadius    = beamg->getTS().getColl51().rOut();
      double coll5HalfLength     = beamg->getTS().getColl51().halfLength();

      TubsParams  vdParamsColl5OutSurf(coll5OuterRadius - 2.*vdHalfLength,
                                       coll5OuterRadius,
                                       coll5HalfLength - 2.*vdHalfLength);

      VolumeInfo vd = nestTubs( VirtualDetector::volumeName(vdId),
                                vdParamsColl5OutSurf, downstreamVacuumMaterial, 0,
                                vdg->getLocal(vdId),
                                parent,
                                vdId, vdIsVisible, G4Color::Red(), vdIsSolid,
                                forceAuxEdgeVisible,
                                placePV,
                                false);

      doSurfaceCheck && checkForOverlaps(vd.physical, _config, verbosityLevel>0);

    }

    /*************************************************/
    /* new virtual detector**************************/


    if ( !_config.getBool("isDumbbell",false) ){
      double Ravr = ds->rIn1();

      if ( _config.getBool("hasTSdA",false) ) {
        Ravr = _config.getDouble("TSdA.rFactorForVDs");
      }

      bool opaflag = false;
      double opaz0, opaz1, opari0, opari1;
      if ( _config.getBool("hasProtonAbsorber", true) ) {
        GeomHandle<MECOStyleProtonAbsorber> pageom;
        if ( pageom->isAvailable(ProtonAbsorberId::opabs1) ) {
          opaflag = true;
          MECOStyleProtonAbsorberPart opa = pageom->part(2);
          opaz0 = opa.center().z()-opa.halfLength();
          opaz1 = opa.center().z()+opa.halfLength();
          opari0 = opa.innerRadiusAtStart();
          opari1 = opa.innerRadiusAtEnd();
        }
      }
      vdId = VirtualDetectorId::STMUpstream;
      if( vdg->exist(vdId) ) {

        if ( verbosityLevel > 0) {
          cout << __func__ << " constructing " << VirtualDetector::volumeName(vdId)  << endl;
        }

        double zvd = vdg->getGlobal(vdId).z();
        if (opaflag) {
          Ravr = (opari1 - opari0)/(opaz1 - opaz0) * (zvd - opaz0) + opari0;
        }
        double rvd = Ravr - 5.0;

        if ( verbosityLevel > 0) {
          cout << __func__ << " " << VirtualDetector::volumeName(vdId) <<
            " z, r : " << zvd << ", " << rvd << endl;
        }

        TubsParams vdParamsSTMUpstream(0.,rvd,vdHalfLength);
        std::string theDS3("DS3Vacuum");
        if ( _config.getBool("inGaragePosition",false) ) theDS3 = "garageFakeDS3Vacuum";
        VolumeInfo const & parent = ( _config.getBool("isDumbbell",false) ) ?
          _helper->locateVolInfo(theDS3) :
          _helper->locateVolInfo("DS2Vacuum"); //DS3Vacuum to move the targets

        CLHEP::Hep3Vector const& parentInMu2e = parent.centerInMu2e();

        if (verbosityLevel >0) {
          cout << __func__ << " " << VirtualDetector::volumeName(vdId) << " Z offset in Mu2e    : " <<
            zvd << endl;
          cout << __func__ << " " << VirtualDetector::volumeName(vdId) << " Z extent in Mu2e    : " <<
            zvd - vdHalfLength << ", " << zvd + vdHalfLength << "\n"
               << " at " << vdg->getGlobal(vdId) <<  " parent: " << parent.centerInMu2e() << endl;
        }

        VolumeInfo vd = nestTubs( VirtualDetector::volumeName(vdId),
                                  vdParamsSTMUpstream, downstreamVacuumMaterial, 0,
                                  vdg->getLocal(vdId),
                                  parent,
                                  vdId,
                                  vdIsVisible,
                                  G4Color::Red(), vdIsSolid,
                                  forceAuxEdgeVisible,
                                  placePV,
                                  false);

        doSurfaceCheck && checkForOverlaps(vd.physical, _config, verbosityLevel>0);

        if ( verbosityLevel > 0) {
          cout << __func__ << " constructing " << VirtualDetector::volumeName(vdId) << endl
               << " at " << vdg->getGlobal(vdId) << endl
               << " at " << vdg->getLocal(vdId) << " w.r.t. parent (DS3Vacuum or DS2Vacuum?) " << endl;
          cout << __func__ << "    VD parameters: " << vdParams << endl;
          cout << __func__ << "    VD rel. posit: " << vdg->getLocal(vdId) << endl;
        }
      }
    }


    /*****************end of new virtual detector************************/


    // Virtual Detectors ST_In, ST_Out are placed inside DS2, just before and after stopping target

    // If there is no neutron absorber, virtual detectors 9 and 10 extend to
    // inner wall of DS2 minus 5 mm.
    // Existence of internal neutron absorber(INA) and/or outer proton absorber(OPA) is checked.
    // Priority on radius determination goes to OPA, and next, INA, and finally DS2.
    // Final radius is the extention to OPA, INA or DS, minus 5 mm.

    if ( !_config.getBool("isDumbbell",false) ){
      double Ravr = ds->rIn1();

      if ( _config.getBool("hasTSdA",false) ) {
        Ravr = _config.getDouble("TSdA.rFactorForVDs");
      }

      bool opaflag = false;
      double opaz0, opaz1, opari0, opari1;
      if ( _config.getBool("hasProtonAbsorber", true) ) {
        GeomHandle<MECOStyleProtonAbsorber> pageom;
        if ( pageom->isAvailable(ProtonAbsorberId::opabs1) ) {
          opaflag = true;
          MECOStyleProtonAbsorberPart opa = pageom->part(2);
          opaz0 = opa.center().z()-opa.halfLength();
          opaz1 = opa.center().z()+opa.halfLength();
          opari0 = opa.innerRadiusAtStart();
          opari1 = opa.innerRadiusAtEnd();
        }
      }

      for( int vdId=VirtualDetectorId::ST_In;
           vdId<=VirtualDetectorId::ST_Out;
           ++vdId) if( vdg->exist(vdId) ) {

          if ( verbosityLevel > 0) {
            cout << __func__ << " constructing " << VirtualDetector::volumeName(vdId)  << endl;
          }

          double zvd = vdg->getGlobal(vdId).z();
          if (opaflag) {
            Ravr = (opari1 - opari0)/(opaz1 - opaz0) * (zvd - opaz0) + opari0;
          }
          double rvd = Ravr - 5.0;

          if ( verbosityLevel > 0) {
            cout << __func__ << " " << VirtualDetector::volumeName(vdId) <<
              " z, r : " << zvd << ", " << rvd << endl;
          }

          TubsParams vdParamsTarget(0.,rvd,vdHalfLength);
          std::string theDS3("DS3Vacuum");
          if ( _config.getBool("inGaragePosition",false) ) theDS3 = "garageFakeDS3Vacuum";

          VolumeInfo const & parent = ( _config.getBool("isDumbbell",false) ) ?
            _helper->locateVolInfo(theDS3) :
            _helper->locateVolInfo("DS2Vacuum"); //DS3Vacuum to move the targets

          if (verbosityLevel >0) {
            cout << __func__ << " " << VirtualDetector::volumeName(vdId) << " Z offset in Mu2e    : " <<
              zvd << endl;
            cout << __func__ << " " << VirtualDetector::volumeName(vdId) << " Z extent in Mu2e    : " <<
              zvd - vdHalfLength << ", " << zvd + vdHalfLength << endl;
          }

          VolumeInfo vd = nestTubs( VirtualDetector::volumeName(vdId),
                                    vdParamsTarget, downstreamVacuumMaterial, 0,
                                    vdg->getLocal(vdId),
                                    parent,
                                    vdId,
                                    vdIsVisible,
                                    G4Color::Red(), vdIsSolid,
                                    forceAuxEdgeVisible,
                                    placePV,
                                    false);

          doSurfaceCheck && checkForOverlaps(vd.physical, _config, verbosityLevel>0);

        }
    }

    if ( _config.getBool("hasTracker",false)  ) {


      // placing virtual detectors in the middle of the tracker

      // check if tracker exists and if the number of planes
      // tracker.numPlanes is even is done in VirtualDetectorMaker

      vdId = VirtualDetectorId::TT_Mid;
      if( vdg->exist(vdId) ) {

        if ( verbosityLevel > 0) {
          cout << __func__ << " constructing " << VirtualDetector::volumeName(vdId)  << endl;
        }

        // the radius of tracker mother
        Tracker const & tracker = *(GeomHandle<Tracker>());
        double orvd = tracker.g4Tracker()->mother().tubsParams().outerRadius();
        double irvd = tracker.g4Tracker()->mother().tubsParams().innerRadius();

        if ( tracker.g4Tracker()->getSupportModel() == SupportModel::detailedv0 ) {
          auto const& beams =  tracker.g4Tracker()->getSupportStructure().beamBody();
          if ( beams.empty() ){
            throw cet::exception("GEOM")
              << "Cannot create virtual detector " << VirtualDetectorId(vdId).name()
              << " unless support beams are defined\n";

          }
          orvd = beams.at(0).tubsParams().innerRadius();
        }


        if ( verbosityLevel > 0) {
          double zvd = vdg->getGlobal(vdId).z();
          cout << __func__  << " " << VirtualDetector::volumeName(vdId) <<
            " z, r : " << zvd << ", " << irvd << " " << orvd << endl;
        }

        TubsParams vdParamsTracker(irvd,orvd,vdHalfLength);

        VolumeInfo const & parent = _helper->locateVolInfo("TrackerMother");

        CLHEP::Hep3Vector vdPos = vdg->getGlobal(vdId)-parent.centerInMu2e();
        //        cout << "foo: TT_Mid: " << vdPos << " " << vdg->getLocal(vdId) << endl;
        //        cout << "foo: TT_Mid: " << vdParamsTracker << endl;

        VolumeInfo vd = nestTubs( VirtualDetector::volumeName(vdId),
                                  vdParamsTracker, downstreamVacuumMaterial, 0,
                                  vdPos,
                                  parent,
                                  vdId, vdIsVisible, G4Color::Red(), vdIsSolid,
                                  forceAuxEdgeVisible,
                                  placePV,
                                  false);

        doSurfaceCheck && checkForOverlaps(vd.physical, _config, verbosityLevel>0);

        vdId = VirtualDetectorId::TT_MidInner;
        if( vdg->exist(vdId) ) {

          if ( verbosityLevel > 0) {
            cout << __func__ << " constructing " << VirtualDetector::volumeName(vdId)  << endl;
          }

          // VD TT_MidInner is placed inside the tracker at the same z position as
          // VD TT_Mid but from radius 0 to the inner radius of the tracker
          // mother volume. However, its mother volume is DS3Vacuum
          // which has a different offset. We will use the global offset
          // here (!) as DS is not in the geometry service yet

          // we need to take into account the "overlap" with the TT_InSurf

          TubsParams vdParamsTrackerInner(0.,irvd-2.*vdHalfLength,vdHalfLength);
          std::string theDS3("DS3Vacuum");
          if ( _config.getBool("inGaragePosition",false) ) theDS3 = "garageFakeDS3Vacuum";
          VolumeInfo const & parent = _helper->locateVolInfo(theDS3);

          G4ThreeVector vdLocalOffset = vdg->getGlobal(vdId) - parent.centerInMu2e();

          if ( verbosityLevel > 0) {
            double zvd = vdg->getGlobal(vdId).z();
            cout << __func__ << " " << VirtualDetector::volumeName(vdId) <<
              " z, r : " << zvd  << ", " << irvd << endl;
          }

          //        cout << "foo: TT_MidInner: " << vdLocalOffset         << endl;
          //        cout << "foo: TT_MidInner: " << vdParamsTrackerInner << endl;

          VolumeInfo vd = nestTubs( VirtualDetector::volumeName(vdId),
                                    vdParamsTrackerInner, downstreamVacuumMaterial, 0,
                                    vdLocalOffset,
                                    parent,
                                    vdId, vdIsVisible, G4Color::Red(), vdIsSolid,
                                    forceAuxEdgeVisible,
                                    placePV,
                                    false);

          doSurfaceCheck && checkForOverlaps(vd.physical, _config, verbosityLevel>0);

        }

      }

      // placing virtual detectors TT_FrontHollow, TT_FrontPA in front
      // of the tracker (in the proton absorber region); check if
      // tracker exist is done in VirtualDetectorMaker

      if (    _config.getBool("hasProtonAbsorber",false)
              && !_config.getBool("protonabsorber.isHelical", false)
              && !_config.getBool("protonabsorber.isShorterCone", false)) {

        // This branch is for the case that the proton absorber penetrates this vd

        vdId = VirtualDetectorId::TT_FrontHollow;
        if( vdg->exist(vdId) ) {

          if ( verbosityLevel > 0) {
            cout << __func__ << " constructing " << VirtualDetector::volumeName(vdId)  << endl;
          }
          if ( !_config.getBool("hasProtonAbsorber",false) ) {
            throw cet::exception("GEOM")
              << "This virtual detector " << VirtualDetectorId(vdId).name()
              << " can only be placed if proton absorber is present\n";
          }

          // the radius of tracker mother
          Tracker const & tracker = *(GeomHandle<Tracker>());
          double orvd = tracker.g4Tracker()->mother().tubsParams().outerRadius();
          double vdZ  = vdg->getGlobal(vdId).z();

          if ( verbosityLevel > 0) {
            cout << __func__ << " " << VirtualDetector::volumeName(vdId) <<
              " z, r : " << vdZ << ", " << orvd << endl;
          }

          // we will create an subtraction solid
          // (we will "subtract" protonAbsorber)
          // and place it (the subtraction solid) in DS3Vacuum

          std::string theDS3("DS3Vacuum");
          if ( _config.getBool("inGaragePosition",false) ) theDS3 = "garageFakeDS3Vacuum";

          VolumeInfo const & parent = _helper->locateVolInfo(theDS3);

          G4ThreeVector vdLocalOffset = vdg->getGlobal(vdId) - parent.centerInMu2e();

          VolumeInfo vdFullInfo;
          vdFullInfo.name = VirtualDetector::volumeName(vdId) + "_FULL";

          TubsParams  vdParamsTrackerFrontFull(0.,orvd,vdHalfLength);

          vdFullInfo.solid = new G4Tubs(vdFullInfo.name,
                                        vdParamsTrackerFrontFull.innerRadius(),
                                        vdParamsTrackerFrontFull.outerRadius(),
                                        vdParamsTrackerFrontFull.zHalfLength(),
                                        vdParamsTrackerFrontFull.phi0(),
                                        vdParamsTrackerFrontFull.phiMax());

          if ( verbosityLevel > 0) {
            cout << __func__ << " constructing " <<  vdFullInfo.name << endl;
          }

          VolumeInfo const & protonabs2Info = _helper->locateVolInfo("protonabs2");

          VolumeInfo vdHollowInfo;
          if ( verbosityLevel > 0) {
            cout << __func__ << " constructing " << VirtualDetector::volumeName(vdId)  << endl;
          }
          vdHollowInfo.name = VirtualDetector::volumeName(vdId);

          // we need to make sure that the vd z is within protonabs2 z
          // in addition the outomatic check if vd is inside ds3Vac is done by G4 itself

          double pabs2Z = protonabs2Info.centerInMu2e()[CLHEP::Hep3Vector::Z];
          double pzhl   = static_cast<G4Cons*>(protonabs2Info.solid)->GetZHalfLength();

          if (verbosityLevel >0) {
            cout << __func__ << " " << protonabs2Info.name << " Z offset in Mu2e    : " <<
              pabs2Z << endl;
            cout << __func__ << " " << protonabs2Info.name << " Z extent in Mu2e    : " <<
              pabs2Z - pzhl  << ", " <<  pabs2Z + pzhl  << endl;
            cout << __func__ << " " << vdFullInfo.name     << " Z offset in Mu2e    : " <<
              vdZ << endl;
            cout << __func__ << " " << vdFullInfo.name     << " Z extent in Mu2e    : " <<
              vdZ - vdHalfLength << ", " << vdZ + vdHalfLength << endl;
          }

          if ( (pabs2Z+pzhl-vdZ-vdHalfLength)<0.0) {
            throw cet::exception("GEOM")
              << "Incorrect positioning of "
              << protonabs2Info.name
              << " and "
              << vdFullInfo.name
              <<"\n";
          }

          // need to find the relative offset of vd & protonAbsorber;
          // we will use the global offsets; they are (in Mu2e)
          // vdg->getGlobal(vdId)
          // protonabs2Info.centerInMu2e()
          // only global offsets can be used for vdet TT_FrontHollow, TT_FrontPA

          if ( verbosityLevel > 0) {
            cout << __func__ << " constructing " <<  vdHollowInfo.name << " name check " << endl;
          }

          vdHollowInfo.solid = new G4SubtractionSolid(vdHollowInfo.name,
                                                      vdFullInfo.solid,
                                                      protonabs2Info.solid,
                                                      0,
                                                      protonabs2Info.centerInMu2e()-
                                                      vdg->getGlobal(vdId));

          vdHollowInfo.centerInParent = vdLocalOffset;
          vdHollowInfo.centerInWorld  = vdHollowInfo.centerInParent + parent.centerInWorld;

          finishNesting(vdHollowInfo,
                        downstreamVacuumMaterial,
                        0,
                        vdLocalOffset,
                        parent.logical,
                        vdId,
                        vdIsVisible,
                        G4Color::Red(),
                        vdIsSolid,
                        forceAuxEdgeVisible,
                        placePV,
                        false);

          doSurfaceCheck && checkForOverlaps(vdHollowInfo.physical, _config, verbosityLevel>0);

          if ( verbosityLevel > 0) {

            // both protonabs2 & vd are placed in DS3Vacuum, do they have proper local offsets?

            double theZ  = vdHollowInfo.centerInMu2e()[CLHEP::Hep3Vector::Z];
            double theHL = static_cast<G4Tubs*>(vdFullInfo.solid)->GetZHalfLength();
            cout << __func__ << " " << vdHollowInfo.name <<
              " Z offset in Mu2e    : " <<
              theZ << endl;
            cout << __func__ << " " << vdHollowInfo.name <<
              " Z extent in Mu2e    : " <<
              theZ - theHL << ", " << theZ + theHL << endl;

            cout << __func__ << " " << vdHollowInfo.name <<
              " local input offset in G4                  : " <<
              vdLocalOffset << endl;
            cout << __func__ << " " << vdHollowInfo.name <<
              " local GetTranslation()       offset in G4 : " <<
              vdHollowInfo.physical->GetTranslation() << endl;

            cout << __func__ << " " << protonabs2Info.name <<
              " local GetTranslation()              offset in G4 : " <<
              protonabs2Info.physical->GetTranslation() << endl;

            cout << __func__ << " " << protonabs2Info.name << " " << vdHollowInfo.name <<
              " local GetTranslation() offset diff in G4 : " <<
              vdHollowInfo.physical->GetTranslation() - protonabs2Info.physical->GetTranslation() <<
              endl;
            cout << __func__ <<
              " protonabs2Info.centerInMu2e() - vdg->getGlobal(vdId) offset           : " <<
              protonabs2Info.centerInMu2e()-vdg->getGlobal(vdId) << endl;
          }

          //  now the complementary solid, it has to be placed in protonabs2

          vdId = VirtualDetectorId::TT_FrontPA;
          if (vdg->exist(vdId)) {
            if ( verbosityLevel > 0) {
              cout << __func__ << " constructing " << VirtualDetector::volumeName(vdId)  << endl;
            }
            VolumeInfo vdIntersectionInfo;
            vdIntersectionInfo.name = VirtualDetector::volumeName(vdId);

            vdIntersectionInfo.solid = new G4IntersectionSolid(vdIntersectionInfo.name + "_INT",
                                                               vdFullInfo.solid,
                                                               protonabs2Info.solid,
                                                               0,
                                                               protonabs2Info.centerInMu2e()-
                                                               vdg->getGlobal(vdId));

            VolumeInfo const & parent = _helper->locateVolInfo("protonabs2");
            vdLocalOffset = vdg->getGlobal(vdId)-protonabs2Info.centerInMu2e();

            vdIntersectionInfo.centerInParent = vdLocalOffset;
            vdIntersectionInfo.centerInWorld  = vdIntersectionInfo.centerInParent +
              parent.centerInWorld;

            finishNesting(vdIntersectionInfo,
                          downstreamVacuumMaterial,
                          0,
                          vdLocalOffset,
                          parent.logical,
                          vdId,
                          vdIsVisible,
                          G4Color::Red(),
                          vdIsSolid,
                          forceAuxEdgeVisible,
                          placePV,
                          false);

            doSurfaceCheck && checkForOverlaps(vdIntersectionInfo.physical, _config, verbosityLevel>0);

            if ( verbosityLevel > 0) {

              // vd is placed in protonabs2

              double theZ  = vdIntersectionInfo.centerInMu2e()[CLHEP::Hep3Vector::Z];
              double theHL = static_cast<G4Tubs*>(vdFullInfo.solid)->GetZHalfLength();
              cout << __func__ << " " << vdIntersectionInfo.name <<
                " Z offset in Mu2e    : " <<
                theZ << endl;
              cout << __func__ << " " << vdIntersectionInfo.name <<
                " Z extent in Mu2e    : " <<
                theZ - theHL << ", " << theZ + theHL << endl;

              cout << __func__ << " " << vdIntersectionInfo.name <<
                " local input offset in G4                  : " <<
                vdLocalOffset << endl;
              cout << __func__ << " " << vdIntersectionInfo.name <<
                " local GetTranslation()       offset in G4 : " <<
                vdIntersectionInfo.physical->GetTranslation() << endl;

              cout << __func__ << " " << protonabs2Info.name <<
                " local GetTranslation()              offset in G4 : " <<
                protonabs2Info.physical->GetTranslation() << endl;

              cout << __func__ << " " << protonabs2Info.name << " " << vdIntersectionInfo.name <<
                " local GetTranslation() offset diff in G4 : " <<
                vdIntersectionInfo.physical->GetTranslation() -
                protonabs2Info.physical->GetTranslation() << endl;
              cout << __func__ <<
                " protonabs2Info.centerInMu2e() - vdg->getGlobal(vdId) offset           : " <<
                protonabs2Info.centerInMu2e()-vdg->getGlobal(vdId) << endl;
            }

          }
        }
      }
      else {

        // If there is no proton absorber that penetrates the front VD, then the VD is a simple disk.
        // Keep the same name even though the "hollow" part of the name is not really appropriate ...

        vdId = VirtualDetectorId::TT_FrontHollow;
        if( vdg->exist(vdId) ) {

          if ( verbosityLevel > 0) {
            cout << __func__ << " constructing " << VirtualDetector::volumeName(vdId)  << endl;
          }

          // the radius of tracker mother
          Tracker const & tracker = *(GeomHandle<Tracker>());
          double orvd = tracker.g4Tracker()->mother().tubsParams().outerRadius();
          double vdZ  = vdg->getGlobal(vdId).z();

          if ( verbosityLevel > 0) {
            cout << __func__ << " " << VirtualDetector::volumeName(vdId) <<
              " z, r : " << vdZ << ", " << orvd << endl;
          }

          std::string theDS3("DS3Vacuum");
          if ( _config.getBool("inGaragePosition",false) ) theDS3 = "garageFakeDS3Vacuum";
          VolumeInfo const & parent = _helper->locateVolInfo(theDS3);

          G4ThreeVector vdLocalOffset = vdg->getGlobal(vdId) - parent.centerInMu2e();

          TubsParams  vdParamsTrackerFrontFull(0.,orvd,vdHalfLength);

          //          cout << "foo: TT_Front: " << vdLocalOffset             << endl;
          //          cout << "foo: TT_Front: " << vdParamsTrackerFrontFull << endl;
          VolumeInfo vdInfo = nestTubs(VirtualDetector::volumeName(vdId),
                                       vdParamsTrackerFrontFull,
                                       downstreamVacuumMaterial,
                                       0,
                                       vdLocalOffset,
                                       parent,
                                       vdId,
                                       vdIsVisible,
                                       G4Color::Red(),
                                       vdIsSolid,
                                       forceAuxEdgeVisible,
                                       placePV,
                                       false);

          doSurfaceCheck && checkForOverlaps(vdInfo.physical, _config, verbosityLevel>0);

        }

      }

      vdId = VirtualDetectorId::TT_Back;
      if( vdg->exist(vdId) ) {

        if ( verbosityLevel > 0) {
          cout << __func__ << " constructing " << VirtualDetector::volumeName(vdId)  << endl;
        }
        // the radius of tracker mother
        Tracker const & tracker = *(GeomHandle<Tracker>());
        double orvd = tracker.g4Tracker()->mother().tubsParams().outerRadius();
        double vdZ  = vdg->getGlobal(vdId).z();

        if ( verbosityLevel > 0) {
          cout << __func__ << " " << VirtualDetector::volumeName(vdId) <<
            " z, r : " << vdZ << ", " << orvd << endl;
        }

        std::string theDS3("DS3Vacuum");
        if ( _config.getBool("inGaragePosition",false) ) theDS3 = "garageFakeDS3Vacuum";
        VolumeInfo const & parent = _helper->locateVolInfo(theDS3);

        G4ThreeVector vdLocalOffset = vdg->getGlobal(vdId) - parent.centerInMu2e();

        TubsParams  vdParamsTrackerBackFull(0.,orvd,vdHalfLength);

        //        cout << "foo: TT_Back: " << vdLocalOffset    << endl;
        //        cout << "foo: TT_Back: " << vdParamsTrackerBackFull << endl;

        VolumeInfo vdInfo = nestTubs(VirtualDetector::volumeName(vdId),
                                     vdParamsTrackerBackFull,
                                     downstreamVacuumMaterial,
                                     0,
                                     vdLocalOffset,
                                     parent,
                                     vdId,
                                     vdIsVisible,
                                     G4Color::Red(),
                                     vdIsSolid,
                                     forceAuxEdgeVisible,
                                     placePV,
                                     false);

        doSurfaceCheck && checkForOverlaps(vdInfo.physical, _config, verbosityLevel>0);

      }

      vdId = VirtualDetectorId::TT_OutSurf;
      if( vdg->exist(vdId) ) {

        if ( verbosityLevel > 0) {
          cout << __func__ << " constructing " << VirtualDetector::volumeName(vdId)  << endl;
        }

        // the radius of tracker mother
        Tracker const & tracker = *(GeomHandle<Tracker>());
        TubsParams const& motherParams = tracker.g4Tracker()->mother().tubsParams();
        double orvd = motherParams.outerRadius();
        double vdZ  = vdg->getGlobal(vdId).z();

        if ( verbosityLevel > 0) {
          cout << __func__ << " " << VirtualDetector::volumeName(vdId) <<
            " z, r : " << vdZ << ", " << orvd << endl;
        }

        std::string theDS3("DS3Vacuum");
        if ( _config.getBool("inGaragePosition",false) ) theDS3 = "garageFakeDS3Vacuum";
        VolumeInfo const & parent = _helper->locateVolInfo(theDS3);

        G4ThreeVector vdLocalOffset = vdg->getGlobal(vdId) - parent.centerInMu2e();

        // the detector is on the outer surface of the tracker envelope
        // it is thin cylinder, NOT a thin disk
        TubsParams  vdParamsTrackerOutSurf(orvd,orvd+2.*vdHalfLength,motherParams.zHalfLength());

        //        cout << "foo: TT_OutSurf: " << vdLocalOffset         << endl;
        //        cout << "foo: TT_OutSurf: " << vdParamsTrackerOutSurf << endl;

        VolumeInfo vdInfo = nestTubs(VirtualDetector::volumeName(vdId),
                                     vdParamsTrackerOutSurf,
                                     downstreamVacuumMaterial,
                                     0,
                                     vdLocalOffset,
                                     parent,
                                     vdId,
                                     vdIsVisible,
                                     G4Color::Red(),
                                     vdIsSolid,
                                     forceAuxEdgeVisible,
                                     placePV,
                                     false);

        doSurfaceCheck && checkForOverlaps(vdInfo.physical, _config, verbosityLevel>0);

      }

      vdId = VirtualDetectorId::TT_InSurf;
      if( vdg->exist(vdId) ) {

        if ( verbosityLevel > 0) {
          cout << __func__ << " constructing " << VirtualDetector::volumeName(vdId)  << endl;
        }

        // the radius of tracker mother
        Tracker const & tracker = *(GeomHandle<Tracker>());
        TubsParams const& motherParams = tracker.g4Tracker()->mother().tubsParams();
        double irvd = motherParams.innerRadius();
        double vdZ  = vdg->getGlobal(vdId).z();

        if ( verbosityLevel > 0) {
          cout << __func__ << " " << VirtualDetector::volumeName(vdId) <<
            " z, r : " << vdZ << ", " << irvd << endl;
        }

        std::string theDS3("DS3Vacuum");
        if ( _config.getBool("inGaragePosition",false) ) theDS3 = "garageFakeDS3Vacuum";
        VolumeInfo const & parent = _helper->locateVolInfo(theDS3);

        G4ThreeVector vdLocalOffset = vdg->getGlobal(vdId) - parent.centerInMu2e();

        // the detector is on the inner surface of the tracker envelope
        // it is thin cylinder, NOT a thin disk
        TubsParams  vdParamsTrackerInSurf(irvd-2.*vdHalfLength,irvd,motherParams.zHalfLength());

        //        cout << "foo: TT_InSurf: " << vdLocalOffset         << endl;
        //        cout << "foo: TT_InSurf: " << vdParamsTrackerInSurf << endl;

        VolumeInfo vdInfo = nestTubs(VirtualDetector::volumeName(vdId),
                                     vdParamsTrackerInSurf,
                                     downstreamVacuumMaterial,
                                     0,
                                     vdLocalOffset,
                                     parent,
                                     vdId,
                                     vdIsVisible,
                                     G4Color::Red(),
                                     vdIsSolid,
                                     forceAuxEdgeVisible,
                                     placePV,
                                     false);

        doSurfaceCheck && checkForOverlaps(vdInfo.physical, _config, verbosityLevel>0);

      }

    } // end hasTracker

    vdId = VirtualDetectorId::EMFC1Entrance;
    if( vdg->exist(vdId) ) {

      if ( verbosityLevel > 0) {
        cout << __func__ << " constructing " << VirtualDetector::volumeName(vdId)  << endl;
      }
      /*
        VolumeInfo const & parent = _helper->locateVolInfo("HallAir");
        GeomHandle<ProtonBeamDump> dump;

        const double vdYmin = dump->frontShieldingCenterInMu2e().y()
        - dump->frontShieldingHalfSize()[1]
        + building->hallFloorThickness()
        ;
        const double vdYmax = std::min(
        dump->frontShieldingCenterInMu2e().y() + dump->frontShieldingHalfSize()[1],
        building->hallInsideYmax()
        );

        std::vector<double> hlen(3);
        hlen[0] = dump->frontShieldingHalfSize()[0];
        hlen[1] = (vdYmax - vdYmin)/2;
        hlen[2] = vdg->getHalfLength();

        // NB: it's not "shielding" center in Y in case the ceiling height is a limitation
        CLHEP::Hep3Vector shieldingFaceCenterInMu2e( (dump->shieldingFaceXmin()+
        dump->shieldingFaceXmax())/2,

        (vdYmax + vdYmin)/2,

        (dump->shieldingFaceZatXmin()+
        dump->shieldingFaceZatXmax())/2
        );

        CLHEP::Hep3Vector vdOffset(dump->coreRotationInMu2e() * CLHEP::Hep3Vector(0, 0, hlen[2]));


        if ( verbosityLevel > 0) {
        std::cout<<"shieldingFaceCenterInMu2e = "<<shieldingFaceCenterInMu2e
        <<", parent.centerInMu2e() = "<<parent.centerInMu2e()
        <<", vdOffset = "<<vdOffset
        <<std::endl;
        }

        VolumeInfo vdInfo = nestBox(VirtualDetector::volumeName(vdId),
        hlen,
        downstreamVacuumMaterial,
        reg.add(dump->coreRotationInMu2e().inverse()),
        shieldingFaceCenterInMu2e + vdOffset - parent.centerInMu2e(),
        parent,
        vdId,
        vdIsVisible,
        G4Color::Red(),
        vdIsSolid,
        forceAuxEdgeVisible,
        placePV,
        false);

        doSurfaceCheck && checkForOverlaps(vdInfo.physical, _config, verbosityLevel>0);

      */
    }


    // placing virtual detector on the exit (beam dump direction) and inside
    // of PS vacuum,  right before the PS enclosure end plate.
    vdId = VirtualDetectorId::PS_FrontExit;
    if ( vdg->exist(vdId) )
      {
        const VolumeInfo& parent = _helper->locateVolInfo("PSVacuum");

        const Tube& psevac = GeomHandle<PSVacuum>()->vacuum();
        GeomHandle<PSEnclosure> pse;
        double zoffset = -psevac.halfLength();
        if(pse->version() > 2) { //move the detector away from the PS endcap
          auto polycone = pse->endPlatePolycone();
          double pseZOffset = pse->getExtraOffset();
          double pseMaxZ = polycone.zPlanes().back() + 0.01; //add an extra 10 um gap
          zoffset += std::max(0., (pseMaxZ + (polycone.originInMu2e().z() + pseZOffset) - (psevac.originInMu2e().z() - psevac.halfLength())));
        }
        TubsParams vdParams(0., psevac.outerRadius(), vdg->getHalfLength());

        G4ThreeVector vdCenterInParent(0., 0., zoffset + vdg->getHalfLength());

        if ( verbosityLevel > 0) {
          cout << __func__ << ": Constructing " << VirtualDetector::volumeName(vdId)
               << " VD center in PS vacuum = " << vdCenterInParent
               << " (zoffset = " << zoffset << ")\n";
        }
        VolumeInfo vdInfo = nestTubs(VirtualDetector::volumeName(vdId),
                                     vdParams,
                                     upstreamVacuumMaterial,
                                     0,
                                     vdCenterInParent,
                                     parent,
                                     vdId,
                                     vdIsVisible,
                                     G4Color::Red(),
                                     vdIsSolid,
                                     forceAuxEdgeVisible,
                                     placePV,
                                     false
                                     );

        doSurfaceCheck && checkForOverlaps(vdInfo.physical, _config, verbosityLevel>0);

      }

    if (_config.getBool("targetPS.hasVD.backward", false)) {
      vdId = VirtualDetectorId::PT_Back;
      if ( vdg->exist(vdId) )
        {
          const VolumeInfo& parent = _helper->locateVolInfo("ProductionTargetMother");
          G4Tubs *PTMoth = static_cast<G4Tubs*>(parent.solid);

          TubsParams vdParams(0., PTMoth->GetOuterRadius(), vdg->getHalfLength());

          G4ThreeVector vdCenterInParent(0., 0., -PTMoth->GetZHalfLength() + vdg->getHalfLength());

          VolumeInfo vdInfo = nestTubs(VirtualDetector::volumeName(vdId),
                                       vdParams,
                                       downstreamVacuumMaterial,
                                       0,
                                       vdCenterInParent,
                                       parent,
                                       vdId,
                                       vdIsVisible,
                                       G4Color::Red(),
                                       vdIsSolid,
                                       forceAuxEdgeVisible,
                                       placePV,
                                       false
                                       );

          doSurfaceCheck && checkForOverlaps(vdInfo.physical, _config, verbosityLevel>0);

        }
    }

    if (_config.getBool("targetPS.hasVD.forward", false)) {
      vdId = VirtualDetectorId::PT_Front;
      if ( vdg->exist(vdId) )
        {
          const VolumeInfo& parent = _helper->locateVolInfo("ProductionTargetMother");
          G4Tubs *PTMoth = static_cast<G4Tubs*>(parent.solid);

          TubsParams vdParams(0., PTMoth->GetOuterRadius(), vdg->getHalfLength());

          G4ThreeVector vdCenterInParent(0., 0., PTMoth->GetZHalfLength() - vdg->getHalfLength());

          VolumeInfo vdInfo = nestTubs(VirtualDetector::volumeName(vdId),
                                       vdParams,
                                       downstreamVacuumMaterial,
                                       0,
                                       vdCenterInParent,
                                       parent,
                                       vdId,
                                       vdIsVisible,
                                       G4Color::Red(),
                                       vdIsSolid,
                                       forceAuxEdgeVisible,
                                       placePV,
                                       false
                                       );

          doSurfaceCheck && checkForOverlaps(vdInfo.physical, _config, verbosityLevel>0);

        }
    }

    // A plane in front of the beam dump and ExtMon collimator entrance
    vdId = VirtualDetectorId::ExtMonCommonPlane;
    if( vdg->exist(vdId) ) {

      if ( verbosityLevel > 0) {
        cout << __func__ << " constructing " << VirtualDetector::volumeName(vdId)  << endl;
      }

      checkForStale( "vd.ExtMonCommonPlane.z", _config);

      const double dz = _config.getDouble("vd.ExtMonCommonPlane.beamDumpFaceDistance");
      const double dxmin = _config.getDouble("vd.ExtMonCommonPlane.dxmin");
      const double dxmax = _config.getDouble("vd.ExtMonCommonPlane.dxmax");
      const double dymin = _config.getDouble("vd.ExtMonCommonPlane.dymin");
      const double dymax = _config.getDouble("vd.ExtMonCommonPlane.dymax");

      VolumeInfo const & parent = _helper->locateVolInfo("HallAir");

      const double halfThick = vdg->getHalfLength();

      GeomHandle<ProtonBeamDump> dump;
      CLHEP::Hep3Vector vdCenterInMu2e =
        dump->mouthCenterInMu2e()
        + dump->coreRotationInMu2e()
        * CLHEP::Hep3Vector((dxmax+dxmin)/2,
                            (dymax+dymin)/2,
                            dump->mouthHalfSize()[2]+dz+halfThick
                            );

      CLHEP::Hep3Vector vdCenterInParent = vdCenterInMu2e - parent.centerInMu2e();


      std::vector<double> hlen(3);
      hlen[0] = (dxmax - dxmin)/2;
      hlen[1] = (dymax - dymin)/2;
      hlen[2] = halfThick;

      VolumeInfo vdInfo = nestBox(VirtualDetector::volumeName(vdId),
                                  hlen,
                                  upstreamVacuumMaterial,
                                  reg.add(dump->coreRotationInMu2e().inverse()),
                                  vdCenterInParent,
                                  parent,
                                  vdId,
                                  vdIsVisible,
                                  G4Color::Red(),
                                  vdIsSolid,
                                  forceAuxEdgeVisible,
                                  placePV,
                                  false
                                  );

      doSurfaceCheck && checkForOverlaps(vdInfo.physical, _config, verbosityLevel>0);
    }

    // placing virtual detector at the dump core face
    vdId = VirtualDetectorId::ProtonBeamDumpCoreFace;
    if( vdg->exist(vdId) ) {
      if ( verbosityLevel > 0) {
        cout << __func__ << " constructing " << VirtualDetector::volumeName(vdId)  << endl;
      }

      VolumeInfo const & parent = _helper->locateVolInfo("ProtonBeamDumpCore");

      GeomHandle<ProtonBeamDump> dump;

      CLHEP::Hep3Vector centerInCore(0, 0, dump->coreHalfSize()[2] - vdg->getHalfLength());

      std::vector<double> hlen(3);
      hlen[0] = dump->coreHalfSize()[0];
      hlen[1] = dump->coreHalfSize()[1];
      hlen[2] = vdg->getHalfLength();

      VolumeInfo vdInfo = nestBox(VirtualDetector::volumeName(vdId),
                                  hlen,
                                  downstreamVacuumMaterial,
                                  0,
                                  centerInCore,
                                  parent,
                                  vdId,
                                  vdIsVisible,
                                  G4Color::Red(),
                                  vdIsSolid,
                                  forceAuxEdgeVisible,
                                  placePV,
                                  false
                                  );

      doSurfaceCheck && checkForOverlaps(vdInfo.physical, _config, verbosityLevel>0);

    }

    // ExtMonFNAL detector VDs - created in constructExtMonFNAL()

    // placing virtual detector at exit of DS neutron shielding
    vdId = VirtualDetectorId::DSNeutronShieldExit;
    if( vdg->exist(vdId) ) {
      if ( verbosityLevel > 0) {
        cout << __func__ << " constructing " << VirtualDetector::volumeName(vdId)  << endl;
      }

      VolumeInfo const & parent = _helper->locateVolInfo("HallAir");
      CLHEP::Hep3Vector const& parentInMu2e = parent.centerInMu2e();

      // The following code was built on the Geometry #13 description of
      // External shielding.  To break that dependency, I have created
      // a work-around that is currently meant to be a temporary solution.
      // David Norvil Brown, December 2014.


      double shieldHoleRadius = _config.getDouble("ExtShieldDownstream.detecHoleRadius")*CLHEP::mm;

      const TubsParams vdParams(0, shieldHoleRadius, vdg->getHalfLength());

      VolumeInfo vdInfo = nestTubs(VirtualDetector::volumeName(vdId),
                                   vdParams,
                                   downstreamVacuumMaterial,
                                   0,
                                   vdg->getGlobal(vdId) - parentInMu2e, //position w.r.t. parent
                                   parent,
                                   vdId,
                                   vdIsVisible,
                                   G4Color::Red(),
                                   vdIsSolid,
                                   forceAuxEdgeVisible,
                                   placePV,
                                   false
                                   );

      doSurfaceCheck && checkForOverlaps(vdInfo.physical, _config, verbosityLevel>0);

      if ( verbosityLevel > 0) {
        cout << __func__ << " constructing " << VirtualDetector::volumeName(vdId) << endl
             << " at " << vdg->getGlobal(vdId) << endl
             << " at " << vdg->getGlobal(vdId) - parentInMu2e <<" or local "<< vdg->getLocal(vdId)<< " w.r.t. parent (HallAir) " << endl;
        cout << __func__ << "    VD parameters: " << vdParams << endl;
        cout << __func__ << "    VD rel. posit: " << vdg->getLocal(vdId) << endl;
      }
    }



    if ( _config.getBool("hasDiskCalorimeter",true)) {

        int vdIdDiskEdge = VirtualDetectorId::EMC_Disk_0_EdgeIn;
        int vdIdDiskSurf = VirtualDetectorId::EMC_Disk_0_SurfIn;
        int vdIdFEBEdge  = VirtualDetectorId::EMC_FEB_0_EdgeIn;
        int vdIdFEBSurf  = VirtualDetectorId::EMC_FEB_0_SurfIn;

        double vdgThick            = 2.0*vdg->getHalfLength();
        DiskCalorimeter const& cal = *(GeomHandle<DiskCalorimeter>());

        for (size_t id = 0; id < cal.nDisks(); id++){

            std::ostringstream diskname; diskname<<"CaloDisk_" <<id;
            const VolumeInfo& caloDisk = _helper->locateVolInfo(diskname.str());
            G4Tubs* disk  = static_cast<G4Tubs*>(caloDisk.logical->GetSolid());

            double diskRadIn  = disk->GetInnerRadius();
            double diskRadOut = disk->GetOuterRadius();
            double diskHalfZ  = disk->GetZHalfLength();

            TubsParams  vdParamsFrontDisk(diskRadIn,           diskRadOut,         vdgThick/2.0);
            TubsParams  vdParamsInnerDisk(diskRadIn,           diskRadIn+vdgThick, diskHalfZ-vdgThick);
            TubsParams  vdParamsOuterDisk(diskRadOut-vdgThick, diskRadOut,         diskHalfZ-vdgThick);

            G4ThreeVector posFrontDisk(0,0,-diskHalfZ+vdgThick/2.0);
            G4ThreeVector posBackDisk(0,0,diskHalfZ-vdgThick/2.0);
            G4ThreeVector posInnerDisk(0,0,0);

            if( vdg->exist(vdIdDiskSurf)){
              VolumeInfo vdInfo = nestTubs(VirtualDetector::volumeName(vdIdDiskSurf),
                                           vdParamsFrontDisk,downstreamVacuumMaterial,0,
                                           posFrontDisk,caloDisk,
                                           vdIdDiskSurf,vdIsVisible,G4Color::Red(),vdIsSolid,forceAuxEdgeVisible,
                                           placePV,false);
              ++vdIdDiskSurf;

              VolumeInfo vdInfo2 = nestTubs(VirtualDetector::volumeName(vdIdDiskSurf),
                                            vdParamsFrontDisk,downstreamVacuumMaterial,0,
                                            posBackDisk,caloDisk,
                                            vdIdDiskSurf,vdIsVisible,G4Color::Red(),vdIsSolid,forceAuxEdgeVisible,
                                            placePV,false);
              ++vdIdDiskSurf;

              doSurfaceCheck && checkForOverlaps(vdInfo.physical, _config, verbosityLevel>0);
              doSurfaceCheck && checkForOverlaps(vdInfo2.physical, _config, verbosityLevel>0);
            }


            if( vdg->exist(vdIdDiskEdge)){
              VolumeInfo vdInfo = nestTubs(VirtualDetector::volumeName(vdIdDiskEdge),
                                           vdParamsInnerDisk,downstreamVacuumMaterial,0,
                                           posInnerDisk,caloDisk,
                                           vdIdDiskSurf,vdIsVisible,G4Color::Red(),vdIsSolid,forceAuxEdgeVisible,
                                           placePV,false);
              ++vdIdDiskEdge;

              //needed to maintain the consistence with the numbering scheme
              VolumeInfo vdInfo2 = nestTubs(VirtualDetector::volumeName(vdIdDiskEdge),
                                            vdParamsOuterDisk,downstreamVacuumMaterial,0,
                                            posInnerDisk,caloDisk,
                                            vdIdDiskSurf,1,G4Color::Red(),vdIsSolid,forceAuxEdgeVisible,
                                            placePV,false);
              ++vdIdDiskEdge;

              doSurfaceCheck && checkForOverlaps(vdInfo.physical, _config, verbosityLevel>0);
              doSurfaceCheck && checkForOverlaps(vdInfo2.physical, _config, verbosityLevel>0);
            }


            std::ostringstream cratename; cratename<<"CaloFEB_" <<id;
            VolumeInfo const& caloFEBParent = _helper->locateVolInfo(cratename.str());
            G4Tubs* crate  = static_cast<G4Tubs*>(caloFEBParent.logical->GetSolid());

            double crateRadIn    = crate->GetInnerRadius();
            double crateRadOut   = crate->GetOuterRadius();
            double crateHalfZ    = crate->GetZHalfLength();
            double crateStartPhi = crate->GetStartPhiAngle();
            double crateDphi     = crate->GetDeltaPhiAngle();

            TubsParams  vdParamsFrontFEB(crateRadIn,           crateRadOut,         vdgThick/2.0,        crateStartPhi, crateDphi);
            TubsParams  vdParamsInnerFEB(crateRadIn,           crateRadIn+vdgThick, crateHalfZ-vdgThick, crateStartPhi, crateDphi);
            TubsParams  vdParamsOuterFEB(crateRadOut-vdgThick, crateRadOut,         crateHalfZ-vdgThick, crateStartPhi, crateDphi);

            G4ThreeVector posFrontFEB(0,0,-crateHalfZ+vdgThick/2.0);
            G4ThreeVector posBackFEB(0,0,crateHalfZ-vdgThick/2.0);
            G4ThreeVector posInnerFEB(0,0,0);
            G4ThreeVector posOuterFEB(0,0,0);


            if( vdg->exist(vdIdFEBSurf) ){
               VolumeInfo vdInfo = nestTubs(VirtualDetector::volumeName(vdIdFEBSurf),
                                            vdParamsFrontFEB,downstreamVacuumMaterial,0,
                                            posFrontFEB,caloFEBParent,
                                            vdIdFEBSurf,vdIsVisible,G4Color::Red(),vdIsSolid,forceAuxEdgeVisible,
                                            placePV,false);
               ++vdIdFEBSurf;

               VolumeInfo vdInfo2 = nestTubs(VirtualDetector::volumeName(vdIdFEBSurf),
                                             vdParamsFrontFEB,downstreamVacuumMaterial,0,
                                             posBackFEB,caloFEBParent,
                                             vdIdFEBSurf,vdIsVisible,G4Color::Red(),vdIsSolid,forceAuxEdgeVisible,
                                             placePV,false);


               ++vdIdFEBSurf;

               doSurfaceCheck && checkForOverlaps(vdInfo.physical, _config, verbosityLevel>0);
               doSurfaceCheck && checkForOverlaps(vdInfo2.physical, _config, verbosityLevel>0);
            }

            if( vdg->exist(vdIdFEBEdge)){
               VolumeInfo vdInfo = nestTubs(VirtualDetector::volumeName(vdIdFEBEdge),
                                            vdParamsInnerFEB,downstreamVacuumMaterial,0,
                                            posInnerFEB,caloFEBParent,
                                            vdIdFEBEdge,vdIsVisible,G4Color::Red(),vdIsSolid,forceAuxEdgeVisible,
                                            placePV,false);
               ++vdIdFEBEdge;

               VolumeInfo vdInfo2 = nestTubs(VirtualDetector::volumeName(vdIdFEBEdge),
                                             vdParamsOuterFEB,downstreamVacuumMaterial,0,
                                             posOuterFEB,caloFEBParent,
                                             vdIdFEBEdge,vdIsVisible,G4Color::Red(),vdIsSolid,forceAuxEdgeVisible,
                                             placePV,false);
               ++vdIdFEBEdge;

               doSurfaceCheck && checkForOverlaps(vdInfo.physical, _config, verbosityLevel>0);
               doSurfaceCheck && checkForOverlaps(vdInfo2.physical, _config, verbosityLevel>0);
            }


          }
      }


    //-----------------------------------------------------------------------------------------------------------------------------
    // placing virtual detector in the MSTM Mother
    // FIXME get MSTM from GeometryService once in there
    vdId = VirtualDetectorId::STM_UpStr;
    if ( vdg->exist(vdId) ) {

      //const VolumeInfo& parent = _helper->locateVolInfo("MSTMMother");
      const VolumeInfo& parent = _helper->locateVolInfo("HallAir");
      GeomHandle<CosmicRayShield> CRS;
      //const double y_crv_max       = CRS->getSectorPosition("D").y() + (CRS->getSectorHalfLengths("D"))[1];
      const double yExtentLow      = std::abs(_config.getDouble("yOfFloorSurface.below.mu2eOrigin") );
      const double x_vd_halflength = (CRS->getSectorHalfLengths("D"))[0];
      //const double y_vd_halflength = (y_crv_max + yExtentLow)/2.0;
      const double y_mother_halflength = yExtentLow;
      const double dimVD[3] = { x_vd_halflength, y_mother_halflength, vdg->getHalfLength() };

      VolumeInfo vdFullInfo;
      vdFullInfo.solid = new G4Box("STM_UpStr_Full",
                                   dimVD[0],dimVD[1],dimVD[2]);

      //Subtracting hole from original VD86
      ///////////////////////////////////////////////////////////////////////////////
      VolumeInfo HoleInfo;
      VolumeInfo vdHollowInfo;
      //double shieldHoleRadius = _config.getDouble("ExtShieldDownstream.detecHoleRadius")*CLHEP::mm;
      double shieldHoleRadius = _config.getDouble("stm.shield.rOut")*CLHEP::mm;
      HoleInfo.solid = new G4Tubs(vdFullInfo.name,
                                  0,
                                  shieldHoleRadius,
                                  2*vdg->getHalfLength(),
                                  0*CLHEP::deg,
                                  360*CLHEP::deg);


      G4ThreeVector vdLocalOffset = vdg->getGlobal(vdId) - parent.centerInMu2e();

      vdHollowInfo.name = VirtualDetector::volumeName(vdId);
      vdHollowInfo.solid = new G4SubtractionSolid(vdHollowInfo.name,
                                                  vdFullInfo.solid,
                                                  HoleInfo.solid,
                                                  0,
                                                  /*HoleInfo.centerInMu2e()-
                                                    vdg->getGlobal(vdId)*/G4ThreeVector(0,0,0));

    vdHollowInfo.centerInParent = vdLocalOffset;
    vdHollowInfo.centerInWorld  = vdHollowInfo.centerInParent + parent.centerInWorld;


    finishNesting(vdHollowInfo,
                  downstreamVacuumMaterial,
                  0,
                  vdLocalOffset,
                  parent.logical,
                  vdId,
                  vdIsVisible,
                  G4Color::Red(),
                  vdIsSolid,
                  forceAuxEdgeVisible,
                  placePV,
                  false);

      if ( verbosityLevel > 0) {
        cout << __func__ << " constructing " << VirtualDetector::volumeName(vdId) << endl
             << " at " << vdg->getGlobal(vdId) << endl
             << " at " << vdg->getLocal(vdId) << " w.r.t. parent (HallAir) " << endl;
        cout << __func__ << "    VD parameters: " << vdParams << endl;
        cout << __func__ << "    VD rel. posit: " << vdg->getLocal(vdId) << endl;
      }

    doSurfaceCheck && checkForOverlaps(vdHollowInfo.physical, _config, verbosityLevel>0);
  }

    vdId = VirtualDetectorId::STM_FieldOfViewCollDnStr;
    if ( vdg->exist(vdId) ) {

      const VolumeInfo& parent = _helper->locateVolInfo("HallAir");
      GeomHandle<CosmicRayShield> CRS;
      //const double y_crv_max       = CRS->getSectorPosition("D").y() + (CRS->getSectorHalfLengths("D"))[1];
      const double yExtentLow      = std::abs(_config.getDouble("yOfFloorSurface.below.mu2eOrigin") );
      const double x_vd_halflength = (CRS->getSectorHalfLengths("D"))[0];
      //const double y_vd_halflength = (y_crv_max + yExtentLow)/2.0;
      const double y_mother_halflength = yExtentLow;
      const double dimVD[3] = { x_vd_halflength, y_mother_halflength, vdg->getHalfLength() };

      VolumeInfo vdInfo = nestBox(VirtualDetector::volumeName(vdId),
                                  dimVD,
                                  downstreamVacuumMaterial,
                                  0,                              //rotation
                                  vdg->getLocal(vdId),
                                  parent,
                                  vdId,
                                  vdIsVisible,
                                  G4Color::White(),
                                  vdIsSolid,
                                  forceAuxEdgeVisible,
                                  placePV,
                                  false
                                  );

      if ( verbosityLevel > 0) {
        cout << __func__ << " constructing " << VirtualDetector::volumeName(vdId) << endl
             << " at " << vdg->getGlobal(vdId) << endl
             << " at " << vdg->getLocal(vdId) << " w.r.t. parent (HallAir) " << endl;
        cout << __func__ << "    VD parameters: " << vdParams << endl;
        cout << __func__ << "    VD rel. posit: " << vdg->getLocal(vdId) << endl;
      }
      doSurfaceCheck && checkForOverlaps(vdInfo.physical, _config, verbosityLevel>0);
    }

    vdId = VirtualDetectorId::STM_MagDnStr;
    if ( vdg->exist(vdId) ) {

      //const VolumeInfo& parent = _helper->locateVolInfo("MSTMMother");
      const VolumeInfo& parent = _helper->locateVolInfo("HallAir");
      GeomHandle<CosmicRayShield> CRS;
      //const double y_crv_max       = CRS->getSectorPosition("D").y() + (CRS->getSectorHalfLengths("D"))[1];
      double yExtentLow      = 0.0;
      double x_vd_halflength = 0.0;
      //double y_vd_halflength = (y_crv_max + yExtentLow)/2.0;

      if (_config.getBool("stm.magnet.build",false)){
        yExtentLow          = _config.getDouble("stm.magnet.holeHalfHeight");
        x_vd_halflength     = _config.getDouble("stm.magnet.holeHalfWidth");
      } else {
        yExtentLow      = std::abs(_config.getDouble("yOfFloorSurface.below.mu2eOrigin") );
        x_vd_halflength = (CRS->getSectorHalfLengths("D"))[0];
      }
      const double y_mother_halflength = yExtentLow;

      const double dimVD[3] = { x_vd_halflength, y_mother_halflength, vdg->getHalfLength() };

      VolumeInfo vdInfo = nestBox(VirtualDetector::volumeName(vdId),
                                  dimVD,
                                  downstreamVacuumMaterial,
                                  0,                              //rotation
                                  vdg->getLocal(vdId),
                                  parent,
                                  vdId,
                                  vdIsVisible,
                                  G4Color::White(),
                                  vdIsSolid,
                                  forceAuxEdgeVisible,
                                  placePV,
                                  false
                                  );

      if ( verbosityLevel > 0) {
        cout << __func__ << " constructing " << VirtualDetector::volumeName(vdId) << endl
             << " at " << vdg->getGlobal(vdId) << endl
             << " at " << vdg->getLocal(vdId) << " w.r.t. parent (HallAir) " << endl;
        cout << __func__ << "    VD parameters: " << vdParams << endl;
        cout << __func__ << "    VD rel. posit: " << vdg->getLocal(vdId) << endl;
      }
      doSurfaceCheck && checkForOverlaps(vdInfo.physical, _config, verbosityLevel>0);
    }

    vdId = VirtualDetectorId::STM_SpotSizeCollUpStr;
    if ( vdg->exist(vdId) ) {
      //const VolumeInfo& parent = _helper->locateVolInfo("MSTMMother");
      const VolumeInfo& parent = _helper->locateVolInfo("stmDownstreamEnvelope");
      const double vdRIn  = 0.0;
      const double vdROut = _config.getDouble("vd.STMSSCollUpStr.r");
      const TubsParams vdParams(vdRIn, vdROut, vdg->getHalfLength());

      VolumeInfo vdInfo = nestTubs(VirtualDetector::volumeName(vdId),
                                   vdParams,
                                   downstreamVacuumMaterial,
                                   0,
                                   vdg->getLocal(vdId), //local position w.r.t. parent
                                   parent,
                                   vdId,
                                   vdIsVisible, //
                                   G4Color::White(),
                                   vdIsSolid,
                                   forceAuxEdgeVisible,
                                   placePV,
                                   false
                                   );

      if ( verbosityLevel > 0) {
        cout << __func__ << " constructing " << VirtualDetector::volumeName(vdId) << endl
             << " at " << vdg->getGlobal(vdId) << endl
             << " at " << vdg->getLocal(vdId) << " w.r.t. parent (HallAir) " << endl;
        cout << __func__ << "    VD parameters: " << vdParams << endl;
        cout << __func__ << "    VD rel. posit: " << vdg->getLocal(vdId) << endl;
      }
      doSurfaceCheck && checkForOverlaps(vdInfo.physical, _config, verbosityLevel>0);
    }

    vdId = VirtualDetectorId::STM_CollDnStr;
    if ( vdg->exist(vdId) ) {
      //const VolumeInfo& parent = _helper->locateVolInfo("MSTMMother");
      const VolumeInfo& parent = _helper->locateVolInfo("stmDownstreamEnvelope");
      const double vdRIn  = 0.0;
      const double vdROut = _config.getDouble("vd.STMCollDnStr.r");
      const TubsParams vdParams(vdRIn, vdROut, vdg->getHalfLength());

      VolumeInfo vdInfo = nestTubs(VirtualDetector::volumeName(vdId),
                                   vdParams,
                                   downstreamVacuumMaterial,
                                   0,
                                   vdg->getLocal(vdId), //local position w.r.t. parent
                                   parent,
                                   vdId,
                                   vdIsVisible, //
                                   G4Color::White(),
                                   vdIsSolid,
                                   forceAuxEdgeVisible,
                                   placePV,
                                   false
                                   );

      if ( verbosityLevel > 0) {
        cout << __func__ << " constructing " << VirtualDetector::volumeName(vdId) << endl
             << " at " << vdg->getGlobal(vdId) << endl
             << " at " << vdg->getLocal(vdId) << " w.r.t. parent (HallAir) " << endl;
        cout << __func__ << "    VD parameters: " << vdParams << endl;
        cout << __func__ << "    VD rel. posit: " << vdg->getLocal(vdId) << endl;
      }
      doSurfaceCheck && checkForOverlaps(vdInfo.physical, _config, verbosityLevel>0);
    }


    vdId = VirtualDetectorId::STM_Det1UpStr;
    if ( vdg->exist(vdId) ) {
      //const VolumeInfo& parent = _helper->locateVolInfo("MSTMMother");
      const VolumeInfo& parent = _helper->locateVolInfo("stmDownstreamEnvelope");
      const double vdRIn  = 0.0;
      const double vdROut = _config.getDouble("stm.det1.rOut");
      const TubsParams vdParams(vdRIn, vdROut, vdg->getHalfLength());

      VolumeInfo vdInfo = nestTubs(VirtualDetector::volumeName(vdId),
                                   vdParams,
                                   downstreamVacuumMaterial,
                                   0,
                                   vdg->getLocal(vdId), //local position w.r.t. parent
                                   parent,
                                   vdId,
                                   vdIsVisible, //
                                   G4Color::White(),
                                   vdIsSolid,
                                   forceAuxEdgeVisible,
                                   placePV,
                                   false
                                   );

      if ( verbosityLevel > 0) {
        cout << __func__ << " constructing " << VirtualDetector::volumeName(vdId) << endl
             << " at " << vdg->getGlobal(vdId) << endl
             << " at " << vdg->getLocal(vdId) << " w.r.t. parent (HallAir) " << endl;
        cout << __func__ << "    VD parameters: " << vdParams << endl;
        cout << __func__ << "    VD rel. posit: " << vdg->getLocal(vdId) << endl;
      }
      doSurfaceCheck && checkForOverlaps(vdInfo.physical, _config, verbosityLevel>0);
    }


    vdId = VirtualDetectorId::STM_Det2UpStr;
    if ( vdg->exist(vdId) ) {
      //const VolumeInfo& parent = _helper->locateVolInfo("MSTMMother");
      const VolumeInfo& parent = _helper->locateVolInfo("stmDownstreamEnvelope");
      const double vdRIn  = 0.0;
      const double vdROut = _config.getDouble("stm.det2.rOut");
      const TubsParams vdParams(vdRIn, vdROut, vdg->getHalfLength());

      VolumeInfo vdInfo = nestTubs(VirtualDetector::volumeName(vdId),
                                   vdParams,
                                   downstreamVacuumMaterial,
                                   0,
                                   vdg->getLocal(vdId), //local position w.r.t. parent
                                   parent,
                                   vdId,
                                   vdIsVisible, //
                                   G4Color::White(),
                                   vdIsSolid,
                                   forceAuxEdgeVisible,
                                   placePV,
                                   false
                                   );

      if ( verbosityLevel > 0) {
        cout << __func__ << " constructing " << VirtualDetector::volumeName(vdId) << endl
             << " at " << vdg->getGlobal(vdId) << endl
             << " at " << vdg->getLocal(vdId) << " w.r.t. parent (HallAir) " << endl;
        cout << __func__ << "    VD parameters: " << vdParams << endl;
        cout << __func__ << "    VD rel. posit: " << vdg->getLocal(vdId) << endl;
      }
      doSurfaceCheck && checkForOverlaps(vdInfo.physical, _config, verbosityLevel>0);
    }


    vdId = VirtualDetectorId::PSPbarIn;
    if ( vdg->exist(vdId) ) {

      const VolumeInfo& parentTS1 = _helper->locateVolInfo("TS1Vacuum");
      const VolumeInfo& parentPS = _helper->locateVolInfo("PSVacuum");

      double pbarTS1InOffset = _config.getDouble("pbar.coll1In.offset", 1.0);
      const double pbarTS1InRecordROut = _config.getDouble("pbar.coll1In.rOutRecord");
      const TubsParams vdParams(0.0, pbarTS1InRecordROut, vdg->getHalfLength());

      VolumeInfo vdInfo;
      if (pbarTS1InOffset < 0.0) {
        vdInfo = nestTubs(VirtualDetector::volumeName(vdId),
                          vdParams,
                          upstreamVacuumMaterial,
                          0,
                          vdg->getLocal(vdId), //local position w.r.t. parent
                          parentPS,
                          vdId,
                          vdIsVisible, //
                          G4Color::White(),
                          vdIsSolid,
                          forceAuxEdgeVisible,
                          placePV,
                          false
                          );
      }
      else {
        vdInfo = nestTubs(VirtualDetector::volumeName(vdId),
                          vdParams,
                          downstreamVacuumMaterial,
                          0,
                          vdg->getLocal(vdId), //local position w.r.t. parent
                          parentTS1,
                          vdId,
                          vdIsVisible, //
                          G4Color::White(),
                          vdIsSolid,
                          forceAuxEdgeVisible,
                          placePV,
                          false
                          );
      }

      if ( verbosityLevel > 0) {
        cout << __func__ << " constructing " << VirtualDetector::volumeName(vdId) << endl
             << " at " << vdg->getGlobal(vdId) << endl
             << " at " << vdg->getLocal(vdId) << " w.r.t. parent (PSVacuum) " << endl;
        cout << __func__ << "    VD parameters: " << vdParams << endl;
        cout << __func__ << "    VD rel. posit: " << vdg->getLocal(vdId) << endl;
        //
      }
      doSurfaceCheck && checkForOverlaps(vdInfo.physical, _config, verbosityLevel>0);

    }

    vdId = VirtualDetectorId::PSPbarOut;
    if ( vdg->exist(vdId) ) {

      const VolumeInfo& parentTS1 = _helper->locateVolInfo("TS1Vacuum");
      const VolumeInfo& parentPS = _helper->locateVolInfo("PSVacuum");

      double pbarTS1InOffset = _config.getDouble("pbar.coll1In.offset", 1.0);
      const double pbarTS1InRecordROut = _config.getDouble("pbar.coll1In.rOutRecord");
      const TubsParams vdParams(0.0, pbarTS1InRecordROut, vdg->getHalfLength());

      VolumeInfo vdInfo;
      if (pbarTS1InOffset < 0.0) {
        vdInfo = nestTubs(VirtualDetector::volumeName(vdId),
                          vdParams,
                          upstreamVacuumMaterial,
                          0,
                          vdg->getLocal(vdId), //local position w.r.t. parent
                          parentPS,
                          vdId,
                          vdIsVisible, //
                          G4Color::White(),
                          vdIsSolid,
                          forceAuxEdgeVisible,
                          placePV,
                          false
                          );
      }
      else {
        vdInfo = nestTubs(VirtualDetector::volumeName(vdId),
                          vdParams,
                          downstreamVacuumMaterial,
                          0,
                          vdg->getLocal(vdId), //local position w.r.t. parent
                          parentTS1,
                          vdId,
                          vdIsVisible, //
                          G4Color::White(),
                          vdIsSolid,
                          forceAuxEdgeVisible,
                          placePV,
                          false
                          );
      }

      if ( verbosityLevel > 0) {
        cout << __func__ << " constructing " << VirtualDetector::volumeName(vdId) << endl
             << " at " << vdg->getGlobal(vdId) << endl
             << " at " << vdg->getLocal(vdId) << " w.r.t. parent (PSVacuum) " << endl;
        cout << __func__ << "    VD parameters: " << vdParams << endl;
        cout << __func__ << "    VD rel. posit: " << vdg->getLocal(vdId) << endl;
      }

      doSurfaceCheck && checkForOverlaps(vdInfo.physical, _config, verbosityLevel>0);

    }

    GeomHandle<CosmicRayShield> CRS;
    for(int vdId=VirtualDetectorId::CRV_R; vdId<=VirtualDetectorId::CRV_U; vdId++)
      {

        if(vdg->exist(vdId))
          {
            std::string vdName;
            vdName = VirtualDetector::volumeName(vdId).back();
            std::vector<double> halfLengths = CRS->getSectorHalfLengths(vdName);
            const CLHEP::Hep3Vector vdDirection = _config.getHep3Vector("crs.vdDirection"+vdName);
            for(int i=0; i<3; i++) {if(vdDirection[i]!=0) halfLengths[i]=vdg->getHalfLength();}

            VolumeInfo const &parent = _helper->locateVolInfo("HallAir");
            G4Material* hallAirMaterial = parent.logical->GetMaterial();

            if(verbosityLevel > 0)
              {
                cout << __func__ << " constructing " << VirtualDetector::volumeName(vdId) << " at " << vdg->getGlobal(vdId) << endl;
                cout << __func__ << "    VD half lengths: (" << halfLengths[0]<<","<<halfLengths[1]<<","<<halfLengths[2]<<")" << endl;
                cout << __func__ << "    VD rel. position: " << vdg->getLocal(vdId) << endl;
              }

            VolumeInfo vdInfo = nestBox(VirtualDetector::volumeName(vdId),
                                        halfLengths,
                                        hallAirMaterial,
                                        0,
                                        vdg->getLocal(vdId),
                                        parent,
                                        vdId,
                                        vdIsVisible,
                                        G4Color::Red(),
                                        vdIsSolid,
                                        forceAuxEdgeVisible,
                                        placePV,
                                        false
                                        );

            doSurfaceCheck && checkForOverlaps(vdInfo.physical, _config, verbosityLevel>0);

          }
      }

    vdId = VirtualDetectorId::STM_Final;
    if ( vdg->exist(vdId) ) {
      const VolumeInfo& parent = _helper->locateVolInfo("HallAir");
      const double vdRIn  = 0.0;
      const double vdROut = _config.getDouble("vd.STMFin.r");
      const TubsParams vdParams(vdRIn, vdROut, vdg->getHalfLength());
      VolumeInfo vdInfo = nestTubs(VirtualDetector::volumeName(vdId),
                                   vdParams,
                                   downstreamVacuumMaterial,
                                   0,
                                   vdg->getLocal(vdId), //local position w.r.t. parent
                                   parent,
                                   vdId,
                                   vdIsVisible, //
                                   G4Color::White(),
                                   vdIsSolid,
                                   forceAuxEdgeVisible,
                                   placePV,
                                   false
                                   );
      if ( verbosityLevel > 0) {
        cout << __func__ << " constructing " << VirtualDetector::volumeName(vdId) << endl
             << " at " << vdg->getGlobal(vdId) << endl
             << " at " << vdg->getLocal(vdId) << " w.r.t. parent (HallAir) " << endl;
        cout << __func__ << "    VD parameters: " << vdParams << endl;
        cout << __func__ << "    VD rel. posit: " << vdg->getLocal(vdId) << endl;
      }
      doSurfaceCheck && checkForOverlaps(vdInfo.physical, _config, verbosityLevel>0);
    }

    // placing virtual detector at exit of DS neutron shielding

    vdId = VirtualDetectorId::STM_UpStrHole;

    VolumeInfo const & parent = _helper->locateVolInfo("HallAir");
    CLHEP::Hep3Vector const& parentInMu2e = parent.centerInMu2e();

    if( vdg->exist(vdId) ) {
      if ( verbosityLevel > 0) {
        cout << __func__ << " constructing " << VirtualDetector::volumeName(vdId)  << endl;
      }


      // The following code is adding the inner part of vd 86 to account for
      // overlapping of old VD boundaries with the STM_CRVShieldPipe

      VolumeInfo vdInfo;

      double shieldHoleRadius = _config.getDouble("stm.shield.rIn")*CLHEP::mm;

      const TubsParams vdParams(0, shieldHoleRadius, vdg->getHalfLength());

      VolumeInfo vdHoleInfo = nestTubs(VirtualDetector::volumeName(vdId),
                                       vdParams,
                                       downstreamVacuumMaterial,
                                       0,
                                       vdg->getGlobal(vdId) - parentInMu2e, //position w.r.t. parent
                                       parent,
                                       vdId,
                                       vdIsVisible,
                                       G4Color::Red(),
                                       vdIsSolid,
                                       forceAuxEdgeVisible,
                                       placePV,
                                       false
                                       );

      doSurfaceCheck && checkForOverlaps(vdHoleInfo.physical, _config, verbosityLevel>0);

      if ( verbosityLevel > 0) {
        cout << __func__ << " constructing " << VirtualDetector::volumeName(vdId) << endl
             << " at " << vdg->getGlobal(vdId) << endl
             << " at " << vdg->getGlobal(vdId) - parentInMu2e <<" or local "<< vdg->getLocal(vdId)<< " w.r.t. parent (HallAir) " << endl;
        cout << __func__ << "    VD parameters: " << vdParams << endl;
        cout << __func__ << "    VD rel. posit: " << vdg->getLocal(vdId) << endl;
      }
    }


  } // constructVirtualDetectors()
}
