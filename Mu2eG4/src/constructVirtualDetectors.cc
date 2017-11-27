//
// Free function to create the virtual detectors
//
// $Id: constructVirtualDetectors.cc,v 1.70 2014/07/29 16:24:44 genser Exp $
// $Author: genser $
// $Date: 2014/07/29 16:24:44 $
//
// Original author KLG based on Mu2eWorld constructVirtualDetectors

// C++ includes
#include <iostream>
#include <string>

// Mu2e includes.
#include "Mu2eG4/inc/constructVirtualDetectors.hh"

#include "BeamlineGeom/inc/Beamline.hh"
#include "CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "DetectorSolenoidGeom/inc/DetectorSolenoid.hh"
#include "G4Helper/inc/G4Helper.hh"
#include "G4Helper/inc/VolumeInfo.hh"
#include "GeomPrimitives/inc/Tube.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/VirtualDetector.hh"
#include "DataProducts/inc/VirtualDetectorId.hh"
#include "MECOStyleProtonAbsorberGeom/inc/MECOStyleProtonAbsorber.hh"
#include "Mu2eG4/inc/SensitiveDetectorName.hh"
#include "Mu2eG4/inc/checkForOverlaps.hh"
#include "Mu2eG4/inc/findMaterialOrThrow.hh"
#include "Mu2eG4/inc/finishNesting.hh"
#include "Mu2eG4/inc/nestTubs.hh"
#include "Mu2eG4/inc/nestBox.hh"
#include "ProductionSolenoidGeom/inc/PSVacuum.hh"
#include "ProtonBeamDumpGeom/inc/ProtonBeamDump.hh"
#include "ProductionTargetGeom/inc/ProductionTarget.hh"
#include "TTrackerGeom/inc/TTracker.hh"

// G4 includes
#include "G4Material.hh"
#include "G4SDManager.hh"
#include "G4VSensitiveDetector.hh"
#include "G4Color.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4SubtractionSolid.hh"
#include "G4IntersectionSolid.hh"

using namespace std;

namespace mu2e {

  // Construct the virtual detectors

  void constructVirtualDetectors( SimpleConfig const & _config ){

    // Place virtual detectors

    int static const verbosityLevel = _config.getInt("vd.verbosityLevel",0);

    bool const vdIsVisible         = _config.getBool("vd.visible",true);
    bool const vdIsSolid           = _config.getBool("vd.solid",true);
    bool const forceAuxEdgeVisible = _config.getBool("g4.forceAuxEdgeVisible",false);
    bool const doSurfaceCheck      = _config.getBool("g4.doSurfaceCheck",false) || _config.getBool("vd.doSurfaceCheck",false);
    bool const placePV             = true;

    GeomHandle<VirtualDetector> vdg;
    if( vdg->nDet()<=0 ) return;

    GeomHandle<Beamline> beamg;

    GeomHandle<DetectorSolenoid> ds;
    TransportSolenoid const&  ts = beamg->getTS();
    G4Material* downstreamVacuumMaterial = findMaterialOrThrow( ds->vacuumMaterial() );
    G4Material* upstreamVacuumMaterial   = findMaterialOrThrow(  ts.upstreamVacuumMaterial() );

    double rCol = ts.getColl51().rOut();
    double rCin = ts.getColl1().rIn1();
    double vdHalfLength = CLHEP::mm * vdg->getHalfLength();
    
    TubsParams vdParams(0,rCol,vdHalfLength);
    TubsParams vdParamsIn(0,rCin,vdHalfLength);

    // Virtual Detectors Coll1_In, COll1_Out are placed inside TS1

    G4VSensitiveDetector* vdSD = G4SDManager::GetSDMpointer()->
      FindSensitiveDetector(SensitiveDetectorName::VirtualDetector());

    G4Helper* _helper = &(*(art::ServiceHandle<G4Helper>()));

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

        vd.logical->SetSensitiveDetector(vdSD);
      }

    // Virtual Detectors Coll1_pBarCollar_In, COll1_pBarCollar_Out are 
    // placed inside Coll1, which is inside TS1

    // Just copy what is done above, with minor edits.  
    // FIXME: one should factorize some the code below; the main
    // things which change: parent and offset
    for( int vdId=VirtualDetectorId::Coll1_pBarCollar_In;
         vdId<=VirtualDetectorId::Coll1_pBarCollar_Out;
         ++vdId) if( vdg->exist(vdId) ) {
        VolumeInfo const & parent = _helper->locateVolInfo("TS1Vacuum");
        if ( verbosityLevel > 0) {
          cout << __func__ << " constructing " << VirtualDetector::volumeName(vdId)
               << " at " << vdg->getGlobal(vdId) << endl;
          cout << __func__ << "    VD parameters: " << vdParamsIn << endl;
          cout << __func__ << "    VD rel. posit: " << vdg->getLocal(vdId) << endl;
        }

        VolumeInfo vd = nestTubs( VirtualDetector::volumeName(vdId),
                                  vdParamsIn, upstreamVacuumMaterial, 0,
                                  vdg->getLocal(vdId),
                                  parent,
                                  vdId, vdIsVisible, G4Color::Red(), vdIsSolid,
                                  forceAuxEdgeVisible,
                                  placePV,
                                  false);

        doSurfaceCheck && checkForOverlaps(vd.physical, _config, verbosityLevel>0);

        vd.logical->SetSensitiveDetector(vdSD);
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

        myvd.logical->SetSensitiveDetector(vdSD);
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

        myvd.logical->SetSensitiveDetector(vdSD);
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

        vd.logical->SetSensitiveDetector(vdSD);
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

        vd.logical->SetSensitiveDetector(vdSD);
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

      vd.logical->SetSensitiveDetector(vdSD);

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

          vd.logical->SetSensitiveDetector(vdSD);
          
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

          vd.logical->SetSensitiveDetector(vdSD);
        }
    }

    if ( _config.getBool("hasTTracker",false)  ) {


      // placing virtual detectors in the middle of the ttracker

      // check if ttracker exists and if the number of planes
      // ttracker.numPlanes is even is done in VirtualDetectorMaker

      vdId = VirtualDetectorId::TT_Mid;
      if( vdg->exist(vdId) ) {

        if ( verbosityLevel > 0) {
          cout << __func__ << " constructing " << VirtualDetector::volumeName(vdId)  << endl;
        }

        // the radius of tracker mother
        TTracker const & ttracker = *(GeomHandle<TTracker>());
        double orvd = ttracker.mother().tubsParams().outerRadius();
        double irvd = ttracker.mother().tubsParams().innerRadius();

        if ( ttracker.getSupportModel() == SupportModel::detailedv0 ) {
          auto const& beams =  ttracker.getSupportStructure().beamBody();
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

        TubsParams vdParamsTTracker(irvd,orvd,vdHalfLength);

        VolumeInfo const & parent = _helper->locateVolInfo("TrackerMother");

        CLHEP::Hep3Vector vdPos = vdg->getGlobal(vdId)-parent.centerInMu2e();
        //        cout << "foo: TT_Mid: " << vdPos << " " << vdg->getLocal(vdId) << endl;
        //        cout << "foo: TT_Mid: " << vdParamsTTracker << endl;

        VolumeInfo vd = nestTubs( VirtualDetector::volumeName(vdId),
                                  vdParamsTTracker, downstreamVacuumMaterial, 0,
                                  vdPos,
                                  parent,
                                  vdId, vdIsVisible, G4Color::Red(), vdIsSolid,
                                  forceAuxEdgeVisible,
                                  placePV,
                                  false);

        doSurfaceCheck && checkForOverlaps(vd.physical, _config, verbosityLevel>0);

        vd.logical->SetSensitiveDetector(vdSD);

        vdId = VirtualDetectorId::TT_MidInner;
        if( vdg->exist(vdId) ) {

          if ( verbosityLevel > 0) {
            cout << __func__ << " constructing " << VirtualDetector::volumeName(vdId)  << endl;
          }

          // VD TT_MidInner is placed inside the ttracker at the same z position as
          // VD TT_Mid but from radius 0 to the inner radius of the ttracker
          // mother volume. However, its mother volume is DS3Vacuum
          // which has a different offset. We will use the global offset
          // here (!) as DS is not in the geometry service yet

          // we need to take into account the "overlap" with the TT_InSurf

          TubsParams vdParamsTTrackerInner(0.,irvd-2.*vdHalfLength,vdHalfLength);
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
          //        cout << "foo: TT_MidInner: " << vdParamsTTrackerInner << endl;

          VolumeInfo vd = nestTubs( VirtualDetector::volumeName(vdId),
                                    vdParamsTTrackerInner, downstreamVacuumMaterial, 0,
                                    vdLocalOffset,
                                    parent,
                                    vdId, vdIsVisible, G4Color::Red(), vdIsSolid,
                                    forceAuxEdgeVisible,
                                    placePV,
                                    false);

          doSurfaceCheck && checkForOverlaps(vd.physical, _config, verbosityLevel>0);

          vd.logical->SetSensitiveDetector(vdSD);
        }

      }

      // placing virtual detectors TT_FrontHollow, TT_FrontPA in front
      // of the ttracker (in the proton absorber region); check if
      // ttracker exist is done in VirtualDetectorMaker

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
          TTracker const & ttracker = *(GeomHandle<TTracker>());
          double orvd = ttracker.mother().tubsParams().outerRadius();
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

          TubsParams  vdParamsTTrackerFrontFull(0.,orvd,vdHalfLength);

          vdFullInfo.solid = new G4Tubs(vdFullInfo.name,
                                        vdParamsTTrackerFrontFull.innerRadius(),
                                        vdParamsTTrackerFrontFull.outerRadius(),
                                        vdParamsTTrackerFrontFull.zHalfLength(),
                                        vdParamsTTrackerFrontFull.phi0(),
                                        vdParamsTTrackerFrontFull.phiMax());

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

          vdHollowInfo.logical->SetSensitiveDetector(vdSD);

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


            vdIntersectionInfo.logical->SetSensitiveDetector(vdSD);

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
          TTracker const & ttracker = *(GeomHandle<TTracker>());
          double orvd = ttracker.mother().tubsParams().outerRadius();
          double vdZ  = vdg->getGlobal(vdId).z();

          if ( verbosityLevel > 0) {
            cout << __func__ << " " << VirtualDetector::volumeName(vdId) <<
              " z, r : " << vdZ << ", " << orvd << endl;
          }

	  std::string theDS3("DS3Vacuum");
	  if ( _config.getBool("inGaragePosition",false) ) theDS3 = "garageFakeDS3Vacuum";
          VolumeInfo const & parent = _helper->locateVolInfo(theDS3);

          G4ThreeVector vdLocalOffset = vdg->getGlobal(vdId) - parent.centerInMu2e();

          TubsParams  vdParamsTTrackerFrontFull(0.,orvd,vdHalfLength);

          //          cout << "foo: TT_Front: " << vdLocalOffset             << endl;
          //          cout << "foo: TT_Front: " << vdParamsTTrackerFrontFull << endl;
          VolumeInfo vdInfo = nestTubs(VirtualDetector::volumeName(vdId),
                                       vdParamsTTrackerFrontFull,
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

          vdInfo.logical->SetSensitiveDetector(vdSD);

        }

      }

      vdId = VirtualDetectorId::TT_Back;
      if( vdg->exist(vdId) ) {

        if ( verbosityLevel > 0) {
          cout << __func__ << " constructing " << VirtualDetector::volumeName(vdId)  << endl;
        }
        // the radius of tracker mother
        TTracker const & ttracker = *(GeomHandle<TTracker>());
        double orvd = ttracker.mother().tubsParams().outerRadius();
        double vdZ  = vdg->getGlobal(vdId).z();

        if ( verbosityLevel > 0) {
          cout << __func__ << " " << VirtualDetector::volumeName(vdId) <<
            " z, r : " << vdZ << ", " << orvd << endl;
        }

	  std::string theDS3("DS3Vacuum");
	  if ( _config.getBool("inGaragePosition",false) ) theDS3 = "garageFakeDS3Vacuum";
        VolumeInfo const & parent = _helper->locateVolInfo(theDS3);

        G4ThreeVector vdLocalOffset = vdg->getGlobal(vdId) - parent.centerInMu2e();

        TubsParams  vdParamsTTrackerBackFull(0.,orvd,vdHalfLength);

        //        cout << "foo: TT_Back: " << vdLocalOffset    << endl;
        //        cout << "foo: TT_Back: " << vdParamsTTrackerBackFull << endl;

        VolumeInfo vdInfo = nestTubs(VirtualDetector::volumeName(vdId),
                                     vdParamsTTrackerBackFull,
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

        vdInfo.logical->SetSensitiveDetector(vdSD);

      }

      vdId = VirtualDetectorId::TT_OutSurf;
      if( vdg->exist(vdId) ) {

        if ( verbosityLevel > 0) {
          cout << __func__ << " constructing " << VirtualDetector::volumeName(vdId)  << endl;
        }

        // the radius of tracker mother
        TTracker const & ttracker = *(GeomHandle<TTracker>());
        TubsParams const& motherParams = ttracker.mother().tubsParams();
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

        // the detector is on the outer surface of the ttracker envelope
        // it is thin cylinder, NOT a thin disk
        TubsParams  vdParamsTTrackerOutSurf(orvd,orvd+2.*vdHalfLength,motherParams.zHalfLength());

        //        cout << "foo: TT_OutSurf: " << vdLocalOffset         << endl;
        //        cout << "foo: TT_OutSurf: " << vdParamsTTrackerOutSurf << endl;

        VolumeInfo vdInfo = nestTubs(VirtualDetector::volumeName(vdId),
                                     vdParamsTTrackerOutSurf,
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

        vdInfo.logical->SetSensitiveDetector(vdSD);

      }

      vdId = VirtualDetectorId::TT_InSurf;
      if( vdg->exist(vdId) ) {

        if ( verbosityLevel > 0) {
          cout << __func__ << " constructing " << VirtualDetector::volumeName(vdId)  << endl;
        }

        // the radius of tracker mother
        TTracker const & ttracker = *(GeomHandle<TTracker>());
        TubsParams const& motherParams = ttracker.mother().tubsParams();
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

        // the detector is on the inner surface of the ttracker envelope
        // it is thin cylinder, NOT a thin disk
        TubsParams  vdParamsTTrackerInSurf(irvd-2.*vdHalfLength,irvd,motherParams.zHalfLength());

        //        cout << "foo: TT_InSurf: " << vdLocalOffset         << endl;
        //        cout << "foo: TT_InSurf: " << vdParamsTTrackerInSurf << endl;

        VolumeInfo vdInfo = nestTubs(VirtualDetector::volumeName(vdId),
                                     vdParamsTTrackerInSurf,
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

        vdInfo.logical->SetSensitiveDetector(vdSD);

      }

    } // end hasTTracker
    else if ( _config.getBool("hasITracker",false) && _config.getBool("itracker.VirtualDetect",false) ) {

      VolumeInfo const & trckrParent = _helper->locateVolInfo("TrackerMother");

      double tModInRd = ((G4Tubs *)trckrParent.solid)->GetInnerRadius();
      double tModOtRd = ((G4Tubs *)trckrParent.solid)->GetOuterRadius();
      double tModDz   = ((G4Tubs *)trckrParent.solid)->GetDz();


      vdId = VirtualDetectorId::IT_VD_InSurf;
      if( vdg->exist(vdId) ) {

        if ( verbosityLevel > 0) {
          cout << __func__ << " constructing " << VirtualDetector::volumeName(vdId)  << endl;
        }

        // the radius of tracker mother
        double irvd = tModInRd; //envelopeParams.innerRadius();
        double vdZ  = vdg->getGlobal(vdId).z();

        if ( verbosityLevel > 0) {
          cout << __func__ << " " << VirtualDetector::volumeName(vdId) <<
            " z, r : " << vdZ << ", " << irvd << endl;
        }

        VolumeInfo const & parent = _helper->locateVolInfo("DS3Vacuum");

        G4ThreeVector vdLocalOffset = vdg->getGlobal(vdId) - parent.centerInMu2e();

        // the detector is on the inner surface of the tracker envelope
        // it is thin cylinder, NOT a thin disk
        TubsParams  vdParamsITrackerInSurf(irvd-2.*vdHalfLength,irvd,tModDz);

        VolumeInfo vdInfo = nestTubs(VirtualDetector::volumeName(vdId),
                                     vdParamsITrackerInSurf,
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

        vdInfo.logical->SetSensitiveDetector(vdSD);
      }

      vdId = VirtualDetectorId::IT_VD_EndCap_Front;
      if( vdg->exist(vdId) ) {

        if ( verbosityLevel > 0) {
          cout << __func__ << " constructing " << VirtualDetector::volumeName(vdId)  << endl;
        }

        // the radius of tracker mother
        double irvd = tModInRd; //envelopeParams.innerRadius();
        double vdZ  = vdg->getGlobal(vdId).z();

        if ( verbosityLevel > 0) {
          cout << __func__ << " " << VirtualDetector::volumeName(vdId) <<
            " z, r : " << vdZ << ", " << irvd << endl;
        }

        VolumeInfo const & parent = _helper->locateVolInfo("DS3Vacuum");

        G4ThreeVector vdLocalOffset = vdg->getGlobal(vdId) - parent.centerInMu2e();// -G4ThreeVector(0.0,0.0,-(tModDz+vdHalfLength));

        // the detector is on the inner surface of the tracker envelope
        // it is thin cylinder, NOT a thin disk
        TubsParams  vdParamsITrackerInSurf(tModInRd,tModOtRd,vdHalfLength);

        VolumeInfo vdInfo = nestTubs(VirtualDetector::volumeName(vdId),
                                     vdParamsITrackerInSurf,
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

        vdInfo.logical->SetSensitiveDetector(vdSD);
      }

      vdId = VirtualDetectorId::IT_VD_EndCap_Back;
      if( vdg->exist(vdId) ) {

        if ( verbosityLevel > 0) {
          cout << __func__ << " constructing " << VirtualDetector::volumeName(vdId)  << endl;
        }

        // the radius of tracker mother
        double irvd = tModInRd; //envelopeParams.innerRadius();
        double vdZ  = vdg->getGlobal(vdId).z();

        if ( verbosityLevel > 0) {
          cout << __func__ << " " << VirtualDetector::volumeName(vdId) <<
            " z, r : " << vdZ << ", " << irvd << endl;
        }

        VolumeInfo const & parent = _helper->locateVolInfo("DS3Vacuum");

        G4ThreeVector vdLocalOffset = vdg->getGlobal(vdId) - parent.centerInMu2e();// -G4ThreeVector(0.0,0.0,(tModDz+vdHalfLength));

        // the detector is on the inner surface of the tracker envelope
        // it is thin cylinder, NOT a thin disk
        TubsParams  vdParamsITrackerInSurf(tModInRd,tModOtRd,vdHalfLength);

        VolumeInfo vdInfo = nestTubs(VirtualDetector::volumeName(vdId),
                                     vdParamsITrackerInSurf,
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

        vdInfo.logical->SetSensitiveDetector(vdSD);
      }

    } // end hasITracker

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

      vdInfo.logical->SetSensitiveDetector(vdSD);
      */
    }


    // placing virtual detector on the exit (beam dump direction) and inside
    // of PS vacuum,  right before the PS enclosure end plate.
    vdId = VirtualDetectorId::PS_FrontExit;
    if ( vdg->exist(vdId) )
      {
        const VolumeInfo& parent = _helper->locateVolInfo("PSVacuum");

        const Tube& psevac = GeomHandle<PSVacuum>()->vacuum();

        TubsParams vdParams(0., psevac.outerRadius(), vdg->getHalfLength());

        G4ThreeVector vdCenterInParent(0., 0., -psevac.halfLength() + vdg->getHalfLength());

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

        vdInfo.logical->SetSensitiveDetector(vdSD);
      }

    if (_config.getBool("targetPS.hasVD.backward", false)) {
            vdId = VirtualDetectorId::PT_Back;
            if ( vdg->exist(vdId) )
            {
                    const VolumeInfo& parent = _helper->locateVolInfo("ProductionTargetMother");
                    G4Tubs *PTMoth = (G4Tubs *) parent.solid;

                    TubsParams vdParams(0., PTMoth->GetOuterRadius(), vdg->getHalfLength());

                    G4ThreeVector vdCenterInParent(0., 0., -PTMoth->GetDz() + vdg->getHalfLength());

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

                    vdInfo.logical->SetSensitiveDetector(vdSD);
            }
    }

    if (_config.getBool("targetPS.hasVD.forward", false)) {
            vdId = VirtualDetectorId::PT_Front;
            if ( vdg->exist(vdId) )
            {
                    const VolumeInfo& parent = _helper->locateVolInfo("ProductionTargetMother");
                    G4Tubs *PTMoth = (G4Tubs *) parent.solid;

                    TubsParams vdParams(0., PTMoth->GetOuterRadius(), vdg->getHalfLength());

                    G4ThreeVector vdCenterInParent(0., 0., PTMoth->GetDz() - vdg->getHalfLength());

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

                    vdInfo.logical->SetSensitiveDetector(vdSD);
            }
    }

    // An XY plane between the PS and anything ExtMon
    vdId = VirtualDetectorId::ExtMonCommonPlane;
    if( vdg->exist(vdId) ) {
      // Not currently supported
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

      vdInfo.logical->SetSensitiveDetector(vdSD);
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

      vdInfo.logical->SetSensitiveDetector(vdSD);
      
      if ( verbosityLevel > 0) {
          cout << __func__ << " constructing " << VirtualDetector::volumeName(vdId) << endl
               << " at " << vdg->getGlobal(vdId) << endl
               << " at " << vdg->getGlobal(vdId) - parentInMu2e <<" or local "<< vdg->getLocal(vdId)<< " w.r.t. parent (HallAir) " << endl;
          cout << __func__ << "    VD parameters: " << vdParams << endl;
          cout << __func__ << "    VD rel. posit: " << vdg->getLocal(vdId) << endl;
      }
    }




    if ( _config.getBool("hasDiskCalorimeter",true) ) {

      int vdIdDiskEdge = VirtualDetectorId::EMC_Disk_0_EdgeIn;
      int vdIdDiskSurf = VirtualDetectorId::EMC_Disk_0_SurfIn;
      int vdIdFEBEdge  = VirtualDetectorId::EMC_FEB_0_EdgeIn;
      int vdIdFEBSurf  = VirtualDetectorId::EMC_FEB_0_SurfIn;

      double delta     = 2*vdg->getHalfLength()+0.02;


      DiskCalorimeter const& cal = *(GeomHandle<DiskCalorimeter>());
      VolumeInfo const& caloParent       = _helper->locateVolInfo("CalorimeterMother");    
      CLHEP::Hep3Vector const& caloParentInMu2e     = caloParent.centerInMu2e();

      for (size_t id = 0; id < cal.nDisk(); id++)
      {
	   std::ostringstream cratename;      cratename<<"CalorimeterFEB_" <<id;

	   VolumeInfo const& caloFEBParent              = _helper->locateVolInfo(cratename.str());
	   CLHEP::Hep3Vector const& caloFEBParentInMu2e = caloFEBParent.centerInMu2e();

           const CLHEP::Hep3Vector & sizeDisk = cal.disk(id).geomInfo().size();
           G4ThreeVector posDisk              = cal.geomInfo().origin() + cal.disk(id).geomInfo().originLocal();
           G4double crateHalfLength           = cal.caloInfo().crateHalfLength();
	   G4double crystalDepth              = 2.0*cal.caloInfo().crystalHalfLength();
	   
	   G4double wrapThickness             = cal.caloInfo().wrapperThickness();
	   G4double wrapHalfDepth             = (crystalDepth + wrapThickness)/2.; 

     	   G4ThreeVector posCrate             = cal.disk(id).geomInfo().origin() + CLHEP::Hep3Vector(0.0,0.0,crateHalfLength-wrapHalfDepth);
	   G4double diskRadIn                 = cal.caloInfo().stepsRadiusIn()  - cal.caloInfo().caseThicknessIn();
	   G4double diskRadOut                = cal.caloInfo().stepsRadiusOut() + cal.caloInfo().caseThicknessOut();
           G4double crateRadIn                = cal.caloInfo().crateRadiusIn();
           G4double crateRadOut               = cal.caloInfo().crateRadiusOut();
      	   G4double phi0Crate                 = (15./360.)*2*CLHEP::pi;                      

	   TubsParams  vdParamsFrontDisk(diskRadIn-delta,   diskRadOut+delta,   vdg->getHalfLength());
           TubsParams  vdParamsInnerDisk(diskRadIn-2*delta, diskRadIn-delta,   sizeDisk[2]/2.0);
           // TubsParams  vdParamsOuterDisk(diskRadOut+delta, diskRadOut+2*delta, sizeDisk[2]/2.0);
 
           TubsParams  vdParamsFrontFEB(crateRadIn + 2*delta, crateRadOut - 2*delta, vdg->getHalfLength()   , -phi0Crate,CLHEP::pi+2*phi0Crate);
           TubsParams  vdParamsInnerFEB(crateRadIn + delta  , crateRadIn +2*delta  , crateHalfLength-2*delta, -phi0Crate,CLHEP::pi+2*phi0Crate);
           TubsParams  vdParamsOuterFEB(crateRadOut - 2*delta , crateRadOut-delta  , crateHalfLength-2*delta, -phi0Crate,CLHEP::pi+2*phi0Crate);

           G4ThreeVector posFrontDisk = posDisk - caloParentInMu2e - G4ThreeVector(0,0,sizeDisk.z()/2.0+delta);
           G4ThreeVector posBackDisk  = posDisk - caloParentInMu2e + G4ThreeVector(0,0,sizeDisk.z()/2.0+delta);
           G4ThreeVector posInnerDisk = posDisk - caloParentInMu2e;

           G4ThreeVector posFrontFEB  = posCrate - caloFEBParentInMu2e - G4ThreeVector(0,0,crateHalfLength-delta);
           G4ThreeVector posBackFEB   = posCrate - caloFEBParentInMu2e + G4ThreeVector(0,0,crateHalfLength-delta);
           G4ThreeVector posInnerFEB  = posCrate - caloFEBParentInMu2e + G4ThreeVector(0,2*delta, 0);
           G4ThreeVector posOuterFEB  = posCrate - caloFEBParentInMu2e;


           if( vdg->exist(vdIdDiskSurf) )
           {
               VolumeInfo vdInfo = nestTubs(VirtualDetector::volumeName(vdIdDiskSurf),
                                            vdParamsFrontDisk,
                                            downstreamVacuumMaterial,
                                            0,
                                            posFrontDisk,
                                            caloParent,
                                            vdIdDiskSurf,
                                            vdIsVisible,
                                            G4Color::Red(),
                                            vdIsSolid,
                                            forceAuxEdgeVisible,
                                            placePV,
                                            false);
               ++vdIdDiskSurf;

               VolumeInfo vdInfo2 = nestTubs(VirtualDetector::volumeName(vdIdDiskSurf),
                                            vdParamsFrontDisk,
                                            downstreamVacuumMaterial,
                                            0,
                                            posBackDisk,
                                            caloParent,
                                            vdIdDiskSurf,
                                            vdIsVisible,
                                            G4Color::Red(),
                                            vdIsSolid,
                                            forceAuxEdgeVisible,
                                            placePV,
                                            false);
               ++vdIdDiskSurf;

               doSurfaceCheck && checkForOverlaps(vdInfo.physical, _config, verbosityLevel>0);
               doSurfaceCheck && checkForOverlaps(vdInfo2.physical, _config, verbosityLevel>0);
               vdInfo.logical->SetSensitiveDetector(vdSD);
               vdInfo2.logical->SetSensitiveDetector(vdSD);
           }

           
           if( vdg->exist(vdIdDiskEdge) )
           {
               VolumeInfo vdInfo = nestTubs(VirtualDetector::volumeName(vdIdDiskEdge),
                                            vdParamsInnerDisk,
                                            downstreamVacuumMaterial,
                                            0,
                                            posInnerDisk,
                                            caloParent,
                                            vdIdDiskEdge,
                                            vdIsVisible,
                                            G4Color::Red(),
                                            vdIsSolid,
                                            forceAuxEdgeVisible,
                                            placePV,
                                            false);
               ++vdIdDiskEdge;

               // VolumeInfo vdInfo2 = nestTubs(VirtualDetector::volumeName(vdIdDiskEdge),
               //                              vdParamsOuterDisk,
               //                              downstreamVacuumMaterial,
               //                              0,
               //                              posInnerDisk,
               //                              caloParent,
               //                              vdIdDiskEdge,
               //                              vdIsVisible,
               //                              G4Color::Red(),
               //                              vdIsSolid,
               //                              forceAuxEdgeVisible,
               //                              placePV,
               //                              false);
               // ++vdIdDiskEdge;

               doSurfaceCheck && checkForOverlaps(vdInfo.physical, _config, verbosityLevel>0);
               // doSurfaceCheck && checkForOverlaps(vdInfo2.physical, _config, verbosityLevel>0);
               vdInfo.logical->SetSensitiveDetector(vdSD);
               // vdInfo2.logical->SetSensitiveDetector(vdSD);
           }
           
           if( vdg->exist(vdIdFEBSurf) )
           {
               VolumeInfo vdInfo = nestTubs(VirtualDetector::volumeName(vdIdFEBSurf),
                                            vdParamsFrontFEB,
                                            downstreamVacuumMaterial,
                                            0,
                                            posFrontFEB,
                                            caloFEBParent,
                                            vdIdFEBSurf,
                                            vdIsVisible,
                                            G4Color::Red(),
                                            vdIsSolid,
                                            forceAuxEdgeVisible,
                                            placePV,
                                            false);
               ++vdIdFEBSurf;

               VolumeInfo vdInfo2 = nestTubs(VirtualDetector::volumeName(vdIdFEBSurf),
                                            vdParamsFrontFEB,
                                            downstreamVacuumMaterial,
                                            0,
                                            posBackFEB,
                                            caloFEBParent,
                                            vdIdFEBSurf,
                                            vdIsVisible,
                                            G4Color::Red(),
                                            vdIsSolid,
                                            forceAuxEdgeVisible,
                                            placePV,
                                            false);
               ++vdIdFEBSurf;

	       doSurfaceCheck && checkForOverlaps(vdInfo.physical, _config, verbosityLevel>0);
	       doSurfaceCheck && checkForOverlaps(vdInfo2.physical, _config, verbosityLevel>0);
	       vdInfo.logical->SetSensitiveDetector(vdSD);
	       vdInfo2.logical->SetSensitiveDetector(vdSD);
           }
           
           if( vdg->exist(vdIdFEBEdge) )
           {
               VolumeInfo vdInfo = nestTubs(VirtualDetector::volumeName(vdIdFEBEdge),
                                            vdParamsInnerFEB,
                                            downstreamVacuumMaterial,
                                            0,
                                            posInnerFEB,
                                            caloFEBParent,
                                            vdIdFEBEdge,
                                            vdIsVisible,
                                            G4Color::Red(),
                                            vdIsSolid,
                                            forceAuxEdgeVisible,
                                            placePV,
                                            false);
               ++vdIdFEBEdge;

               VolumeInfo vdInfo2 = nestTubs(VirtualDetector::volumeName(vdIdFEBEdge),
                                            vdParamsOuterFEB,
                                            downstreamVacuumMaterial,
                                            0,
                                            posOuterFEB,
                                            caloFEBParent,
                                            vdIdFEBEdge,
                                            vdIsVisible,
                                            G4Color::Red(),
                                            vdIsSolid,
                                            forceAuxEdgeVisible,
                                            placePV,
                                            false);
               ++vdIdFEBEdge;

               doSurfaceCheck && checkForOverlaps(vdInfo.physical, _config, verbosityLevel>0);
	       doSurfaceCheck && checkForOverlaps(vdInfo2.physical, _config, verbosityLevel>0);
               vdInfo.logical->SetSensitiveDetector(vdSD);
	       vdInfo2.logical->SetSensitiveDetector(vdSD);
           }

       }// end loop on disks
    }//hasDiskCalorimeter

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
      vdInfo.logical->SetSensitiveDetector(vdSD);
    }

//     vdId = VirtualDetectorId::STM_CRVShieldDnStr;
//     if ( vdg->exist(vdId) ) {
// 
//       const VolumeInfo& parent = _helper->locateVolInfo("HallAir");
//       GeomHandle<CosmicRayShield> CRS;
//       //const double y_crv_max       = CRS->getSectorPosition("D").y() + (CRS->getSectorHalfLengths("D"))[1];
//       const double yExtentLow      = std::abs(_config.getDouble("yOfFloorSurface.below.mu2eOrigin") );
//       const double x_vd_halflength = (CRS->getSectorHalfLengths("D"))[0];
//       //const double y_vd_halflength = (y_crv_max + yExtentLow)/2.0;
//       const double y_mother_halflength = yExtentLow;
//       const double dimVD[3] = { x_vd_halflength, y_mother_halflength, vdg->getHalfLength() };
// 
//       VolumeInfo vdInfo = nestBox(VirtualDetector::volumeName(vdId),
//                                   dimVD,
//                                   downstreamVacuumMaterial,
//                                   0,                              //rotation
//                                   vdg->getLocal(vdId),
//                                   parent,
//                                   vdId,
//                                   vdIsVisible,
//                                   G4Color::White(),
//                                   vdIsSolid,
//                                   forceAuxEdgeVisible,
//                                   placePV,
//                                   false
//                                   );
// 
//       if ( verbosityLevel > 0) {
//           cout << __func__ << " constructing " << VirtualDetector::volumeName(vdId) << endl
//                << " at " << vdg->getGlobal(vdId) << endl
//                << " at " << vdg->getLocal(vdId) << " w.r.t. parent (HallAir) " << endl;
//           cout << __func__ << "    VD parameters: " << vdParams << endl;
//           cout << __func__ << "    VD rel. posit: " << vdg->getLocal(vdId) << endl;
//       }
//       doSurfaceCheck && checkForOverlaps(vdInfo.physical, _config, verbosityLevel>0);
//       vdInfo.logical->SetSensitiveDetector(vdSD);
//     }    

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
      vdInfo.logical->SetSensitiveDetector(vdSD);
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
        yExtentLow          = _config.getDouble("stm.magnet.halfHeight");
        x_vd_halflength     = _config.getDouble("stm.magnet.halfWidth");
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
      vdInfo.logical->SetSensitiveDetector(vdSD);
    }

    vdId = VirtualDetectorId::STM_SpotSizeCollUpStr;
    if ( vdg->exist(vdId) ) {

      //const VolumeInfo& parent = _helper->locateVolInfo("MSTMMother");
      const VolumeInfo& parent = _helper->locateVolInfo("HallAir");
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
      vdInfo.logical->SetSensitiveDetector(vdSD);
    }    
    
    vdId = VirtualDetectorId::STM_CollDnStr;
    if ( vdg->exist(vdId) ) {

      //const VolumeInfo& parent = _helper->locateVolInfo("MSTMMother");
      const VolumeInfo& parent = _helper->locateVolInfo("HallAir");
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
      vdInfo.logical->SetSensitiveDetector(vdSD);
    }

    
    vdId = VirtualDetectorId::STM_Det1UpStr;
    if ( vdg->exist(vdId) ) {

      //const VolumeInfo& parent = _helper->locateVolInfo("MSTMMother");
      const VolumeInfo& parent = _helper->locateVolInfo("HallAir");
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
      vdInfo.logical->SetSensitiveDetector(vdSD);
    }

    
    vdId = VirtualDetectorId::STM_Det2UpStr;
    if ( vdg->exist(vdId) ) {

      //const VolumeInfo& parent = _helper->locateVolInfo("MSTMMother");
      const VolumeInfo& parent = _helper->locateVolInfo("HallAir");
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
      vdInfo.logical->SetSensitiveDetector(vdSD);
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
      }

      doSurfaceCheck && checkForOverlaps(vdInfo.physical, _config, verbosityLevel>0);

      vdInfo.logical->SetSensitiveDetector(vdSD);
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

      vdInfo.logical->SetSensitiveDetector(vdSD);
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

        vdInfo.logical->SetSensitiveDetector(vdSD);
      }
    }
   
  } // constructVirtualDetectors()
}
