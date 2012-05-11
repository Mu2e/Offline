//
// Free function to create the virtual detectors
//
// $Id: constructVirtualDetectors.cc,v 1.35 2012/05/11 23:28:27 youzy Exp $
// $Author: youzy $
// $Date: 2012/05/11 23:28:27 $
//
// Original author KLG based on Mu2eWorld constructVirtualDetectors
//
// C++ includes

#include <iostream>
#include <string>

// Mu2e includes.

#include "Mu2eG4/inc/constructVirtualDetectors.hh"

#include "BeamlineGeom/inc/Beamline.hh"
#include "G4Helper/inc/G4Helper.hh"
#include "G4Helper/inc/VolumeInfo.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "ProductionSolenoidGeom/inc/PSEnclosure.hh"
#include "GeomPrimitives/inc/Tube.hh"
#include "Mu2eBuildingGeom/inc/Mu2eBuilding.hh"
#include "ProtonBeamDumpGeom/inc/ProtonBeamDump.hh"
#include "ProductionTargetGeom/inc/ProductionTarget.hh"
#include "ExtinctionMonitorUCIGeom/inc/ExtMonUCI.hh"
#include "Mu2eG4/inc/MaterialFinder.hh"
#include "Mu2eG4/inc/SensitiveDetectorName.hh"
#include "Mu2eG4/inc/VirtualDetectorSD.hh"
#include "Mu2eG4/inc/nestTubs.hh"
#include "Mu2eG4/inc/nestBox.hh"
#include "Mu2eG4/inc/finishNesting.hh"
#include "GeometryService/inc/VirtualDetector.hh"
#include "TTrackerGeom/inc/TTracker.hh"
#include "MCDataProducts/inc/VirtualDetectorId.hh"
#include "G4Helper/inc/AntiLeakRegistry.hh"

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

  void constructVirtualDetectors( SimpleConfig const * const _config ){

    // Place virtual detectors

    int static const verbosityLevel = _config->getInt("vd.verbosityLevel",0);

    bool vdIsVisible         = _config->getBool("vd.visible",true);
    bool vdIsSolid           = _config->getBool("vd.solid",true);
    bool forceAuxEdgeVisible = _config->getBool("g4.forceAuxEdgeVisible",false);
    bool doSurfaceCheck      = _config->getBool("g4.doSurfaceCheck",false);
    bool const placePV       = true;

    int const nSurfaceCheckPoints = 100000; // for a more thorrow check due to the vd shape

    AntiLeakRegistry& reg = art::ServiceHandle<G4Helper>()->antiLeakRegistry();
    GeomHandle<VirtualDetector> vdg;
    if( vdg->nDet()<=0 ) return;

    GeomHandle<Mu2eBuilding> building;

    GeomHandle<Beamline> beamg;
    double rVac         = CLHEP::mm * beamg->getTS().innerRadius();

    double vdHalfLength = CLHEP::mm * vdg->getHalfLength();

    MaterialFinder materialFinder(*_config);
    G4Material* vacuumMaterial     = materialFinder.get("toyDS.insideMaterialName");

    TubsParams vdParams(0,rVac,vdHalfLength);

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
        VolumeInfo const & parent = _helper->locateVolInfo("ToyTS1Vacuum");
        if ( verbosityLevel > 0) {
          cout << __func__ << " constructing " << VirtualDetector::volumeName(vdId)  << endl;
        }
        VolumeInfo vd = nestTubs( VirtualDetector::volumeName(vdId),
                                  vdParams, vacuumMaterial, 0,
                                  vdg->getLocal(vdId),
                                  parent,
                                  vdId, vdIsVisible, G4Color::Red(), vdIsSolid,
                                  forceAuxEdgeVisible,
                                  placePV,
                                  false);
        // vd are very thin, a more thorough check is needed
        doSurfaceCheck && vd.physical->CheckOverlaps(nSurfaceCheckPoints,0.0,true);
        vd.logical->SetSensitiveDetector(vdSD);
      }

    // Virtual Detectors Coll31_In, Coll31_Out, Coll32_In, Coll32_Out are placed inside TS3

    for( int vdId=VirtualDetectorId::Coll31_In;
         vdId<=VirtualDetectorId::Coll32_Out;
         ++vdId) if( vdg->exist(vdId) ) {
        VolumeInfo const & parent = _helper->locateVolInfo("ToyTS3Vacuum");
        if ( verbosityLevel > 0) {
          cout << __func__ << " constructing " << VirtualDetector::volumeName(vdId)  << endl;
        }
        VolumeInfo vd = nestTubs( VirtualDetector::volumeName(vdId),
                                  vdParams, vacuumMaterial, 0,
                                  vdg->getLocal(vdId),
                                  parent,
                                  vdId, vdIsVisible, G4Color::Red(), vdIsSolid,
                                  forceAuxEdgeVisible,
                                  placePV,
                                  false);
        // vd are very thin, a more thorough check is needed
        doSurfaceCheck && vd.physical->CheckOverlaps(nSurfaceCheckPoints,0.0,true);
        vd.logical->SetSensitiveDetector(vdSD);
      }

    // Virtual Detectors Coll5_In, Coll5_Out are placed inside TS5

    for( int vdId=VirtualDetectorId::Coll5_In; 
         vdId<=VirtualDetectorId::Coll5_Out; 
         ++vdId) if( vdg->exist(vdId) ) {
        VolumeInfo const & parent = _helper->locateVolInfo("ToyTS5Vacuum");
        if ( verbosityLevel > 0) {
          cout << __func__ << " constructing " << VirtualDetector::volumeName(vdId)  << endl;
        }
        VolumeInfo vd = nestTubs( VirtualDetector::volumeName(vdId),
                                  vdParams, vacuumMaterial, 0,
                                  vdg->getLocal(vdId),
                                  parent,
                                  vdId, vdIsVisible, G4Color::Red(), vdIsSolid,
                                  forceAuxEdgeVisible,
                                  placePV,
                                  false);
        // vd are very thin, a more thorough check is needed
        doSurfaceCheck && vd.physical->CheckOverlaps(nSurfaceCheckPoints,0.0,true);
        vd.logical->SetSensitiveDetector(vdSD);
      }

    // Virtual Detectors ST_In, ST_Out are placed inside DS2, just before and after stopping target

    // If there is no neutron absorber, virtual detectors 9 and 10 extend to
    // inner wall of DS2 minus 5 mm. If neutron absorber is defined, these
    // detectors extend to neutron absorber minus 5 mm.
    if ( !_config->getBool("isDumbbell",false) ){
    double Ravr = _config->getDouble("toyDS.rIn");
    double deltaR = 0;
    double Z0 = 0;
    double deltaZ = 1.0;

    if ( _config->getBool("hasNeutronAbsorber",false) ) {
      double NAIInnerRadius0     = _config->getDouble("neutronabsorber.internalInnerRadius0");
      double NAIInnerRadius1     = _config->getDouble("neutronabsorber.internalInnerRadius1");
      Ravr   = (NAIInnerRadius0+NAIInnerRadius1)/2;
      deltaR = (NAIInnerRadius1-NAIInnerRadius0);
      Z0     = _config->getDouble("neutronabsorber.internalZ01");
      deltaZ = 2.0 *_config->getDouble("neutronabsorber.internalHalfLengthZ01");
    }

    for( int vdId=VirtualDetectorId::ST_In; 
         vdId<=VirtualDetectorId::ST_Out; 
         ++vdId) if( vdg->exist(vdId) ) {
      
        if ( verbosityLevel > 0) {
          cout << __func__ << " constructing " << VirtualDetector::volumeName(vdId)  << endl;
        }

        double zvd = vdg->getGlobal(vdId).z();
        double rvd = Ravr + deltaR/deltaZ*(zvd-Z0) - 5.0;

        if ( verbosityLevel > 0) {
          cout << __func__ << " " << VirtualDetector::volumeName(vdId) <<
            " z, r : " << zvd << ", " << rvd << endl;
        }

        TubsParams vdParamsTarget(0.,rvd,vdHalfLength);

        VolumeInfo const & parent = ( _config->getBool("isDumbbell",false) ) ? _helper->locateVolInfo("ToyDS3Vacuum") : _helper->locateVolInfo("ToyDS2Vacuum"); //ToyDS3Vacuum to move the targets

        if (verbosityLevel >0) {
          cout << __func__ << " " << VirtualDetector::volumeName(vdId) << " Z offset in Mu2e    : " <<
            zvd << endl;      
          cout << __func__ << " " << VirtualDetector::volumeName(vdId) << " Z extent in Mu2e    : " <<
            zvd - vdHalfLength << ", " << zvd + vdHalfLength << endl;
        }

        VolumeInfo vd = nestTubs( VirtualDetector::volumeName(vdId),
                                  vdParamsTarget, vacuumMaterial, 0,
                                  vdg->getLocal(vdId),
                                  parent,
                                  vdId,
                                  vdIsVisible, 
                                  G4Color::Red(), vdIsSolid,
                                  forceAuxEdgeVisible,
                                  placePV,
                                  false);
        // vd are very thin, a more thorough check is needed
        doSurfaceCheck && vd.physical->CheckOverlaps(nSurfaceCheckPoints,0.0,true);

        vd.logical->SetSensitiveDetector(vdSD);
      }
    }

    // placing virtual detectors in the middle of the ttracker

    // check if ttracker exists and if the number of devices
    // ttracker.numDevices is even is done in VirtualDetectorMaker

    int vdId = VirtualDetectorId::TT_Mid;
    if( vdg->exist(vdId) ) {
      
      if ( verbosityLevel > 0) {
        cout << __func__ << " constructing " << VirtualDetector::volumeName(vdId)  << endl;
      }

      // the radius of tracker mother
      TTracker const & ttracker = *(GeomHandle<TTracker>());
      double orvd = ttracker.getTrackerEnvelopeParams().outerRadius();
      double irvd = ttracker.getTrackerEnvelopeParams().innerRadius();

      if ( verbosityLevel > 0) {
        double zvd = vdg->getGlobal(vdId).z();
        cout << __func__  << " " << VirtualDetector::volumeName(vdId) <<
          " z, r : " << zvd << ", " << irvd << " " << orvd << endl;
      }

      TubsParams vdParamsTTracker(irvd,orvd,vdHalfLength);

      VolumeInfo const & parent = _helper->locateVolInfo("TrackerMother");

      VolumeInfo vd = nestTubs( VirtualDetector::volumeName(vdId),
                                vdParamsTTracker, vacuumMaterial, 0,
                                vdg->getLocal(vdId),
                                parent,
                                vdId, vdIsVisible, G4Color::Red(), vdIsSolid,
                                forceAuxEdgeVisible,
                                placePV,
                                false);
      // vd are very thin, a more thorough check is needed
      doSurfaceCheck && vd.physical->CheckOverlaps(nSurfaceCheckPoints,0.0,true);
      vd.logical->SetSensitiveDetector(vdSD);

      vdId = VirtualDetectorId::TT_MidInner;
      if( vdg->exist(vdId) ) {

        if ( verbosityLevel > 0) {
          cout << __func__ << " constructing " << VirtualDetector::volumeName(vdId)  << endl;
        }

        // VD TT_MidInner is placed inside the ttracker at the same z position as
        // VD TT_Mid but from radius 0 to the inner radius of the ttracker
        // mother volume. However, its mother volume is ToyDS3Vacuum
        // which has a different offset. We will use the global offset
        // here (!) as DS is not in the geometry service yet

        // we need to take into account the "overlap" with the TT_InSurf

        TubsParams vdParamsTTrackerInner(0.,irvd-2.*vdHalfLength,vdHalfLength);

        VolumeInfo const & parent = _helper->locateVolInfo("ToyDS3Vacuum");

        G4ThreeVector vdLocalOffset = vdg->getGlobal(vdId) - parent.centerInMu2e();

        if ( verbosityLevel > 0) {
          double zvd = vdg->getGlobal(vdId).z();
          cout << __func__ << " " << VirtualDetector::volumeName(vdId) <<
            " z, r : " << zvd  << ", " << irvd << endl;
        }

        VolumeInfo vd = nestTubs( VirtualDetector::volumeName(vdId),
                                  vdParamsTTrackerInner, vacuumMaterial, 0,
                                  vdLocalOffset,
                                  parent,
                                  vdId, vdIsVisible, G4Color::Red(), vdIsSolid,
                                  forceAuxEdgeVisible,
                                  placePV,
                                  false);
        // vd are very thin, a more thorough check is needed
        doSurfaceCheck && vd.physical->CheckOverlaps(nSurfaceCheckPoints,0.0,true);
        vd.logical->SetSensitiveDetector(vdSD);
      }

    }

    // placing virtual detectors TT_FrontHollow, TT_FrontPA in front
    // of the ttracker (in the proton absorber region); check if
    // ttracker exist is done in VirtualDetectorMaker

    if ( _config->getBool("hasProtonAbsorber",false) ) {

      vdId = VirtualDetectorId::TT_FrontHollow;
      if( vdg->exist(vdId) ) {

        if ( verbosityLevel > 0) {
          cout << __func__ << " constructing " << VirtualDetector::volumeName(vdId)  << endl;
        }
        if ( !_config->getBool("hasProtonAbsorber",false) ) {
          throw cet::exception("GEOM")
            << "This virtual detector " << VirtualDetectorId(vdId).name()
            << " can only be placed if proton absorber is present\n";
        }

        // the radius of tracker mother
        TTracker const & ttracker = *(GeomHandle<TTracker>());
        double orvd = ttracker.getTrackerEnvelopeParams().outerRadius();
        double vdZ  = vdg->getGlobal(vdId).z();

        if ( verbosityLevel > 0) {
          cout << __func__ << " " << VirtualDetector::volumeName(vdId) <<
            " z, r : " << vdZ << ", " << orvd << endl;
        }

        // we will create an subtraction solid 
        // (we will "subtract" protonAbsorber) 
        // and place it (the subtraction solid) in ToyDS3Vacuum

        VolumeInfo const & parent = _helper->locateVolInfo("ToyDS3Vacuum");
      
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
                      vacuumMaterial,
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
        // vd are very thin, a more thorough check is needed
        doSurfaceCheck && vdHollowInfo.physical->
          CheckOverlaps(nSurfaceCheckPoints,0.0,true);

        if ( verbosityLevel > 0) {

          // both protonabs2 & vd are placed in ToyDS3Vacuum, do they have proper local offsets?

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
                        vacuumMaterial,
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
          // vd are very thin, a more thorough check is needed
          doSurfaceCheck && vdIntersectionInfo.physical->
            CheckOverlaps(nSurfaceCheckPoints,0.0,true);

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

      // if there is no proton absorber only one simple vd is placed
      // the hollow part of the name may not be quite right here...

      vdId = VirtualDetectorId::TT_FrontHollow;
      if( vdg->exist(vdId) ) {

        if ( verbosityLevel > 0) {
          cout << __func__ << " constructing " << VirtualDetector::volumeName(vdId)  << endl;
        }

        // the radius of tracker mother
        TTracker const & ttracker = *(GeomHandle<TTracker>());
        double orvd = ttracker.getTrackerEnvelopeParams().outerRadius();
        double vdZ  = vdg->getGlobal(vdId).z();

        if ( verbosityLevel > 0) {
          cout << __func__ << " " << VirtualDetector::volumeName(vdId) <<
            " z, r : " << vdZ << ", " << orvd << endl;
        }

        VolumeInfo const & parent = _helper->locateVolInfo("ToyDS3Vacuum");
      
        G4ThreeVector vdLocalOffset = vdg->getGlobal(vdId) - parent.centerInMu2e();

        TubsParams  vdParamsTTrackerFrontFull(0.,orvd,vdHalfLength);

        VolumeInfo vdInfo = nestTubs(VirtualDetector::volumeName(vdId), 
                                     vdParamsTTrackerFrontFull, 
                                     vacuumMaterial, 
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
        // vd are very thin, a more thorough check is needed
        doSurfaceCheck && vdInfo.physical->CheckOverlaps(nSurfaceCheckPoints,0.0,true);

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
      double orvd = ttracker.getTrackerEnvelopeParams().outerRadius();
      double vdZ  = vdg->getGlobal(vdId).z();

      if ( verbosityLevel > 0) {
        cout << __func__ << " " << VirtualDetector::volumeName(vdId) <<
          " z, r : " << vdZ << ", " << orvd << endl;
      }

      VolumeInfo const & parent = _helper->locateVolInfo("ToyDS3Vacuum");
      
      G4ThreeVector vdLocalOffset = vdg->getGlobal(vdId) - parent.centerInMu2e();

      TubsParams  vdParamsTTrackerBackFull(0.,orvd,vdHalfLength);

      VolumeInfo vdInfo = nestTubs(VirtualDetector::volumeName(vdId), 
                                   vdParamsTTrackerBackFull, 
                                   vacuumMaterial, 
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
      // vd are very thin, a more thorough check is needed
      doSurfaceCheck && vdInfo.physical->CheckOverlaps(nSurfaceCheckPoints,0.0,true);

      vdInfo.logical->SetSensitiveDetector(vdSD);

    }

    vdId = VirtualDetectorId::TT_OutSurf;
    if( vdg->exist(vdId) ) {

      if ( verbosityLevel > 0) {
        cout << __func__ << " constructing " << VirtualDetector::volumeName(vdId)  << endl;
      }

      // the radius of tracker mother
      TTracker const & ttracker = *(GeomHandle<TTracker>());
      TubsParams const envelopeParams = ttracker.getTrackerEnvelopeParams();
      double orvd = envelopeParams.outerRadius();
      double vdZ  = vdg->getGlobal(vdId).z();

      if ( verbosityLevel > 0) {
        cout << __func__ << " " << VirtualDetector::volumeName(vdId) <<
          " z, r : " << vdZ << ", " << orvd << endl;
      }

      VolumeInfo const & parent = _helper->locateVolInfo("ToyDS3Vacuum");
      
      G4ThreeVector vdLocalOffset = vdg->getGlobal(vdId) - parent.centerInMu2e();

      // the detector is on the outer surface of the ttracker envelope
      // it is thin cylinder, NOT a thin disk
      TubsParams  vdParamsTTrackerOutSurf(orvd,orvd+2.*vdHalfLength,envelopeParams.zHalfLength());

      VolumeInfo vdInfo = nestTubs(VirtualDetector::volumeName(vdId), 
                                   vdParamsTTrackerOutSurf, 
                                   vacuumMaterial, 
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
      // vd are very thin, a more thorough check is needed
      doSurfaceCheck && vdInfo.physical->CheckOverlaps(nSurfaceCheckPoints,0.0,true);

      vdInfo.logical->SetSensitiveDetector(vdSD);

    }

    vdId = VirtualDetectorId::TT_InSurf;
    if( vdg->exist(vdId) ) {

      if ( verbosityLevel > 0) {
        cout << __func__ << " constructing " << VirtualDetector::volumeName(vdId)  << endl;
      }

      // the radius of tracker mother
      TTracker const & ttracker = *(GeomHandle<TTracker>());
      TubsParams const envelopeParams = ttracker.getTrackerEnvelopeParams();
      double irvd = envelopeParams.innerRadius();
      double vdZ  = vdg->getGlobal(vdId).z();

      if ( verbosityLevel > 0) {
        cout << __func__ << " " << VirtualDetector::volumeName(vdId) <<
          " z, r : " << vdZ << ", " << irvd << endl;
      }

      VolumeInfo const & parent = _helper->locateVolInfo("ToyDS3Vacuum");
      
      G4ThreeVector vdLocalOffset = vdg->getGlobal(vdId) - parent.centerInMu2e();

      // the detector is on the inner surface of the ttracker envelope
      // it is thin cylinder, NOT a thin disk
      TubsParams  vdParamsTTrackerInSurf(irvd-2.*vdHalfLength,irvd,envelopeParams.zHalfLength());

      VolumeInfo vdInfo = nestTubs(VirtualDetector::volumeName(vdId), 
                                   vdParamsTTrackerInSurf, 
                                   vacuumMaterial, 
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
      // vd are very thin, a more thorough check is needed
      doSurfaceCheck && vdInfo.physical->CheckOverlaps(nSurfaceCheckPoints,0.0,true);

      vdInfo.logical->SetSensitiveDetector(vdSD);

    }

    vdId = VirtualDetectorId::EMFC1Entrance;
    if( vdg->exist(vdId) ) {

      if ( verbosityLevel > 0) {
        cout << __func__ << " constructing " << VirtualDetector::volumeName(vdId)  << endl;
      }

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
                                  vacuumMaterial, 
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

      // vd are very thin, a more thorough check is needed
      doSurfaceCheck && vdInfo.physical->CheckOverlaps(nSurfaceCheckPoints,0.0,true);

      vdInfo.logical->SetSensitiveDetector(vdSD);
    }


    // placing virtual detector on the exit (beam dump direction) and inside
    // of PS vacuum,  right before the PS enclosure end plate.
    vdId = VirtualDetectorId::PS_FrontExit;
    if ( vdg->exist(vdId) )
    {
      const VolumeInfo& parent = _helper->locateVolInfo("PSEnclosureVacuum");

      GeomHandle<PSEnclosure> pse;

      const Tube& psevac = GeomHandle<PSEnclosure>()->vacuum();

      TubsParams vdParams(0., psevac.outerRadius(), vdg->getHalfLength());

      G4ThreeVector vdCenterInParent(0., 0., -psevac.halfLength() + vdg->getHalfLength());

      VolumeInfo vdInfo = nestTubs(VirtualDetector::volumeName(vdId),
                                   vdParams,
                                   vacuumMaterial,
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
      // vd are very thin, a more thorough check is needed
      doSurfaceCheck && vdInfo.physical->CheckOverlaps(nSurfaceCheckPoints,0.0,true);

      vdInfo.logical->SetSensitiveDetector(vdSD);
    }

    // placing virtual detector infront of ExtMonUCI removable shielding
    vdId = VirtualDetectorId::EMIEntrance1;
    if( vdg->exist(vdId) ) {
      VolumeInfo const & parent = _helper->locateVolInfo("HallAir");
      GeomHandle<ExtMonUCI::ExtMon> extmon;

      std::vector<double> hlen(3);
      hlen[0] = _config->getDouble("vd.ExtMonUCIEntrance1.HalfLengths", 100.0);
      hlen[1] = _config->getDouble("vd.ExtMonUCIEntrance1.HalfLengths", 100.0);
      hlen[2] = vdg->getHalfLength();

      std::vector<double> pos;
      _config->getVectorDouble("vd.ExtMonUCIEntrance1.Position", pos, 3);
      G4ThreeVector vdLocalOffset = G4ThreeVector( pos[0], pos[1], pos[2] ) - parent.centerInMu2e();

      VolumeInfo vdInfo = nestBox(VirtualDetector::volumeName(vdId),
                                  hlen,
                                  vacuumMaterial,
                                  0,
                                  vdLocalOffset,
                                  parent,
                                  vdId,
                                  vdIsVisible,
                                  G4Color::Red(),
                                  vdIsSolid,
                                  forceAuxEdgeVisible,
                                  placePV,
                                  false
                                 );

      // vd are very thin, a more thorough check is needed
      doSurfaceCheck && vdInfo.physical->CheckOverlaps(nSurfaceCheckPoints,0.0,true);

      vdInfo.logical->SetSensitiveDetector(vdSD);
    }

    // placing virtual detector between ExtMonUCI removable shielding and front shielding
    vdId = VirtualDetectorId::EMIEntrance2;
    if( vdg->exist(vdId) ) {
      VolumeInfo const & parent = _helper->locateVolInfo("ExtMonUCI");
      GeomHandle<ExtMonUCI::ExtMon> extmon;

      std::vector<double> hlen(3);
      hlen[0] = extmon->envelopeParams()[0];
      hlen[1] = extmon->envelopeParams()[1];
      hlen[2] = vdg->getHalfLength();

      VolumeInfo vdInfo = nestBox(VirtualDetector::volumeName(vdId),
                                  hlen,
                                  vacuumMaterial,
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

      // vd are very thin, a more thorough check is needed
      doSurfaceCheck && vdInfo.physical->CheckOverlaps(nSurfaceCheckPoints,0.0,true);

      vdInfo.logical->SetSensitiveDetector(vdSD);
    }

    // placing virtual detector at the entrance and exit collimator0-3
    vdId = VirtualDetectorId::EMIC0Entrance;
    for (int iCol = 0; iCol < 4; iCol++) {
      if( vdg->exist(vdId) ) {
        VolumeInfo const & parent = _helper->locateVolInfo("ExtMonUCI");
        GeomHandle<ExtMonUCI::ExtMon> extmon;

        std::vector<double> hlen(3);
        hlen[0] = extmon->col(iCol)->paramsOuter()[0];
        hlen[1] = extmon->col(iCol)->paramsOuter()[1];
        hlen[2] = vdg->getHalfLength();

        for (int iFrontBack = 0; iFrontBack < 2; iFrontBack++) {
          VolumeInfo vdInfo = nestBox(VirtualDetector::volumeName(vdId),
                                      hlen,
                                      vacuumMaterial,
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

          // vd are very thin, a more thorough check is needed
          doSurfaceCheck && vdInfo.physical->CheckOverlaps(nSurfaceCheckPoints,0.0,true);

          vdInfo.logical->SetSensitiveDetector(vdSD);

          vdId++;
        }
      }
    }

    // An XY plane between the PS and anything ExtMon
    vdId = VirtualDetectorId::ExtMonCommonPlane;
    if( vdg->exist(vdId) ) {
      if ( verbosityLevel > 0) {
        cout << __func__ << " constructing " << VirtualDetector::volumeName(vdId)  << endl;
      }

      VolumeInfo const & parent = _helper->locateVolInfo("HallAir");
      GeomHandle<ProtonBeamDump> dump;

      const double requested_z = _config->getDouble("vd.ExtMonCommonPlane.z");

      const double xmin = std::min(building->hallInsideXPSCorner(),
                                   dump->frontShieldingCenterInMu2e()[0]
                                   + (requested_z - dump->frontShieldingCenterInMu2e()[2])*tan(dump->coreRotY())
                                   - dump->frontShieldingHalfSize()[0]/cos(dump->coreRotY())
                                   );

      CLHEP::Hep3Vector centerInMu2e(
                                     (building->hallInsideXmax() + xmin)/2,
                                     (building->hallInsideYmax() + building->hallInsideYmin())/2,
                                     requested_z
                                     );

      std::vector<double> hlen(3);
      hlen[0] = (building->hallInsideXmax() - xmin)/2;
      hlen[1] = (building->hallInsideYmax() - building->hallInsideYmin())/2;
      hlen[2] = vdg->getHalfLength();

      VolumeInfo vdInfo = nestBox(VirtualDetector::volumeName(vdId),
                                  hlen,
                                  vacuumMaterial,
                                  0,
                                  centerInMu2e - parent.centerInMu2e(),
                                  parent,
                                  vdId,
                                  vdIsVisible,
                                  G4Color::Red(),
                                  vdIsSolid,
                                  forceAuxEdgeVisible,
                                  placePV,
                                  false
                                 );

      // vd are very thin, a more thorough check is needed
      doSurfaceCheck && vdInfo.physical->CheckOverlaps(nSurfaceCheckPoints,0.0,true);

      vdInfo.logical->SetSensitiveDetector(vdSD);
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
                                  vacuumMaterial,
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

      // vd are very thin, a more thorough check is needed
      doSurfaceCheck && vdInfo.physical->CheckOverlaps(nSurfaceCheckPoints,0.0,true);

      vdInfo.logical->SetSensitiveDetector(vdSD);
    }

    // ExtMonFNAL detector VDs - created in constructExtMonFNAL()

  } // constructVirtualDetectors()
}
