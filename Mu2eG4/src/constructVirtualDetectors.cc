//
// Free function to create the virtual detectors
//
// $Id: constructVirtualDetectors.cc,v 1.8 2011/06/10 17:48:48 genser Exp $
// $Author: genser $
// $Date: 2011/06/10 17:48:48 $
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
#include "Mu2eG4/inc/MaterialFinder.hh"
#include "Mu2eG4/inc/SensitiveDetectorName.hh"
#include "Mu2eG4/inc/VirtualDetectorSD.hh"
#include "Mu2eG4/inc/nestTubs.hh"
#include "Mu2eG4/inc/finishNesting.hh"
#include "VirtualDetectorGeom/inc/VirtualDetector.hh"
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

  void constructVirtualDetectors( SimpleConfig const * const _config ){

    // Place virtual detectors

    int static const verbosityLevel = _config->getInt("vd.verbosityLevel",0);

    bool vdIsVisible         = _config->getBool("vd.visible",true);
    bool vdIsSolid           = _config->getBool("vd.solid",true);
    bool forceAuxEdgeVisible = _config->getBool("g4.forceAuxEdgeVisible",false);
    bool doSurfaceCheck      = _config->getBool("g4.doSurfaceCheck",false);
    bool const placePV       = true;

    int const nSurfaceCheckPoints = 100000;

    GeomHandle<VirtualDetector> vdg;
    if( vdg->nDet()<=0 ) return;

    GeomHandle<Beamline> beamg;
    double rVac         = CLHEP::mm * beamg->getTS().innerRadius();

    double vdHalfLength = CLHEP::mm * vdg->getHalfLength();

    MaterialFinder materialFinder(*_config);
    G4Material* vacuumMaterial     = materialFinder.get("toyDS.insideMaterialName");

    TubsParams vdParams(0,rVac,vdHalfLength);

    // Virtual Detectors 1 and 2 are placed inside TS1

    G4VSensitiveDetector* vdSD = G4SDManager::GetSDMpointer()->
      FindSensitiveDetector(SensitiveDetectorName::VirtualDetector());

    G4Helper* _helper = &(*(art::ServiceHandle<G4Helper>()));

    string static const vdBaseName = "VirtualDetector";
    ostringstream name;

    for( int vdId=1; vdId<=2; ++vdId) if( vdg->exist(vdId) ) {
      VolumeInfo const & parent = _helper->locateVolInfo("ToyTS1Vacuum");
      name.str("");
      name << vdBaseName << vdId;
      VolumeInfo vd = nestTubs( name.str(), vdParams, vacuumMaterial, 0,
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

    // Virtual Detectors 3-6 are placed inside TS3

    for( int vdId=3; vdId<=6; ++vdId) if( vdg->exist(vdId) ) {
      VolumeInfo const & parent = _helper->locateVolInfo("ToyTS3Vacuum");
      name.str("");
      name << vdBaseName << vdId;
      VolumeInfo vd = nestTubs( name.str(), vdParams, vacuumMaterial, 0,
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

    // Virtual Detectors 7-8 are placed inside TS5

    for( int vdId=7; vdId<=8; ++vdId) if( vdg->exist(vdId) ) {
      VolumeInfo const & parent = _helper->locateVolInfo("ToyTS5Vacuum");
      name.str("");
      name << vdBaseName << vdId;
      VolumeInfo vd = nestTubs( name.str(), vdParams, vacuumMaterial, 0,
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

    // Virtual Detectors 9-10 are placed inside DS2, just before and after stopping target

    // If there is no neutron absorber, virtual detectors 9 and 10 extend to
    // inner wall of DS2 minus 5 mm. If neutron absorber is defined, these
    // detectors extend to neutron absorber minus 5 mm.

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

    for( int vdId=9; vdId<=10; ++vdId) if( vdg->exist(vdId) ) {
      
      double zvd = vdg->getGlobal(vdId).z();
      double rvd = Ravr + deltaR/deltaZ*(zvd-Z0) - 5.0;

      if ( verbosityLevel > 0) {
        cout << __func__ << " VD " << vdId << 
	  " z, r : " << zvd << ", " << rvd << endl;
      }

      TubsParams vdParamsTarget(0.,rvd,vdHalfLength);

      VolumeInfo const & parent = _helper->locateVolInfo("ToyDS2Vacuum");
      name.str("");
      name << vdBaseName << vdId;
      VolumeInfo vd = nestTubs( name.str(), vdParamsTarget, vacuumMaterial, 0,
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

      if (verbosityLevel >0) {
        cout << __func__ << " " << name.str() << " Z offset in Mu2e    : " <<
          zvd << endl;      
        cout << __func__ << " " << name.str() << " Z extent in Mu2e    : " <<
          zvd - vdHalfLength << ", " << zvd + vdHalfLength << endl;
      }

      vd.logical->SetSensitiveDetector(vdSD);
    }

    // placing virtual detectors in the middle of the ttracker

    // check if ttracker exists and if the number of devices
    // ttracker.numDevices is even is done in VirtualDetectorMaker

    int vdId = 11;
    if( vdg->exist(vdId) ) {
      
      // the radius of tracker mother
      TTracker const & ttracker = *(GeomHandle<TTracker>());
      double orvd = ttracker.getTrackerEnvelopeParams().outerRadius();
      double irvd = ttracker.getTrackerEnvelopeParams().innerRadius();

      if ( verbosityLevel > 0) {
        double zvd = vdg->getGlobal(vdId).z();
        cout << __func__  << " VD " << vdId << 
	  " z, r : " << zvd << ", " << irvd << " " << orvd << endl;
      }

      TubsParams vdParamsTTracker(irvd,orvd,vdHalfLength);

      VolumeInfo const & parent = _helper->locateVolInfo("TrackerMother");
      name.str("");
      name << vdBaseName << vdId;
      VolumeInfo vd = nestTubs( name.str(), vdParamsTTracker, vacuumMaterial, 0,
                                vdg->getLocal(vdId),
                                parent,
                                vdId, vdIsVisible, G4Color::Red(), vdIsSolid,
                                forceAuxEdgeVisible,
                                placePV,
                                false);
      // vd are very thin, a more thorough check is needed
      doSurfaceCheck && vd.physical->CheckOverlaps(nSurfaceCheckPoints,0.0,true);
      vd.logical->SetSensitiveDetector(vdSD);

      vdId = 12;
      if( vdg->exist(vdId) ) {

        TubsParams vdParamsTTrackerInner(0.,irvd,vdHalfLength);

        // VD 12 is placed inside the ttracker at the same z position as
        // VD 11 but from radius 0 to the inner radius of the ttracker
        // mother volume. However, its mother volume is ToyDS3Vacuum
        // which has a different offset. We will use the global offset
        // here (!) as DS is not in the geometry service yet

        VolumeInfo const & parent = _helper->locateVolInfo("ToyDS3Vacuum");

        G4ThreeVector vdLocalOffset = vdg->getGlobal(vdId) - parent.centerInMu2e();

        if ( verbosityLevel > 0) {
          double zvd = vdg->getGlobal(vdId).z();
          cout << __func__ << " VD " << vdId << 
	    " z, r : " << zvd  << ", " << irvd << endl;
        }

        name.str("");
        name << vdBaseName << vdId;
        VolumeInfo vd = nestTubs( name.str(), vdParamsTTrackerInner, vacuumMaterial, 0,
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

    // placing virtual detectors in front of the ttracker (in the proton absorber region)
    // check if ttracker exist is done in VirtualDetectorMaker

    if ( _config->getBool("hasProtonAbsorber",false) ) {

      vdId = 13;
      if( vdg->exist(vdId) ) {

        if ( !_config->getBool("hasProtonAbsorber",false) ) {
          throw cet::exception("GEOM")
            << "This virtual detector " << vdId
            << " can only be placed if proton absorber is present\n";
        }

        // the radius of tracker mother
        TTracker const & ttracker = *(GeomHandle<TTracker>());
        double orvd = ttracker.getTrackerEnvelopeParams().outerRadius();
        double vdZ  = vdg->getGlobal(vdId).z();

        if ( verbosityLevel > 0) {
          cout << __func__ << " VD " << vdId << 
            " z, r : " << vdZ << ", " << orvd << endl;
        }

        // we will create an subtraction solid 
        // (we will "subtract" protonAbsorber) 
        // and place it in ToyDS3Vacuum

        VolumeInfo const & parent = _helper->locateVolInfo("ToyDS3Vacuum");
      
        G4ThreeVector vdLocalOffset = vdg->getGlobal(vdId) - parent.centerInMu2e();

        VolumeInfo vdFullInfo;
        name.str("");
        name << vdBaseName << vdId;
        vdFullInfo.name = name.str();

        TubsParams  vdParamsTTrackerFrontFull(0.,orvd,vdHalfLength);

        G4Tubs* vdFullSolid = new G4Tubs(vdFullInfo.name,
                                         vdParamsTTrackerFrontFull.innerRadius(),
                                         vdParamsTTrackerFrontFull.outerRadius(),
                                         vdParamsTTrackerFrontFull.zHalfLength(),
                                         vdParamsTTrackerFrontFull.phi0(),
                                         vdParamsTTrackerFrontFull.phiMax());

        VolumeInfo const & protonabs2Info = _helper->locateVolInfo("protonabs2");

        VolumeInfo vdHollowInfo;
        vdHollowInfo.name = name.str();

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

        vdHollowInfo.solid = new G4SubtractionSolid(vdHollowInfo.name,
                                                    vdFullSolid, 
                                                    protonabs2Info.solid,
                                                    0,
                                                    protonabs2Info.centerInMu2e()-vdg->getGlobal(vdId));

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
          double theHL = vdFullSolid->GetZHalfLength();
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
            vdHollowInfo.physical->GetTranslation() - protonabs2Info.physical->GetTranslation() << endl;
          cout << __func__ << " protonabs2Info.centerInMu2e() - vdg->getGlobal(vdId) offset           : " <<
            protonabs2Info.centerInMu2e()-vdg->getGlobal(vdId) << endl;
        }

        vdHollowInfo.logical->SetSensitiveDetector(vdSD);

        //  now the complementary solid, it has to be placed in protonabs2

        vdId = 14;
        if (vdg->exist(vdId)) {
          name.str("");
          name << vdBaseName << vdId;
          VolumeInfo vdIntersectionInfo;
          vdIntersectionInfo.name = name.str();

          vdIntersectionInfo.solid = new G4IntersectionSolid(vdIntersectionInfo.name,
                                                             vdFullSolid, 
                                                             protonabs2Info.solid,
                                                             0,
                                                             protonabs2Info.centerInMu2e()-vdg->getGlobal(vdId));

          VolumeInfo const & parent = _helper->locateVolInfo("protonabs2");
          vdLocalOffset = vdg->getGlobal(vdId)-protonabs2Info.centerInMu2e();

          vdIntersectionInfo.centerInParent = vdLocalOffset;
          vdIntersectionInfo.centerInWorld  = vdIntersectionInfo.centerInParent + parent.centerInWorld;
 
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
            double theHL = vdFullSolid->GetZHalfLength();
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
              vdIntersectionInfo.physical->GetTranslation() - protonabs2Info.physical->GetTranslation() << endl;
            cout << __func__ << " protonabs2Info.centerInMu2e() - vdg->getGlobal(vdId) offset           : " <<
              protonabs2Info.centerInMu2e()-vdg->getGlobal(vdId) << endl;
          }


          vdIntersectionInfo.logical->SetSensitiveDetector(vdSD);

        }
      }
    }
    else {

      // if there is no proton absorber only one simple vd is placed

      vdId = 13;
      if( vdg->exist(vdId) ) {

        // the radius of tracker mother
        TTracker const & ttracker = *(GeomHandle<TTracker>());
        double orvd = ttracker.getTrackerEnvelopeParams().outerRadius();
        double vdZ  = vdg->getGlobal(vdId).z();

        if ( verbosityLevel > 0) {
          cout << __func__ << " VD " << vdId << 
            " z, r : " << vdZ << ", " << orvd << endl;
        }

        VolumeInfo const & parent = _helper->locateVolInfo("ToyDS3Vacuum");
      
        G4ThreeVector vdLocalOffset = vdg->getGlobal(vdId) - parent.centerInMu2e();

        VolumeInfo vdFullInfo;
        name.str("");
        name << vdBaseName << vdId;
        vdFullInfo.name = name.str();

        TubsParams  vdParamsTTrackerFrontFull(0.,orvd,vdHalfLength);

        VolumeInfo vdInfo = nestTubs(name.str(), 
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

  }
}
