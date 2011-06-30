//
// Define a sensitive detector for Straws.
// This version does not use G4HCofThisEvent etc...
// Framwork DataProducts are used instead
//
// $Id: StrawSD.cc,v 1.33 2011/06/30 04:55:13 kutschke Exp $
// $Author: kutschke $
// $Date: 2011/06/30 04:55:13 $
//
// Original author Rob Kutschke
//

#include <cstdio>

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/exception.h"

// Mu2e includes
#include "Mu2eG4/inc/StrawSD.hh"
#include "Mu2eG4/inc/EventNumberList.hh"
#include "Mu2eG4/inc/PhysicsProcessInfo.hh"
#include "LTrackerGeom/inc/LTracker.hh"
#include "TTrackerGeom/inc/TTracker.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "Mu2eUtilities/inc/TwoLinePCA.hh"
#include "Mu2eUtilities/inc/LinePointPCA.hh"
#include "Mu2eUtilities/inc/SimpleConfig.hh"

// G4 includes
#include "G4RunManager.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4ios.hh"

//
// Outstanding questions:
//
// 1) Why is diffAngle so big?
//

using namespace std;

namespace mu2e {

  StrawSD::StrawSD(G4String name, const SimpleConfig& config ):
    G4VSensitiveDetector(name),
    _collection(0),
    _processInfo(0),
    _nStrawsPerDevice(0),
    _nStrawsPerSector(0),
    _TrackerVersion(0),
    _debugList(0),
    _sizeLimit(config.getInt("g4.stepsSizeLimit",0)),
    _currentSize(0),
    _simID(0),
    _productGetter(0){

    // Get list of events for which to make debug printout.
    string key("g4.strawSDEventList");
    if ( config.hasName(key) ){
      vector<int> list;
      config.getVectorInt(key,list);
      _debugList.add(list);
    }

    art::ServiceHandle<GeometryService> geom;

    if ( !geom->hasElement<TTracker>() && !geom->hasElement<LTracker>() ) {
      throw cet::exception("GEOM")
        << "Expected one of L or T Trackers but found neither.\n";
    }

    if ( geom->hasElement<TTracker>() ) {

      GeomHandle<TTracker> ttracker;

      const Device& device = ttracker->getDevice(0);
      const Sector& sector = device.getSector(0);
      const Layer&  layer  = sector.getLayer(0);

      _nStrawsPerSector = sector.nLayers()  * layer.nStraws();
      _nStrawsPerDevice = device.nSectors() * _nStrawsPerSector;

      _TrackerVersion = config.getInt("TTrackerVersion",3);

      if ( _TrackerVersion != 3) {
        throw cet::exception("StrawSD")
          << "Expected TTrackerVersion of 3 but found " << _TrackerVersion <<endl;
        // esp take a look at the detectorOrigin calculation
      }

    }

    if ( geom->hasElement<LTracker>() ) {

      GeomHandle<LTracker> ltracker;

      _TrackerVersion = config.getInt("LTrackerVersion",3); // also see Mu2eWorld.cc

      if ( _TrackerVersion != 3) {
        throw cet::exception("StrawSD")
          << "Expected LTrackerVersion of 3 but found " << _TrackerVersion <<endl;
        // esp take a look at the detectorOrigin calculation
      }

    }

  }

  StrawSD::~StrawSD(){ }

  void StrawSD::Initialize(G4HCofThisEvent* HCE){
    _currentSize=0;
  }

  G4bool StrawSD::ProcessHits(G4Step* aStep,G4TouchableHistory*){

    _currentSize += 1;

    if( _sizeLimit>0 && _currentSize>_sizeLimit ) {
      if( (_currentSize - _sizeLimit)==1 ) {
        mf::LogWarning("G4") << "Maximum number of particles reached in StrawSD: "
                              << _currentSize << endl;
      }
      return false;
    }

    G4double edep = aStep->GetTotalEnergyDeposit();
    G4double step = aStep->GetStepLength();

    // I am not sure why we get these cases but we do.  Skip them.
    if ( edep == 0. && step == 0. ) return false;

    const G4TouchableHandle & touchableHandle = aStep->GetPreStepPoint()->GetTouchableHandle();

    // Origin of the LTracker.  Need to get this from G4.
    //    static G4ThreeVector detectorOrigin( -3904., -7350., 6200.);
    static G4ThreeVector detectorOrigin = GetTrackerOrigin(touchableHandle);

//     cout << "Debugging detectorOrigin   " << detectorOrigin  << endl;
//     cout << "Debugging det name         " << touchableHandle->GetVolume(2)->GetName() << endl;
//     cout << "Debugging det G4 origin    " << cdo  << endl;

    // Position at start of step point, in world system and in
    // a system in which the center of the tracking detector is the origin.
    G4ThreeVector prePosWorld   = aStep->GetPreStepPoint()->GetPosition();
    G4ThreeVector prePosTracker = prePosWorld - detectorOrigin;

    G4ThreeVector preMomWorld = aStep->GetPreStepPoint()->GetMomentum();

    G4Event const* event = G4RunManager::GetRunManager()->GetCurrentEvent();

//     G4int en = event->GetEventID();
//     G4int ti = aStep->GetTrack()->GetTrackID();
//     G4int cn = touchableHandle->GetCopyNumber();
//     G4int rn = touchableHandle->GetReplicaNumber();

//     G4TouchableHistory* theTouchable =
//       (G4TouchableHistory*)( aStep->GetPreStepPoint()->GetTouchable() );

//     cout << "Debugging history depth " <<
//       setw(4) << theTouchable->GetHistoryDepth() << endl;

//     cout << "Debugging replica 0 1 2 " <<
//       setw(4) << theTouchable->GetReplicaNumber(0) <<
//       setw(4) << theTouchable->GetReplicaNumber(1) <<
//       setw(4) << theTouchable->GetReplicaNumber(2) << endl;

//     cout << "Debugging replica 0 1 2 3 4" <<
//       setw(4) << touchableHandle->GetReplicaNumber(0) <<
//       setw(4) << touchableHandle->GetReplicaNumber(1) <<
//       setw(4) << touchableHandle->GetReplicaNumber(2) <<
//       setw(4) << touchableHandle->GetReplicaNumber(3) <<
//       setw(4) << touchableHandle->GetReplicaNumber(4) << endl;

//     cout << "Debugging PV Name Mother Name" <<
//       theTouchable->GetVolume(0)->GetName() << " " <<
//       theTouchable->GetVolume(1)->GetName() << " " <<
//       theTouchable->GetVolume(2)->GetName() << " " <<
//       theTouchable->GetVolume(3)->GetName() << " " <<
//       theTouchable->GetVolume(4)->GetName() << endl;

//     cout << "Debugging PV Name Mother Name" <<
//       touchableHandle->GetVolume(0)->GetName() << " " <<
//       touchableHandle->GetVolume(1)->GetName() << " " <<
//       touchableHandle->GetVolume(2)->GetName() << " " <<
//       touchableHandle->GetVolume(3)->GetName() << " " <<
//       touchableHandle->GetVolume(4)->GetName() << endl;

//     cout << "Debugging hit info event track copyn replican: " <<
//       setw(4) << en << " " <<
//       setw(4) << ti << " " <<
//       setw(4) << cn << " " <<
//       setw(4) << rn << endl;

    // getting the sector/device number

    G4int sdcn = 0;

    if ( _TrackerVersion == 3) {

//       cout << "Debugging _nStrawsPerDevice " << _nStrawsPerDevice << endl;
//       cout << "Debugging _nStrawsPerSector " << _nStrawsPerSector << endl;

      sdcn = touchableHandle->GetCopyNumber(1) +
        _nStrawsPerSector*(touchableHandle->GetReplicaNumber(2)) +
        _nStrawsPerDevice*(touchableHandle->GetReplicaNumber(3));

    } else {

      throw cet::exception("GEOM")
        << "Expected TrackerVersion of 3 but found " << _TrackerVersion << endl;

    }

    //    cout << "Debugging sdcn " << sdcn << endl;

    // We add the hit object to the framework strawHit collection created in produce

    G4String const& pname  = aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
    ProcessCode endCode(_processInfo->findAndCount(pname));

    _collection->push_back( StepPointMC(art::Ptr<SimParticle>( *_simID, aStep->GetTrack()->GetTrackID(), _productGetter ),
                                        sdcn,
                                        edep,
                                        aStep->GetNonIonizingEnergyDeposit(),
                                        aStep->GetPreStepPoint()->GetGlobalTime(),
                                        aStep->GetPreStepPoint()->GetProperTime(),
                                        prePosTracker,
                                        preMomWorld,
                                        step,
                                        endCode
                                        ));

    // Some debugging tests.
    if ( !_debugList.inList() ) return true;

    // Transformations between world to local coordinate systems.
    G4AffineTransform const& toLocal = touchableHandle->GetHistory()->GetTopTransform();
    G4AffineTransform        toWorld = toLocal.Inverse();

    G4ThreeVector postPosWorld = aStep->GetPostStepPoint()->GetPosition();
    G4ThreeVector postPosLocal = toLocal.TransformPoint(postPosWorld);

    G4ThreeVector prePosLocal  = toLocal.TransformPoint(prePosWorld);

    G4ThreeVector preMomLocal = toLocal.TransformAxis(preMomWorld);

    // Compute the directed chord in both local and world coordinates.
    G4ThreeVector deltaWorld = postPosWorld - prePosWorld;
    G4ThreeVector deltaLocal = postPosLocal - prePosLocal;

    // Angle between the directed chord and the momentum.
    G4ThreeVector dT( deltaWorld.x(), deltaWorld.y(), 0.);
    G4ThreeVector pT( preMomWorld.x(), preMomWorld.y(), 0. );
    //G4double angle = dT.angle(pT);

    G4ThreeVector dTLocal( deltaLocal.x(), deltaLocal.y(), 0.);
    G4ThreeVector pTLocal( preMomLocal.x(), preMomLocal.y(), 0. );
    //G4double angleLocal = dTLocal.angle(pTLocal);

    // This is took big. O(1.e-5 CLHEP::radians) or about 1% of the value. Why?
    //G4double diffAngle = angle-angleLocal;

    G4ThreeVector localOrigin(0.,0.,0.);
    G4ThreeVector worldOrigin = toWorld.TransformPoint(localOrigin) - detectorOrigin;

    G4ThreeVector localZUnit(0.,0.,1.);
    G4ThreeVector worldZUnit = toWorld.TransformAxis(localZUnit);

    // make sure it works with the constructTTrackerv3
    //    int copy = touchableHandle->GetCopyNumber();
    int copy = sdcn;
    int eventNo = event->GetEventID();


    /*
      // Works for both TTracker and LTracker.
    printf ( "Addhit: %4d %4d %6d %3d %3d | %10.2f %10.2f %10.2f | %10.2f %10.2f %10.2f | %10.7f %10.7f\n",
             eventNo,  _collection->size(), copy,
             aStep->IsFirstStepInVolume(), aStep->IsLastStepInVolume(),
             prePosTracker.x(), prePosTracker.y(), prePosTracker.z(),
             preMomWorld.x(),   preMomWorld.y(),   preMomWorld.z(),
             prePosLocal.perp(),  postPosLocal.perp()  );
    fflush(stdout);
    */

    // Reconstruction Geometry for the LTracker.
    // Need to make this work for the TTracker too.
    art::ServiceHandle<GeometryService> geom;
    if ( geom->hasElement<LTracker>() ) {

      GeomHandle<LTracker> ltracker;
      Straw const& straw = ltracker->getStraw( StrawIndex(copy) );
      G4ThreeVector mid  = straw.getMidPoint();
      G4ThreeVector w    = straw.getDirection();

      // Point of closest approach of track to wire.
      // Straight line approximation.
      TwoLinePCA pca( mid, w, prePosTracker, preMomWorld);

      // Point on the wire that is closest to the step point.
      LinePointPCA lppca( mid, w, prePosTracker);
      double ddd = lppca.unit().cosTheta(preMomWorld);

      double ttt = lppca.unit().cosTheta(w);

      printf ( "Addhit: %4d %4d %6d %3d %3d | %10.2f %10.2f %10.2f | %10.2f %10.2f %10.2f | %6.3f %10.7f | %10.7f %10.7f\n",
               eventNo,  int(_collection->size()), copy,
               aStep->IsFirstStepInVolume(), aStep->IsLastStepInVolume(),
               prePosTracker.x(), prePosTracker.y(), prePosTracker.z(),
               preMomWorld.x(),   preMomWorld.y(),   preMomWorld.z(), ddd, ttt,
               prePosLocal.perp(),  postPosLocal.perp()  );
      fflush(stdout);

      // This works.  uvw.perp() is always 2.500000000xxxx or 2.499999999xxxx
      // Note that uhat and vhat are not the same as the local (xhat, yhat)
      // from G4.  They differ by a rotation about the local zhat.
      G4ThreeVector z(0.,0.,1.);
      G4ThreeVector v = worldZUnit.cross(z).unit();
      G4ThreeVector u = v.cross(worldZUnit);
      G4ThreeVector detLocal( prePosTracker-mid);
      double u0 = detLocal.dot(u);
      double v0 = detLocal.dot(v);
      double w0 = detLocal.dot(w);
      G4ThreeVector uvw(u0, v0, w0);

    } else if ( geom->hasElement<TTracker>() ) {

      GeomHandle<TTracker> ttracker;
      Straw const& straw = ttracker->getStraw( StrawIndex(copy) );
      G4ThreeVector mid  = straw.getMidPoint();
      G4ThreeVector w    = straw.getDirection();

      // Point of closest approach of track to wire.
      // Straight line approximation.
      TwoLinePCA pca( mid, w, prePosTracker, preMomWorld);

      // Check that the radius of the reference point in the local
      // coordinates of the straw.  Should be 2.5 mm.
      double s = w.dot(prePosTracker-mid);
      CLHEP::Hep3Vector point = prePosTracker - (mid + s*w);

//      this works with transporOnly.py
//      StrawDetail const& strawDetail = straw.getDetail();
//       if (pca.dca()>strawDetail.outerRadius() || abs(point.mag()-strawDetail.outerRadius())>1.e-6 ) {

//         cerr << "*** Bad hit?: eid "
//              << event->GetEventID()    << ", cs "
//              << _collection->GetSize() << " tid "
//              << newHit->trackId()      << " vid "
//              << newHit->volumeId()     << " "
//              << straw.id()         << " | "
//              << pca.dca()          << " "
//              << prePosTracker      << " "
//              << preMomWorld        << " "
//              << point.mag()        << " "
//              << newHit->eDep()     << " "
//              << s                  << " | "
//              << mid
//              << endl;
//       }

    }

    // End of debug section.

    //    newHit->Print();
    //    newHit->Draw();

    return true;

  }

  void StrawSD::EndOfEvent(G4HCofThisEvent*){

    if( _sizeLimit>0 && _currentSize>=_sizeLimit ) {
      mf::LogWarning("G4") << "Total of " << _currentSize
                            << " straw hits were generated in the event."
                            << endl
                            << "Only " << _sizeLimit << " are saved in output collection."
                            << endl;
      cout << "Total of " << _currentSize
           << " straw hits were generated in the event."
           << endl
           << "Only " << _sizeLimit << " are saved in output collection."
           << endl;
    }

    if (verboseLevel>0) {
      G4int NbHits = _collection->size();
      G4cout << "\n-------->Hits Collection: in this event they are " << NbHits
             << " hits in the straw chambers: " << G4endl;
      for (G4int i=0;i<NbHits;i++) (*_collection)[i].print(G4cout);
    }
  }

  G4ThreeVector StrawSD::GetTrackerOrigin(const G4TouchableHandle & touchableHandle) {

    // how deep in the hierachy is the tracker a.k.a tracker depth
    // (depends on the tracker version)
    // for LTracker 1,2,3 it is 3
    // the tracker version is set in the constructor

    size_t td = 3;

    art::ServiceHandle<GeometryService> geom;
    if ( geom->hasElement<TTracker>() ) {
      td =_TrackerVersion +1;
    }

    //    cout << "Debugging: tracker depth/version: " << td << "/" << _TrackerVersion << endl;

    size_t hdepth = touchableHandle->GetHistoryDepth();

    G4ThreeVector cdo;

    for (size_t dd=0; dd!=hdepth; ++dd) {

      if (dd>=td) cdo += touchableHandle->GetVolume(dd)->GetTranslation();

//       cout << "Debugging: det depth name copy#: " << dd << " " <<
//         touchableHandle->GetVolume(dd)->GetName() << " " <<
//         touchableHandle->GetVolume(dd)->GetCopyNo() << endl;

//       G4LogicalVolume* lvp = touchableHandle->GetVolume(dd)->GetLogicalVolume();
//       G4int nd = lvp->GetNoDaughters();
//       for (G4int d = 0;d!=nd; ++d) {
//         cout << "Debugging: daughter: " << lvp->GetDaughter(d)->GetName() << " " <<
//           lvp->GetDaughter(d)->GetCopyNo() << endl;
//       }

//       cout << "Debugging det origin: " << touchableHandle->GetVolume(dd)->GetTranslation() <<
//         " " << cdo << endl;

    }

    return cdo;

  }

  void StrawSD::beforeG4Event(StepPointMCCollection& outputHits, 
                              PhysicsProcessInfo& processInfo,
                              art::ProductID const& simID, 
                              art::EDProductGetter const* productGetter ){
    _collection    = &outputHits;
    _processInfo   = &processInfo;
    _simID         = &simID;
    _productGetter = productGetter;
    return;
  } // end of beforeG4Event

} //namespace mu2e
