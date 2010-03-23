//
// Define a sensitive detector for Straws.
// ( Not sure yet if I can use this for both LTracker and TTracker?)
// 
// $Id: StrawSD.cc,v 1.4 2010/03/23 20:29:48 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2010/03/23 20:29:48 $
//
// Original author Rob Kutschke
//

#include <cstdio>

// Mu2e incldues
#include "Mu2eG4/inc/StrawSD.hh"
#include "Mu2eG4/inc/EventNumberList.hh"
#include "LTrackerGeom/inc/LTracker.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "Mu2eUtilities/inc/TwoLinePCA.hh"
#include "Mu2eUtilities/inc/LinePointPCA.hh"

// G4 includes
#include "G4RunManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

//
// Outstanding questions:
//
// 1) Why is diffAngle so big?
//

using namespace std;

namespace mu2e {

  StrawSD::StrawSD(G4String name) :G4VSensitiveDetector(name){
    G4String HCname;
    collectionName.insert(HCname="StepPointG4Collection");
  }


  StrawSD::~StrawSD(){ }

  void StrawSD::Initialize(G4HCofThisEvent* HCE){

    _collection = new StepPointG4Collection
      (SensitiveDetectorName,collectionName[0]); 
    static G4int HCID = -1;
    if(HCID<0){ 
      HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]); 
    }
    HCE->AddHitsCollection( HCID, _collection ); 

  }
  

  G4bool StrawSD::ProcessHits(G4Step* aStep,G4TouchableHistory*){

    G4Event const* event = G4RunManager::GetRunManager()->GetCurrentEvent();

    G4double edep = aStep->GetTotalEnergyDeposit();

    // Eventually we will want this but not now.
    //if(edep==0.) return false;

    // Origin of the LTracker.  Need to get this from G4.
    static G4ThreeVector detectorOrigin( -3904., -7350., 6200.);

    // Position at start of step point, in world system and in
    // a system in which the center of the tracking detector is the origin.
    G4ThreeVector prePosWorld   = aStep->GetPreStepPoint()->GetPosition();
    G4ThreeVector prePosTracker = prePosWorld - detectorOrigin;

    G4ThreeVector preMomWorld = aStep->GetPreStepPoint()->GetMomentum();

    StepPointG4* newHit = 
      new StepPointG4(aStep->GetTrack()->GetTrackID()-1,
		      aStep->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(),
		      edep,
		      prePosTracker,
		      preMomWorld,
		      aStep->GetPreStepPoint()->GetGlobalTime()
		      );
    
    // The collection takes ownership of the hit. 
    _collection->insert( newHit );

    // Some debugging tests.

    //
    EventNumberList list;
    if ( !list.inList() ) return true;

    // Transformations between world to local coordinate systems.
    G4AffineTransform const& toLocal = aStep->GetPreStepPoint()->GetTouchableHandle()->GetHistory()->GetTopTransform();
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
    G4double dot = dT.unit().dot(pT.unit());
    G4double angle = dT.angle(pT);

    G4ThreeVector dTLocal( deltaLocal.x(), deltaLocal.y(), 0.);
    G4ThreeVector pTLocal( preMomLocal.x(), preMomLocal.y(), 0. );
    G4double dotLocal = dTLocal.unit().dot(pTLocal.unit());
    G4double angleLocal = dTLocal.angle(pTLocal);

    // This is took big. O(1.e-5 radians) or about 1% of the value. Why?
    G4double diffAngle = angle-angleLocal;

    G4ThreeVector localOrigin(0.,0.,0.);
    G4ThreeVector worldOrigin = toWorld.TransformPoint(localOrigin) - detectorOrigin;

    G4ThreeVector localZUnit(0.,0.,1.);
    G4ThreeVector worldZUnit = toWorld.TransformAxis(localZUnit);

    int copy = aStep->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber();

    int eventNo = event->GetEventID();

    /*
    G4cout << "Addhit: "
	   << setw(4) << eventNo << " "
	   << setw(4) << _collection->entries() << " " 
	   << setw(6) << copy <<  "  | "
	   << std::setprecision(6)
	   << prePosTracker << " | "
	   << preMomWorld 
	   << std::setprecision(6)
	   << endl;
    */
 
    /*
    G4cout << "Add hit: "
	   << std::setprecision(8)
	   << _collection->entries() << " "
	   << aStep->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber() << " "
	   << aStep->GetTrack()->GetTrackID() << " "
	   << prePosWorld  << " "
	   << preMomWorld   << " | "
	   << angle << "   "
	   << angleLocal << "   "
	   << diffAngle << "  | "
	   << endl; 

    G4cout << "Line: " 
	   << prePosLocal.x()  << " "
	   << prePosLocal.y()  << " "
	   << postPosLocal.x()     << " "
	   << postPosLocal.y()     << " "
	   << preMomLocal.unit().x()    << " "
	   << preMomLocal.unit().y()    << " "
	   << worldOrigin.x() << " "
	   << worldOrigin.y() << " "
	   << std::setprecision(6)
	   << endl;
    */



    // Reco Geometry for the LTracker.
    GeomHandle<LTracker> ltracker;
    Straw const& straw = ltracker->getStraw( StrawIndex(copy) );
    G4ThreeVector mid = straw.getMidPoint();
    G4ThreeVector w   = straw.getDirection();

    // Point of closest approach of track to wire.
    // Straight line approximation.
    TwoLinePCA pca( mid, w, prePosTracker, preMomWorld);

    
    // Point on the wire that is closest to the step point.
    LinePointPCA lppca( mid, w, prePosTracker);
    double ddd = lppca.unit().cosTheta(preMomWorld);

    double ttt = lppca.unit().cosTheta(w);

    printf ( "Addhit: %4d %4d %6d %3d %3d | %10.2f %10.2f %10.2f | %10.2f %10.2f %10.2f | %6.3f %10.7f | %10.7f %10.7f\n",
	     eventNo,  _collection->entries(), copy,
	     aStep->IsFirstStepInVolume(), aStep->IsLastStepInVolume(),
	     prePosTracker.x(), prePosTracker.y(), prePosTracker.z(), 
	     preMomWorld.x(),   preMomWorld.y(),   preMomWorld.z(), ddd, ttt,
             prePosLocal.perp(),  postPosLocal.perp()  );
    fflush(stdout);

    /*
    G4cout << "Measurement: "
	   << aStep->GetTrack()->GetTrackID() << " "
	   << pca.dca()  << " "
	   << pca.dca2d()  << " "
	   << endl;
    */

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

    

    // End of debug section.

    //newHit->Print();
    //newHit->Draw();
    
    return true;
  }

  void StrawSD::EndOfEvent(G4HCofThisEvent*){

    if (verboseLevel>0) { 
      G4int NbHits = _collection->entries();
      G4cout << "\n-------->Hits Collection: in this event they are " << NbHits 
	     << " hits in the straw chambers: " << G4endl;
      for (G4int i=0;i<NbHits;i++) (*_collection)[i]->Print();
    } 
  }
  
} //namespace mu2e

