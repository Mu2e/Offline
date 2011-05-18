#ifndef Mu2eG4_ITGasLayerSD_hh
#define Mu2eG4_ITGasLayerSD_hh

//// Mu2e includes
#include "ITrackerGeom/inc/ITracker.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
//#include "Mu2eUtilities/inc/TwoLinePCA.hh"
//#include "Mu2eG4/inc/StepPointG4.hh"
//
//// G4 includes
//#include "G4VSensitiveDetector.hh"
//#include "G4HCofThisEvent.hh"
//#include "G4Step.hh"
//#include "G4ThreeVector.hh"
//#include "G4SDManager.hh"
//#include "G4ios.hh"
//
//#include "CLHEP/Units/PhysicalConstants.h"
// Mu2e includes
#include "Mu2eG4/inc/EventNumberList.hh"
#include "ToyDP/inc/StepPointMCCollection.hh"
//#include "GeometryService/inc/GeometryService.hh"
//#include "GeometryService/inc/GeomHandle.hh"

// G4 includes
#include "G4VSensitiveDetector.hh"

class G4Step;
class G4HCofThisEvent;

namespace mu2e {

  // Forward declarations in mu2e namespace
  class SimpleConfig;

  class ITGasLayerSD : public G4VSensitiveDetector{

  public:
    ITGasLayerSD(G4String, const SimpleConfig& config);
    ~ITGasLayerSD();

    void Initialize(G4HCofThisEvent*);
    virtual G4bool ProcessHits(G4Step*, G4TouchableHistory*) { return false; }
//    G4bool ProcessHits(G4Step*, G4TouchableHistory*) { return false; }
    void EndOfEvent(G4HCofThisEvent*);

    void beforeG4Event(StepPointMCCollection& outputHits);

    static void setMu2eDetCenterInWorld(const G4ThreeVector &origin) {
            _mu2eDetCenter = origin;
    }

  protected:
//    StepPointG4Collection* _collection;
    int _superlayer;
    int _ring;
    int _nwires;
    double _Dphi;

    StepPointMCCollection* _collection;

    // List of events for which to enable debug printout.
    EventNumberList _debugList;

    GeomHandle<ITracker> itracker;
    //ITracker::GeomType _ittype;

    // Limit maximum size of the steps collection
    int _sizeLimit;
    int _currentSize;

    // Mu2e point of origin
    static G4ThreeVector _mu2eDetCenter;

  };

} // namespace mu2e

#endif /* Mu2eG4_ITGasLayerSD_hh */
