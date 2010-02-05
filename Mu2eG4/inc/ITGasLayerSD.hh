#ifndef ITGasLayerSD_h
#define ITGasLayerSD_h 1

// Mu2e includes
#include "ITrackerGeom/inc/ITracker.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "Mu2eUtilities/inc/TwoLinePCA.hh"
#include "Mu2eG4/inc/StepPointG4.hh"

// G4 includes
#include "G4VSensitiveDetector.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

#include "CLHEP/config/CLHEP.h"

namespace mu2e {

  class ITGasLayerSD : public G4VSensitiveDetector{

  public:
    ITGasLayerSD(G4String);
    ~ITGasLayerSD();
    
    void Initialize(G4HCofThisEvent*);
    virtual G4bool ProcessHits(G4Step*, G4TouchableHistory*) { return false; }
    void EndOfEvent(G4HCofThisEvent*);

  protected:
    StepPointG4Collection* _collection;
    int _superlayer;
    int _ring;
    int _nwires;
    double _Dphi;
    
  };

} // namespace mu2e

#endif /*ITGasLayerSD_h*/
