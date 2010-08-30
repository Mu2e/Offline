#ifndef SteppingAction_h
#define SteppingAction_h 1
//
// Called at every G4 step.
//
// $Id: SteppingAction.hh,v 1.3 2010/08/30 22:21:13 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2010/08/30 22:21:13 $
//
// Original author Rob Kutschke
//

// Mu2e includes
#include "Mu2eG4/inc/EventNumberList.hh"

// G4 includes
#include "G4UserSteppingAction.hh"
#include "G4ThreeVector.hh"


namespace mu2e {

  // Forward declarations in mu2e namespace
  class SimpleConfig;

  class SteppingAction : public G4UserSteppingAction{

  public:
    SteppingAction( const SimpleConfig& config );
    ~SteppingAction(){};
    
    void UserSteppingAction(const G4Step*);

    G4ThreeVector const& lastPosition() const { return _lastPosition; }
    G4ThreeVector const& lastMomentum() const { return _lastMomentum; }

    void setZRef( G4double zref){ 
      _zref=zref;
    }
  
  private:
    G4ThreeVector _lastPosition;
    G4ThreeVector _lastMomentum;
    
    G4double _zref;
    
    // Lists of events and tradcks for which to enable debug printout.
    EventNumberList _debugEventList;
    EventNumberList _debugTrackList;

  };
  
} // end namespace mu2e
#endif
