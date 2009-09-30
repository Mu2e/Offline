#ifndef SteppingAction_h
#define SteppingAction_h 1
//
// Called at every G4 step.
//
// $Id: SteppingAction.hh,v 1.1 2009/09/30 22:57:47 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2009/09/30 22:57:47 $
//
// Original author Rob Kutschke
//

#include "G4UserSteppingAction.hh"
#include "G4ThreeVector.hh"

namespace mu2e {
  class SteppingAction : public G4UserSteppingAction{
  public:
    SteppingAction();
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
    
  };
  
} // end namespace mu2e
#endif
