// This interface is used to share the same TrackingAction between
// G4_module and Mu2eG4_module, that use different stepping actions.
// Once G4_module is removed, this interface can be removed as well.
//
// Andrei Gaponenko, 2015

#ifndef Mu2eG4_IMu2eG4SteppingAction_hh
#define Mu2eG4_IMu2eG4SteppingAction_hh

#include <vector>
#include "CLHEP/Vector/LorentzVector.h"

namespace mu2e {

  class IMu2eG4SteppingAction {

  public:
    virtual ~IMu2eG4SteppingAction() {}

    virtual void BeginOfTrack() = 0;
    virtual void EndOfTrack() = 0;

    virtual std::vector<CLHEP::HepLorentzVector> const&  trajectory()  = 0;

    // swap the stored trajectory object with the given one
    virtual void swapTrajectory(std::vector<CLHEP::HepLorentzVector>& trajectory) = 0;
  };

} // end namespace mu2e

#endif /* Mu2eG4_IMu2eG4SteppingAction_hh */
