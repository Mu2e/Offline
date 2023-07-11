// Andrei Gaponenko, 2015

#ifndef Mu2eG4_Mu2eG4StackingAction_hh
#define Mu2eG4_Mu2eG4StackingAction_hh

#include "Geant4/G4UserStackingAction.hh"

#include "Offline/Mu2eG4/inc/IMu2eG4Cut.hh"

namespace mu2e {

  class Mu2eG4StackingAction: public G4UserStackingAction{
  public:
    Mu2eG4StackingAction(IMu2eG4Cut& stackingCuts,
                         IMu2eG4Cut& commonCuts);

    G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track* aTrack) override;

  private:
    // owned by Mu2eG4 module.
    IMu2eG4Cut* stackingCuts_ = nullptr;
    IMu2eG4Cut* commonCuts_ = nullptr;
  };

} // end namespace mu2e

#endif /* Mu2eG4_Mu2eG4StackingAction_hh */
