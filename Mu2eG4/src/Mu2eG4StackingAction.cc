// Andrei Gaponenko, 2015

#include "Mu2eG4/inc/Mu2eG4StackingAction.hh"

#include <vector>
#include <string>
#include <iostream>

namespace mu2e {

  Mu2eG4StackingAction::Mu2eG4StackingAction(const fhicl::ParameterSet&,
                                             IMu2eG4Cut& stackingCuts,
                                             IMu2eG4Cut& commonCuts)
    : stackingCuts_(&stackingCuts)
    , commonCuts_(&commonCuts)
  {}
    
    G4ClassificationOfNewTrack Mu2eG4StackingAction::ClassifyNewTrack(const G4Track* trk){
        if(stackingCuts_->stackingActionCut(trk)) {
            return fKill;
        }
        if(commonCuts_->stackingActionCut(trk)) {
            return fKill;
        }
        return fUrgent;
    }



} // end namespace mu2e
