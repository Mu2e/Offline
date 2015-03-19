// Andrei Gaponenko, 2015

#include "Mu2eG4/inc/Mu2eG4StackingAction.hh"

#include <vector>
#include <string>

#include "ConfigTools/inc/SimpleConfig.hh"
#include "cetlib/exception.h"

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

  void Mu2eG4StackingAction::NewStage(){ }

  void Mu2eG4StackingAction::PrepareNewEvent(){ }

  //================================================================
  void Mu2eG4StackingAction::checkConfigRelics(const SimpleConfig& config) {
    static const std::vector<std::string> keys = {
      "g4.doCosmicKiller",
      "g4.cosmicKillLevel",
      "g4.cosmicVerbose",
      "g4.cosmicPcut",
      "g4.yaboveDirtYmin",
      "g4.stackPrimaryOnly",
      "g4.killLowEKine",
      "g4.killPitchToLowToStore",
      "g4.minPitch",
      "g4.stackingActionDropPDG",
      "g4.stackingActionKeepPDG",
      "g4.eKineMin",
      "g4.killLowEKinePDG",
      "g4.eKineMinPDG"
    };

    std::string present;
    for(const auto k: keys) {
      if(config.hasName(k)) {
        present += k+" ";
      }
    }
    if(!present.empty()) {
      throw cet::exception("CONFIG")<<"Mu2eG4StackingAction: Please use fcl to configure Mu2eG4_module. "
                                    <<"Detected obsolete SimpleConfig parameters: "<<present;
    }
  }



} // end namespace mu2e
