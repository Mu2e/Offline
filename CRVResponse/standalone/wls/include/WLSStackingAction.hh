#ifndef WLSStackingAction_h
#define WLSStackingAction_h 1

#include "G4UserStackingAction.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ClassificationOfNewTrack.hh"
#include "G4Track.hh"
#include "G4VProcess.hh"


//This class is only used for testing purpose.
//Don't use for productioon!!!
//It can be used to disable/enable one or more of the produced "photon types".

class WLSStackingAction : public G4UserStackingAction
{
  public:

  virtual G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track* track)
  {
    if(track->GetCreatorProcess()!=NULL)
    {
      if(track->GetCreatorProcess()->GetProcessName()=="Scintillation")
      {
        return(G4ClassificationOfNewTrack::fKill);
      }
      if(track->GetCreatorProcess()->GetProcessName()=="Cerenkov" && track->GetOriginTouchable()->GetVolume()->GetName()=="Scintillator")
      {
        return(G4ClassificationOfNewTrack::fKill);
      }
      if(track->GetCreatorProcess()->GetProcessName()=="Cerenkov" && track->GetOriginTouchable()->GetVolume()->GetName()=="WLSFiber")
      {
        return(G4ClassificationOfNewTrack::fKill);
      }
    }
    return(G4UserStackingAction::ClassifyNewTrack(track));
  }
};

#endif
