#ifndef WLSStackingAction_h
#define WLSStackingAction_h 1

#include "G4UserStackingAction.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ClassificationOfNewTrack.hh"
#include "G4Track.hh"
#include "G4VProcess.hh"



class WLSStackingAction : public G4UserStackingAction
{
  public:

  virtual G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track* track)
  {
//These sections are used for testing purposes to remove certain types of photons.
//Don't use for production!
/*
    if(track->GetCreatorProcess()!=NULL && track->GetOriginTouchable()!=NULL)
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
*/

//This section can be left for production
    if(track->GetOriginTouchable()!=NULL)
    {
      if(track->GetOriginTouchable()->GetVolume()->GetName()=="World")
      {
        return(G4ClassificationOfNewTrack::fKill);
      }
    }

    return(G4UserStackingAction::ClassifyNewTrack(track));
  }
};

#endif
