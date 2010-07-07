//
// Override the G4RunManager class so that the Mu2e framework can drive
// the event loop. 
//
// $Id: Mu2eG4RunManager.cc,v 1.3 2010/07/07 16:58:35 genser Exp $
// $Author: genser $ 
// $Date: 2010/07/07 16:58:35 $
//
// Original author Rob Kutschke
//
// This version is derived from Geant4's G4RunManager class with the following cvs info:
//    .cc: 2007/11/09, verson 1.108
//    .hh: 2007/11/13, version 1.51
//


// Mu2e includes.
#include "Mu2eG4/inc/Mu2eG4RunManager.hh"

// Includes from G4.
#include "G4UImanager.hh"
#include "G4ScoringManager.hh"
#include "G4Timer.hh"

using namespace std;

namespace mu2e {

  // If the c'tor is called a second time, the c'tor of base will
  // generate an exception.
  Mu2eG4RunManager::Mu2eG4RunManager():
    G4RunManager(),
    _macroFile(0),
    _n_select(-1),
    _i_event(0),
    _realElapsed(0.),
    _systemElapsed(0.),
    _userElapsed(0.),
    _msg(""){
  }

  // Destructor of base is called automatically.  No need to do anything.
  Mu2eG4RunManager::~Mu2eG4RunManager(){
  }

  // Do the "begin run" parts of BeamOn.
  void Mu2eG4RunManager::BeamOnBeginRun(const char* macroFile, G4int n_select){
  
    bool cond = ConfirmBeamOnCondition();
    if(!cond){
      // throw here
      return;
    }

    // Initialize member data for the new run.
    _macroFile     = macroFile;
    _n_select      = n_select;
    _msg           = "";
    _i_event       = 0;
    _realElapsed   = 0;
    _systemElapsed = 0;
    _userElapsed   = 0;

    // Construct command to execute the macro.
    if(_macroFile!=0){ 
      _msg = "/control/execute ";
      _msg += _macroFile;
    } else{ 
      _n_select = -1; 
    }

    numberOfEventToBeProcessed = 1;
    ConstructScoringWorlds();
    RunInitialization();

  }

  // Do the "per event" part of DoEventLoop.
  void Mu2eG4RunManager::BeamOnDoOneEvent(){

    timer->Start();

    // This is the body of the event loop from DoEventLoop().
    currentEvent = GenerateEvent(_i_event);
    eventManager->ProcessOneEvent(currentEvent);
    AnalyzeEvent(currentEvent);
    UpdateScoring();

    if(_i_event<_n_select) G4UImanager::GetUIpointer()->ApplyCommand(_msg);

    // Should pause, not stop, if I can do that.
    timer->Stop();

    // Accumulate time spent in G4 for all events in this run.
    _realElapsed   += timer->GetRealElapsed();
    _systemElapsed += timer->GetSystemElapsed();
    _userElapsed   += timer->GetUserElapsed();

  }

  void Mu2eG4RunManager::BeamOnEndEvent(){

    StackPreviousEvent(currentEvent);
    currentEvent = 0;

    ++_i_event;
  }

  // Do the "end of run" parts of DoEventLoop and BeamOn.
  void Mu2eG4RunManager::BeamOnEndRun(){

    // From G4RunManager::DoEventLoop
    if(verboseLevel>0){

      G4cout << "G4Run terminated." << G4endl;
      G4cout << "G4Run Summary" << G4endl;
      if(runAborted)
        { G4cout << "  G4Run Aborted after " << _i_event << " events processed." << G4endl; }
      else
        { G4cout << "  Number of events processed : " << _i_event << G4endl; }
      G4cout << "  User="  << _userElapsed 
             << "s Real="  << _realElapsed 
             << "s Sys="   << _systemElapsed 
             << G4endl;
    }

    // From G4RunManager::BeamOn.
    RunTermination();
  }

//   // Copied verbatim from G4RunManager::UpdateScoring().
//   // Does all of its work via singletons and public/protected objects in the 
//   // base.  So it should work identically.  No idea why it is private in the base.
//   void Mu2eG4RunManager::UpdateScoring()
//   {
//     G4ScoringManager* ScM = ::G4ScoringManager::GetScoringManagerIfExist();
//     if(!ScM) return;
//     G4int nPar = ScM->GetNumberOfMesh();
//     if(nPar<1) return;

//     G4HCofThisEvent* HCE = currentEvent->GetHCofThisEvent();
//     if(!HCE) return;
//     G4int nColl = HCE->GetCapacity();
//     for(G4int i=0;i<nColl;i++)
//       {
//         G4VHitsCollection* HC = HCE->GetHC(i);
//         if(HC) ScM->Accumulate(HC);
//       }
//   }

} // end namespace mu2e
