//
// Override the G4RunManager class so that the Mu2e framework can drive
// the event loop.
//
// $Id: Mu2eG4RunManager.cc,v 1.6 2011/05/18 02:27:18 wb Exp $
// $Author: wb $
// $Date: 2011/05/18 02:27:18 $
//
// Original author Rob Kutschke
//
// This version is derived from Geant4's G4RunManager class with the following cvs info:
//    .cc: 2007/11/09, verson 1.108
//    .hh: 2009/11/13, version 1.52
//
// Notes:
// 1) In G4RunManager the counter i_event is used as the event number.
//    In this code it is taken from the event number of the art::event.
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
    _nProcessed(0),
    _realElapsed(0.),
    _systemElapsed(0.),
    _userElapsed(0.),
    _msg(""){
  }

  // Destructor of base is called automatically.  No need to do anything.
  Mu2eG4RunManager::~Mu2eG4RunManager(){
  }

  // Do the "begin run" parts of BeamOn.
  void Mu2eG4RunManager::BeamOnBeginRun( unsigned int runNumber, const char* macroFile, G4int n_select){

    SetRunIDCounter(runNumber);

    bool cond = ConfirmBeamOnCondition();
    if(!cond){
      // throw here
      return;
    }

    // Initialize member data for the new run.
    _macroFile     = macroFile;
    _n_select      = n_select;
    _msg           = "";
    _nProcessed    = 0;
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
  void Mu2eG4RunManager::BeamOnDoOneEvent( int eventNumber){

    timer->Start();

    // This is the body of the event loop from DoEventLoop().
    currentEvent = GenerateEvent(eventNumber);
    eventManager->ProcessOneEvent(currentEvent);
    AnalyzeEvent(currentEvent);
    UpdateScoring();

    if(_nProcessed<_n_select) G4UImanager::GetUIpointer()->ApplyCommand(_msg);

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

    ++_nProcessed;
  }

  // Do the "end of run" parts of DoEventLoop and BeamOn.
  void Mu2eG4RunManager::BeamOnEndRun(){

    // From G4RunManager::DoEventLoop
    if(verboseLevel>0){

      G4cout << "G4Run terminated." << G4endl;
      G4cout << "G4Run Summary" << G4endl;
      if(runAborted)
        { G4cout << "  G4Run Aborted after " << _nProcessed << " events processed." << G4endl; }
      else
        { G4cout << "  Number of events processed : " << _nProcessed << G4endl; }
      G4cout << "  User="  << _userElapsed
             << "s Real="  << _realElapsed
             << "s Sys="   << _systemElapsed
             << G4endl;
    }

    // From G4RunManager::BeamOn.
    RunTermination();
  }

} // end namespace mu2e
