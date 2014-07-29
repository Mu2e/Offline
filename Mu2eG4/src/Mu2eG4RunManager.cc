//
// Override the G4RunManager class so that the Mu2e framework can drive
// the event loop.
//
// $Id: Mu2eG4RunManager.cc,v 1.9 2014/07/29 21:39:21 genser Exp $
// $Author: genser $
// $Date: 2014/07/29 21:39:21 $
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
#include "G4Run.hh"

// C++
#include <limits>

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

    //numberOfEventToBeProcessed should be the total number of events to be processed 
    // or a large number and NOT 1 for G4 to work properly
    numberOfEventToBeProcessed = std::numeric_limits<int>::max(); // largest int for now
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


    if(verboseLevel>1){

      G4int oldPrecision = G4cout.precision(3);
      std::ios::fmtflags oldFlags = G4cout.flags();
      G4cout.setf(std::ios::fixed,std::ios::floatfield); 

      G4cout << "TimeEvent> "
             << eventNumber << " "
             << G4RunManager::GetRunManager()->GetCurrentRun()->GetRunID() << " "
             << timer->GetRealElapsed() << " "
             << timer->GetUserElapsed() + timer->GetSystemElapsed() << " "
             << _userElapsed+_systemElapsed
             << G4endl;

      G4cout.setf(oldFlags);
      G4cout.precision(oldPrecision);

    }

  }

  void Mu2eG4RunManager::BeamOnEndEvent(){

    StackPreviousEvent(currentEvent);
    currentEvent = 0;

    ++_nProcessed;
  }

  // Do the "end of run" parts of DoEventLoop and BeamOn.
  void Mu2eG4RunManager::BeamOnEndRun(){

    // From G4RunManager::DoEventLoop
    if(verboseLevel>1){

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

      G4int oldPrecision = G4cout.precision(3);
      std::ios::fmtflags oldFlags = G4cout.flags();
      G4cout.setf(std::ios::fixed,std::ios::floatfield); 

      G4cout << "TimeReport> Time report complete in ";
      G4cout << _realElapsed;
      G4cout << " seconds"
             << G4endl;

      G4cout.setf(oldFlags);
      G4cout.precision(oldPrecision);

    }

    // From G4RunManager::BeamOn.
    RunTermination();
  }

} // end namespace mu2e
