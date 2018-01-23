//
// Override the G4RunManager class so that the Mu2e framework can drive
// the event loop.
//
// Original author Lisa Goodenough
//
//
// Notes:
//
// Implementation file for Mu2eG4RunManager



#include "Mu2eG4/inc/Mu2eG4RunManager.hh"


using namespace std;

namespace mu2e {
  
    // If the c'tor is called a second time, the c'tor of base will
    // generate an exception.
    Mu2eG4RunManager::Mu2eG4RunManager():
        G4RunManager()
        {}
  
    // Destructor of base is called automatically.  No need to do anything.
    Mu2eG4RunManager::~Mu2eG4RunManager()
        {}
    
    void Mu2eG4RunManager::Mu2eG4WaitForEndEventLoopWorkers()
    {
        //nothing to do in sequential mode
    }
    
    void Mu2eG4RunManager::Mu2eG4TerminateWorkers()
    {
        //nothing to do in sequential mode
    
    }
    
    void Mu2eG4RunManager::SetNumberOfThreads(G4int num_threads)
    {
        
        //nothing to do in sequential mode
    }
  
  
  
} // end namespace artg4
