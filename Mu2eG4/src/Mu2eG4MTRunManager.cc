//
// Override the G4RunManager class so that the Mu2e framework can drive
// the event loop.
//
// Original author Lisa Goodenough
//
//
// Notes:
//
// Implementation file for Mu2eG4MTRunManager



#include "Mu2eG4/inc/Mu2eG4MTRunManager.hh"


using namespace std;

namespace mu2e {
  
    // If the c'tor is called a second time, the c'tor of base will
    // generate an exception.
    Mu2eG4MTRunManager::Mu2eG4MTRunManager():
        G4MTRunManager()
        {}
  
    // Destructor of base is called automatically.  No need to do anything.
    Mu2eG4MTRunManager::~Mu2eG4MTRunManager()
        {}
    
    //this function is a protected member of G4MTRunManager but we need to access it
    //from Mu2eG4_module, so we must make it piblic here
    void Mu2eG4MTRunManager::Mu2eG4WaitForEndEventLoopWorkers()
    {
        WaitForEndEventLoopWorkers();
    }
  
} // end namespace artg4
