#ifndef Mu2eG4_RunManager_hh
#define Mu2eG4_RunManager_hh
//
// Override the G4RunManager class so that the Mu2eG4 framework can drive
// the event loop.
//
// Original author Lisa Goodenough
//


// Included from Geant4
#include "G4RunManager.hh"

namespace mu2e {
  
    class Mu2eG4RunManager : public G4RunManager{
    
  public:
    
    Mu2eG4RunManager();
    virtual ~Mu2eG4RunManager();
    
    void Mu2eG4WaitForEndEventLoopWorkers();
    void Mu2eG4TerminateWorkers();
    void SetNumberOfThreads(G4int);
    
    
  private:
    
    // Private and unimplemented to prevent copying.
    Mu2eG4RunManager( Mu2eG4RunManager const & );
    Mu2eG4RunManager& operator=( Mu2eG4RunManager const & );
    
  };
  
} // end namespace mu2e

#endif /* Mu2eG4_RunManager_hh */
