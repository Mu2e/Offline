#ifndef Mu2eG4_MTRunManager_hh
#define Mu2eG4_MTRunManager_hh
//
// Override the G4MTRunManager class so that the Mu2eG4 framework can drive
// the event loop.
//
// Original author Lisa Goodenough
//


// Included from Geant4
#include "G4MTRunManager.hh"

namespace mu2e {
  
    class Mu2eG4MTRunManager : public G4MTRunManager{
    
  public:
    
    Mu2eG4MTRunManager();
    virtual ~Mu2eG4MTRunManager();
    
    //we need our own versions of these functions in order to correctly control the event loop
    void Mu2eG4Initialize(G4int);
    void Mu2eG4InitializeEventLoop(G4int);
    void Mu2eG4WaitForEndEventLoopWorkers();
    void Mu2eG4RunTermination();
        
    static Mu2eG4MTRunManager* GetMasterRunManager();
        
    
    void TestFunc();
    
        
        
  private:
    
    // Private and unimplemented to prevent copying.
    Mu2eG4MTRunManager( Mu2eG4MTRunManager const & );
    Mu2eG4MTRunManager& operator=( Mu2eG4MTRunManager const & );
        
    //Singleton implementing master thread behavior
    static Mu2eG4MTRunManager* fMu2eMasterRM;
        
        
    //G4MTRunManagerKernel* MTkernel;
    G4int my_nworkers;
    
  };
  
} // end namespace mu2e

#endif /* Mu2eG4_MTRunManager_hh */
