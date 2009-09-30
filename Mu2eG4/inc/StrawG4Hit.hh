#ifndef StrawG4Hit_h
#define StrawG4Hit_h 1
//
// A class to hold hits created by G4 in straw detectors.
// These hits are made by the sensitive detector class for straws, StrawSD.
// 
// $Id: StrawG4Hit.hh,v 1.1 2009/09/30 22:57:47 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2009/09/30 22:57:47 $
//
// Original author Rob Kutschke
//
// The core of the class is the persisent hit object of type StrawMCHit.
// The rest of the class provides G4 related functionality.
//

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

#include "ToyDP/inc/StrawMCHit.hh"

namespace mu2e { 

  class StrawG4Hit : public G4VHit{
  public:
  
    StrawG4Hit():
      _hit(){
    }
    
    StrawG4Hit( G4int trackId,
		G4int strawIndex,
		G4double edep,
		G4ThreeVector const& position,
		G4ThreeVector const& momentum,
		G4double time
		):
      _hit(trackId,strawIndex,edep,time,position,momentum){
    }
    
    ~StrawG4Hit(){
    }
    
    // Compiler generated copy c'tor and assignment operator will be OK.

    G4int operator==(const StrawG4Hit&) const;
    
    // G4 requires these.
    inline void* operator new(size_t);
    inline void  operator delete(void*);
  
    void Draw();
    void Print();
  public:
  
    mu2e::StrawMCHit const& hit() const { return _hit; }

    G4int trackId()          const { return _hit.trackId();    }
    G4int strawIndex()       const { return _hit.strawIndex(); }
    G4double eDep()          const { return _hit.eDep();       }
    G4ThreeVector position() const { return _hit.position();   }
    G4ThreeVector momentum() const { return _hit.momentum();   }
    G4double time()          const { return _hit.time();       }

  private:
  
    // This is the persistent hit object.
    mu2e::StrawMCHit _hit;
    
  };


  typedef G4THitsCollection<StrawG4Hit> StrawG4HitsCollection;

  extern G4Allocator<StrawG4Hit> StrawG4HitAllocator;

  inline void* StrawG4Hit::operator new(size_t){
    void *aHit;
    aHit = (void *) StrawG4HitAllocator.MallocSingle();
    return aHit;
  }


  inline void StrawG4Hit::operator delete(void *aHit){
    StrawG4HitAllocator.FreeSingle((StrawG4Hit*) aHit);
  }

} // namespace mu2e

#endif
