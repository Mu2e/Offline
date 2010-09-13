#ifndef StepPointG4_h
#define StepPointG4_h 1
//
// A class to hold hits created by G4 in most types of sensitive detectors.
// A different class may be needed for calorimeter objects.  This wil
// be Ok for trackers, scintillators etc.
//
// For details, see: ToyDP/inc/StepPointMC.hh .
// 
// $Id: StepPointG4.hh,v 1.3 2010/09/13 23:43:58 logash Exp $
// $Author: logash $
// $Date: 2010/09/13 23:43:58 $
//
// Original author Rob Kutschke
//
// The core of the class is the persisent hit object of type StepPointMC.
// The rest of the class provides G4 related functionality.
//

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

#include "ToyDP/inc/StepPointMC.hh"

namespace mu2e { 
  
  class StepPointG4 : public G4VHit{
  public:
    
    inline
    StepPointG4():
      _hit(){
    }
    
    inline
    StepPointG4( G4int trackId,
                 G4int volumeId,
                 G4double edep,
                 G4ThreeVector const& position,
                 G4ThreeVector const& momentum,
                 G4double time,
                 G4double proper,
                 G4double stepLength
                 ):
      _hit(trackId,volumeId,edep,time,proper,position,momentum,stepLength){
    }
    
    // Accept compiler generated versions of:
    //  d'tor
    //  copy c'tor 
    //  assignment operator

    G4int operator==(const StepPointG4&) const;
    
    // G4 requires these.
    inline void* operator new(size_t);
    inline void  operator delete(void*);
  
    void Draw();
    void Print();

  public:
  
    inline mu2e::StepPointMC const&   hit()      const { return _hit; }
    inline G4int                      trackId()  const { return _hit.trackId();  }
    inline G4double                   eDep()     const { return _hit.eDep();     }
    inline G4ThreeVector const&       position() const { return _hit.position(); }
    inline G4ThreeVector const&       momentum() const { return _hit.momentum(); }
    inline G4double                   time()     const { return _hit.time();     }
    inline G4double                   properTime()     const { return _hit.properTime();     }
    inline StepPointMC::VolumeId_type volumeId() const { return _hit.volumeId(); }

  private:
  
    // This is the persistent hit object.
    mu2e::StepPointMC _hit;
    
  };


  typedef G4THitsCollection<StepPointG4> StepPointG4Collection;

  extern G4Allocator<StepPointG4> StepPointG4Allocator;

  inline void* StepPointG4::operator new(size_t){
    void *aHit;
    aHit = (void *) StepPointG4Allocator.MallocSingle();
    return aHit;
  }


  inline void StepPointG4::operator delete(void *aHit){
    StepPointG4Allocator.FreeSingle((StepPointG4*) aHit);
  }

} // namespace mu2e

#endif
