#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1
//
// Give generated tracks to G4. This implements two algorithms:
//
// 1) testTrack - a trivial 1 track generator for debugging geometries.
// 2) fromEvent - copies generated tracks from the event.
//
// $Id: PrimaryGeneratorAction.hh,v 1.2 2010/03/13 00:09:16 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2010/03/13 00:09:16 $
//
// Original author Rob Kutschke
//

// C++ includes
#include <string>
#include <memory>

// Framework includes
#include "FWCore/Framework/interface/Event.h"

// G4 includes
#include "globals.hh"
#include "G4VUserPrimaryGeneratorAction.hh"

// Mu2e inclues
#include "Mu2eUtilities/inc/RandomUnitSphere.hh"

class G4ParticleDefinition;
class G4ParticleGun;
class G4Event;
class TH1D;

namespace mu2e {

  class SteppingAction;
  class WorldInfo;
  class Mu2eWorld;

  class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction{
  public:
    PrimaryGeneratorAction( const std::string& generatorModuleLabel);
    ~PrimaryGeneratorAction();
    
  public:
    
    // This is the interface specified by G4.
    void GeneratePrimaries(G4Event*);

    // Should change the interface for Primary 
    void setEvent( edm::Event const& event) {_event = &event;}

    void setWorld( Mu2eWorld const* world ){
      _world=world;
    }

  private:

    void fromEvent( G4Event* );
    void testTrack( G4Event* );

    // The event we are working on;
    // Must be set before the call to GeneratePrimaries.
    // Should change to a pull method, rather than a push.
    edm::Event const* _event;

    // The particle I am generating.
    G4ParticleDefinition* _particleDefinition;

    // Module label used to find the event generator input.
    std::string _generatorModuleLabel;

    // Generate random directions tracks on a unit sphere.
    std::auto_ptr<RandomUnitSphere> _randomUnitSphere;
    
    // Non-owing pointer to the detector information.
    Mu2eWorld const* _world;
    
    
    TH1D* _totalMultiplicity;
    
  };
  
}  // end namespace mu2e
#endif
