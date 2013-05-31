#ifndef Mu2eG4_PrimaryGeneratorAction_hh
#define Mu2eG4_PrimaryGeneratorAction_hh
//
// Give generated tracks to G4 by copying information from a GenParticleCollection.
//
// $Id: PrimaryGeneratorAction.hh,v 1.10 2013/05/31 22:50:47 gandr Exp $
// $Author: gandr $
// $Date: 2013/05/31 22:50:47 $
//
// Original author Rob Kutschke
//

// C++ includes
#include <string>

// Framework includes
#include "art/Framework/Principal/Event.h"

// G4 includes
#include "globals.hh"
#include "G4VUserPrimaryGeneratorAction.hh"

class G4ParticleDefinition;
class G4ParticleGun;
class G4Event;
class TH1D;

namespace mu2e {

  class SteppingAction;

  class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction{
  public:
    PrimaryGeneratorAction( const std::string& generatorModuleLabel);
    ~PrimaryGeneratorAction();

  public:

    // This is the interface specified by G4.
    void GeneratePrimaries(G4Event*);

    // Should change the interface for Primary
    void setEvent( art::Event const& event) {_event = &event;}

  private:

    void fromEvent( G4Event* );

    // The event we are working on;
    // Must be set before the call to GeneratePrimaries.
    // Should change to a pull method, rather than a push.
    art::Event const* _event;

    // Module label used to find the event generator input.
    std::string _generatorModuleLabel;

    TH1D* _totalMultiplicity;

  };

}  // end namespace mu2e
#endif /* Mu2eG4_PrimaryGeneratorAction_hh */
