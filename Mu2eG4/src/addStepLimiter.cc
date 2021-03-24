//
// Free functions to add step limiters to some particles.
//
//
// Original author Rob Kutschke
//
// Notes
// 1) For notes about the iterator class:
//      G4ParticleTable::G4PTblDicIterator
//    See http://mu2e.fnal.gov/atwork/computing/G4Notes.shtml .

// C++ includes
#include <vector>
#include <string>

// Mu2e includes
#include "Mu2eG4/inc/addStepLimiter.hh"

// G4 includes
#include "Geant4/G4ParticleTable.hh"
#include "Geant4/G4ProcessManager.hh"
#include "Geant4/G4StepLimiter.hh"

using namespace std;

namespace mu2e{

  // Add step limiters for a standard list of particles:
  // e, mu, CLHEP::pi, K, p (particles and anti-particles) plus chargedgeantino.
  void addStepLimiter (){

    // Create the standard list.
    vector<G4String> list;
    list.push_back( "e+"  );
    list.push_back( "e-"  );
    list.push_back( "mu+" );
    list.push_back( "mu-" );
    list.push_back( "pi+" );
    list.push_back( "pi-" );
    list.push_back( "kaon+"  );
    list.push_back( "kaon-"  );
    list.push_back( "proton" );
    list.push_back( "anti_proton"     );
    list.push_back( "chargedgeantino" );
    list.push_back( "deuteron" );
    list.push_back( "anti_deuteron" );
    list.push_back( "triton" );
    list.push_back( "anti_triton" );
    list.push_back( "He3" );
    list.push_back( "anti_He3" );
    list.push_back( "alpha" );
    list.push_back( "anti_alpha" );

    // Do the work for particles in this list.
    addStepLimiter( list );
  }

  // Add step limiters for a user specified list of particles.
  void addStepLimiter (const vector<G4String>& list ){

    // Get an iterator over existing particles.
    G4ParticleTable* ptable = G4ParticleTable::GetParticleTable();
    G4ParticleTable::G4PTblDicIterator* iter = ptable->GetIterator();

    // See note 1.
    iter->reset();

    // Check each existing particle to see if it is in the list.  See note 1.
    while( (*iter)() ){
      G4ParticleDefinition* particle = iter->value();
      G4ProcessManager* pmanager     = particle->GetProcessManager();
      G4String particleName          = particle->GetParticleName();

      // Is this particle in the list?
      if ( find( list.begin(), list.end(), particleName ) != list.end() ){

        // The process manager takes ownership of the G4StepLimiter object.
        pmanager->AddDiscreteProcess(new G4StepLimiter);
      }

    }
  }


}  // end namespace mu2e
