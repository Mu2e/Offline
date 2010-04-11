//
// Free functions to add step limiters to some particles.
//
// $Id: addStepLimiter.cc,v 1.1 2010/04/11 15:15:12 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/04/11 15:15:12 $
//
// Original author Rob Kutschke
//

// C++ includes
#include <vector>
#include <string>

// Mu2e includes
#include "Mu2eG4/inc/addStepLimiter.hh"

// G4 includes
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"
#include "G4StepLimiter.hh"

using namespace std;

namespace mu2e{

  // Add step limiters for a standard list of particles:
  // e, mu, pi, K, p (particles and anti-particles) plus chargedgeantino.
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

    // Do the work for particles in this list.
    addStepLimiter( list );
  }

  // Add step limiters for a user specified list of particles.
  void addStepLimiter (const vector<G4String>& list ){

    // Get an iterator over existing particles.
    G4ParticleTable* ptable = G4ParticleTable::GetParticleTable();
    G4ParticleTable::G4PTblDicIterator* iter = ptable->GetIterator();

    // This step is required.
    iter->reset();

    // Check each existing particle to see if it is in the list.
    while( (*iter)() ){
      G4ParticleDefinition* particle = iter->value();
      G4ProcessManager* pmanager     = particle->GetProcessManager();
      G4String particleName            = particle->GetParticleName();

      // Is this particle in the list?
      if ( find( list.begin(), list.end(), particleName ) != list.end() ){

        // The process manager takes ownership of the G4StepLimiter object.
        pmanager->AddDiscreteProcess(new G4StepLimiter);
      }

    }
  }
  

}  // end namespace mu2e
