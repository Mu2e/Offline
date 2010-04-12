//
// Free functions to add step limiters to some particles.
//
// $Id: addStepLimiter.cc,v 1.2 2010/04/12 13:11:44 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/04/12 13:11:44 $
//
// Original author Rob Kutschke
//
// Notes
// 1) About the iterator class:
//      G4ParticleTable::G4PTblDicIterator
//    The G4ParticleTable class internally contains a 
//      std::map<G4String,ParticleDefinition*>
//    Users of the class need to be able to loop over the map.  To provide 
//    this functionality the designers of the class chose to provide a 
//    custom written iterator.  ( Current best practice would be to provide
//    a typedef to the underlying stl iterator type.  It is possible that
//    this design was frozen before std::map was robust and that this class
//    provides a work around that once had been necessary. )
//
//    The rules for using this iterator are:
//     1) You must call reset() before starting the loop.
//        If you forget you will usually, but not always, start the iteration
//        past the end of the map, which produces undefined behaviour.
//     2) The call to operator ():
//          a) on the first call after reset, it will initialize an internal 
//             iterator to point at the first element of the map.
//          b) on subsequent calls it will increment the internal iterator.
//          c) On all calls it will return true(false) if the interal iterator
//             is !=(==) to the .end() of the map.
//     3) At any time, a call to key() or to value() will return the information
//        pointed to by the interally held iterator.
//
//     A particularly bizarre feature is that the G4ParticleTable class holds
//     its own iterator.  When you ask for an iterator, you actually get a pointer
//     to this one iterator.  Therefore nested loops over the map have totally
//     unexpected behaviour.

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
