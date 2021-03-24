//
// Free function wrapper around
//   G4NistManager::FindOrBuildMaterial
// The wrapper does the job of throwing if the pointer comes
// back null.
//
//
// Original author Rob Kutschke
//
//

#include "Mu2eG4/inc/findMaterialOrThrow.hh"

// Framework includes
#include "cetlib_except/exception.h"

// G4 includes
#include "Geant4/G4String.hh"
#include "Geant4/G4NistManager.hh"

namespace mu2e {
  G4Material* findMaterialOrThrow( G4String const& name){

    // Look up the material.
    G4Material* m = G4NistManager::Instance()->
      FindOrBuildMaterial(name,true,true);

    // Throw if necessary.
    if ( !m ){
      throw cet::exception("GEOM")
        << "Could not find a material with the name: "
        << name
        << "\n";
    }
    return m;
  }
}
