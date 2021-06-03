//
// Umbrela for the the Mu2e G4 world classes
//
//
// Original author Rob Kutschke
//
//

// C++ includes
#include <iostream>
#include <vector>
#include <iomanip>

// Framework includes
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "cetlib_except/exception.h"
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"

// Mu2e includes
#include "Mu2eG4Helper/inc/Mu2eG4Helper.hh"
#include "Mu2eG4/inc/Mu2eUniverse.hh"

// G4 includes
#include "Geant4/G4PhysicalVolumeStore.hh"

using namespace std;

namespace mu2e {

  Mu2eUniverse::Mu2eUniverse(const Mu2eG4Config::Debug& debug):
    _geom(*(art::ServiceHandle<GeometryService>())),
    _config(_geom.config()),
    _helper(&(*(art::ServiceHandle<Mu2eG4Helper>()))),
    _verbosityLevel(debug.worldVerbosityLevel()),
    _g4VerbosityLevel(debug.diagLevel())
  {} // beware of the order of initialization/declarations

  Mu2eUniverse::~Mu2eUniverse(){
  }

  // Convert to base units for all of the items in the vector.
  void Mu2eUniverse::setUnits( vector<double>& V, G4double unit ){
    for ( vector<double>::iterator b=V.begin(), e=V.end();
          b!=e; ++b){
      *b *= unit;
    }
  }

  // A helper function for debugging.  Print a subset of the physical volume store.
  void Mu2eUniverse::printPhys() {
    G4PhysicalVolumeStore* pstore = G4PhysicalVolumeStore::GetInstance();
    int n(0);
    for ( std::vector<G4VPhysicalVolume*>::const_iterator i=pstore->begin(); i!=pstore->end(); i++){
      cout << "Physical Volume: "
           << setw(5) << n++
           << (*i)->GetName()
           << endl;
      if ( n > 25 ) break;
    }

  }

} // end namespace mu2e
