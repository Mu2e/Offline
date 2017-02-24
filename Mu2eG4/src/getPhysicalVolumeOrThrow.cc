//
// Free function wrapper around:
//
//  G4PhysicalVolumeStore::GetInstance ()->GetVolume("HallAir");
//
// The wrapper does the job of throwing if the pointer comes back null.
//
// $Id: getPhysicalVolumeOrThrow.cc,v 1.3 2011/05/18 02:27:18 wb Exp $
// $Author: wb $
// $Date: 2011/05/18 02:27:18 $
//
// Original author Rob Kutschke
//
//

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"

// Mu2e includes
#include "Mu2eG4/inc/getPhysicalVolumeOrThrow.hh"

// G4 includes
#include "G4String.hh"
#include "G4PhysicalVolumeStore.hh"

namespace mu2e {

  G4VPhysicalVolume* getPhysicalVolumeOrThrow( G4String const& name, bool mustHave ){

    G4VPhysicalVolume* pVol = G4PhysicalVolumeStore::GetInstance ()->GetVolume(name);

    if ( pVol == 0 ){
      throw cet::exception("GEOM")
        << "Could not find requested Physical Volume: "
        << name
        << "\n";
    }

    return pVol;

  }

} // end namespace mu2e
