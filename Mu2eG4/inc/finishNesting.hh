#ifndef Mu2eG4_finishNesting_hh
#define Mu2eG4_finishNesting_hh
//
// Free function to be used by the nest... functions
//
//
// Original author KLG
//

#include <string>

#include "Mu2eG4Helper/inc/VolumeInfo.hh"
#include "GeometryService/inc/G4GeometryOptions.hh"
#include "GeometryService/inc/GeometryService.hh"

#include "Geant4/G4Colour.hh"
#include "Geant4/G4ThreeVector.hh"
#include "Geant4/G4RotationMatrix.hh"

class G4Material;
class G4LogicalVolume;

namespace mu2e {

  void finishNesting(VolumeInfo& info,
                     G4Material* material,
                     G4RotationMatrix const* rot,
                     G4ThreeVector const& offset,
                     G4LogicalVolume* parent,
                     int copyNo,
                     bool const isVisible,
                     G4Colour const& color,
                     bool const forceSolid,
                     bool const forceAuxEdgeVisible,
                     bool const placePV,
                     bool const doSurfaceCheck,
                     bool const verbose = false
                     );

  // Version that accesses G4GeometryOptions
  // ...inline to avoid violating ODR
  inline void finishNesting(VolumeInfo& info,
                            G4Material* material,
                            G4RotationMatrix const* rot,
                            G4ThreeVector const& offset,
                            G4LogicalVolume* parent,
                            int copyNo,
                            G4Colour const& color,
                            std::string const& lookupToken = "",
                            bool const verbose = false
                            ) 
  {

    const std::string& lookupString = lookupToken.empty() ? info.name : lookupToken;
      
    G4GeometryOptions* geomOptions = art::ServiceHandle<GeometryService>()->geomOptions();

    finishNesting( info,
                   material,
                   rot,
                   offset,
                   parent,
                   copyNo,
		   geomOptions->isVisible( lookupString ),
		   color,
		   geomOptions->isSolid( lookupString ),
		   geomOptions->forceAuxEdgeVisible( lookupString ),
		   geomOptions->placePV( lookupString ),
		   geomOptions->doSurfaceCheck( lookupString ),
		   verbose );
  }

}

#endif /* Mu2eG4_finishNesting_hh */
