#ifndef Mu2eG4_nestPolycone_hh
#define Mu2eG4_nestPolycone_hh
//
// Free function to create and place a new G4Polycone, place inside a logical volume.
//
//
// Original author Rob Kutschke
//

#include <string>
#include <vector>

#include "Mu2eG4Helper/inc/VolumeInfo.hh"

class G4CSGSolid;
class G4LogicalVolume;
class G4Material;
class G4VPhysicalVolume;

// G4 includes
#include "Geant4/G4Polycone.hh"
#include "Geant4/G4Colour.hh"
#include "Geant4/G4RotationMatrix.hh"
#include "Geant4/G4ThreeVector.hh"


namespace mu2e {

  class Polycone;
  class PolyconsParams;

  VolumeInfo nestPolycone ( std::string const& name,
                            PolyconsParams const & polyObj,
                            G4Material* material,
                            G4RotationMatrix const* rot,
                            G4ThreeVector const & offset,
                            VolumeInfo const & parent,
                            int copyNo,
                            bool const isVisible,
                            G4Colour const color,
                            bool const forceSolid,
                            bool const forceAuxEdgeVisible,
                            bool const placePV,
                            bool const doSurfaceCheck
                            );

  VolumeInfo nestPolycone ( std::string const& name,
                            PolyconsParams const & polyObj,
                            G4Material* material,
                            G4RotationMatrix const* rot,
                            G4ThreeVector const & offset,
                            VolumeInfo const & parent,
                            int copyNo,
                            G4Colour const color,
			    std::string const& lookupToken
                            );

}

#endif /* Mu2eG4_nestPolycone_hh */
