#ifndef Mu2eG4_nestPolyhedra_hh
#define Mu2eG4_nestPolyhedra_hh
//
// Free function to create and place a new G4Polyhedra, place inside a logical volume.
//
//
// Original author Rob Kutschke
//

#include <string>
#include <vector>

#include "G4Helper/inc/VolumeInfo.hh"

class G4CSGSolid;
class G4LogicalVolume;
class G4Material;
class G4VPhysicalVolume;

// G4 includes
#include "G4Polyhedra.hh"
#include "G4Colour.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"


namespace mu2e {

  class Polyhedra;
  class PolyhedraParams;

  VolumeInfo nestPolyhedra ( std::string const& name,
                            PolyhedraParams const & polyObj,
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

  VolumeInfo nestPolyhedra ( std::string const& name,
                            PolyhedraParams const & polyObj,
                            G4Material* material,
                            G4RotationMatrix const* rot,
                            G4ThreeVector const & offset,
                            VolumeInfo const & parent,
                            int copyNo,
                            G4Colour const color,
			    std::string const& lookupToken
                            );

}

#endif /* Mu2eG4_nestPolyhedra_hh */
