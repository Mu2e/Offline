#ifndef Mu2eG4_nestTorus_hh
#define Mu2eG4_nestTorus_hh
//
// Free function to create and place a new G4Torus inside a logical volume.
//
//

#include <array>
#include <string>
#include <vector>

#include "G4Helper/inc/VolumeInfo.hh"
#include "GeomPrimitives/inc/TorusParams.hh"

class G4Material;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4CSGSolid;

#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4Colour.hh"


namespace mu2e {

  const unsigned long torusDim = 5;

  VolumeInfo nestTorus ( std::string const& name,
                         std::array<double,torusDim> const halfDim,
                         G4Material* material,
                         G4RotationMatrix const* rot,
                         G4ThreeVector const& offset,
                         G4LogicalVolume* parent,
                         int copyNo,
                         bool const isVisible,
                         G4Colour const color,
                         bool const forceSolid,
                         bool const forceAuxEdgeVisible,
                         bool const placePV,
                         bool const doSurfaceCheck
                         );

  // Alternate argument list (and different behavior)
  // using  VolumeInfo object
  VolumeInfo nestTorus ( std::string const& name,
                         std::array<double,torusDim> const halfDim,
                         G4Material* material,
                         G4RotationMatrix const* rot,
                         G4ThreeVector const& offset,
                         const VolumeInfo& parent,
                         int copyNo,
                         bool const isVisible,
                         G4Colour const color,
                         bool const forceSolid,
                         bool const forceAuxEdgeVisible,
                         bool const placePV,
                         bool const doSurfaceCheck
                         );

  // Alternate argument list (and different behavior)
  // using  VolumeInfo object
  VolumeInfo nestTorus ( std::string const& name,
                         std::array<double,torusDim> const halfDim,
                         G4Material* material,
                         G4RotationMatrix const* rot,
                         G4ThreeVector const& offset,
                         const VolumeInfo& parent,
                         int copyNo,
                         G4Colour const color,
			 std::string const& lookupToken = ""
                         );


  VolumeInfo nestTorus ( std::string const& name,
			 TorusParams const& halfDim,
			 G4Material* material,
			 G4RotationMatrix const* rot,
			 G4ThreeVector const& offset,
			 const VolumeInfo& parent,
			 int copyno,
                         bool const isVisible,
                         G4Colour const color,
                         bool const forceSolid,
                         bool const forceAuxEdgeVisible,
                         bool const placePV,
                         bool const doSurfaceCheck
                         );
			 
}

#endif /* Mu2eG4_nestTorus_hh */
