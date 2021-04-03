#ifndef Mu2eG4_nestExtrudedSolid_hh
#define Mu2eG4_nestExtrudedSolid_hh
//
// Free function to create and place a new G4ExtrudedSolid inside a
// logical volume.
//

#include <string>
#include <vector>

#include "Mu2eG4Helper/inc/VolumeInfo.hh"

class G4Material;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4CSGSolid;

// G4 includes
#include "Geant4/G4ThreeVector.hh"
#include "Geant4/G4TwoVector.hh"
#include "Geant4/G4RotationMatrix.hh"
#include "Geant4/G4Colour.hh"
#include "Geant4/G4ExtrudedSolid.hh"


namespace mu2e {

  VolumeInfo nestExtrudedSolid( std::string const& name,
                                double hz,
                                std::vector<double> const& x,
                                std::vector<double> const& y,
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
                                bool const doSurfaceCheck
                                );

  // Alternate argument list (and different behavior)
  // using VolumeInfo object for the parameters.
  VolumeInfo nestExtrudedSolid (std::string const& name,
                                double hz,
                                std::vector<double> const& x,
                                std::vector<double> const& y,
                                G4Material* material,
                                G4RotationMatrix const* rot,
                                G4ThreeVector const& offset,
                                VolumeInfo const& parent,
                                int copyNo,
                                bool const isVisible,
                                G4Colour const& color,
                                bool const forceSolid,
                                bool const forceAuxEdgeVisible,
                                bool const placePV,
                                bool const doSurfaceCheck
                                );
  // same with lookupToken
  VolumeInfo nestExtrudedSolid (std::string const& name,
                                double hz,
                                std::vector<double> const& x,
                                std::vector<double> const& y,
                                G4Material* material,
                                G4RotationMatrix const* rot,
                                G4ThreeVector const& offset,
                                VolumeInfo const & parent,
                                int copyNo,
                                G4Colour const& color,
                                std::string const& lookupToken
                                );

  // Using VolumeInfo object and the zsections version of the solid constructor
  VolumeInfo nestExtrudedSolid (std::string const& name,
                                std::vector<G4TwoVector> const& polygon,
                                std::vector<G4ExtrudedSolid::ZSection> const& zsections,
                                G4Material* material,
                                G4RotationMatrix const* rot,
                                G4ThreeVector const& offset,
                                VolumeInfo const& parent,
                                int copyNo,
                                bool const isVisible,
                                G4Colour const& color,
                                bool const forceSolid,
                                bool const forceAuxEdgeVisible,
                                bool const placePV,
                                bool const doSurfaceCheck
                                );

  // same with lookupToken
  VolumeInfo nestExtrudedSolid (std::string const& name,
                                std::vector<G4TwoVector> const& polygon,
                                std::vector<G4ExtrudedSolid::ZSection> const& zsections,
                                G4Material* material,
                                G4RotationMatrix const* rot,
                                G4ThreeVector const& offset,
                                VolumeInfo const & parent,
                                int copyNo,
                                G4Colour const& color,
                                std::string const& lookupToken
                                );

  // using VolumeInfo object for the parameters and full set of G4ExtrudedSolid params
  VolumeInfo nestExtrudedSolid (std::string const& name,
                                std::vector<G4TwoVector> const& polygon,
                                G4double hz,
                                G4TwoVector const& offset1, G4double scale1,
                                G4TwoVector const& offset2, G4double scale2,
                                G4Material* material,
                                G4RotationMatrix const* rot,
                                G4ThreeVector const& offset,
                                VolumeInfo const & parent,
                                int copyNo,
                                bool const isVisible,
                                G4Colour const& color,
                                bool const forceSolid,
                                bool const forceAuxEdgeVisible,
                                bool const placePV,
                                bool const doSurfaceCheck
                                );

  // using VolumeInfo object for the parameters and full set of G4ExtrudedSolid params
  VolumeInfo nestExtrudedSolid (std::string const& name,
                                std::vector<G4TwoVector> const& polygon,
                                G4double hz,
                                G4TwoVector const& offset1, G4double scale1,
                                G4TwoVector const& offset2, G4double scale2,
                                G4Material* material,
                                G4RotationMatrix const* rot,
                                G4ThreeVector const& offset,
                                VolumeInfo const & parent,
                                int copyNo,
                                G4Colour const& color,
                                std::string const& lookupToken
                                );

}

#endif /* Mu2eG4_nestExtrudedSolid_hh */
