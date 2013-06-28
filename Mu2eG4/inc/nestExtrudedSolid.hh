#ifndef Mu2eG4_nestExtrudedSolid_hh
#define Mu2eG4_nestExtrudedSolid_hh
//
// Free function to create and place a new G4ExtrudedSolid inside a logical volume.
//
// $Id: nestExtrudedSolid.hh,v 1.6 2013/06/28 19:51:12 wieschie Exp $
// $Author: wieschie $
// $Date: 2013/06/28 19:51:12 $
//
// Original author Rob Kutschke
//

#include <string>
#include <vector>

#include "G4Helper/inc/VolumeInfo.hh"

class G4Material;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4CSGSolid;

// G4 includes
#include "G4ThreeVector.hh"
#include "G4TwoVector.hh"
#include "G4RotationMatrix.hh"
#include "G4Colour.hh"
#include "G4ExtrudedSolid.hh"


namespace mu2e {

  VolumeInfo nestExtrudedSolid( std::string const& name,
                                double hz,
                                std::vector<double> &x,
                                std::vector<double> &y,
                                G4Material* material,
                                G4RotationMatrix const* rot,
                                const G4ThreeVector& offset,
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
  // using VolumeInfo object for the parameters.
  VolumeInfo nestExtrudedSolid (std::string const& name,
                                double hz,
                                std::vector<double> &x,
                                std::vector<double> &y,
                                G4Material* material,
                                G4RotationMatrix const* rot,
                                const G4ThreeVector& offset,
                                VolumeInfo const & parent,
                                int copyNo,
                                bool const isVisible,
                                G4Colour const color,
                                bool const forceSolid,
                                bool const forceAuxEdgeVisible,
                                bool const placePV,
                                bool const doSurfaceCheck
                                );

  // Using VolumeInfo object and the zsections version of the solid constructor
  VolumeInfo nestExtrudedSolid (std::string const& name,
                                const std::vector<G4TwoVector>& polygon,
                                const std::vector<G4ExtrudedSolid::ZSection>& zsections,
                                G4Material* material,
                                G4RotationMatrix const* rot,
                                const G4ThreeVector& offset,
                                VolumeInfo const & parent,
                                int copyNo,
                                bool const isVisible,
                                G4Colour const color,
                                bool const forceSolid,
                                bool const forceAuxEdgeVisible,
                                bool const placePV,
                                bool const doSurfaceCheck
                                );

}

#endif /* Mu2eG4_nestExtrudedSolid_hh */
