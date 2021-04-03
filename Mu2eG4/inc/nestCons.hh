#ifndef Mu2eG4_nestCons_hh
#define Mu2eG4_nestCons_hh
//
// Free function to create and place a new G4Cons inside a logical volume.
//
//
// Original author Rob Kutschke
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
#include "Geant4/G4RotationMatrix.hh"
#include "Geant4/G4Colour.hh"


namespace mu2e {

  const unsigned long consDim = 7;

  VolumeInfo nestCons ( std::string const& name,
                        double const params[consDim],
                        G4Material* material,
                        G4RotationMatrix const* rot,
                        G4ThreeVector const & offset,
                        G4LogicalVolume* parent,
                        int copyNo,
                        bool const isVisible,
                        G4Colour const color,
                        bool const forceSolid,
                        bool const forceAuxEdgeVisible,
                        bool const placePV,
                        bool const doSurfaceCheck
                        );



  // Alternate argument list, using a vector for the parameters.
  inline VolumeInfo nestCons ( std::string const& name,
                               std::vector<double>&  params,
                               G4Material* material,
                               G4RotationMatrix const* rot,
                               G4ThreeVector const & offset,
                               G4LogicalVolume* parent,
                               int copyNo,
                               bool const isVisible,
                               G4Colour const color,
                               bool const forceSolid,
                               bool const forceAuxEdgeVisible,
                               bool const placePV,
                               bool const doSurfaceCheck
                               ){
    return nestCons( name,
                     &params[0],
                     material,
                     rot,
                     offset,
                     parent,
                     copyNo,
                     isVisible,
                     color,
                     forceSolid,
                     forceAuxEdgeVisible,
                     placePV,
                     doSurfaceCheck
                     );
  }

  // Alternate argument list (and different behavior)
  // using VolumeInfo object for the parameters.
  VolumeInfo nestCons ( std::string const& name,
                        double const params[consDim],
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

 
  // Alternate argument list (and different behavior)
  // using VolumeInfo object for the parameters.
  VolumeInfo nestCons ( std::string const& name,
                        double const params[consDim],
                        G4Material* material,
                        G4RotationMatrix const* rot,
                        G4ThreeVector const & offset,
                        VolumeInfo const & parent,
                        int copyNo,
                        G4Colour const color,
			std::string const& lookupToken
                        ); 

  inline VolumeInfo nestCons ( std::string const& name,
                        std::array<double,consDim> const& params,
                        G4Material* material,
                        G4RotationMatrix const* rot,
                        G4ThreeVector const & offset,
                        VolumeInfo const & parent,
                        int copyNo,
                        G4Colour const color,
			std::string const& lookupToken
                        ) {
    return nestCons ( name,
		      &params[0],
		      material,
		      rot,
		      offset,
		      parent,
		      copyNo,
		      color,
		      lookupToken
		      ); 
  }

  inline VolumeInfo nestCons ( std::string const& name,
			       std::vector<double> const& params,
			       G4Material* material,
			       G4RotationMatrix const* rot,
			       G4ThreeVector const & offset,
			       VolumeInfo const & parent,
			       int copyNo,
			       G4Colour const color,
			       std::string const& lookupToken
			       ) {
    return nestCons ( name,
		      &params[0],
		      material,
		      rot,
		      offset,
		      parent,
		      copyNo,
		      color,
		      lookupToken
		      ); 
  }

}

#endif /* Mu2eG4_nestCons_hh */
