#ifndef Mu2eG4_nestTubs_hh
#define Mu2eG4_nestTubs_hh
//
// Free function to create and place a new G4Tubs, place inside a logical volume.
//
// $Id: nestTubs.hh,v 1.11 2011/05/20 19:18:44 wb Exp $
// $Author: wb $
// $Date: 2011/05/20 19:18:44 $
//
// Original author Rob Kutschke
//

#include <string>
#include <vector>

#include "G4Helper/inc/VolumeInfo.hh"
#include "TrackerGeom/inc/TubsParams.hh"

class G4CSGSolid;
class G4LogicalVolume;
class G4Material;
class G4VPhysicalVolume;

// G4 includes
#include "G4Colour.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"


namespace mu2e {

  VolumeInfo nestTubs ( std::string const& name,
                        double const params[5],
                        G4Material* material,
                        G4RotationMatrix* rot,
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
  inline VolumeInfo nestTubs ( std::string const& name,
                               std::vector<double>&  params,
                               G4Material* material,
                               G4RotationMatrix* rot,
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
    return nestTubs( name,
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

  // Alternate argument list, using a TubsParams object for the parameters.
  inline VolumeInfo nestTubs ( std::string const& name,
                               TubsParams const & params,
                               G4Material* material,
                               G4RotationMatrix* rot,
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
    return nestTubs( name,
                     params.data(),
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
  // using a TubsParams & VolumeInfo object for the parameters.
  VolumeInfo nestTubs ( std::string const& name,
                        double const params[5],
                        G4Material* material,
                        G4RotationMatrix* rot,
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
  // using a TubsParams & VolumeInfo object for the parameters.
  inline VolumeInfo nestTubs ( std::string const& name,
                               TubsParams const & params,
                               G4Material* material,
                               G4RotationMatrix* rot,
                               G4ThreeVector const & offset,
                               VolumeInfo const & parent,
                               int copyNo,
                               bool const isVisible,
                               G4Colour const color,
                               bool const forceSolid,
                               bool const forceAuxEdgeVisible,
                               bool const placePV,
                               bool const doSurfaceCheck
                               ){
    return nestTubs( name,
                     params.data(),
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


}

#endif /* Mu2eG4_nestTubs_hh */
