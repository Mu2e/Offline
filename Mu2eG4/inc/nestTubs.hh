#ifndef NESTTUBS_HH
#define NESTTUBS_HH
//
// Free function to create and place a new G4Tubs, place inside a logical volume.
// 
// $Id: nestTubs.hh,v 1.7 2010/12/02 17:47:05 genser Exp $
// $Author: genser $ 
// $Date: 2010/12/02 17:47:05 $
//
// Original author Rob Kutschke
//

#include <string>
#include <vector>

#include "Mu2eG4/inc/VolumeInfo.hh"
#include "TrackerGeom/inc/TubsParams.hh"

class G4Material;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4CSGSolid;

// G4 includes
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4Colour.hh"


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
                     &params.innerRadius,
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
                     &params.innerRadius,
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

#endif
