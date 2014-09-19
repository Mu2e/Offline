#ifndef Mu2eG4_nestBox_hh
#define Mu2eG4_nestBox_hh
//
// Free function to create a new G4 Box, placed inside a logical volume.
//
// $Id: nestBox.hh,v 1.11 2014/09/19 19:14:53 knoepfel Exp $
// $Author: knoepfel $
// $Date: 2014/09/19 19:14:53 $
//
// Original author Rob Kutschke
//

#include <string>
#include <vector>

#include "G4Helper/inc/VolumeInfo.hh"

class Box;

class G4Material;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4CSGSolid;

#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4Colour.hh"

namespace mu2e {

  VolumeInfo nestBox ( std::string const& name,
                       double const halfDim[3],
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

  // Alternate argument list, using a vector for the half dimensions.
  inline VolumeInfo nestBox ( std::string const& name,
                              std::vector<double> const&  halfDim,
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
    return nestBox( name,
                    &halfDim[0],
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
  // using VolumeInfo object
  VolumeInfo nestBox ( std::string const& name,
                       double const halfDim[3],
                       G4Material* material,
                       G4RotationMatrix const* rot,
                       G4ThreeVector const& offset,
                       VolumeInfo const & parent,
                       int copyNo,
                       bool const isVisible,
                       G4Colour const color,
                       bool const forceSolid,
                       bool const forceAuxEdgeVisible,
                       bool const placePV,
                       bool const doSurfaceCheck
                       );

  // Alternate argument list, (and different behavior)
  // using VolumeInfo object and using a vector for the half dimensions.
  inline VolumeInfo nestBox ( std::string const& name,
                              std::vector<double> const&  halfDim,
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
                              ){
    return nestBox( name,
                    &halfDim[0],
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
  // using VolumeInfo object
  VolumeInfo nestBox ( std::string const& name,
                       Box const & box,
                       G4Material* material,
                       G4RotationMatrix const* rot,
                       G4ThreeVector const& offset,
                       VolumeInfo const & parent,
                       int copyNo,
                       bool const isVisible,
                       G4Colour const color,
                       bool const forceSolid,
                       bool const forceAuxEdgeVisible,
                       bool const placePV,
                       bool const doSurfaceCheck
                       );

  // Alternate argument list 
  // using VolumeInfo object
  VolumeInfo nestBox ( std::string const& name,
		       std::vector<double> const&  halfDim,
		       G4Material* material,
		       G4RotationMatrix const* rot,
		       G4ThreeVector const& offset,
		       G4LogicalVolume* parent,
		       int copyNo,
		       G4Colour const color,
		       std::string const& lookupToken = ""
		       );

  // Alternate argument list 
  // using VolumeInfo object
  VolumeInfo nestBox ( std::string const& name,
		       double const halfDim[3],
		       G4Material* material,
		       G4RotationMatrix const* rot,
		       G4ThreeVector const& offset,
		       VolumeInfo const & parent,
		       int copyNo,
		       G4Colour const color,
		       std::string const& lookupToken = ""
		       );

  // Alternate argument list 
  // using VolumeInfo object
  inline VolumeInfo nestBox ( std::string const& name,
			      std::vector<double> const&  halfDim,
			      G4Material* material,
			      G4RotationMatrix const* rot,
			      G4ThreeVector const& offset,
			      VolumeInfo const & parent,
			      int copyNo,
			      G4Colour const color,
			      std::string const& lookupToken = ""
			      ) {
    return nestBox( name,
		    halfDim,
		    material,
		    rot,
		    offset,
		    parent.logical,
		    copyNo,
		    color,
		    lookupToken
		    );

  }
  

}

#endif /* Mu2eG4_nestBox_hh */
