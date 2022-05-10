//
// Free function to create a new G4 Box, placed inside a logical volume.
//
//
// Original author Rob Kutschke
//

#include <string>

// Mu2e includes
#include "Offline/GeomPrimitives/inc/Box.hh"
#include "Offline/Mu2eG4/inc/nestBox.hh"
#include "Offline/Mu2eG4/inc/finishNesting.hh"

// G4 includes
#include "Geant4/G4Box.hh"
#include "Geant4/G4LogicalVolume.hh"
#include "Geant4/G4Material.hh"
#include "Geant4/G4ThreeVector.hh"
#include "Geant4/G4VisAttributes.hh"

using namespace std;

namespace mu2e {

  //
  // Place a box inside a logical volume.
  //
  VolumeInfo nestBox ( string const& name,
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
                       ){

    VolumeInfo info;

    info.name    = name;

    info.solid   = new G4Box( name, halfDim[0], halfDim[1], halfDim[2] );

    finishNesting(info,
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

    return info;
  }

  VolumeInfo nestBox ( string const& name,
                       double const halfDim[3],
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
                       ){

    VolumeInfo info(name,offset,parent.centerInWorld);

    info.solid   = new G4Box( name, halfDim[0], halfDim[1], halfDim[2] );

    finishNesting(info,
                  material,
                  rot,
                  offset,
                  parent.logical,
                  copyNo,
                  isVisible,
                  color,
                  forceSolid,
                  forceAuxEdgeVisible,
                  placePV,
                  doSurfaceCheck
                  );

    return info;

  }

  VolumeInfo nestBox ( string const& name,
                       Box const& box,
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
                       ){

    VolumeInfo info(name,offset,parent.centerInWorld);

    info.solid   = new G4Box( name, box.getXhalfLength(), box.getYhalfLength(), box.getZhalfLength() );

    finishNesting(info,
                  material,
                  rot,
                  offset,
                  parent.logical,
                  copyNo,
                  isVisible,
                  color,
                  forceSolid,
                  forceAuxEdgeVisible,
                  placePV,
                  doSurfaceCheck
                  );

    return info;

  }

  VolumeInfo nestBox ( string const& name,
                       std::vector<double> const& halfDim,
                       G4Material* material,
                       G4RotationMatrix const* rot,
                       G4ThreeVector const& offset,
                       G4LogicalVolume* parent,
                       int copyNo,
                       G4Colour const color,
                       string const& lookupToken
                       ){

    VolumeInfo info;

    info.name    = name;

    info.solid   = new G4Box( name, halfDim[0], halfDim[1], halfDim[2] );

    finishNesting(info,
                  material,
                  rot,
                  offset,
                  parent,
                  copyNo,
                  color,
                  lookupToken
                  );

    return info;

  }


  VolumeInfo nestBox ( string const& name,
                       double const halfDim[3],
                       G4Material* material,
                       G4RotationMatrix const* rot,
                       G4ThreeVector const& offset,
                       const VolumeInfo& parent,
                       int copyNo,
                       G4Colour const color,
                       string const& lookupToken
                       ){

    VolumeInfo info(name,offset,parent.centerInWorld);

    info.solid   = new G4Box( name, halfDim[0], halfDim[1], halfDim[2] );

    finishNesting(info,
                  material,
                  rot,
                  offset,
                  parent.logical,
                  copyNo,
                  color,
                  lookupToken
                  );

    return info;

  }

}
