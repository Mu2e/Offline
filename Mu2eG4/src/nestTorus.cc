//
// Free function to create and place a new G4Torus, place inside a logical volume.
//
// $Id: nestTorus.cc,v 1.5 2011/05/18 02:27:18 wb Exp $
// $Author: wb $
// $Date: 2010/03/15
//

#include <string>

#include "Mu2eG4/inc/nestTorus.hh"
#include "Mu2eG4/inc/finishNesting.hh"

#include "G4Torus.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4ThreeVector.hh"

using namespace std;

namespace mu2e {

  //
  // Create and place a G4Torus inside a logical volume.
  //
  VolumeInfo nestTorus ( string const& name,
                         double const halfDim[5],
                         G4Material* material,
                         G4RotationMatrix* rot,
                         G4ThreeVector const& offset,
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

    info.name     = name;

    info.solid   = new G4Torus( name, halfDim[0], halfDim[1], halfDim[2], halfDim[3], halfDim[4]  );

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


  VolumeInfo nestTorus ( string const & name,
                         double const halfDim[5],
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


    VolumeInfo info(name,offset,parent.centerInWorld);

    info.solid    = new G4Torus( name, halfDim[0], halfDim[1], halfDim[2], halfDim[3], halfDim[4] );

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

}
