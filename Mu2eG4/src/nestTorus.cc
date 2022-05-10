//
// Free function to create and place a new G4Torus, place inside a logical volume.
//
//

#include <string>

#include "Offline/Mu2eG4/inc/nestTorus.hh"
#include "Offline/Mu2eG4/inc/finishNesting.hh"
#include "Offline/GeomPrimitives/inc/TorusParams.hh"

#include "Geant4/G4Torus.hh"
#include "Geant4/G4LogicalVolume.hh"
#include "Geant4/G4Material.hh"
#include "Geant4/G4ThreeVector.hh"

using namespace std;

namespace mu2e {

  //
  // Create and place a G4Torus inside a logical volume.
  //
  VolumeInfo nestTorus ( string const& name,
                         array<double,torusDim> const halfDim,
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
                         array<double,torusDim> const halfDim,
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

  VolumeInfo nestTorus ( string const & name,
                         array<double,torusDim> const halfDim,
                         G4Material* material,
                         G4RotationMatrix const* rot,
                         G4ThreeVector const & offset,
                         VolumeInfo const & parent,
                         int copyNo,
                         G4Colour const color,
                         string const & lookupToken
                         ){


    VolumeInfo info(name,offset,parent.centerInWorld);

    info.solid    = new G4Torus( name, halfDim[0], halfDim[1], halfDim[2], halfDim[3], halfDim[4] );

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

  VolumeInfo nestTorus( string const& name,
                        TorusParams const& halfDim,
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
                        bool const doSurfaceCheck )
  {

    VolumeInfo info(name,offset,parent.centerInWorld);

    info.solid    = new G4Torus( name, halfDim.innerRadius(),
                                 halfDim.outerRadius(), halfDim.torusRadius(),
                                 halfDim.phi0(), halfDim.phiMax() );

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
