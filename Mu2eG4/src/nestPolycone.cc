//
// Free function to create and place a new G4Tubs, place inside a logical volume.
//
// $Id: nestPolycone.cc,v 1.2 2014/09/19 19:15:23 knoepfel Exp $
// $Author: knoepfel $
// $Date: 2014/09/19 19:15:23 $
//
// Original author Rob Kutschke
//

#include <string>

// Mu2e includes
#include "GeomPrimitives/inc/PolyconsParams.hh"
#include "Mu2eG4/inc/findMaterialOrThrow.hh"
#include "Mu2eG4/inc/nestPolycone.hh"
#include "Mu2eG4/inc/finishNesting.hh"

// G4 includes
#include "G4Polycone.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4ThreeVector.hh"

using namespace std;

namespace mu2e {

  //
  // Create and place a G4Polycone inside a logical volume.
  //
  VolumeInfo nestPolycone ( string const & name,
			    PolyconsParams const & polyParams,
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
    
    info.solid    = new G4Polycone( name, 
				    polyParams.phi0(),
				    polyParams.phiTotal(),
				    polyParams.numZPlanes(),
				    &polyParams.zPlanes()[0],
				    &polyParams.rInner()[0],
				    &polyParams.rOuter()[0] );

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


  //
  // Create and place a G4Polycone inside a logical volume.
  //
  VolumeInfo nestPolycone ( string const & name,
			    PolyconsParams const & polyParams,
			    G4Material* material,
			    G4RotationMatrix const* rot,
			    G4ThreeVector const & offset,
			    VolumeInfo const & parent,
			    int copyNo,
			    G4Colour const color,
			    string const & lookupToken
			    ){
    
    VolumeInfo info(name,offset,parent.centerInWorld);
    
    info.solid    = new G4Polycone( name, 
				    polyParams.phi0(),
				    polyParams.phiTotal(),
				    polyParams.numZPlanes(),
				    &polyParams.zPlanes()[0],
				    &polyParams.rInner()[0],
				    &polyParams.rOuter()[0] );

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
