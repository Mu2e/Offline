//
// Free function to create and place a new G4Cons inside a logical volume.
//
// $Id: nestCons.cc,v 1.4 2014/09/19 19:15:21 knoepfel Exp $
// $Author: knoepfel $
// $Date: 2014/09/19 19:15:21 $
//
// Original author Rob Kutschke
//

#include <string>

// Mu2e includes
#include "Mu2eG4/inc/nestCons.hh"
#include "Mu2eG4/inc/finishNesting.hh"

// G4 includes
#include "G4Cons.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4ThreeVector.hh"

using namespace std;

namespace mu2e {

  //
  // Create and place a G4Cons inside a logical volume.
  //
  VolumeInfo nestCons ( string const & name,
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
                        ){


    VolumeInfo info;

    info.name    = name;

    info.solid   = new G4Cons( name, params[0], params[1], params[2], params[3],
                               params[4], params[5], params[6]  );

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

  VolumeInfo nestCons ( string const & name,
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
                        ){


    VolumeInfo info(name,offset,parent.centerInWorld);

    info.solid   = new G4Cons( name, params[0], params[1], params[2], params[3],
                               params[4], params[5], params[6]  );

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

  VolumeInfo nestCons ( string const & name,
                        double const params[consDim],
                        G4Material* material,
                        G4RotationMatrix const* rot,
                        G4ThreeVector const & offset,
                        VolumeInfo const & parent,
                        int copyNo,
                        G4Colour const color,
			string const & lookupToken
                        ){


    VolumeInfo info(name,offset,parent.centerInWorld);

    info.solid   = new G4Cons( name, params[0], params[1], params[2], params[3],
                               params[4], params[5], params[6]  );

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
