#ifndef Mu2eG4_nestPolycone_hh
#define Mu2eG4_nestPolycone_hh
//
// Free function to create and place a new G4Polycone, place inside a logical volume.
//
// $Id: nestPolycone.hh,v 1.2 2013/12/02 20:12:16 genser Exp $
// $Author: genser $
// $Date: 2013/12/02 20:12:16 $
//
// Original author Rob Kutschke
//

#include <string>
#include <vector>

#include "G4Helper/inc/VolumeInfo.hh"

class G4CSGSolid;
class G4LogicalVolume;
class G4Material;
class G4VPhysicalVolume;

// G4 includes
#include "G4Polycone.hh"
#include "G4Colour.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"


namespace mu2e {

  class Polycone;
  class PolyconsParams;

  VolumeInfo nestPolycone ( std::string const& name,
                            PolyconsParams const & polyObj,
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

}

#endif /* Mu2eG4_nestPolycone_hh */
