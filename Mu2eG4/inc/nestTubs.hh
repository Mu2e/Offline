#ifndef NESTTUBS_HH
#define NESTTUBS_HH
//
// Free function to create and place a new G4Tubs, place inside a logical volume.
// 
// $Id: nestTubs.hh,v 1.5 2010/08/31 16:54:52 genser Exp $
// $Author: genser $ 
// $Date: 2010/08/31 16:54:52 $
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

#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4Colour.hh"


namespace mu2e {

  VolumeInfo nestTubs ( std::string const& name,
                        double params[5],
                        G4Material* material,
                        G4RotationMatrix* rot,
                        const G4ThreeVector& offset,
                        G4LogicalVolume* parent,
                        int copyNo,
                        G4Colour color,
                        bool forceSolid,
                        bool doSurfaceCheck
                        );
  


  // Alternate argument list, using a vector for the parameters.
  inline VolumeInfo nestTubs ( std::string const& name,
                               std::vector<double>&  params,
                               G4Material* material,
                               G4RotationMatrix* rot,
                               const G4ThreeVector& offset,
                               G4LogicalVolume* parent,
                               int copyNo,
                               G4Colour color,
                               bool forceSolid,
                               bool doSurfaceCheck
                               ){
    return nestTubs( name, 
                     &params[0],
                     material,
                     rot,
                     offset,
                     parent,
                     copyNo,
                     color,
                     forceSolid,
                     doSurfaceCheck
                     );
  }

  // Alternate argument list, using a TubsParams object for the parameters.
  inline VolumeInfo nestTubs ( std::string const& name,
                               TubsParams& params,
                               G4Material* material,
                               G4RotationMatrix* rot,
                               const G4ThreeVector& offset,
                               G4LogicalVolume* parent,
                               int copyNo,
                               G4Colour color,
                               bool forceSolid,
                               bool doSurfaceCheck
                               ){
    return nestTubs( name, 
                     &params.innerRadius,
                     material,
                     rot,
                     offset,
                     parent,
                     copyNo,
                     color,
                     forceSolid,
                     doSurfaceCheck
                     );
  }
}

#endif
