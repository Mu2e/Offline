#ifndef NESTBOX_HH
#define NESTBOX_HH
//
// Free function to create a new G4 Box, placed inside a logical volume.
// 
// $Id: nestBox.hh,v 1.1 2009/09/30 22:57:47 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2009/09/30 22:57:47 $
//
// Original author Rob Kutschke
//

#include <string>
#include <vector>

#include "Mu2eG4/inc/VolumeInfo.hh"

class G4Material;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4CSGSolid;

#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4Colour.hh"

namespace mu2e {

  VolumeInfo nestBox ( std::string const& name,
		       double halfDim[3],
		       G4Material* material,
		       G4RotationMatrix* rot,
		       G4ThreeVector const& offset,
		       G4LogicalVolume* parent,
		       int copyNo,
		       G4Colour color = G4Colour::Black(),
		       bool forceSolid = false
		       );
  


  // Alternate argument list, using a vector for the half dimensions.
  //
  inline VolumeInfo nestBox ( std::string const& name,
			      std::vector<double>&  halfDim,
			      G4Material* material,
			      G4RotationMatrix* rot,
			      G4ThreeVector& offset,
			      G4LogicalVolume* parent,
			      int copyNo,
			      G4Colour color = G4Colour::Black(),
			      bool forceSolid = false
			      ){
    return nestBox( name, 
		    &halfDim[0],
		    material,
		    rot,
		    offset,
		    parent,
		    copyNo,
		    color,
		    forceSolid
		    );
  }
}

#endif
