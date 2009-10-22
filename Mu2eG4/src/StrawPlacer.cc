//
// Class to place one straw within the tracker mother volume.
//
// $Id: StrawPlacer.cc,v 1.2 2009/10/22 16:27:59 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2009/10/22 16:27:59 $
//
// Original author Rob Kutschke
//
// This does not know about segmentation of the tracker mother volume.
// Segmenation should be added to improve navigation speed.
// 

// C++ includes
#include <iostream>
#include <sstream>

// Mu2e includes
#include "Mu2eG4/inc/StrawPlacer.hh"
#include "LTrackerGeom/inc/Straw.hh"

// G4 includes
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"

using namespace std;

namespace mu2e {

  StrawPlacer::StrawPlacer( std::string basename
			    , G4LogicalVolume* logical
			    , G4LogicalVolume* motherLogical
			    ):
    _basename(basename),
    _logical(logical),
    _motherLogical(motherLogical){
  }
  
  StrawPlacer::~StrawPlacer(){}
  
  void StrawPlacer::operator() ( const Straw& s ) {

    // Copy number of this volume.
    int copyNo = s.Index().asInt();

    // Name sfor the volume.
    ostringstream os;
    os << _basename << "_" << copyNo;

    // Midpoint of the wire.
    const G4ThreeVector& mid = s.getMidPoint();

    // Unit vector in the direction of the wire.
    const G4ThreeVector& w   = s.getDirection();

    // Compute rotation of straw volume.
    // I think that this leaks the rotation matrix.
    // Also do not need separate matrices for each
    // straw, only one per wedge/vane sector.
    // Need to:
    //  1) One rotation matrix per wedge/vane.
    //  2) Delete them at end of run.

    double theta = w.theta()*radian;
    double phi   = w.phi()*radian;
    double alpha = M_PI/2.-phi;

    // This choice of Euler angles ensures that the z axis fixed to the straw
    // object points in the direction of w, as measured in the world frame
    G4RotationMatrix* rot = new G4RotationMatrix( -alpha, -theta, alpha );

    // Place the volume.
    G4VPhysicalVolume* physical =  new G4PVPlacement( rot 
						      ,mid
						      ,_logical
						      ,os.str()
						      ,_motherLogical
						      ,0
						      ,copyNo
						      );

  }

} // namespace mu2e
