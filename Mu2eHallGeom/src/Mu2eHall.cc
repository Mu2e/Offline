// Implementation file for an accessor function for Mu2eHall.
// Allows users to know inner and outer wall placement for 
// reference.
// Original author: David No. Brown, Louisville.  August 2015

// This class's header file
#include "Mu2eHallGeom/inc/Mu2eHall.hh"

// C++ header files
#include <iostream>
#include <string>
#include <vector>

// other Mu2e include files
#include "GeomPrimitives/inc/ExtrudedSolid.hh"

// CLHEP include files
#include "CLHEP/Units/SystemOfUnits.h"
#include "CLHEP/Vector/TwoVector.h"
#include "CLHEP/Vector/ThreeVector.h"


namespace mu2e {

  double Mu2eHall::getWallExtentz( const std::string& bldgSolName, const int side) const {
    
    const ExtrudedSolid& solid = bldgSolids_.find( bldgSolName )->second;

    const std::vector<CLHEP::Hep2Vector> verts = solid.getVertices();

    double posHigh = -99999.0, pos2nd = -99999.0, negHigh = 99999.0, 
      neg2nd = 99999.0;

    for ( unsigned int i = 0; i < verts.size(); i++ ) {
      double z = verts[i].x();// x of vertices get rotated into z
      if ( z > posHigh ) {
	pos2nd = posHigh;
	posHigh = z;
      }
      if ( z > pos2nd && z < posHigh ) {
	pos2nd = z;
      }

      if ( z < negHigh ) {
	neg2nd = negHigh;
	negHigh = z;
      }
      if ( z < neg2nd && z > negHigh ) {
	neg2nd = z;
      }
    }
    if ( side == 2 ) return posHigh  + solid.getOffsetFromMu2eOrigin().z();
    if ( side == 1 ) return pos2nd   + solid.getOffsetFromMu2eOrigin().z();
    if ( side == -1 ) return neg2nd  + solid.getOffsetFromMu2eOrigin().z();
    if ( side == -2 ) return negHigh + solid.getOffsetFromMu2eOrigin().z();
    return 0.0;
  }
}
