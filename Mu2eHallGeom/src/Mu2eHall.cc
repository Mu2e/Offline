// Implementation file for an accessor function for Mu2eHall.
// Allows users to know inner and outer wall placement for
// reference.
// Original author: David No. Brown, Louisville.  August 2015

// This class's header file
#include "Offline/Mu2eHallGeom/inc/Mu2eHall.hh"

// C++ header files
#include <iostream>
#include <string>
#include <vector>

// other Mu2e include files
#include "Offline/GeomPrimitives/inc/ExtrudedSolid.hh"

// CLHEP include files
#include "CLHEP/Units/SystemOfUnits.h"
#include "CLHEP/Vector/TwoVector.h"
#include "CLHEP/Vector/ThreeVector.h"

#include "cetlib_except/exception.h"

namespace mu2e {
  const ExtrudedSolid& Mu2eHall::getBldgSolid(const std::string& name) const {
    auto it = bldgSolids_.find(name);
    if(it == bldgSolids_.end()) {
      throw cet::exception("GEOM")
        <<"Mu2eHall::getBldgSolid(): unknown volume \""<< name <<"\"\n";
    }
    return it->second;
  }

  const ExtrudedSolid& Mu2eHall::getDirtSolid(const std::string& name) const {
    auto it = dirtSolids_.find(name);
    if(it == dirtSolids_.end()) {
      throw cet::exception("GEOM")
        <<"Mu2eHall::getDirtSolid(): unknown volume \""<< name <<"\"\n";
    }
    return it->second;
  }

  const RotExtrudedSolid& Mu2eHall::getRotSolid(const std::string& name) const {
    auto it = rotatedSolids_.find(name);
    if(it == rotatedSolids_.end()) {
      throw cet::exception("GEOM")
        <<"Mu2eHall::getRotSolid(): unknown volume \""<< name <<"\"\n";
    }
    return it->second;
  }

  const GenericTrap& Mu2eHall::getDirtTrapSolid(const std::string& name) const {
    auto it = dirtTrapSolids_.find(name);
    if(it == dirtTrapSolids_.end()) {
      throw cet::exception("GEOM")
        <<"Mu2eHall::getDirtTrapSolid(): unknown volume \""<< name <<"\"\n";
    }
    return it->second;
  }

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
