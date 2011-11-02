#include "GeometryService/inc/WorldG4Maker.hh"

#include <iostream>

#include "cetlib/exception.h"

#include "CLHEP/Units/SystemOfUnits.h"

#include "Mu2eUtilities/inc/SimpleConfig.hh"

#include "GeometryService/inc/WorldG4.hh"

namespace mu2e {

  // The body of the function is modified from Mu2eWorld::defineMu2eOrigin().
  WorldG4Maker::WorldG4Maker(const SimpleConfig& c) 
    : _wg4(new WorldG4())
  {
    const int diagLevel = c.getInt("world.verbosityLevel", 0);

    c.getVectorDouble("world.halfLengths", _wg4->_halfLengths, 3);

    // Floor thickness.
    _wg4->_hallFloorThickness = c.getDouble("hall.floorThick");

    // Top of the floor in G4 world coordinates.
    double yFloor = -_wg4->_halfLengths[1] + _wg4->_hallFloorThickness;

    if ( diagLevel > 0) {
      std::cout << __func__ << " yFloor : " <<  yFloor  << std::endl;
    }

    // The height above the floor of the y origin of the Mu2e coordinate system.
    double yOriginHeight = c.getDouble("world.mu2eOrigin.height" )*CLHEP::mm;

    // Position of the origin of the mu2e coordinate system
    _wg4->_mu2eOriginInWorld = CLHEP::Hep3Vector(
						 c.getDouble("world.mu2eOrigin.xoffset")*CLHEP::mm,
						 yFloor + yOriginHeight,
						 c.getDouble("world.mu2eOrigin.zoffset")*CLHEP::mm
						 );

    if ( diagLevel > 0) {
      std::cout << __func__ << " mu2eOrigin : " <<  _wg4->_mu2eOriginInWorld  << std::endl;
    }

    // Origin used to construct the MECO detector.
    // Magic number to fix:
    _wg4->_trackerOrigin = _wg4->_mu2eOriginInWorld + CLHEP::Hep3Vector( -3904., 0., 12000.);
    
    if ( diagLevel > 0) {
      std::cout << __func__ << " mu2eDetectorOrigin : " <<  _wg4->_trackerOrigin  << std::endl;
    }

    _wg4->_hallCeilingThickness = c.getDouble("hall.ceilingThick");
    _wg4->_dirtOverburdenDepth  = c.getDouble("dirt.overburdenDepth");
    _wg4->_dirtCapHalfHeight    = c.getDouble("dirt.capHalfHeight");

    c.getVectorDouble("hall.insideHalfLengths",_wg4->_hallInsideHalfLenghts,3);

    // Bottom of the ceiling in G4 world coordinates.
    double yCeilingInSide = yFloor + 2.*_wg4->_hallInsideHalfLenghts[1];

    // Top of the ceiling in G4 world coordinates.
    double yCeilingOutside  = yCeilingInSide + _wg4->_hallCeilingThickness;

    // Surface of the earth in G4 world coordinates.
    double ySurface  = yCeilingOutside + _wg4->_dirtOverburdenDepth;

    // Top of the world.
    double yEverest = ySurface + 2.*_wg4->_dirtCapHalfHeight;

    // Build the reference points that others will use.
    _wg4->_cosmicReferencePoint = CLHEP::Hep3Vector(0., yEverest, 0.);

    _wg4->_dirtG4Ymax = ySurface;
    _wg4->_dirtG4Ymin = yCeilingOutside;

    if ( diagLevel > 0) {
      std::cout << __func__ << " yEverest : " <<  yEverest  << std::endl;
    }

    // Selfconsistency check.
    if ( yEverest > 2.*_wg4->_halfLengths[1] ){
      throw cet::exception("GEOM")
        << "Top of the world is outside of the world volume! \n";
    }
  }
}
