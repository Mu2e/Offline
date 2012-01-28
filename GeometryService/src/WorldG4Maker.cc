#include "GeometryService/inc/WorldG4Maker.hh"

#include <iostream>
#include <algorithm>
#include <iterator>
#include <cmath>

#include "cetlib/exception.h"

#include "CLHEP/Units/SystemOfUnits.h"

#include "Mu2eUtilities/inc/SimpleConfig.hh"

#include "GeometryService/inc/WorldG4.hh"
#include "GeometryService/inc/DetectorSystem.hh"

#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/Mu2eBuilding.hh"
#include "GeometryService/inc/ProtonBeamDump.hh"
#include "ExtinctionMonitorFNAL/inc/ExtMonFNAL.hh"

namespace mu2e {

  // The body of the function is modified from Mu2eWorld::defineMu2eOrigin().
  WorldG4Maker::WorldG4Maker(const SimpleConfig& c)
    : _wg4(new WorldG4())
  {
    const int diagLevel = c.getInt("world.verbosityLevel", 0);

    // The WorldG4Maker is special: it's called after all other detector objects
    // are available in geometry service, therefore it can access their data.
    GeomHandle<Mu2eBuilding> building;
    GeomHandle<ProtonBeamDump> dump;
    GeomHandle<ExtMonFNAL::ExtMon> emf;

    const double worldBottomInMu2e = building->hallInsideYmin()
      - building->hallFloorThickness()
      - c.getDouble("world.margin.bottom");

    const double worldTopInMu2e = building->hallInsideYmax()
      + building->hallCeilingThickness()
      + building->dirtOverburdenDepth()
      + 2*building->dirtCapHalfHeight()
      + c.getDouble("world.margin.top");

    const double worldXminInMu2e = building->hallInsideXmin()
      - building->hallWallThickness()
      - c.getDouble("world.margin.xmin");

    const double worldXmaxInMu2e = building->hallInsideXmax()
      + building->hallWallThickness()
      + c.getDouble("world.margin.xmax");

    _wg4->_hallFormalZminInMu2e =
      std::min(building->hallInsideZExtMonUCIWall(),
               emf->roomCenterInMu2e()[2]
               - emf->roomHalfSize()[2]*std::abs(cos(dump->coreRotY()))
               - emf->roomHalfSize()[0]*std::abs(sin(dump->coreRotY()))
               );

    const double worldZminInMu2e = _wg4->_hallFormalZminInMu2e
      - c.getDouble("world.margin.zmin");

    const double worldZmaxInMu2e = building->hallInsideZmax()
      + building->hallWallThickness()
      + c.getDouble("world.margin.zmax");

    // Dimensions of the world box
    _wg4->_halfLengths = std::vector<double>(3);
    _wg4->_halfLengths[0] = (worldXmaxInMu2e - worldXminInMu2e)/2;
    _wg4->_halfLengths[1] = (worldTopInMu2e - worldBottomInMu2e)/2;
    _wg4->_halfLengths[2] = (worldZmaxInMu2e - worldZminInMu2e)/2;

    // Position of the origin of the mu2e coordinate system
    _wg4->_mu2eOriginInWorld = CLHEP::Hep3Vector(
                                                 -(worldXminInMu2e + worldXmaxInMu2e)/2,
                                                 -(worldBottomInMu2e + worldTopInMu2e)/2,
                                                 -(worldZminInMu2e + worldZmaxInMu2e)/2
                                                 );

    if ( diagLevel > 0) {
      std::cout << __func__ << " world halfLengths = (";
      std::copy(_wg4->_halfLengths.begin(), _wg4->_halfLengths.end(), std::ostream_iterator<double>(std::cout, ", "));
      std::cout << ")"<<std::endl;

      std::cout << __func__ << " mu2eOrigin : " <<  _wg4->_mu2eOriginInWorld  << std::endl;
    }

    // Top of the ceiling in G4 world coordinates.
    double yCeilingOutside  = _wg4->_mu2eOriginInWorld[1]
      + building->hallInsideYmax()
      + building->hallCeilingThickness();

    // Surface of the earth in G4 world coordinates.
    double ySurface  = yCeilingOutside + building->dirtOverburdenDepth();

    // Top of the dirt cap.
    double yEverest = ySurface + 2.*building->dirtCapHalfHeight();

    // Build the center for the cosmic ray production (in world coordinates)
    // at the xz origin of the detector system and on top of the dirt body (incl. the dirt cap).
    CLHEP::Hep3Vector const& detectorSystemOriginInWorld = GeomHandle<DetectorSystem>()->getOrigin() + _wg4->_mu2eOriginInWorld;
    _wg4->_cosmicReferencePoint = CLHEP::Hep3Vector(detectorSystemOriginInWorld.x(), yEverest, detectorSystemOriginInWorld.z());
    _wg4->_dirtG4Ymax = ySurface;
    _wg4->_dirtG4Ymin = yCeilingOutside;

    if ( diagLevel > 0) {
      std::cout << __func__ << " cosmic reference point : " <<  _wg4->_cosmicReferencePoint  << std::endl;
    }

    // Selfconsistency check.
    CLHEP::Hep3Vector const& cosmicReferenceInMu2e = _wg4->_cosmicReferencePoint - _wg4->_mu2eOriginInWorld;
    _wg4->inWorldOrThrow(cosmicReferenceInMu2e);
  }
}
