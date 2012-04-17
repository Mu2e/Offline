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
#include "Mu2eBuildingGeom/inc/Mu2eBuilding.hh"
#include "ProtonBeamDumpGeom/inc/ProtonBeamDump.hh"
#include "ExtinctionMonitorFNAL/inc/ExtMonFNAL.hh"

namespace mu2e {

  std::auto_ptr<WorldG4> WorldG4Maker::make(const SimpleConfig& c) {

    std::auto_ptr<WorldG4> res(new WorldG4());

    const int diagLevel = c.getInt("world.verbosityLevel", 0);

    // The WorldG4Maker is special among geometry makers:
    // it is guaranteed it will be called after all other detector objects
    // are available in geometry service, therefore it can access their data.

    GeomHandle<Mu2eBuilding> building;
    GeomHandle<ProtonBeamDump> dump;
    GeomHandle<ExtMonFNAL::ExtMon> emf;

    const double hallFormalZminInMu2e =
      std::min(
               building->hallInsideZExtMonUCIWall() - building->hallWallThickness(),

               // NB: no walls for this room at the moment
               emf->roomCenterInMu2e()[2]
               - emf->roomHalfSize()[2]*std::abs(cos(dump->coreRotY()))
               - emf->roomHalfSize()[0]*std::abs(sin(dump->coreRotY()))
               );

    const double hallFormalZmaxInMu2e = building->hallInsideZmax() + building->hallWallThickness();

    const double hallFormalXminInMu2e = building->hallInsideXmin() - building->hallWallThickness();

    const double hallFormalXmaxInMu2e = building->hallInsideXmax() + building->hallWallThickness();

    const double hallFormalYminInMu2e = building->hallInsideYmin() - building->hallFloorThickness();

    const double dirtFormalYmax  = std::max(
                                            emf->roomCenterInMu2e()[1] + emf->roomHalfSize()[1]
                                            ,
                                            dump->enclosureCenterInMu2e()[1] + dump->enclosureHalfSize()[1]
                                            );

    const double hallFormalYmaxInMu2e =
      std::max(
               dirtFormalYmax
               ,
               building->hallInsideYmax() + building->hallCeilingThickness()
               // + building->dirtOverburdenDepth() + 2*building->dirtCapHalfHeight()
               );

    res->_dumpDirtFormalYminInMu2e = hallFormalYminInMu2e;
    res->_dumpDirtFormalYmaxInMu2e = dirtFormalYmax;

    //----------------------------------------------------------------
    const double worldBottomInMu2e = hallFormalYminInMu2e - c.getDouble("world.margin.bottom");
    const double worldTopInMu2e    = hallFormalYmaxInMu2e + c.getDouble("world.margin.top");
    const double worldXminInMu2e   = hallFormalXminInMu2e - c.getDouble("world.margin.xmin");
    const double worldXmaxInMu2e   = hallFormalXmaxInMu2e + c.getDouble("world.margin.xmax");
    const double worldZminInMu2e   = hallFormalZminInMu2e - c.getDouble("world.margin.zmin");
    const double worldZmaxInMu2e   = hallFormalZmaxInMu2e + c.getDouble("world.margin.zmax");

    // Dimensions of the world box
    res->_halfLengths = std::vector<double>(3);
    res->_halfLengths[0] = (worldXmaxInMu2e - worldXminInMu2e)/2;
    res->_halfLengths[1] = (worldTopInMu2e - worldBottomInMu2e)/2;
    res->_halfLengths[2] = (worldZmaxInMu2e - worldZminInMu2e)/2;

    // Position of the origin of the mu2e coordinate system
    res->_mu2eOriginInWorld = CLHEP::Hep3Vector(
                                                 -(worldXminInMu2e + worldXmaxInMu2e)/2,
                                                 -(worldBottomInMu2e + worldTopInMu2e)/2,
                                                 -(worldZminInMu2e + worldZmaxInMu2e)/2
                                                 );

    res->_hallFormalHalfSize.resize(3);
    res->_hallFormalHalfSize[0] = (hallFormalXmaxInMu2e - hallFormalXminInMu2e)/2;
    res->_hallFormalHalfSize[1] = (hallFormalYmaxInMu2e - hallFormalYminInMu2e)/2;
    res->_hallFormalHalfSize[2] = (hallFormalZmaxInMu2e - hallFormalZminInMu2e)/2;

    const CLHEP::Hep3Vector hallFormalCenterInMu2e(
                                                   (hallFormalXmaxInMu2e + hallFormalXminInMu2e)/2,
                                                   (hallFormalYmaxInMu2e + hallFormalYminInMu2e)/2,
                                                   (hallFormalZmaxInMu2e + hallFormalZminInMu2e)/2
                                                   );

    res->_hallFormalCenterInWorld = hallFormalCenterInMu2e + res->_mu2eOriginInWorld;

    // By construction everything fits inside formalHallBox.
    // Therefore it is save to start cosmic rays at the Word coordinate
    double yEverest = res->hallFormalCenterInWorld()[1] + res->hallFormalHalfSize()[1];

    // Build the center for the cosmic ray production (in world coordinates)
    // at the xz origin of the detector system and on top of the dirt body (incl. the dirt cap).
    CLHEP::Hep3Vector const& detectorSystemOriginInWorld = GeomHandle<DetectorSystem>()->getOrigin() + res->_mu2eOriginInWorld;
    res->_cosmicReferencePoint = CLHEP::Hep3Vector(detectorSystemOriginInWorld.x(), yEverest, detectorSystemOriginInWorld.z());
    //res->_dirtG4Ymax = ySurface;
    //res->_dirtG4Ymin = yCeilingOutside;

    if ( diagLevel > 0) {
      std::cout<<*res<<std::endl;
    }

    // Selfconsistency check.
    CLHEP::Hep3Vector const& cosmicReferenceInMu2e = res->_cosmicReferencePoint - res->_mu2eOriginInWorld;
    res->inWorldOrThrow(cosmicReferenceInMu2e);

    return res;

  } // make()

} // namespace mu2e
