#include "GeometryService/inc/WorldG4Maker.hh"

#include <iostream>

#include "CLHEP/Units/SystemOfUnits.h"

#include "ConfigTools/inc/SimpleConfig.hh"

#include "Mu2eHallGeom/inc/Mu2eHall.hh"
#include "GeometryService/inc/WorldG4.hh"

#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/Mu2eEnvelope.hh"

namespace mu2e {

  //----------------------------------------------------------------
  std::unique_ptr<WorldG4> WorldG4Maker::make(const Mu2eHall& hall,
                                              const SimpleConfig& c) {

    std::unique_ptr<WorldG4> res(new WorldG4());

    const int diagLevel = c.getInt("world.verbosityLevel", 0);

    // This is where we want mu2eOriginInWorld
    const CLHEP::Hep3Vector requestedWorldCenterInMu2e( -c.getHep3Vector("world.mu2eOriginInWorld"));
    
    // The WorldG4Maker is special among geometry makers:
    // it is guaranteed it will be called after all other detector objects
    // are available in geometry service, therefore it can access their data.

    GeomHandle<Mu2eEnvelope> env;

    const double prelimWorldBottomInMu2e = env->ymin() - c.getDouble("world.minimalMargin.bottom");
    const double prelimWorldTopInMu2e    = env->ymax() + c.getDouble("world.minimalMargin.top");
    const double prelimWorldXminInMu2e   = env->xmin() - c.getDouble("world.minimalMargin.xmin");
    const double prelimWorldXmaxInMu2e   = env->xmax() + c.getDouble("world.minimalMargin.xmax");
    const double prelimWorldZminInMu2e   = env->zmin() - c.getDouble("world.minimalMargin.zmin");
    const double prelimWorldZmaxInMu2e   = env->zmax() + c.getDouble("world.minimalMargin.zmax");

    // Using the above we'd get the following position for the origin of the mu2e coordinate system
    CLHEP::Hep3Vector prelimWorldCenterInMu2e(
                                              (prelimWorldXminInMu2e + prelimWorldXmaxInMu2e)/2,
                                              (prelimWorldBottomInMu2e + prelimWorldTopInMu2e)/2,
                                              (prelimWorldZminInMu2e + prelimWorldZmaxInMu2e)/2
                                              );

    // the requestedWorldCenterInMu2e shifts we need to correct for
    const double dx(requestedWorldCenterInMu2e[0] - prelimWorldCenterInMu2e[0]);
    const double dy(requestedWorldCenterInMu2e[1] - prelimWorldCenterInMu2e[1]);
    const double dz(requestedWorldCenterInMu2e[2] - prelimWorldCenterInMu2e[2]);

    // Compute the final world boundaries
    const double worldBottomInMu2e = prelimWorldBottomInMu2e + (dy < 0 ? 2*dy :    0);
    const double worldTopInMu2e    = prelimWorldTopInMu2e    + (dy < 0 ?    0 : 2*dy);
    const double worldXminInMu2e   = prelimWorldXminInMu2e   + (dx < 0 ? 2*dx :    0);
    const double worldXmaxInMu2e   = prelimWorldXmaxInMu2e   + (dx < 0 ?    0 : 2*dx);
    const double worldZminInMu2e   = prelimWorldZminInMu2e   + (dz < 0 ? 2*dz :    0);
    const double worldZmaxInMu2e   = prelimWorldZmaxInMu2e   + (dz < 0 ?    0 : 2*dz);

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

    //
    res->_hallFormalHalfSize.resize(3);
    res->_hallFormalHalfSize[0] = (env->xmax() - env->xmin())/2;
    res->_hallFormalHalfSize[1] = (env->ymax() - env->ymin())/2;
    res->_hallFormalHalfSize[2] = (env->zmax() - env->zmin())/2;

    const CLHEP::Hep3Vector hallFormalCenterInMu2e(
                                                   (env->xmax() + env->xmin())/2,
                                                   (env->ymax() + env->ymin())/2,
                                                   (env->zmax() + env->zmin())/2
                                                   );

    res->_hallFormalCenterInWorld = hallFormalCenterInMu2e + res->_mu2eOriginInWorld;

    const auto& vol  = hall.getDirtSolid("dirtDsAreaFirstFloorS");
    res->_dirtG4Ymax = vol.getOffsetFromMu2eOrigin().y()+vol.getYhalfThickness();

    if ( diagLevel > 0) {
      std::cout<<*env<<std::endl;
      std::cout<<*res<<std::endl;
    }

    return res;

  } // make()

} // namespace mu2e
