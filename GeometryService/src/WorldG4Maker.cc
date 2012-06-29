#include "GeometryService/inc/WorldG4Maker.hh"

#include <iostream>

#include "CLHEP/Units/SystemOfUnits.h"

#include "Mu2eUtilities/inc/SimpleConfig.hh"

#include "GeometryService/inc/WorldG4.hh"

#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/Mu2eEnvelope.hh"

namespace mu2e {

  //----------------------------------------------------------------
  std::auto_ptr<WorldG4> WorldG4Maker::make(const SimpleConfig& c) {

    std::auto_ptr<WorldG4> res(new WorldG4());

    const int diagLevel = c.getInt("world.verbosityLevel", 0);

    // The WorldG4Maker is special among geometry makers:
    // it is guaranteed it will be called after all other detector objects
    // are available in geometry service, therefore it can access their data.

    GeomHandle<Mu2eEnvelope> env;

    const double worldBottomInMu2e = env->ymin() - c.getDouble("world.margin.bottom");
    const double worldTopInMu2e    = env->ymax() + c.getDouble("world.margin.top");
    const double worldXminInMu2e   = env->xmin() - c.getDouble("world.margin.xmin");
    const double worldXmaxInMu2e   = env->xmax() + c.getDouble("world.margin.xmax");
    const double worldZminInMu2e   = env->zmin() - c.getDouble("world.margin.zmin");
    const double worldZmaxInMu2e   = env->zmax() + c.getDouble("world.margin.zmax");

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
    res->_hallFormalHalfSize[0] = (env->xmax() - env->xmin())/2;
    res->_hallFormalHalfSize[1] = (env->ymax() - env->ymin())/2;
    res->_hallFormalHalfSize[2] = (env->zmax() - env->zmin())/2;

    const CLHEP::Hep3Vector hallFormalCenterInMu2e(
                                                   (env->xmax() + env->xmin())/2,
                                                   (env->ymax() + env->ymin())/2,
                                                   (env->zmax() + env->zmin())/2
                                                   );

    res->_hallFormalCenterInWorld = hallFormalCenterInMu2e + res->_mu2eOriginInWorld;

    if ( diagLevel > 0) {
      std::cout<<*env<<std::endl;
      std::cout<<*res<<std::endl;
    }

    return res;

  } // make()

} // namespace mu2e
