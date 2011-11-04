#include "GeometryService/inc/Mu2eBuildingMaker.hh"

#include <iostream>

#include "cetlib/exception.h"

#include "CLHEP/Units/SystemOfUnits.h"

#include "Mu2eUtilities/inc/SimpleConfig.hh"

#include "GeometryService/inc/Mu2eBuilding.hh"

namespace mu2e {

  Mu2eBuildingMaker::Mu2eBuildingMaker(const SimpleConfig& c) 
    : _b(new Mu2eBuilding())
  {
    _b->_dirtOverburdenDepth  = c.getDouble("dirt.overburdenDepth");
    _b->_dirtCapHalfHeight    = c.getDouble("dirt.capHalfHeight");
    _b->_dirtCapBottomRadius    = c.getDouble("dirt.capBottomRadius");
    _b->_dirtCapTopRadius    = c.getDouble("dirt.capTopRadius");

    c.getVectorDouble("hall.insideHalfLengths",_b->_hallInsideHalfLenghts,3);    
    _b->_hallFloorThickness = c.getDouble("hall.floorThick");
    _b->_hallCeilingThickness = c.getDouble("hall.ceilingThick");
    _b->_hallWallThickness = c.getDouble("hall.wallThick");

    double cx = c.getDouble("hall.centerInMu2e.x");
    double cz = c.getDouble("hall.centerInMu2e.z");
    double h =  c.getDouble("mu2e.origin.heightAboveHallFloor");
    _b->_hallCenterInMu2e = CLHEP::Hep3Vector(cx, _b->_hallInsideHalfLenghts[1] - h, cz);
    
    const int diagLevel = c.getInt("world.verbosityLevel", 0);
    if(diagLevel > 0) {
      std::cout << __func__ << " hallCenterInMu2e = " <<  _b->_hallCenterInMu2e  << std::endl;
    }
  }
}
