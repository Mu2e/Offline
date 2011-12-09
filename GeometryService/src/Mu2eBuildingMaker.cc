#include "GeometryService/inc/Mu2eBuildingMaker.hh"

#include <iostream>

#include "cetlib/exception.h"

#include "CLHEP/Units/SystemOfUnits.h"

#include "Mu2eUtilities/inc/SimpleConfig.hh"

#include "GeometryService/inc/Mu2eBuilding.hh"

#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/ProtonBeamDump.hh"


namespace mu2e {

  Mu2eBuildingMaker::Mu2eBuildingMaker(const SimpleConfig& c) 
    : _b(new Mu2eBuilding())
  {
    using CLHEP::Hep2Vector;

    _b->_hallInsideXmin = c.getDouble("hall.insideXmin");
    _b->_hallInsideXmax = c.getDouble("hall.insideXmax");
    _b->_hallInsideZmax = c.getDouble("hall.insideZmax");

    _b->_hallInsideZBeamDumpWall = c.getDouble("hall.insideZBeamDumpWall");
    _b->_hallInsideZExtMonUCIWall = c.getDouble("hall.insideZExtMonUCIWall");

    _b->_hallInsideYmin = -c.getDouble("mu2e.origin.heightAboveHallFloor");
    _b->_hallInsideYmax = _b->_hallInsideYmin + c.getDouble("hall.insideFullHeight");

    _b->_dirtOverburdenDepth  = c.getDouble("dirt.overburdenDepth");
    _b->_dirtCapHalfHeight    = c.getDouble("dirt.capHalfHeight");
    _b->_dirtCapBottomRadius    = c.getDouble("dirt.capBottomRadius");
    _b->_dirtCapTopRadius    = c.getDouble("dirt.capTopRadius");

    _b->_hallFloorThickness = c.getDouble("hall.floorThick");
    _b->_hallCeilingThickness = c.getDouble("hall.ceilingThick");
    _b->_hallWallThickness = c.getDouble("hall.wallThick");

    // Origin used to construct the MECO detector.
    // Magic number to fix:
    _b->_trackerOriginInMu2e = CLHEP::Hep3Vector( -3904., 0., 12000.);

    //----------------
    //  Computed stuff
    GeomHandle<ProtonBeamDump> dump;
    
    _b->_hallInsideXmaxAtBeamDumpWall = dump->enclosureCenterInMu2e()[0]
      + (_b->_hallInsideZBeamDumpWall - dump->enclosureCenterInMu2e()[2]) * tan(dump->coreRotY())
      - dump->enclosureHalfSize()[0] / cos(dump->coreRotY())
      ;

    //----------------
    //  Check the assumptions used to construct G4 geometry downstream
    
    if(dump->shieldingFaceZatXmin() > _b->_hallInsideZBeamDumpWall) {
      throw cet::exception("GEOM")<<"Mu2eBuildingMaker: hallInsideZBeamDumpWall is too small - conflicts with ProtonBeamDump\n";
    }
    if(dump->shieldingFaceZatXmax() <= _b->_hallInsideZExtMonUCIWall) {
      throw cet::exception("GEOM")<<"Mu2eBuildingMaker: hallInsideZExtMonUCIWall is too large - conflicts with ProtonBeamDump\n";
    }

    //----------------------------------------------------------------
    // Fragment 1
    _b->_concreteOuterOutline1.push_back(Hep2Vector(dump->shieldingFaceXmax(), dump->shieldingFaceZatXmax()));

    const double concreteZmin = _b->hallInsideZExtMonUCIWall() - _b->hallWallThickness();
    if(_b->hallWallThickness() < (dump->shieldingFaceZatXmax() - concreteZmin)*tan(dump->coreRotY())) {
      // two-point case
      _b->_concreteOuterOutline1.push_back(Hep2Vector(dump->shieldingFaceXmax() - _b->hallWallThickness(),
						      dump->shieldingFaceZatXmax() - _b->hallWallThickness()/tan(dump->coreRotY())));

      _b->_concreteOuterOutline1.push_back(Hep2Vector(dump->shieldingFaceXmax() - _b->hallWallThickness(),
						      concreteZmin));
    }
    else { // one-point case
      const double dz = dump->shieldingFaceZatXmax() - concreteZmin;
      _b->_concreteOuterOutline1.push_back(Hep2Vector(dump->shieldingFaceXmax() - dz*tan(dump->coreRotY()),
						      concreteZmin));
    }

    //----------------------------------------------------------------
    // Fragment 2
    _b->_concreteOuterOutline2.push_back(Hep2Vector(_b->hallInsideXmax() + _b->hallWallThickness(),
						    concreteZmin));
    
    _b->_concreteOuterOutline2.push_back(Hep2Vector(_b->hallInsideXmax() + _b->hallWallThickness(),
						    _b->hallInsideZmax() + _b->hallWallThickness()));
    
    
    _b->_concreteOuterOutline2.push_back(Hep2Vector(_b->hallInsideXmin() - _b->hallWallThickness(),
						    _b->hallInsideZmax() + _b->hallWallThickness()));
    
    _b->_concreteOuterOutline2.push_back(Hep2Vector(_b->hallInsideXmin() - _b->hallWallThickness(),
						    _b->hallInsideZBeamDumpWall() - _b->hallWallThickness()));
    
    //----------------------------------------------------------------
    // Fragment 3

    _b->_concreteOuterOutline3.push_back(Hep2Vector(_b->hallInsideXmaxAtBeamDumpWall() 
						    - _b->hallWallThickness()/cos(dump->coreRotY())
						    - _b->hallWallThickness()*tan(dump->coreRotY())
						    ,
						    _b->hallInsideZBeamDumpWall() - _b->hallWallThickness()
						    ));
    
    _b->_concreteOuterOutline3.push_back(Hep2Vector(dump->shieldingFaceXmin()
						    - _b->hallWallThickness()*cos(dump->coreRotY())
						    ,
						    dump->shieldingFaceZatXmin()
						    + _b->hallWallThickness()*sin(dump->coreRotY())
						    ));

    //----------------------------------------------------------------
    const int diagLevel = c.getInt("world.verbosityLevel", 0);
    if(diagLevel > 0) {
      std::cout << __func__ << " trackerOriginInMu2e : " <<  _b->_trackerOriginInMu2e  << std::endl;
    }
  }
}
