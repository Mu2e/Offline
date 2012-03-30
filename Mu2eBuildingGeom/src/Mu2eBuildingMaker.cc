#include "Mu2eBuildingGeom/inc/Mu2eBuildingMaker.hh"

#include <iostream>

#include "cetlib/exception.h"

#include "CLHEP/Units/SystemOfUnits.h"

#include "Mu2eUtilities/inc/SimpleConfig.hh"

#include "Mu2eBuildingGeom/inc/Mu2eBuilding.hh"

#include "ProtonBeamDumpGeom/inc/ProtonBeamDump.hh"


namespace mu2e {

  std::auto_ptr<Mu2eBuilding> Mu2eBuildingMaker::make(const SimpleConfig& c, const ProtonBeamDump& dump)
  {
    std::auto_ptr<Mu2eBuilding> _b(new Mu2eBuilding());

    using CLHEP::Hep2Vector;

    _b->_hallInsideXmin = c.getDouble("hall.insideXmin");
    _b->_hallInsideXmax = c.getDouble("hall.insideXmax");
    _b->_hallInsideZmax = c.getDouble("hall.insideZmax");

    _b->_hallInsideXPSCorner = c.getDouble("hall.insideXPSCorner");
    _b->_hallInsideZPSCorner = c.getDouble("hall.insideZPSCorner");

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
    _b->_hallWallExtMonUCIThickness = c.getDouble("hall.wallExtMonUCIThick");

    // Origin used to construct the MECO detector.
    // Magic number to fix:
    _b->_trackerOriginInMu2e = CLHEP::Hep3Vector( -3904., 0., 12000.);

    //----------------
    //  Computed stuff

    _b->_hallInsideZPStoBeamDumpCorner = dump.enclosureCenterInMu2e()[2]
      + (_b->_hallInsideXPSCorner - dump.enclosureCenterInMu2e()[0])/tan(dump.coreRotY())
      + dump.enclosureHalfSize()[0]/sin(dump.coreRotY())
      ;

    //----------------
    //  Check the assumptions used to construct G4 geometry downstream

    if(_b->_hallInsideZPStoBeamDumpCorner >= _b->_hallInsideZPSCorner) {
      throw cet::exception("GEOM")<<"Mu2eBuildingMaker: can not satisfy (hall.insideXPSCorner, hall.insideZPSCorner)\n";
    }

    if(dump.shieldingFaceZatXmax() <= _b->_hallInsideZExtMonUCIWall) {
      throw cet::exception("GEOM")<<"Mu2eBuildingMaker: hallInsideZExtMonUCIWall is too large - conflicts with ProtonBeamDump\n";
    }

    //----------------------------------------------------------------
    // Fragment 1
    _b->_concreteOuterOutline1.push_back(Hep2Vector(dump.shieldingFaceXmax(), dump.shieldingFaceZatXmax()));
    _b->_hallInsideOutline.push_back(_b->_concreteOuterOutline1.back());

    const double concreteZmin = _b->hallInsideZExtMonUCIWall() - _b->hallWallThickness();
    if(_b->hallWallExtMonUCIThickness() < (dump.shieldingFaceZatXmax() - concreteZmin)*tan(dump.coreRotY())) {
      // two-point case
      _b->_concreteOuterOutline1.push_back(Hep2Vector(dump.shieldingFaceXmax() - _b->hallWallExtMonUCIThickness(),
                                                      dump.shieldingFaceZatXmax() - _b->hallWallExtMonUCIThickness()
                                                      /tan(dump.coreRotY())));

      _b->_concreteOuterOutline1.push_back(Hep2Vector(dump.shieldingFaceXmax() - _b->hallWallExtMonUCIThickness(),
                                                      concreteZmin));
    }
    else { // one-point case
      const double dz = dump.shieldingFaceZatXmax() - concreteZmin;
      _b->_concreteOuterOutline1.push_back(Hep2Vector(dump.shieldingFaceXmax() - dz*tan(dump.coreRotY()),
                                                      concreteZmin));
    }
    _b->_hallInsideOutline.push_back(Hep2Vector(dump.shieldingFaceXmax(), concreteZmin + _b->hallWallThickness()));

    _b->_concreteOuterOutline1.push_back(Hep2Vector(_b->hallInsideXmax() + _b->hallWallThickness(), concreteZmin));
    _b->_hallInsideOutline.push_back(Hep2Vector(_b->hallInsideXmax(), concreteZmin + _b->hallWallThickness()));

    //----------------------------------------------------------------
    // Fragment 2
    _b->_concreteOuterOutline2.push_back(Hep2Vector(_b->hallInsideXmax() + _b->hallWallThickness(),
                                                    _b->hallInsideZmax() + _b->hallWallThickness()));

    _b->_hallInsideOutline.push_back(Hep2Vector(_b->hallInsideXmax(), _b->hallInsideZmax()));


    _b->_concreteOuterOutline2.push_back(Hep2Vector(_b->hallInsideXmin() - _b->hallWallThickness(),
                                                    _b->hallInsideZmax() + _b->hallWallThickness()));

    _b->_hallInsideOutline.push_back(Hep2Vector(_b->hallInsideXmin(), _b->hallInsideZmax()));

    //----------------------------------------------------------------
    // Fragment 3

    _b->_concreteOuterOutline3.push_back(Hep2Vector(_b->hallInsideXmin() - _b->hallWallThickness(),
                                                    _b->hallInsideZPSCorner() - _b->hallWallThickness()));

    _b->_hallInsideOutline.push_back(Hep2Vector(_b->hallInsideXmin(), _b->hallInsideZPSCorner()));

    //----------------
    _b->_concreteOuterOutline3.push_back(Hep2Vector(_b->hallInsideXPSCorner() - _b->hallWallThickness(),
                                                    _b->hallInsideZPSCorner() - _b->hallWallThickness()));

    _b->_hallInsideOutline.push_back(Hep2Vector(_b->hallInsideXPSCorner(), _b->hallInsideZPSCorner()));

    //----------------
    _b->_concreteOuterOutline3.push_back(Hep2Vector(_b->hallInsideXPSCorner() - _b->hallWallThickness(),
                                                    _b->hallInsideZPStoBeamDumpCorner() + _b->hallWallThickness()*tan(dump.coreRotY()/2)
                                                    ));

    _b->_hallInsideOutline.push_back(Hep2Vector(_b->hallInsideXPSCorner(), _b->hallInsideZPStoBeamDumpCorner()));

    //----------------
    _b->_concreteOuterOutline3.push_back(Hep2Vector(dump.shieldingFaceXmin()
                                                    - _b->hallWallThickness()*cos(dump.coreRotY())
                                                    ,
                                                    dump.shieldingFaceZatXmin()
                                                    + _b->hallWallThickness()*sin(dump.coreRotY())
                                                    ));

    _b->_hallInsideOutline.push_back(Hep2Vector(dump.shieldingFaceXmin(), dump.shieldingFaceZatXmin()));

    //----------------------------------------------------------------
    const int diagLevel = c.getInt("world.verbosityLevel", 0);
    if(diagLevel > 0) {
      std::cout << __func__ << " trackerOriginInMu2e : " <<  _b->_trackerOriginInMu2e  << std::endl;
    }

    return _b;
  }
}
