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
    std::auto_ptr<Mu2eBuilding> b(new Mu2eBuilding());

    using CLHEP::Hep2Vector;

    b->_hallInsideXmin = c.getDouble("hall.insideXmin");
    b->_hallInsideXmax = c.getDouble("hall.insideXmax");
    b->_hallInsideZmax = c.getDouble("hall.insideZmax");

    b->_hallInsideXDSCorner = c.getDouble("hall.insideXDSCorner");
    b->_hallInsideZDSCorner = c.getDouble("hall.insideZDSCorner");

    b->_hallInsideXPSCorner = c.getDouble("hall.insideXPSCorner");
    b->_hallInsideZPSCorner = c.getDouble("hall.insideZPSCorner");

    b->_hallInsideZExtMonUCIWall = c.getDouble("hall.insideZExtMonUCIWall");

    b->_hallInsideYmin = -c.getDouble("mu2e.origin.heightAboveHallFloor");
    b->_hallInsideYmax = b->_hallInsideYmin + c.getDouble("hall.insideFullHeight");

    b->_dirtOverburdenDepth  = c.getDouble("dirt.overburdenDepth");
    b->_dirtCapHalfHeight    = c.getDouble("dirt.capHalfHeight");
    b->_dirtCapBottomRadius    = c.getDouble("dirt.capBottomRadius");
    b->_dirtCapTopRadius    = c.getDouble("dirt.capTopRadius");

    b->_hallFloorThickness = c.getDouble("hall.floorThick");
    b->_hallCeilingThickness = c.getDouble("hall.ceilingThick");
    b->_hallWallThickness = c.getDouble("hall.wallThick");
    b->_hallWallExtMonUCIThickness = c.getDouble("hall.wallExtMonUCIThick");

    // Origin used to construct the MECO detector.
    // Magic number to fix:
    b->_trackerOriginInMu2e = CLHEP::Hep3Vector( -3904., 0., 12000.);

    //----------------
    //  Computed stuff

    b->_hallInsideZPStoBeamDumpCorner = dump.enclosureCenterInMu2e()[2]
      + (b->_hallInsideXPSCorner - dump.enclosureCenterInMu2e()[0])/tan(dump.coreRotY())
      + dump.enclosureHalfSize()[0]/sin(dump.coreRotY())
      ;

    //----------------
    //  Check the assumptions used to construct G4 geometry downstream

    if(b->_hallInsideZPStoBeamDumpCorner >= b->_hallInsideZPSCorner) {
      throw cet::exception("GEOM")<<"Mu2eBuildingMaker: can not satisfy (hall.insideXPSCorner, hall.insideZPSCorner)\n";
    }

    if(dump.shieldingFaceZatXmax() <= b->_hallInsideZExtMonUCIWall) {
      throw cet::exception("GEOM")<<"Mu2eBuildingMaker: hallInsideZExtMonUCIWall is too large - conflicts with ProtonBeamDump\n";
    }

    //----------------------------------------------------------------
    // Fragment 1
    b->_concreteOuterOutline1.push_back(Hep2Vector(dump.shieldingFaceXmax(), dump.shieldingFaceZatXmax()));
    b->_hallInsideOutline.push_back(b->_concreteOuterOutline1.back());

    const double concreteZmin = b->hallInsideZExtMonUCIWall() - b->hallWallThickness();
    if(b->hallWallExtMonUCIThickness() < (dump.shieldingFaceZatXmax() - concreteZmin)*tan(dump.coreRotY())) {
      // two-point case
      b->_concreteOuterOutline1.push_back(Hep2Vector(dump.shieldingFaceXmax() - b->hallWallExtMonUCIThickness(),
                                                      dump.shieldingFaceZatXmax() - b->hallWallExtMonUCIThickness()
                                                      /tan(dump.coreRotY())));

      b->_concreteOuterOutline1.push_back(Hep2Vector(dump.shieldingFaceXmax() - b->hallWallExtMonUCIThickness(),
                                                      concreteZmin));
    }
    else { // one-point case
      const double dz = dump.shieldingFaceZatXmax() - concreteZmin;
      b->_concreteOuterOutline1.push_back(Hep2Vector(dump.shieldingFaceXmax() - dz*tan(dump.coreRotY()),
                                                      concreteZmin));
    }
    b->_hallInsideOutline.push_back(Hep2Vector(dump.shieldingFaceXmax(), concreteZmin + b->hallWallThickness()));

    b->_concreteOuterOutline1.push_back(Hep2Vector(b->hallInsideXmax() + b->hallWallThickness(), concreteZmin));
    b->_hallInsideOutline.push_back(Hep2Vector(b->hallInsideXmax(), concreteZmin + b->hallWallThickness()));

    //----------------------------------------------------------------
    // Fragment 2
    b->_concreteOuterOutline2.push_back(Hep2Vector(b->hallInsideXmax() + b->hallWallThickness(),
                                                    b->hallInsideZDSCorner() + b->hallWallThickness()));

    b->_hallInsideOutline.push_back(Hep2Vector(b->hallInsideXmax(), b->hallInsideZDSCorner()));

    //----------------
    b->_concreteOuterOutline2.push_back(Hep2Vector(b->hallInsideXDSCorner() + b->hallWallThickness(),
                                                    b->hallInsideZDSCorner() + b->hallWallThickness()));

    b->_hallInsideOutline.push_back(Hep2Vector(b->hallInsideXDSCorner(), b->hallInsideZDSCorner()));

    //----------------
    b->_concreteOuterOutline2.push_back(Hep2Vector(b->hallInsideXDSCorner() + b->hallWallThickness(),
                                                    b->hallInsideZmax() + b->hallWallThickness()));

    b->_hallInsideOutline.push_back(Hep2Vector(b->hallInsideXDSCorner(), b->hallInsideZmax()));

    //----------------

    b->_concreteOuterOutline2.push_back(Hep2Vector(b->hallInsideXmin() - b->hallWallThickness(),
                                                    b->hallInsideZmax() + b->hallWallThickness()));

    b->_hallInsideOutline.push_back(Hep2Vector(b->hallInsideXmin(), b->hallInsideZmax()));

    //----------------------------------------------------------------
    // Fragment 3

    b->_concreteOuterOutline3.push_back(Hep2Vector(b->hallInsideXmin() - b->hallWallThickness(),
                                                    b->hallInsideZPSCorner() - b->hallWallThickness()));

    b->_hallInsideOutline.push_back(Hep2Vector(b->hallInsideXmin(), b->hallInsideZPSCorner()));

    //----------------
    b->_concreteOuterOutline3.push_back(Hep2Vector(b->hallInsideXPSCorner() - b->hallWallThickness(),
                                                    b->hallInsideZPSCorner() - b->hallWallThickness()));

    b->_hallInsideOutline.push_back(Hep2Vector(b->hallInsideXPSCorner(), b->hallInsideZPSCorner()));

    //----------------
    b->_concreteOuterOutline3.push_back(Hep2Vector(b->hallInsideXPSCorner() - b->hallWallThickness(),
                                                    b->hallInsideZPStoBeamDumpCorner() + b->hallWallThickness()*tan(dump.coreRotY()/2)
                                                    ));

    b->_hallInsideOutline.push_back(Hep2Vector(b->hallInsideXPSCorner(), b->hallInsideZPStoBeamDumpCorner()));

    //----------------
    b->_concreteOuterOutline3.push_back(Hep2Vector(dump.shieldingFaceXmin()
                                                    - b->hallWallThickness()*cos(dump.coreRotY())
                                                    ,
                                                    dump.shieldingFaceZatXmin()
                                                    + b->hallWallThickness()*sin(dump.coreRotY())
                                                    ));

    b->_hallInsideOutline.push_back(Hep2Vector(dump.shieldingFaceXmin(), dump.shieldingFaceZatXmin()));

    //----------------------------------------------------------------
    const int diagLevel = c.getInt("world.verbosityLevel", 0);
    if(diagLevel > 0) {
      std::cout << __func__ << " trackerOriginInMu2e : " <<  b->_trackerOriginInMu2e  << std::endl;
    }

    return b;
  }
}
