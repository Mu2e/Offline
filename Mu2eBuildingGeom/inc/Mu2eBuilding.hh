// Geometry of the hall, dirt, etc.
//
// Andrei Gaponenko, 2011

#ifndef MU2EBUILDING_HH
#define MU2EBUILDING_HH

#include <vector>

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/TwoVector.h"

#include "Mu2eInterfaces/inc/Detector.hh"
#include "Mu2eBuildingGeom/inc/BuildingBasics.hh"

#include "art/Persistency/Common/Wrapper.h"

namespace mu2e {

  class Mu2eBuildingMaker;

  class Mu2eBuilding : virtual public Detector {
  public:

    double hallInsideXmin() const { return _hallInsideXmin; }
    double hallInsideXmax() const { return _hallInsideXmax; }

    double hallInsideYmin() const { return basics_.detectorHallFloorTopY(); }
    double hallInsideYmax() const { return basics_.detectorHallFloorTopY() + basics_.detectorHallInsideFullHeight(); }

    double hallInsideZmax() const { return _hallInsideZmax; }

    double hallInsideXPSCorner() const { return _hallInsideXPSCorner; }
    double hallInsideZPSCorner() const { return _hallInsideZPSCorner; }
    double hallInsideZPStoBeamDumpCorner() const { return _hallInsideZPStoBeamDumpCorner; }

    double hallInsideZExtMonUCIWall() const { return _hallInsideZExtMonUCIWall; }

    double hallFloorThickness() const { return basics_.detectorHallFloorThickness(); }
    double hallCeilingThickness() const { return basics_.detectorHallCeilingThickness(); }
    double hallWallThickness() const { return _hallWallThickness; }
    double hallWallExtMonUCIThickness() const { return _hallWallExtMonUCIThickness; }

    const CLHEP::Hep3Vector& trackerOriginInMu2e() const { return _trackerOriginInMu2e; }

    //----------------------------------------------------------------
    // Outlines used to create hall floor, ceiling, and walls.
    // Points go clockwise, as requred by G4 extruded solid.

    // Starts at the Xmax dump shielding face corner and goes to
    // (hallInsideXmax+wallThick, hallInsideZExtMonUCIWall-wallThick)
    const std::vector<CLHEP::Hep2Vector>& concreteOuterOutline1() const { return _concreteOuterOutline1; }

    // Just two points at the outer Zmax line, from Xmax to Xmin
    const std::vector<CLHEP::Hep2Vector>& concreteOuterOutline2() const { return _concreteOuterOutline2; }

    // Outline continuation that starts starts at (Xmin-wallThick) and
    // goes to the Xmin dump face shielding corner.
    const std::vector<CLHEP::Hep2Vector>& concreteOuterOutline3() const { return _concreteOuterOutline3; }

    // The inner outline of hall walls, all in one piece, clockwise
    // from the Xmax dump shielding face core around the hall and to
    // the Xmin dump shielding face corner.
    const std::vector<CLHEP::Hep2Vector>& hallInsideOutline() const { return _hallInsideOutline; }

    //----------------------------------------------------------------
  private:
    friend class Mu2eBuildingMaker;

    // Private ctr: the class should be only obtained via the maker
    Mu2eBuilding(const BuildingBasics& basics = BuildingBasics() /* need default for I/O */);

    // Or read back from persistent storage
    template<class T> friend class art::Wrapper;

    BuildingBasics basics_;

    double _hallInsideXmin;
    double _hallInsideXmax;

    double _hallInsideZmax;

    double _hallInsideXPSCorner;
    double _hallInsideZPSCorner;
    double _hallInsideZPStoBeamDumpCorner;

    double _hallInsideZExtMonUCIWall;

    CLHEP::Hep3Vector _trackerOriginInMu2e;

    double _hallWallThickness;
    double _hallWallExtMonUCIThickness;

    std::vector<CLHEP::Hep2Vector> _concreteOuterOutline1;
    std::vector<CLHEP::Hep2Vector> _concreteOuterOutline2;
    std::vector<CLHEP::Hep2Vector> _concreteOuterOutline3;
    std::vector<CLHEP::Hep2Vector> _hallInsideOutline;
  };

}

#endif/*MU2EBUILDING_HH*/
