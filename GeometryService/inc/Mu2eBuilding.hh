// Geometry of the hall, dirt, etc.
// 
// Andrei Gaponenko, 2011

#ifndef MU2EBUILDING_HH
#define MU2EBUILDING_HH

#include <vector>

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/TwoVector.h"

#include "Mu2eInterfaces/inc/Detector.hh"

namespace mu2e {

  class Mu2eBuildingMaker;

  class Mu2eBuilding : public Detector {
  public:

    //    const CLHEP::Hep3Vector& hallCenterInMu2e() const { return _hallCenterInMu2e; }
    //    const std::vector<double>& hallInsideHalfLengths() const { return _hallInsideHalfLenghts; }

    double hallInsideXmin() const { return _hallInsideXmin; }
    double hallInsideXmax() const { return _hallInsideXmax; }

    double hallInsideYmin() const { return _hallInsideYmin; }
    double hallInsideYmax() const { return _hallInsideYmax; }

    double hallInsideZmax() const { return _hallInsideZmax; }

    double hallInsideZBeamDumpWall() const { return _hallInsideZBeamDumpWall; }
    double hallInsideXmaxAtBeamDumpWall() const { return _hallInsideXmaxAtBeamDumpWall; }

    double hallInsideZExtMonUCIWall() const { return _hallInsideZExtMonUCIWall; }

    double hallFloorThickness() const { return _hallFloorThickness; }
    double hallCeilingThickness() const { return _hallCeilingThickness; }
    double hallWallThickness() const { return _hallWallThickness; }
    double hallWallExtMonUCIThickness() const { return _hallWallExtMonUCIThickness; }

    double dirtOverburdenDepth() const { return _dirtOverburdenDepth; }
    double dirtCapHalfHeight() const { return _dirtCapHalfHeight; }
    double dirtCapBottomRadius() const { return _dirtCapBottomRadius; }
    double dirtCapTopRadius() const { return _dirtCapTopRadius; }

    const CLHEP::Hep3Vector& trackerOriginInMu2e() const { return _trackerOriginInMu2e; }    

    //----------------------------------------------------------------
    // These coordinates are used in more than one place

    const std::vector<CLHEP::Hep2Vector>& concreteOuterOutline1() const { return _concreteOuterOutline1; }
    const std::vector<CLHEP::Hep2Vector>& concreteOuterOutline2() const { return _concreteOuterOutline2; }
    const std::vector<CLHEP::Hep2Vector>& concreteOuterOutline3() const { return _concreteOuterOutline3; }

    //----------------------------------------------------------------
  private: 
    friend class Mu2eBuildingMaker;

    // Private ctr: the class should be only obtained via the maker
    Mu2eBuilding() {}

    double _hallInsideXmin;
    double _hallInsideXmax;
    double _hallInsideYmin;
    double _hallInsideYmax;

    double _hallInsideZmax;

    double _hallInsideZBeamDumpWall;
    double _hallInsideXmaxAtBeamDumpWall;

    double _hallInsideZExtMonUCIWall;

    CLHEP::Hep3Vector _trackerOriginInMu2e;

    double _hallFloorThickness;
    double _hallCeilingThickness;
    double _hallWallThickness;
    double _hallWallExtMonUCIThickness;

    double _dirtOverburdenDepth;
    double _dirtCapHalfHeight;
    double _dirtCapBottomRadius;
    double _dirtCapTopRadius;

    std::vector<CLHEP::Hep2Vector> _concreteOuterOutline1;
    std::vector<CLHEP::Hep2Vector> _concreteOuterOutline2;
    std::vector<CLHEP::Hep2Vector> _concreteOuterOutline3;
  };

}

#endif/*MU2EBUILDING_HH*/
