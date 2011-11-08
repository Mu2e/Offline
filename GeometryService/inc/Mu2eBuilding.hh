// Geometry of the hall, dirt, etc.
// 
// Andrei Gaponenko, 2011

#ifndef MU2EBUILDING_HH
#define MU2EBUILDING_HH

#include <vector>

#include "CLHEP/Vector/ThreeVector.h"

#include "GeometryService/inc/Detector.hh"

namespace mu2e {

  class Mu2eBuildingMaker;

  class Mu2eBuilding : public Detector {
  public:

    const CLHEP::Hep3Vector& hallCenterInMu2e() const { return _hallCenterInMu2e; }

    const std::vector<double>& hallInsideHalfLengths() const { return _hallInsideHalfLenghts; }
    double hallFloorThickness() const { return _hallFloorThickness; }
    double hallCeilingThickness() const { return _hallCeilingThickness; }
    double hallWallThickness() const { return _hallWallThickness; }

    double dirtOverburdenDepth() const { return _dirtOverburdenDepth; }
    double dirtCapHalfHeight() const { return _dirtCapHalfHeight; }
    double dirtCapBottomRadius() const { return _dirtCapBottomRadius; }
    double dirtCapTopRadius() const { return _dirtCapTopRadius; }

    const CLHEP::Hep3Vector& trackerOriginInMu2e() const { return _trackerOriginInMu2e; }    

    // implement Detector's method
    virtual std::string name() const { return "Mu2eBuilding"; }

    //----------------------------------------------------------------
  private: 
    friend class Mu2eBuildingMaker;

    // Private ctr: the class should be only obtained via the maker
    Mu2eBuilding() {}

    CLHEP::Hep3Vector _hallCenterInMu2e;
    CLHEP::Hep3Vector _trackerOriginInMu2e;
    std::vector<double> _hallInsideHalfLenghts;
    double _hallFloorThickness;
    double _hallCeilingThickness;
    double _hallWallThickness;

    double _dirtOverburdenDepth;
    double _dirtCapHalfHeight;
    double _dirtCapBottomRadius;
    double _dirtCapTopRadius;
  };

}

#endif/*MU2EBUILDING_HH*/
