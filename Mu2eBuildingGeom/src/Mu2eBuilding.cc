#include "Mu2eBuildingGeom/inc/Mu2eBuilding.hh"

namespace mu2e {
  Mu2eBuilding::Mu2eBuilding(const BuildingBasics& basics)
    : basics_(basics)
    , _hallInsideXmin(0.)
    , _hallInsideXmax(0.)
    , _hallInsideZmax(0.)
    , _hallInsideXPSCorner(0.)
    , _hallInsideZPSCorner(0.)
    , _hallInsideZPStoBeamDumpCorner(0.)
    , _hallInsideZExtMonUCIWall(0.)
    , _hallWallThickness(0.)
    , _hallWallExtMonUCIThickness(0.)
  {}
}
