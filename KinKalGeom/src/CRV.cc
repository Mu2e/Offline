#include "Offline/KinKalGeom/inc/CRV.hh"

/*
  // 04/04/25 : made with newest geometry. Rectangles are set within CRV modules so that they are at the depth of the center of the first layer. In offline they are treated as volumes.
  // norm here = layerDirection in offline
  // udir = gapDirection in offline
  // simplified all sectors into single rectangle with additional rectangles representing holes.
*/

namespace mu2e {
  namespace KinKalGeom {
    using KinKal::VEC3;
    using KinKal::Rectangle;
    // currently use hard-coded geometry
    CRV::CRV() :
      rightSectorBlock_{ std::make_shared<Rectangle>(VEC3(-1.0,0.0,0.0),VEC3(0.0,0.0,1.0), VEC3(-2537.19,275.3,-2008.77),10331.25,2275.)},
      rightSectorHole_{ std::make_shared<Rectangle>(VEC3(-1.0,0.0,0.0),VEC3(0.0,0.0,1.0), VEC3(-2537.19,-1324.47,6669.48),1653.,675.)},
      cryoSectorHole_{ std::make_shared<Rectangle>(VEC3(-1.0,0.0,0.0),VEC3(0.0,0.0,1.0), VEC3(-2537.19,1273.03,2950.23),413.25,232.5)},
      leftSectorBlock_{ std::make_shared<Rectangle>(VEC3(1.0,0.0,0.0),VEC3(0.0,0.0,1.0), VEC3(2537.19,275.53,470.73),7851.75,2275.)},
      leftSectorHole_{ std::make_shared<Rectangle>(VEC3(1.0,0.0,0.0),VEC3(0.0,0.0,1.0), VEC3(2537.19,-1324.47,6669.48),1653.,675.)},
      topSectorBlock_{ std::make_shared<Rectangle>(VEC3(0.0,1.0,0.0),VEC3(0.0,0.0,1.0), VEC3(0.,2663.12,-2008.77),10331.25,3000.)},
      extSectorBlock_{ std::make_shared<Rectangle>(VEC3(0.0,1.0,0.0),VEC3(-1.0,0.0,0.0), VEC3(3595.25,2917.12,-9985.2),826.5,2500.)},
      upstreamSectorBlock_{ std::make_shared<Rectangle>(VEC3(0.0,0.0,1.0),VEC3(0.0,-1.0,0.0), VEC3(650.,1597.92,-12599.76),1653.,3450.)},
      downstreamSectorBlock_{ std::make_shared<Rectangle>(VEC3(0.0,0.0,-1.0),VEC3(0.0,-1.0,0.0), VEC3(0.,933.17,8472.26),2066.25,2850.)},
      downstreamSectorHole_{ std::make_shared<Rectangle>(VEC3(0.0,0.0,-1.0),VEC3(0.0,-1.0,0.0), VEC3(0.,106.67,8472.26),413.25,480.)},
      cryoSector1_{ std::make_shared<Rectangle>(VEC3(-1.0,0.0,0.0),VEC3(0.0,-1.0,0.0), VEC3(-3226.34,2006.65,3340.7),516.5625,850.)},
      cryoSector2_{ std::make_shared<Rectangle>(VEC3(-1.0,0.0,0.0),VEC3(0.0,-1.0,0.0), VEC3(-3226.34,1079.,3618.2),103.3125,572.5)}
    {}
  }
}
