#include "Offline/KinKalGeom/inc/CRV.hh"

/*
  // 04/04/25 : made with newest geometry. Rectangles are set within CRV modules so that they are at the depth of the center of the first layer. In offline they are treated as volumes.
  // norm here = layerDirection in offline
  // udir = gapDirection in offline
*/

namespace mu2e {
  namespace KinKalGeom {
    using KinKal::VEC3;
    using KinKal::Rectangle;
    // currently use hard-coded geometry
    CRV::CRV() :
      rightSector1_{ std::make_shared<Rectangle>(VEC3(-1.0,0.0,0.0),VEC3(0.0,0.0,1.0), VEC3(-2537.19,275.3,-10273.77),2066.25,2275.)},
      rightSector2_{ std::make_shared<Rectangle>(VEC3(-1.0,0.0,0.0),VEC3(0.0,0.0,1.0), VEC3(-2537.19,275.3,-2835.27),5372.25,2275.)},
      rightSector3_{ std::make_shared<Rectangle>(VEC3(-1.0,0.0,0.0),VEC3(0.0,0.0,1.0), VEC3(-2537.19,2058.03,2950.23),413.25,522.5)},
      rightSector4_{ std::make_shared<Rectangle>(VEC3(-1.0,0.0,0.0),VEC3(0.0,0.0,1.0), VEC3(-2537.19,-471.97,2950.23),413.25,1520.)},
      rightSector5_{ std::make_shared<Rectangle>(VEC3(-1.0,0.0,0.0),VEC3(0.0,0.0,1.0), VEC3(-2537.19,275.3,4189.98),826.5,2275.)},
      rightSector6_{ std::make_shared<Rectangle>(VEC3(-1.0,0.0,0.0),VEC3(0.0,0.0,1.0), VEC3(-2537.19,950.53,6669.48),1653.,1600.)},
      leftSector1_{ std::make_shared<Rectangle>(VEC3(1.0,0.0,0.0),VEC3(0.0,0.0,1.0), VEC3(2537.19,275.53,-5728.02),1653.,2275.)},
      leftSector2_{ std::make_shared<Rectangle>(VEC3(1.0,0.0,0.0),VEC3(0.0,0.0,1.0), VEC3(2537.19,275.53,470.73),4545.75,2275.)},
      leftSector3_{ std::make_shared<Rectangle>(VEC3(1.0,0.0,0.0),VEC3(0.0,0.0,1.0), VEC3(2537.19,950.53,6669.48),1653.,1600.)},
      topSector1_{ std::make_shared<Rectangle>(VEC3(0.0,1.0,0.0),VEC3(0.0,0.0,1.0), VEC3(0.,2663.12,-11100.27),1239.5,3000.)},
      topSector2_{ std::make_shared<Rectangle>(VEC3(0.0,1.0,0.0),VEC3(0.0,0.0,1.0), VEC3(0.,2663.12,-9034.02),826.5,3000.)},
      topSector3_{ std::make_shared<Rectangle>(VEC3(0.0,1.0,0.0),VEC3(0.0,0.0,1.0), VEC3(0.,2663.12,-6554.52),1653.,3000.)},
      topSector4_{ std::make_shared<Rectangle>(VEC3(0.0,1.0,0.0),VEC3(0.0,0.0,1.0), VEC3(0.,2663.12,470.73),5372.25,3000.)},
      topSector5_{ std::make_shared<Rectangle>(VEC3(0.0,1.0,0.0),VEC3(0.0,0.0,1.0), VEC3(0.,2663.12,7082.73),1239.75,3000.)},
      extSector1_{ std::make_shared<Rectangle>(VEC3(0.0,1.0,0.0),VEC3(-1.0,0.0,0.0), VEC3(4034.17,2917.12,-9985.2),413.25,2500.)},
      extSector2_{ std::make_shared<Rectangle>(VEC3(0.0,1.0,0.0),VEC3(-1.0,0.0,0.0), VEC3(3207.67,2917.12,-9985.2),413.25,2500.)},
      upstreamSector1_{ std::make_shared<Rectangle>(VEC3(0.0,0.0,1.0),VEC3(0.0,-1.0,0.0), VEC3(650.,1597.12,-12599.76),1653.,3450.)},
      downstreamSector1_{ std::make_shared<Rectangle>(VEC3(0.0,0.0,-1.0),VEC3(0.0,-1.0,0.0), VEC3(0.,1759.67,8472.26),1239.75,2850.)},
      downstreamSector2_{ std::make_shared<Rectangle>(VEC3(0.0,0.0,-1.0),VEC3(0.0,-1.0,0.0), VEC3(-1665.,106.67,8472.26),413.25,1185.)},
      downstreamSector3_{ std::make_shared<Rectangle>(VEC3(0.0,0.0,-1.0),VEC3(0.0,-1.0,0.0), VEC3(1665.,106.67,8472.26),413.25,1185.)},
      downstreamSector4_{ std::make_shared<Rectangle>(VEC3(0.0,0.0,-1.0),VEC3(0.0,-1.0,0.0), VEC3(0.,-719.83,8472.26),413.25,2850.)},
      cryoSector1_{ std::make_shared<Rectangle>(VEC3(-1.0,0.0,0.0),VEC3(0.0,-1.0,0.0), VEC3(-7130.34,2006.65,3340.7),413.25,850.)},
      cryoSector2_{ std::make_shared<Rectangle>(VEC3(-1.0,0.0,0.0),VEC3(0.0,-1.0,0.0), VEC3(-7130.34,1079.,3618.2),413.25,572.5)}
    {}
  }
}
