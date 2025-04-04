#include "Offline/KinKalGeom/inc/CRV.hh"
namespace mu2e {
  namespace KinKalGeom {
    using KinKal::VEC3;
    using KinKal::Rectangle;
    // currently use hard-coded geometry
    CRV::CRV() :
      r1_{ std::make_shared<Rectangle>(VEC3(-1.0,0.0,0.0),VEC3(0.0,0.0,1.0), VEC3(-2537.28,275.53,-10281.1),2070,2275)},
      r2_{ std::make_shared<Rectangle>(VEC3(-1.0,0.0,0.0),VEC3(0.0,0.0,1.0), VEC3(-2537.28,275.53,-2829.1),5382,2275)},
      r3_{ std::make_shared<Rectangle>(VEC3(-1.0,0.0,0.0),VEC3(0.0,0.0,1.0), VEC3(-2537.28,2028.03,2966.9),414,522.5)},
      r4_{ std::make_shared<Rectangle>(VEC3(-1.0,0.0,0.0),VEC3(0.0,0.0,1.0), VEC3(-2537.28,-479.47,2966.9),414,1520)},
      r5_{ std::make_shared<Rectangle>(VEC3(-1.0,0.0,0.0),VEC3(0.0,0.0,1.0), VEC3(-2537.28,275.53,4213.9),828,522.5)},
      r6_{ std::make_shared<Rectangle>(VEC3(-1.0,0.0,0.0),VEC3(0.0,0.0,1.0), VEC3(-2537.28,950.53,5460.9),1656,522.5)},
      l1_{ std::make_shared<Rectangle>(VEC3(-1.0,0.0,0.0),VEC3(0.0,0.0,1.0), VEC3(2537.28,275.53,-5727.1),1656,2275)},
      l2_{ std::make_shared<Rectangle>(VEC3(-1.0,0.0,0.0),VEC3(0.0,0.0,1.0), VEC3(2537.28,275.53,482.9),4554,2275)},
      l3_{ std::make_shared<Rectangle>(VEC3(-1.0,0.0,0.0),VEC3(0.0,0.0,1.0), VEC3(2537.28,950.53,6692.9),1656,1600)},
      t1_{ std::make_shared<Rectangle>(VEC3(0.0,1.0,0.0),VEC3(0.0,0.0,1.0), VEC3(0,2663.21,-11109.1),1242,3000)},
      t2_{ std::make_shared<Rectangle>(VEC3(0.0,1.0,0.0),VEC3(0.0,0.0,1.0), VEC3(0,2663.21,-9039.1),828,3000)},
      t3_{ std::make_shared<Rectangle>(VEC3(0.0,1.0,0.0),VEC3(0.0,0.0,1.0), VEC3(0,2663.21,-6555.1),1656,3000)},
      t4_{ std::make_shared<Rectangle>(VEC3(0.0,1.0,0.0),VEC3(0.0,0.0,1.0), VEC3(0,2663.21,482.9),5382,3000)},
      t5_{ std::make_shared<Rectangle>(VEC3(0.0,1.0,0.0),VEC3(0.0,0.0,1.0), VEC3(0,2663.21,7106.9),1242,3000)},
      e1_{ std::make_shared<Rectangle>(VEC3(0.0,-1.0,0.0),VEC3(1.0,0.0,0.0), VEC3(4033.5,2837.58,-9986.8),414,2500)},
      e2_{ std::make_shared<Rectangle>(VEC3(0.0,-1.0,0.0),VEC3(1.0,0.0,0.0), VEC3(6398,2837.58,-9986.8),414,2500)},
      u1_{ std::make_shared<Rectangle>(VEC3(0.0,0.0,-1.0),VEC3(1.0,0.0,0.0), VEC3(1453,636.2,8465.76),3450,1656)},
      d1_{ std::make_shared<Rectangle>(VEC3(0.0,0.0,-1.0),VEC3(1.0,0.0,0.0), VEC3(0,1719,8465.76),2850,1242)},
      d2_{ std::make_shared<Rectangle>(VEC3(0.0,0.0,-1.0),VEC3(1.0,0.0,0.0), VEC3(-1665,63,8465.76),1185,424)},
      d3_{ std::make_shared<Rectangle>(VEC3(0.0,0.0,-1.0),VEC3(1.0,0.0,0.0), VEC3(1665,63,8465.76),1185,424)},
      d4_{ std::make_shared<Rectangle>(VEC3(0.0,0.0,-1.0),VEC3(1.0,0.0,0.0), VEC3(0,-765,8465.76),2850,424)},
      c1_{ std::make_shared<Rectangle>(VEC3(-1.0,0.0,0.0),VEC3(0.0,0.0,1.0), VEC3(-2687.07,1329.88,2616.75),414,450)},
      c2_{ std::make_shared<Rectangle>(VEC3(-1.0,0.0,0.0),VEC3(0.0,0.0,1.0), VEC3(-2687.07,1329.88,3483.4),414,450)},
      c3_{ std::make_shared<Rectangle>(VEC3(-1.0,0.0,0.0),VEC3(0.0,0.0,1.0), VEC3(-3216.21,1534.75,3947.4),414,1050)},
      c4_{ std::make_shared<Rectangle>(VEC3(-1.0,0.0,0.0),VEC3(0.0,0.0,1.0), VEC3(-3216.21,697.15,3808.4),414,1050)} {}
  }
}
