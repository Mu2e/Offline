#include "Offline/KinKalGeom/inc/TestCRV.hh"
namespace mu2e {
  namespace KinKalGeom {
    using KinKal::VEC3;
    using KinKal::Rectangle;
    // currently use hard-coded geometry
    TestCRV::TestCRV() :
      ex1_{ std::make_shared<Rectangle>(VEC3(0.0,1.0,0.0),VEC3(1.0,0.0,0.0), VEC3(0.0,4775,-438),3000,1675)}, // layer widths are approximate FIXME
      t1_{ std::make_shared<Rectangle>(VEC3(0.0,1.0,0.0),VEC3(0.0,0.0,1.0), VEC3(0.0,4625,-438),1185,850)},
      t2_{ std::make_shared<Rectangle>(VEC3(0.0,1.0,0.0),VEC3(1.0,0.0,0.0), VEC3(0.0,4925,-438),1600,850)} {}
  }
}
