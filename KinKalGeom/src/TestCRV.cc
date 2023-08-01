#include "Offline/KinKalGeom/inc/TestCRV.hh"
namespace mu2e {
  namespace KinKalGeom {
    using KinKal::VEC3;
    using KinKal::Rectangle;
    // currently use hard-coded geometry
    TestCRV::TestCRV() :
      ex1_(VEC3(0.0,1.0,0.0),VEC3(1.0,0.0,0.0), VEC3(0.0,4730,-438),3000,1675), // layer widths are approximate FIXME
      t1_(VEC3(0.0,1.0,0.0),VEC3(0.0,0.0,1.0), VEC3(0.0,4580,-438),1185,850),
      t2_(VEC3(0.0,1.0,0.0),VEC3(1.0,0.0,0.0), VEC3(0.0,4880,-438),1600,850) {}
  }
}
