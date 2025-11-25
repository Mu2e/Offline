#include "Offline/KinKalGeom/inc/Calo.hh"
namespace mu2e {
  namespace KinKalGeom {
    using KinKal::VEC3;
    using KinKal::Cylinder;
    using KinKal::Disk;

    Calo::Calo() :
      EMC_Disk_0_SurfIn_ { std::make_shared<Cylinder>(VEC3(0.0,0.0,1.0),VEC3(0.0,0.0,4.0),850.11,1635.11)},

    {}
  }
}
