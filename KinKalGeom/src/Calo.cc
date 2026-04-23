#include "Offline/KinKalGeom/inc/Calo.hh"
namespace mu2e {
  namespace KKGeom {
    using KinKal::VEC3;
    using KinKal::Cylinder;
    using KinKal::Disk;

    Calo::Calo(double z0, double z1, double r0_inner, double r0_outer, double r1_inner, double r1_outer,
               double z0_front, double z0_back, double z1_front, double z1_back) :
      EMC_Disk_0_Inner_ { std::make_shared<Cylinder>(VEC3(0.0,0.0,1.0),VEC3(0.0,0.0,z0),r0_inner, 0.5*(z0_back-z0_front))},
      EMC_Disk_0_Outer_ { std::make_shared<Cylinder>(VEC3(0.0,0.0,1.0),VEC3(0.0,0.0,z0),r0_outer, 0.5*(z0_back-z0_front))},
      EMC_Disk_1_Inner_ { std::make_shared<Cylinder>(VEC3(0.0,0.0,1.0),VEC3(0.0,0.0,z1),r1_inner, 0.5*(z1_back-z1_front))},
      EMC_Disk_1_Outer_ { std::make_shared<Cylinder>(VEC3(0.0,0.0,1.0),VEC3(0.0,0.0,z1),r1_outer, 0.5*(z1_back-z1_front))},

      EMC_Disk_0_Front_ { std::make_shared<Disk>(VEC3(0.0,0.0,1.0),VEC3(1.0,0.0,0.0),VEC3(0.0,0.0,z0_front),r0_outer)},
      EMC_Disk_0_Back_ { std::make_shared<Disk>(VEC3(0.0,0.0,1.0),VEC3(1.0,0.0,0.0),VEC3(0.0,0.0,z0_back),r0_outer)},
      EMC_Disk_1_Front_ { std::make_shared<Disk>(VEC3(0.0,0.0,1.0),VEC3(1.0,0.0,0.0),VEC3(0.0,0.0,z1_front),r1_outer)},
      EMC_Disk_1_Back_ { std::make_shared<Disk>(VEC3(0.0,0.0,1.0),VEC3(1.0,0.0,0.0),VEC3(0.0,0.0,z1_back),r1_outer)},
      z0_front_(z0_front), z0_back_(z0_back), z1_front_(z1_front), z1_back_(z1_back)
    {}
  }
}
