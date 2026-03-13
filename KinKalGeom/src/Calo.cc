#include "Offline/KinKalGeom/inc/Calo.hh"
namespace mu2e {
  namespace KinKalGeom {
    using KinKal::VEC3;
    using KinKal::Cylinder;
    using KinKal::Disk;

    // D0 = 11842 -10175.0 = 1667
    // D1 = 13220 -10175.0 = 3045
    /*
    double calorimeter.caloDiskRadiusIn           = 335;
    double calorimeter.caloDiskRadiusOut          = 719;
    half length = 192.295


    Cylinder(VEC3 const& axis, VEC3 const& center, double radius, double halflen );
    Disk(VEC3 const& norm,VEC3 const& udir, VEC3 const& center, double radius)
    */
    Calo::Calo() :
      EMC_Disk_0_Inner_ { std::make_shared<Cylinder>(VEC3(0.0,0.0,1.0),VEC3(0.0,0.0,1667),335,192.295)},
      EMC_Disk_0_Outer_ { std::make_shared<Cylinder>(VEC3(0.0,0.0,1.0),VEC3(0.0,0.0,1667),719,192.295)},
      EMC_Disk_1_Inner_ { std::make_shared<Cylinder>(VEC3(0.0,0.0,1.0),VEC3(0.0,0.0,3045),335,192.295)},
      EMC_Disk_1_Outer_ { std::make_shared<Cylinder>(VEC3(0.0,0.0,1.0),VEC3(0.0,0.0,3045),719,192.295)},

      EMC_Disk_0_Front_ { std::make_shared<Disk>(VEC3(0.0,0.0,1.0),VEC3(1.0,0.0,0.0),VEC3(0.0,0.0,1667-192.295),719.)},
      EMC_Disk_0_Back_ { std::make_shared<Disk>(VEC3(0.0,0.0,1.0),VEC3(1.0,0.0,0.0),VEC3(0.0,0.0,1667+192.295),719.)},
      EMC_Disk_1_Front_ { std::make_shared<Disk>(VEC3(0.0,0.0,1.0),VEC3(1.0,0.0,0.0),VEC3(0.0,0.0,3045-192.295),719.)},
      EMC_Disk_1_Back_ { std::make_shared<Disk>(VEC3(0.0,0.0,1.0),VEC3(1.0,0.0,0.0),VEC3(0.0,0.0,3045+192.295 ),719.)}
    {}
  }
}
