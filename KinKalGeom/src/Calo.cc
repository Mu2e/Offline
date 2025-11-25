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

    D0 front = 1667
    D0 back =
    D1 front =
    D1 back =

    Cylinder(VEC3 const& axis, VEC3 const& center, double radius, double halflen );
    Disk(VEC3 const& norm,VEC3 const& udir, VEC3 const& center, double radius)
    */
    Calo::Calo() :
      EMC_Disk_0_SurfIn_ { std::make_shared<Cylinder>(VEC3(0.0,0.0,1.0),VEC3(0.0,0.0,1667),719,192.295)},
      EMC_Disk_0_SurfOut_ { std::make_shared<Cylinder>(VEC3(0.0,0.0,1.0),VEC3(0.0,0.0,1667),719,192.295)},
      EMC_Disk_1_SurfIn_ { std::make_shared<Cylinder>(VEC3(0.0,0.0,1.0),VEC3(0.0,0.0,3045),719,192.295)},
      EMC_Disk_1_SurfOut_ { std::make_shared<Cylinder>(VEC3(0.0,0.0,1.0),VEC3(0.0,0.0,3045),719,192.295)},

      EMC_Disk_0_EdgeIn_ { std::make_shared<Cylinder>(VEC3(0.0,1.0,0.0),VEC3(0.0,0.0,1667),719,192.295)},
      EMC_Disk_0_EdgeOut_ { std::make_shared<Cylinder>(VEC3(0.0,1.0,0.0),VEC3(0.0,0.0,1667),719,192.295)},
      EMC_Disk_1_EdgeIn_ { std::make_shared<Cylinder>(VEC3(0.0,1.0,0.0),VEC3(0.0,0.0,3045),719,192.295)},
      EMC_Disk_1_EdgeOut_ { std::make_shared<Cylinder>(VEC3(0.0,1.0,0.0),VEC3(0.0,0.0,3045),719,192.295)},

      EMC_0_FrontIn_ { std::make_shared<Disk>(VEC3(0.0,0.0,1.0),VEC3(1.0,0.0,0.0),VEC3(0.0,0.0,-1631.12),950.)},
      EMC_0_FrontOut_ { std::make_shared<Disk>(VEC3(0.0,0.0,1.0),VEC3(1.0,0.0,0.0),VEC3(0.0,0.0,-1631.12),950.)},
      EMC_1_FrontIn_ { std::make_shared<Disk>(VEC3(0.0,0.0,1.0),VEC3(1.0,0.0,0.0),VEC3(0.0,0.0,-1631.12),950.)},
      EMC_1_FrontOut_ { std::make_shared<Disk>(VEC3(0.0,0.0,1.0),VEC3(1.0,0.0,0.0),VEC3(0.0,0.0,-1631.12),950.)},

      EMC_2_FrontIn_ { std::make_shared<Disk>(VEC3(0.0,0.0,1.0),VEC3(1.0,0.0,0.0),VEC3(0.0,0.0,-1631.12),950.)},
      EMC_2_FrontOut_ { std::make_shared<Disk>(VEC3(0.0,0.0,1.0),VEC3(1.0,0.0,0.0),VEC3(0.0,0.0,-1631.12),950.)},
      EMC_3_FrontIn_ { std::make_shared<Disk>(VEC3(0.0,0.0,1.0),VEC3(1.0,0.0,0.0),VEC3(0.0,0.0,-1631.12),950.)},
      EMC_3_FrontOut_ { std::make_shared<Disk>(VEC3(0.0,0.0,1.0),VEC3(1.0,0.0,0.0),VEC3(0.0,0.0,-1631.12),950.)},

    {}
  }
}
