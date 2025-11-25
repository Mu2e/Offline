#include "Offline/KinKalGeom/inc/Calo.hh"
namespace mu2e {
  namespace KinKalGeom {
    using KinKal::VEC3;
    using KinKal::Cylinder;
    using KinKal::Disk;
    // currently use hard-coded geometry.  Note: these are only comparable to the MC virtual detector positions
    // at the ~10 um level, as the G4 virtual detector steppoints are random by roughly that amount (step tolerance)
    //FIXME need to update values for the Calorimeter
    Calo::Calo() :
      d0_outer_ { std::make_shared<Cylinder>(VEC3(0.0,0.0,1.0),VEC3(0.0,0.0,4.0),850.11,1635.11)},
      d0_inner_{ std::make_shared<Cylinder>(VEC3(0.0,0.0,1.0),VEC3(0.0,0.0,4.0),376.9,1635.11)},
      d0_front_{ std::make_shared<Disk>(VEC3(0.0,0.0,1.0),VEC3(1.0,0.0,0.0),VEC3(0.0,0.0,-1631.12),950.)},
      d0_back_{ std::make_shared<Disk>(VEC3(0.0,0.0,1.0),VEC3(1.0,0.0,0.0),VEC3(0.0,0.0,1639.10),950.)}

      d1_outer_ { std::make_shared<Cylinder>(VEC3(0.0,0.0,1.0),VEC3(0.0,0.0,4.0),850.11,1635.11)},
      d1_inner_{ std::make_shared<Cylinder>(VEC3(0.0,0.0,1.0),VEC3(0.0,0.0,4.0),376.9,1635.11)},
      d1_front_{ std::make_shared<Disk>(VEC3(0.0,0.0,1.0),VEC3(1.0,0.0,0.0),VEC3(0.0,0.0,-1631.12),950.)},
      d1_back_{ std::make_shared<Disk>(VEC3(0.0,0.0,1.0),VEC3(1.0,0.0,0.0),VEC3(0.0,0.0,1639.10),950.)}
    {}
  }
}
