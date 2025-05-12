#include "Offline/KinKalGeom/inc/Tracker.hh"
namespace mu2e {
  namespace KinKalGeom {
    using KinKal::VEC3;
    using KinKal::Cylinder;
    using KinKal::Disk;
    // currently use hard-coded geometry.  Note: these are only comparable to the MC virtual detector positions
    // at the ~10 um level, as the G4 virtual detector steppoints are random by roughly that amount (step tolerance)
    Tracker::Tracker() :
      // cylinders are defined by TT_outer (_inner) virtual detectors
      // Disks are defined to match TT_front (mid, back) virtual detectors
      outer_ { std::make_shared<Cylinder>(VEC3(0.0,0.0,1.0),VEC3(0.0,0.0,4.0),850.11,1635.11)},
      inner_{ std::make_shared<Cylinder>(VEC3(0.0,0.0,1.0),VEC3(0.0,0.0,4.0),376.9,1635.11)},
      // expand the disk radii to meet the DS
      front_{ std::make_shared<Disk>(VEC3(0.0,0.0,1.0),VEC3(1.0,0.0,0.0),VEC3(0.0,0.0,-1631.12),950.)},
      mid_{ std::make_shared<Disk>(VEC3(0.0,0.0,1.0),VEC3(1.0,0.0,0.0),VEC3(0.0,0.0,10.09),950.)},
      back_{ std::make_shared<Disk>(VEC3(0.0,0.0,1.0),VEC3(1.0,0.0,0.0),VEC3(0.0,0.0,1639.10),950.)}
    {}
  }
}
