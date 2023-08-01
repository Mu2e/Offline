#include "Offline/KinKalGeom/inc/Tracker.hh"
namespace mu2e {
  namespace KinKalGeom {
    using KinKal::VEC3;
    // currently use hard-coded geometry.  Note: these are only comparable to the MC virtual detector positions
    // at the ~10 um level, as the G4 virtual detector steppoints are random by roughly that amount (step tolerance)
    Tracker::Tracker() :
      // cylinders are defined by TT_outer (_inner) virtual detectors
      // Disks are defined to match TT_ent (mid, exit) virtual detectors
      outercyl_(VEC3(0.0,0.0,1.0),VEC3(0.0,0.0,4.0),850.1,1635.11),
      innercyl_(VEC3(0.0,0.0,1.0),VEC3(0.0,0.0,4.0),376.9,1635.11),
      ent_(VEC3(0.0,0.0,1.0),VEC3(1.0,0.0,0.0),VEC3(0.0,0.0,-1631.11),850.1),
      mid_(VEC3(0.0,0.0,1.0),VEC3(1.0,0.0,0.0),VEC3(0.0,0.0,10.1),850.1),
      exit_(VEC3(0.0,0.0,1.0),VEC3(1.0,0.0,0.0),VEC3(0.0,0.0,1639.11),850.1)
    {}
  }
}
