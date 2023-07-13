//
//  Data payload for the intersection point of a PiecewiseTrajectory with a surface
//  original author: David Brown (LBNL) 2023
//
#ifndef RecoGeom_InterData_hh
#define RecoGeom_InterData_hh
#include "Offline/DataProducts/inc/GenVector.hh"
#include "Offline/RecoGeom/inc/IntersectFlag.hh"

namespace mu2e {
  namespace RecoGeom {
    struct InterData {
      IntersectFlag flag_; // intersection status
      XYZVectorD pos_; // intersection position
      XYZVectorD norm_; // surface normal at intersection
      XYZVectorD pdir_; // particle direction at intersection
      double time_; // time at intersection (from particle)
    };
  }
}
#endif
