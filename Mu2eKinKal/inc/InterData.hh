//
//  Data payload for the intersection point of a PiecewiseTrajectory with a surface
//  original author: David Brown (LBNL) 2023
//
#ifndef Mu2eKinKal_InterData_hh
#define Mu2eKinKal_InterData_hh
#include "Offline/DataProducts/inc/GenVector.hh"
namespace mu2e {
  namespace Mu2eKinKal {
    struct InterData {
      enum Status{unknown=-1,nointersection=0,good,outsidetime,outsideboundary};
      Status status_ = unknown;
      XYZVectorD inter_; // intersection point
      XYZVectorD norm_; // surface normal at intersection
      XYZVectorD pdir_; // particle direction at intersection
      double time_; // particle time at intersection
    };
  }
}
#endif
