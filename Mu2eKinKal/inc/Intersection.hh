//
//  Calculate the intersection point of a simple Trajectory with a surface
//  This must be specialized for every case (every pair of trajectory and surface)
//  original author: David Brown (LBNL) 2023
//
#ifndef Mu2eKinKal_Intersection_hh
#define Mu2eKinKal_Intersection_hh
#include "Offline/Mu2eKinKal/inc/InterData.hh"

namespace mu2e {
  namespace Mu2eKinKal {
    template <class TRAJ, class SURF> class Intersection {
      public:
        // constructor performs the intersction
        explicit Intersection(TRAJ const& ktraj, SURF const& surface,double timehint);
        auto const& interData() const { return data_; }
        auto const& trajectory() const { return ktraj_; }
        auto const& surface() const { return surf_; }
      private:
        InterData data_; // result of this intersction
        TRAJ const& ktraj_; // trajectory
        SURF const& surf_; // surface
    };
  }
}
#endif
