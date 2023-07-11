//
//  Specialization of the intersection template for a KKLine with a cylinder
//  original author: David Brown (LBNL) 2023
//
#include "Offline/Mu2eKinKal/inc/Intersection.hh"
#include "Offline/Mu2eKinKal/inc/Cylinder.hh"
#include "KinKal/Trajectory/LoopHelix.hh"
#include "KinKal/General/Parameters.hh"
#include "KinKal/General/Vectors.hh"
using KTRAJ= KinKal::LoopHelix;
using SURFACE = mu2e::Mu2eKinKal::Cylinder;
namespace mu2e {
  namespace Mu2eKinKal {
    template<> Intersection<KTRAJ,SURFACE> {
      Intersection(KTRAJ const& ktraj, SURFACE const& surface) : ktraj_(ktraj), surf_(surface) {
      // temporary: do nothing
      std::cout << "trajectory " << ktraj << std::endl;
    }
    };


//    KTRAJ lhelix;
//    XYZVectorD axis(0.0,0.0,1.0), center;
//    SURFACE cyl(axis, center, 1.0, 1.0);
//    Intersection inter(lhelix,cyl);

  }
}
