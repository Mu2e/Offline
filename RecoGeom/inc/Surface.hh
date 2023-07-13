//
//  Abstract interface to a 2-dimensional surface in 3-space
//  These are used in reconstruction to describe materials and
//  common reference locations
//  original author: David Brown (LBN) 2023
//
#ifndef RecoGeom_Surface_hh
#define RecoGeom_Surface_hh
#include "Offline/RecoGeom/inc/Ray.hh"
#include "Offline/RecoGeom/inc/IntersectFlag.hh"
#include "Math/Vector3D.h"
using ROOT::Math::XYZVectorD;
namespace mu2e {
  namespace RecoGeom {
    class Surface {
      public:
        // determine if a point is on the surface
        virtual bool onSurface(XYZVectorD const& point, double tol=1e-8) const = 0;
        // find the distance along a ray where it would intersect this surface; Returned flag describes what happened
        virtual IntersectFlag intersect(Ray const& ray,double& dist, double tol=1e-12) const = 0;
        // find the normal to the surface at the given point.  Direction convention is surface-dependent
        virtual XYZVectorD normal(XYZVectorD const& point) const  = 0;
    };
  }
}
#endif
