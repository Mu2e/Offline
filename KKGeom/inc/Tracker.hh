//
//  Define the nominal tracker boundary and reference surfaces.  These are used to extrapolate and sample track fits, and to build
//  passive material models
//  original author: David Brown (LBN) 2023
//
#ifndef RecoGeom_Tracker_hh
#define RecoGeom_Tracker_hh
#include "KinKal/Geometry/inc/Cylinder.hh"
#include "KinKal/Geometry/inc/Disk.hh"
#include "KinKal/Geometry/inc/Annulus.hh"
#include <exception>
namespace mu2e {
  namespace RecoGeom {
    class Tracker {
      public:
        // construction parameters should come from a database or the geometry service
        explicit Tracker(XYZVectorD const& axis, XYZVectorD const& center, double halflen,
            double activeradius, double outradius) :
        activevol_(axis,center,halflen,activeradius),
        totalvol_(axis,center,halflen,outradius),
        ent_(axis,XYZvectorD(center.X(),center.Y(),center.Z()-halflen),activeradius),
        mid_(axis,center,activeradius),
        xit_(axis,XYZvectorD(center.X(),center.Y(),center.Z()+halflen),activeradius),
      {}
        // accessors
        // 'entrance' refers to upstream (negative z) end, 'middle' to the middle of the tracker, 'exit' refers to downstream (positive z) end
        auto const& activeVolume() const { return activevol_; }
        auto const& totalVolume() const { return totalvol_; }
        auto const& entrancePlane() const { return ent_; }
        auto const& mid() const { return mid_; }
        auto const& exitActiveDisk() const { return xit_; }

      private:
        Cylinder activevol_; // active volume cylinder
        Cylinder totalvol_; // total volume cylinder (including manifolds)
        Disk ent_, mid_, xit_; // standard reference point planes
        std::vector<Annulus> manifolds_; // manifolds
    };
  }
}

#endif
