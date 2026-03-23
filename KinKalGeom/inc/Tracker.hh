//
//  Define the nominal tracker boundary and reference surfaces, used to extrapolate and sample KinKal track fits, and to build
//  the passive materials in the fit
//  original author: David Brown (LBNL) 2023
//
#ifndef KinKalGeom_Tracker_hh
#define KinKalGeom_Tracker_hh
#include "KinKal/Geometry/Cylinder.hh"
#include "KinKal/Geometry/Disk.hh"
#include "KinKal/Geometry/Annulus.hh"
#include <memory>
namespace mu2e {
  namespace KinKalGeom {
    class Tracker {
      public:
        using CylPtr = std::shared_ptr<KinKal::Cylinder>;
        using DiskPtr = std::shared_ptr<KinKal::Disk>;
        // default constructor with nominal geometry
        Tracker();
        // accessors
        // return by reference
        auto const& outer() const { check_init(); return *outer_; }
        auto const& inner() const { check_init(); return *inner_; }
        auto const& front() const { check_init(); return *front_; }
        auto const& middle() const { check_init(); return *mid_; }
        auto const& back() const { check_init(); return *back_; }
        auto const& outerPtr() const { check_init(); return outer_; }
        auto const& innerPtr() const { check_init(); return inner_; }
        auto const& frontPtr() const { check_init(); return front_; }
        auto const& middlePtr() const { check_init(); return mid_; }
        auto const& backPtr() const { check_init(); return back_; }

      private:
        mutable bool initialized_ = false; // defer construction to allow services to be established
        void check_init() const;
        void initialize();
        CylPtr outer_, inner_; // active volume boundary
        DiskPtr front_, mid_, back_; // standard reference planes
        // TODO: add passive material for manifolds, absorber rings, beams, ...
    };
  }
}

#endif
