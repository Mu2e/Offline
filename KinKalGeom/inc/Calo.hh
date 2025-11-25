//
//  Define the nominal calorimeter boundary and reference surfaces, used to extrapolate and sample KinKal track fits, and to build
//  the passive materials in the fit
//  original author: Sophie Middleton (2025)
//
#ifndef KinKalGeom_Calorimeter_hh
#define KinKalGeom_Calorimeter_hh
#include "KinKal/Geometry/Cylinder.hh"
#include "KinKal/Geometry/Disk.hh"
#include "KinKal/Geometry/Annulus.hh"
#include <memory>
namespace mu2e {
  namespace KinKalGeom {
    class Calorimeter {
      public:
        using CylPtr = std::shared_ptr<KinKal::Cylinder>;
        using DiskPtr = std::shared_ptr<KinKal::Disk>;
        // default constructor with nominal geometry
        Calorimeter();
        // accessors - d0 - disk 0, d1 - disk 1
        // return by reference
        auto const& d0_outer() const { return *d0_outer_; }
        auto const& d0_inner() const { return *d0_inner_; }
        auto const& d0_front() const { return *d0_front_; }
        auto const& d0_back() const  { return *d0_back_;  }

        auto const& d1_outer() const { return *d1_outer_; }
        auto const& d1_inner() const { return *d1_inner_; }
        auto const& d1_front() const { return *d1_front_; }
        auto const& d1_back() const  { return *d1_back_;  }

        // return by ptr
        auto const& d0outerPtr() const { return d0_outer_; }
        auto const& d0innerPtr() const { return d0_inner_; }
        auto const& d0frontPtr() const { return d0_front_; }
        auto const& d0backPtr() const  { return d0_back_;  }

        auto const& d1outerPtr() const { return d1_outer_; }
        auto const& d1innerPtr() const { return d1_inner_; }
        auto const& d1frontPtr() const { return d1_front_; }
        auto const& d1backPtr() const  { return d1_back_;  }

      private:
        CylPtr d0_outer_, d0_inner_, d1_outer_, d1_inner_; // active volume boundary
        DiskPtr d0_front_, d0_back_, d1_front_, d1_back_; // disk front and back
    };
  }
}

#endif
