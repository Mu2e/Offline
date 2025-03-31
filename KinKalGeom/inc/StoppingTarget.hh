//
//  Define the foils and bounding surfaces of the target
//  original author: David Brown (LBNL) 2023
//
#ifndef KinKalGeom_StoppingTarget_hh
#define KinKalGeom_StoppingTarget_hh
#include "KinKal/Geometry/Cylinder.hh"
#include "KinKal/Geometry/Disk.hh"
#include "KinKal/Geometry/Annulus.hh"
#include <memory>
#include <vector>
namespace mu2e {
  namespace KinKalGeom {
    class StoppingTarget {
      public:
        using CylPtr = std::shared_ptr<KinKal::Cylinder>;
        using DiskPtr = std::shared_ptr<KinKal::Disk>;
        using AnnPtr = std::shared_ptr<KinKal::Annulus>;
        using FoilCol = std::vector<AnnPtr>;
        // default constructor with nominal geometry
        StoppingTarget();
        // accessors
        // return by reference
        auto const& outer() const { return *outer_; }
        auto const& inner() const { return *inner_; }
        auto const& front() const { return *front_; }
        auto const& back() const { return *back_; }
        auto const& foil(size_t ifoil) const { return *foils_[ifoil]; }
        // return by ptr
        auto const& outerPtr() const { return outer_; }
        auto const& innerPtr() const { return inner_; }
        auto const& frontPtr() const { return front_; }
        auto const& backPtr() const { return back_; }
        auto const& foils() const { return foils_; }
        auto const& foilPtr(size_t ifoil) const { return foils_[ifoil]; }
      private:
        CylPtr outer_, inner_; // boundaries
        DiskPtr front_, back_;
        FoilCol foils_; // target foils
    };
  }
}

#endif
