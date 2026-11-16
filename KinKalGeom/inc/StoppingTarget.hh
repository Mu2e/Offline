//
//  Define the foils and bounding surfaces of the target
//  original author: David Brown (LBNL) 2023
//
#ifndef KinKalGeom_StoppingTarget_hh
#define KinKalGeom_StoppingTarget_hh
#include "KinKal/Geometry/Cylinder.hh"
#include "KinKal/Geometry/Disk.hh"
#include "KinKal/Geometry/Annulus.hh"
#include "cetlib_except/exception.h"
#include <memory>
#include <vector>
namespace mu2e {
  namespace KKGeom {
    class StoppingTarget {
      public:
        using CylPtr = std::shared_ptr<KinKal::Cylinder>;
        using DiskPtr = std::shared_ptr<KinKal::Disk>;
        using AnnPtr = std::shared_ptr<KinKal::Annulus>;
        using FoilCol = std::vector<AnnPtr>;
        // built by KinKalGeomMaker from the stopping target served by GeometryService.
        // inner is null when the hole of some foil reaches the axis (a foil without a hole): the target then has no inner boundary.
        // thicknesses are the full foil thicknesses (mm), one per foil
        StoppingTarget(CylPtr outer, CylPtr inner, DiskPtr front, DiskPtr back, FoilCol foils, std::vector<double> thicknesses) :
        outer_(outer), inner_(inner), front_(front), back_(back), foils_(foils), thicknesses_(thicknesses) {
          if(thicknesses_.size() != foils_.size())
            throw cet::exception("GEOM") << "KKGeom::StoppingTarget: " << foils_.size() << " foils but " << thicknesses_.size() << " thicknesses";
        }

        // accessors
        // return by reference
        auto const& outer() const { return *outer_; }
        bool hasInner() const { return static_cast<bool>(inner_); }
        auto const& inner() const {
          if(!inner_) throw cet::exception("GEOM") << "KKGeom::StoppingTarget: no inner boundary, a foil has no hole";
          return *inner_;
        }
        auto const& front() const { return *front_; }
        auto const& back() const { return *back_; }
        auto const& foil(size_t ifoil) const { return *foils_[ifoil]; }
        double foilThickness(size_t ifoil) const { return thicknesses_.at(ifoil); }
        // return by ptr
        auto const& outerPtr() const { return outer_; }
        auto const& innerPtr() const { return inner_; } // null when hasInner() is false
        auto const& frontPtr() const { return front_; }
        auto const& backPtr() const { return back_; }
        auto const& foils() const { return foils_; }
        auto const& foilPtr(size_t ifoil) const { return foils_[ifoil]; }
        auto const& foilThicknesses() const { return thicknesses_; }
      private:
        CylPtr outer_, inner_; // boundaries
        DiskPtr front_, back_;
        FoilCol foils_; // target foils
        std::vector<double> thicknesses_; // full thickness of each foil (mm)
    };
  }
}

#endif
