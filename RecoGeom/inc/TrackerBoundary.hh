//
//  Define the surfaces which define the nominal tracker boundary.  It also includes the midplane.  These Surfaces are used to extrapolate and sample track fits
//  original author: David Brown (LBN) 2023
//
#ifndef RecoGeom_TrackerBoundary_hh
#define RecoGeom_TrackerBoundary_hh
#include "Offline/RecoGeom/inc/Cylinder.hh"
#include "Offline/RecoGeom/inc/Annulus.hh"
#include <exception>
namespace mu2e {
  namespace RecoGeom {
    class TrackerBoundary {
      public:
        // construction parameters should come from a database or the geometry service
        explicit TrackerBoundary(XYZVectorD const& axis, XYZVectorD const& center, double halflen,
            double inradius, double outradius) :
        incyl_(axis,center,halflen,inradius),
        outcyl_(axis,center,halflen,outradius),
        entdisk_(axis,XYZvectorD(center.X(),center.Y(),center.Z()-halflen),0.0,inradius),
        middisk_(axis,center,0.0,inradius),
        xitdisk_(axis,XYZvectorD(center.X(),center.Y(),center.Z()+halflen),0.0,inradius),
        entann_(axis,XYZvectorD(center.X(),center.Y(),center.Z()-halflen),inradius,outradius),
        midann_(axis,center,inradius,outradius),
        xitann_(axis,XYZvectorD(center.X(),center.Y(),center.Z()+halflen),inradius,outradius)
      {}
        // accessors
        // 'entrance' refers to upstream (negative z) end, 'middle' to the middle of the tracker, 'exit' refers to downstream (positive z) end
        auto const& activeCylinder() const { return incyl_; }
        auto const& manifoldCylinder() const { return outcyl_; }
        auto const& entranceActiveDisk() const { return entdisk_; }
        auto const& middleActiveDisk() const { return middisk_; }
        auto const& exitActiveDisk() const { return xitdisk_; }
        auto const& entranceManifoldAnnulus() const { return entann_; }
        auto const& middleManifoldAnnulus() const { return midann_; }
        auto const& exitManifoldAnnulus() const { return xitann_; }

      private:
        Cylinder incyl_; // active volume cylinder
        Cylinder outcyl_; // outer manifold volume cylinder
        Annulus entdisk_,middisk_,xitdisk_; // active volume disks
        Annulus entann_,midann_,xitann_; // manifold volume annuli
    };
  }
}

#endif
