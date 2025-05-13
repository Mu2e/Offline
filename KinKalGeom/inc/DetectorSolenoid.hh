//
//  Define the bounding surfaces and passive materials in the detector solenoid
//  original author: David Brown (LBNL) 2023
//
#ifndef KinKalGeom_DetectorSolenoid_hh
#define KinKalGeom_DetectorSolenoid_hh
#include "KinKal/Geometry/Cylinder.hh"
#include "KinKal/Geometry/Frustrum.hh"
#include "KinKal/Geometry/Annulus.hh"
#include "KinKal/Geometry/Disk.hh"
#include <memory>
#include <vector>
namespace mu2e {
  namespace KinKalGeom {
    class DetectorSolenoid {
      public:
        using CylPtr = std::shared_ptr<KinKal::Cylinder>;
        using FruPtr = std::shared_ptr<KinKal::Frustrum>;
        using DiskPtr = std::shared_ptr<KinKal::Disk>;
        using AnnPtr = std::shared_ptr<KinKal::Annulus>;
        // default constructor with nominal geometry
        DetectorSolenoid();
        // accessors
        // return by reference
        auto const& outer() const { return *outer_; }
        auto const& inner() const { return *inner_; }
        auto const& front() const { return *front_; }
        auto const& back() const { return *back_; }
        auto const& innerProtonAbsorber() const { return *ipa_; }
        auto const& innerProtonAbsorberFront() const { return *ipa_front_; }
        auto const& innerProtonAbsorberBack() const { return *ipa_back_; }
        auto const& outerProtonAbsorber() const { return *opa_; }
        auto const& upstreamAbsorber() const { return *tsda_; }
        // return by ptr
        auto const& outerPtr() const { return outer_; }
        auto const& innerPtr() const { return inner_; }
        auto const& frontPtr() const { return front_; }
        auto const& backPtr() const { return back_; }
        auto const& innerProtonAbsorberPtr() const { return ipa_; }
        auto const& innerProtonAbsorberFrontPtr() const { return ipa_front_; }
        auto const& innerProtonAbsorberBackPtr() const { return ipa_back_; }
        auto const& outerProtonAbsorberPtr() const { return opa_; }
        auto const& upstreamAbsorberPtr() const { return tsda_; }
      private:
        CylPtr inner_; //  inner cryostat cylinder
        CylPtr outer_; // outer cryostat cylinder
        DiskPtr front_; // front (upstream) and back (downstream) of DS
        DiskPtr back_;
        CylPtr ipa_; // inner proton absorber
        DiskPtr ipa_front_; // front (upstream) and back (downstream) of IPA
        DiskPtr ipa_back_;
        FruPtr opa_; // outer proton absorber
        AnnPtr tsda_; // TS downstream absorber
    };
  }
}

#endif
