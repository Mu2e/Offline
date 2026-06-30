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
#include "Offline/DataProducts/inc/SurfaceId.hh"
#include <memory>
#include <string>
#include <vector>
namespace mu2e {
  namespace KKGeom {
    class DetectorSolenoid {
      public:
        using CylPtr = std::shared_ptr<KinKal::Cylinder>;
        using FruPtr = std::shared_ptr<KinKal::Frustrum>;
        using DiskPtr = std::shared_ptr<KinKal::Disk>;
        using AnnPtr = std::shared_ptr<KinKal::Annulus>;
        struct MaterialCylinder {
          SurfaceId sid_;
          CylPtr surface_;
          std::string material_;
          double thickness_;
          MaterialCylinder(SurfaceId const& sid, CylPtr const& surface, std::string const& material, double thickness) :
              sid_(sid), surface_(surface), material_(material), thickness_(thickness) {}
        };
        using MaterialCylinderCollection = std::vector<MaterialCylinder>;
        // default constructor with nominal geometry
        DetectorSolenoid( CylPtr inner, CylPtr outer, DiskPtr front, DiskPtr back,
            CylPtr ipa, DiskPtr ipafront, DiskPtr ipaback,
            FruPtr opa, AnnPtr tsda, MaterialCylinderCollection const& materialCylinders = MaterialCylinderCollection()) :
          inner_(inner) , outer_(outer), front_(front), back_(back),
          ipa_(ipa), ipa_front_(ipafront), ipa_back_(ipaback), opa_(opa), tsda_(tsda),
          materialCylinders_(materialCylinders)
      {}

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
        auto const& materialCylinders() const { return materialCylinders_; }
      private:
        CylPtr inner_; //  inner cryostat cylinder boundary
        CylPtr outer_; // outer cryostat cylinder boundary
        DiskPtr front_; // front (upstream) and back (downstream) of DS
        DiskPtr back_;
        CylPtr ipa_; // inner proton absorber
        DiskPtr ipa_front_; // front (upstream) and back (downstream) of IPA
        DiskPtr ipa_back_;
        FruPtr opa_; // outer proton absorber
        AnnPtr tsda_; // TS downstream absorber
        MaterialCylinderCollection materialCylinders_; // passive cylindrical material shells
    };
  }
}

#endif
