//
//  Define the bounding surfaces and passive materials in the detector solenoid
//  original author: David Brown (LBN) 2023
//
#ifndef KinKalGeom_DetectorSolenoid_hh
#define KinKalGeom_DetectorSolenoid_hh
#include "KinKal/Geometry/Cylinder.hh"
#include "KinKal/Geometry/Frustrum.hh"
#include "KinKal/Geometry/Annulus.hh"
#include <vector>
namespace mu2e {
  namespace KinKalGeom {
    class DetectorSolenoid {
      public:
        // default constructor with nominal geometry
        DetectorSolenoid();
        // accessors
        auto const& outerCylinder() const { return outercyl_; }
        auto const& innerCylinder() const { return innercyl_; }
        auto const& innerProtonAbsorber() const { return ipa_; }
        auto const& outerProtonAbsorber() const { return ipa_; }
        auto const& upstreamAbsorber() const { return tsda_; }
      private:
        KinKal::Cylinder outercyl_; // outer cryostat cylinder
        KinKal::Cylinder innercyl_; //  inner cryostat cylinder
        KinKal::Cylinder ipa_; // inner proton absorber
        KinKal::Frustrum opa_; // outer proton absorber
        KinKal::Annulus tsda_; // TS downstream absorber
    };
  }
}

#endif
