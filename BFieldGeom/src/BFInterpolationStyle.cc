//
// Enum-matched-to-String class for defining how to do interpolation
// within a magnetic field map.
//
//
// Contact person, Rob Kutschke
//

#include "BFieldGeom/inc/BFInterpolationStyle.hh"

namespace mu2e {

    std::string const& BFInterpolationStyleDetail::typeName() {
        static std::string type("BFInterpolationStyle");
        return type;
    }

    std::map<BFInterpolationStyleDetail::enum_type, std::string> const&
    BFInterpolationStyleDetail::names() {
        static std::map<enum_type, std::string> nam;

        if (nam.empty()) {
            nam[unknown] = "unknown";
            nam[meco] = "meco";
            nam[trilinear] = "trilinear";
            nam[fit] = "fit";
        }

        return nam;
    }

}  // namespace mu2e
