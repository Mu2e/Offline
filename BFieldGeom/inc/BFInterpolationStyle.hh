#ifndef BFieldGeom_BFInterpolationStyle_hh
#define BFieldGeom_BFInterpolationStyle_hh
//
// Enum-matched-to-String class for defining how to do interpolation
// within a magnetic field map.
//
//
// Contact person, Rob Kutschke
//

#include "Offline/GeneralUtilities/inc/EnumToStringSparse.hh"

namespace mu2e {

    class BFInterpolationStyleDetail {
       public:
        enum enum_type { unknown, unused, trilinear, fit };

        static std::string const& typeName();

        static std::map<enum_type, std::string> const& names();
    };

    typedef EnumToStringSparse<BFInterpolationStyleDetail> BFInterpolationStyle;
}  // namespace mu2e

#endif
