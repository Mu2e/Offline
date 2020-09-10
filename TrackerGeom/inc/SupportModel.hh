#ifndef TrackerGeom_SupportModel_hh
#define TrackerGeom_SupportModel_hh
//
// Enum-matched-to-String class for defining which version of the Tracker
// support structure is to be used.
//
//
// Contact person, Rob Kutschke
//

#include "GeneralUtilities/inc/EnumToStringSparse.hh"

namespace mu2e {

  class SupportModelDetail{
  public:

    enum enum_type { unknown, simple, detailedv0 };

    static std::string const& typeName();

    static std::map<enum_type,std::string> const& names();

  };

  typedef EnumToStringSparse<SupportModelDetail> SupportModel;
}

#endif
