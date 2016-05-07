#ifndef Mu2eUtilities_TrackPatRecType_hh
#define Mu2eUtilities_TrackPatRecType_hh
//
// Enum-matched-to-String class for defining which
// track fit code produced a track.
//
// Contact person, Rob Kutschke
//

#include "GeneralUtilities/inc/EnumToStringSparse.hh"

namespace mu2e {

  class TrackPatRecTypeDetail{
  public:

    enum enum_type { unknown, TrkRecFit, CalPatRec };

    static std::string const& typeName();

    static std::map<enum_type,std::string> const& names();

  };

  typedef EnumToStringSparse<TrackPatRecTypeDetail> TrackPatRecType;
}

#endif
