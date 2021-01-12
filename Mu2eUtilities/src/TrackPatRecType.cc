//
// Enum-matched-to-String class for defining which
// pattern recognition code found a track.
//
// Per Chris Jones
//  1) static initializing nam inside the static member function
//     ensures that there will not be order of initialization problems.
//  2) as written this, improves multi-threading performance since
//     initialization of a const object is done at .so load time.
//     The symbol nam is not visible outside of this compilation unit
//     so it can only be accessed via that the static member function name().
//  3) An even better solution would be to replace std::map with
//     two std::array<>s; then the static initialization can be
//     made constexpr. It's likely that this can be done inside
//     the static member fucntion and still retain the good thread
//     performance.  So we have the best of both worlds.  But he
//     has not looked at any generated code to be sure.
//     There is another problem with the two std::array solution:
//     ensuring order so that a binary search can be used for lookup.
//     The map ensures this.
//  4) For enum to string classes, we only invoke the names() function
//     in conjunction with IO. So thread efficiency is probably not
//     a critical issue.  It is also very unlikely to trip any order of
//     initialization problems - this is an attribute of how the
//     class is used - but it cannot be guaranteed.
//
// Contact person, Rob Kutschke
//

#include <type_traits>
#include <utility>

#include "Mu2eUtilities/inc/TrackPatRecType.hh"

namespace mu2e {

  std::string const& TrackPatRecTypeDetail::typeName() {
    static const std::string type("TrackPatRecType");
    return type;
  }

  static const std::map<TrackPatRecType::enum_type,std::string> nam{
    std::make_pair(TrackPatRecType::unknown,   "unknown"  ),
    std::make_pair(TrackPatRecType::TrkRecFit, "TrkRecFit"),
    std::make_pair(TrackPatRecType::CalPatRec, "CalPatRec")
  };

  std::map<TrackPatRecType::enum_type,std::string> const& TrackPatRecTypeDetail::names(){
    return nam;
  }

}
