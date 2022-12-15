#ifndef Mu2eUtilities_fromStrings_hh
#define Mu2eUtilities_fromStrings_hh
//
// A free function template that translates from std::vector<string>
// to a std::vector of an enum-matched-to-string type specified as the template type argument.
//
// An additional function for PDGCode needed temporarily for backward compatibility.
//
// Original Author Rob Kutschke
//

#include "Offline/DataProducts/inc/PDGCode.hh"

#include <string>
#include <vector>

namespace mu2e {

  // For any type that is constructable from an std::string,
  // take an std::vector of strings and return a std::vector of that type.
  template<typename T>
  std::vector<T> fromStrings( std::vector<std::string> const& vs ){
    std::vector<T> r;
    r.reserve(vs.size());
    for ( auto const& s : vs ){
      r.emplace_back(s);
    }
    return r;
  }

  // A special case for the class PDGCode, return the enum_type, not the class tupe.
  // Needed for backwards compatibiltiy.
  std::vector<mu2e::PDGCode::type> fromStrings_type( std::vector<std::string> const& v );


} // end namespace mu2e

#endif
