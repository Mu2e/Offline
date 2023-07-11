//
// A special case for the class PDGCode, return a vector of the enum_type,
// of not a vector of the class type. Needed for backwards compatibiltiy.
//
// Original Author Rob Kutschke
//

#include "Offline/Mu2eUtilities/inc/fromStrings.hh"

#include <string>
#include <vector>

namespace mu2e {

  std::vector<mu2e::PDGCode::type> fromStrings_type( std::vector<std::string> const& v ){
    std::vector<mu2e::PDGCode::type> r;
    for ( auto const& s : v){
      r.emplace_back( PDGCode(s) ); // PDGCode has implicit conversion to enum_type
    }
    return r;
  }

} // end namespace mu2e
