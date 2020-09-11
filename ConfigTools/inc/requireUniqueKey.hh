#ifndef ConfigTools_requireUniqueKey_hh
#define ConfigTools_requireUniqueKey_hh
//
// Given a list of keys and a SimpleConfig object, count how
// many of the keys have a value of true. Throw if more than
// one is true.  Optionally, throw if none are true.
//
//
// Contact person Rob Kutschke
//

#include <vector>
#include <string>

namespace mu2e {

  class SimpleConfig;

  int requireUniqueKey ( const std::vector<std::string>& keys,
                         const SimpleConfig&             config,
                         bool  throwOnZero=false );

} // namespace mu2e

#endif /* ConfigTools_requireUniqueKey_hh */
