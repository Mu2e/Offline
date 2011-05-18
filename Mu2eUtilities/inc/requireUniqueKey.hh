#ifndef Mu2eUtilities_requireUniqueKey_hh
#define Mu2eUtilities_requireUniqueKey_hh
//
// Given a list of keys and a SimpleConfig object, count how
// many of the keys have a value of true. Throw if more than
// one is true.  Optionally, throw if none are true.
//
// $Id: requireUniqueKey.hh,v 1.4 2011/05/18 02:27:18 wb Exp $
// $Author: wb $
// $Date: 2011/05/18 02:27:18 $
//
// Original author Rob Kutschke
//

#include <vector>
#include <string>

namespace mu2e {

  class SimpleConfig;

  int requireUniqueKey ( const std::vector<std::string>& keys,
                         const SimpleConfig&             config,
                         bool  throwOnZero=false );

} // namespace mu2e

#endif /* Mu2eUtilities_requireUniqueKey_hh */
