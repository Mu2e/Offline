#ifndef Mu2eUtilities_requireUniqueKey_hh
#define Mu2eUtilities_requireUniqueKey_hh
//
// Given a list of keys and a SimpleConfig object, count how
// many of the keys have a value of true. Throw if more than 
// one is true.  Optionally, throw if none are true.
//
// $Id: requireUniqueKey.hh,v 1.3 2011/05/17 15:41:36 greenc Exp $
// $Author: greenc $
// $Date: 2011/05/17 15:41:36 $
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
