#ifndef REQUIREUNIQUEKEY_HH
#define REQUIREUNIQUEKEY_HH
//
// Given a list of keys and a SimpleConfig object, count how
// many of the keys have a value of true. Throw if more than 
// one is true.  Optionally, throw if none are true.
//
// $Id: requireUniqueKey.hh,v 1.1 2010/03/15 21:27:24 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/03/15 21:27:24 $
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

#endif
