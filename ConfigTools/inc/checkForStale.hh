#ifndef ConfigTools_checkForStale_hh
#define ConfigTools_checkForStale_hh
//
// Given a list of stale names and a SimpleConfig object, count how
// many of the names have a value of true. Throw if more than
// one is true.  
//
// $Id: checkForStale.hh,v 1.1 2013/06/14 16:15:04 knoepfel Exp $
// $Author: knoepfel $
// $Date: 2013/06/14 16:15:04 $
//
// Contact person Rob Kutschke
//

#include <vector>
#include <string>

namespace mu2e {

  class SimpleConfig;

  // Check for stale variables using vector of names
  // - looks for exact match between stale list and config-file list
  void checkForStale ( const std::vector<std::string>& staleNames,
                       const SimpleConfig& config,
                       bool  throwOnMatch=true );

  // Check for stale variables using string
  // - uses string::find method, effectively allowing wildcards
  void checkForStale ( std::string staleString,
                       const SimpleConfig& config,
                       bool throwOnMatch=true );

} // namespace mu2e

#endif /* ConfigTools_checkForStale_hh */
