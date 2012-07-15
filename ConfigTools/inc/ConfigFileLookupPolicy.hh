#ifndef ConfigTools_ConfigFileLookupPolicy_hh
#define ConfigTools_ConfigFileLookupPolicy_hh
//
// Complete a relative path by grafting on a path-prefix from
// one of the elements in the environment variable MU2E_SEARCH_PATH,
// which contains a colon separated list of directory names.
//
// $Id: ConfigFileLookupPolicy.hh,v 1.1 2012/07/15 22:00:35 kutschke Exp $
// $Author: kutschke $
// $Date: 2012/07/15 22:00:35 $
//
// Contact person Rob Kutschke

#include "cetlib/filepath_maker.h"

namespace mu2e {
  class ConfigFileLookupPolicy;
}

class mu2e::ConfigFileLookupPolicy :
  public cet::filepath_maker {
public:
  virtual std::string operator() (std::string const &filename);
  virtual ~ConfigFileLookupPolicy() noexcept {}

private:
  inline cet::search_path const &path() const;
};

inline
cet::search_path const &
mu2e::ConfigFileLookupPolicy::
path() const {
  // path_s is function-local-static rather than class-scope so that all
  // instances within an execution have the same search path even if
  // MU2E_SEARCH_PATH gets tweaked within a job -- unlikely, but
  // possible.
  static cet::search_path const path_s("MU2E_SEARCH_PATH");
  return path_s;
}

inline
std::string
mu2e::ConfigFileLookupPolicy::
operator()(std::string const &filename) {
  return path().find_file(filename);
}

#endif /* ConfigTools_ConfigFileLookupPolicy_hh */
