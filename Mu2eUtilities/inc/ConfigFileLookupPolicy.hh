#ifndef Mu2eUtilities_ConfigFileLookupPolicy_hh
#define Mu2eUtilities_ConfigFileLookupPolicy_hh

#include "cetlib/filepath_maker.h"

// $MU2E_SEARCH_PATH contains the colon-separated search path

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

#endif /* Mu2eUtilities_ConfigFileLookupPolicy_hh */
