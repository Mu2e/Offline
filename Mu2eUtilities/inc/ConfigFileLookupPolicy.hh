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
  ConfigFileLookupPolicy();
  virtual std::string operator() (std::string const &filename);
  virtual ~ConfigFileLookupPolicy() {}

private:
  cet::search_path path_;
};

inline
std::string
mu2e::ConfigFileLookupPolicy::
operator()(std::string const &filename) {
  return path_.find_file(filename);
}

#endif /* Mu2eUtilities_ConfigFileLookupPolicy_hh */
