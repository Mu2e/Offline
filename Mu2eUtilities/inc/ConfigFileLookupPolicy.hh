#ifndef x
#define x

#include "cetlib/filepath_maker.h"

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

#endif
