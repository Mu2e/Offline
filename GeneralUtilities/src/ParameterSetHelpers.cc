// Andrei Gaponenko, 2012

#include "GeneralUtilities/inc/ParameterSetHelpers.hh"

#include <sstream>
#include <stdexcept>

#include <vector>
#include "CLHEP/Vector/ThreeVector.h"

template<>
bool fhicl::ParameterSet::get_if_present<CLHEP::Hep3Vector>(std::string const & key, CLHEP::Hep3Vector& value) const {
  std::vector<double> val;
  const bool present = get_if_present<std::vector<double> >(key,val);
  if(present) {
    if(val.size() != 3) {
      std::ostringstream os;
      os<<"Error converting std::vector<double> to CLHEP::Hep3Vector for key \""<<key<<"\": wrong input size = "<<val.size();
      throw std::runtime_error(os.str());
    }
    value = CLHEP::Hep3Vector(val[0], val[1], val[2]);
  }
  return present;
}
