#include "DataProducts/inc/FilterEfficiency.hh"
#include <ostream>

namespace mu2e {
  std::ostream& operator<<(std::ostream& os, const FilterEfficiency& eff) {
    return os<<"FilterEfficiency(accepted="<<eff.accepted()<<", total="<<eff.total()<<")";
  }
}
