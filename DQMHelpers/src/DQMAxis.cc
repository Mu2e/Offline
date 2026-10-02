// One histogram axis, as a constant.
//
// Original Author: R. Mina

#include "Offline/DQMHelpers/inc/DQMAxis.hh"

#include <sstream>

namespace mu2e {

std::string DQMAxis::describe() const
{
  std::ostringstream out;
  out << n << "," << lo << "," << hi;
  return out.str();
}

} // namespace mu2e
