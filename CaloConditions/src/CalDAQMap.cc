#include "Offline/CaloConditions/inc/CalDAQMap.hh"
#include "cetlib_except/exception.h"
#include <iomanip>

using namespace std;
namespace mu2e {

CaloSiPMId CalDAQMap::offlineId(CaloRawSiPMId rawId) const {
  if (!rawId.isValid()) {
    throw cet::exception("CALODAQMP_RANGE")
        << "CalDAQMap::offlineId invalid input rawId" << rawId << std::endl;
  }
  return _raw2Offline[rawId.id()];
}

CaloRawSiPMId CalDAQMap::rawId(CaloSiPMId offId) const {
  if (!offId.isValid()) {
    throw cet::exception("CALODAQMP_RANGE")
        << "CalDAQMap::rawId invalid input offlineId" << offId << std::endl;
  }
  return _offline2Raw[offId.id()];
}

void CalDAQMap::print(std::ostream& os) const {
  os << endl;
  os << endl << "CalDAQMap : " << endl;
  os << endl << " offline  raw " << endl;
  for (uint16_t i = 0; i < CaloConst::_nRawChannel; i++) {
    os << setw(6) << _raw2Offline[i].id() << setw(6) << i << endl;
  }
  os << endl;
}
} // namespace mu2e
