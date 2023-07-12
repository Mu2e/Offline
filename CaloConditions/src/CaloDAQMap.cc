#include "Offline/CaloConditions/inc/CaloDAQMap.hh"
#include "cetlib_except/exception.h"
#include <iomanip>

using namespace std;
namespace mu2e {

  CaloSiPMId CaloDAQMap::offlineId(CaloRawSiPMId rawId) const {
    if(!rawId.isValid()) {
      throw cet::exception("CALODAQMP_RANGE")
        << "CaloDAQMap::offlineId invalid input rawId" << rawId << std::endl;
    }
    CaloSiPMId offId = _raw2Offline[rawId.id()];
    if(!offId.isValid()) {
      throw cet::exception("CALODAQMP_RANGE")
        << "CaloDAQMap::offlineId no offline id for rawId " << rawId << std::endl;
    }
    return offId;
  }

  CaloRawSiPMId CaloDAQMap::rawId(CaloSiPMId offId) const {
    if(!offId.isValid()) {
      throw cet::exception("CALODAQMP_RANGE") << "CaloDAQMap::rawId invalid input offlineId" << offId << std::endl;
    }
    return _offline2Raw[offId.id()];
  }

  void CaloDAQMap::print(std::ostream& os) const {
    os << endl;
    os << endl << "CaloDAQMap : "  << endl;
    os << endl << " offline  raw "  << endl;
    for(uint16_t i=0; i<CaloConst::_nRawChannel;i++) {
      os << setw(6) << _raw2Offline[i].id() << setw(6) << i << endl;
    }
    os << endl;
  }
}
