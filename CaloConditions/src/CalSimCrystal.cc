#include "Offline/CaloConditions/inc/CalSimCrystal.hh"
#include "cetlib_except/exception.h"
#include <iomanip>

namespace mu2e {

  float CalSimCrystal::LRU(const CrystalId& Id) const {
    if (!Id.isValid()) {
      throw cet::exception("CalSimCrystal_RANGE")
          << "CalSimCrystal::id invalid input" << Id <<  "\n";
    }
    return _LRUs[Id.id()];
  }

  const std::vector<float> CalSimCrystal::pePerMeVs(const CrystalId& Id) const {
    if (!Id.isValid()) {
      throw cet::exception("CalSimCrystal_RANGE")
          << "CalSimCrystal::id invalid input" << Id << "\n";
    }
    return _pePerMeVs[Id.id()];
  }

  void CalSimCrystal::print(std::ostream& os) const {
    os << "\n";
    os << "\n" << "CalSimCrystal : \n" ;
    os << "\n" << " Crystalid  LRU   pePerMeV_0   pePerMeV_1" ;
    for (uint16_t i = 0; i < CaloConst::_nCrystal; i++) {
      os << std::setw(6) << i << std::setw(6) << _LRUs[i] << std::setw(6) << _pePerMeVs[i].at(0) << std::setw(6) << _pePerMeVs[i].at(1) <<"\n";
    }
    os << "\n";
  }
}
