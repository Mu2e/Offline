#include "Offline/CaloConditions/inc/CalSimParams.hh"
#include "cetlib_except/exception.h"
#include <iomanip>

namespace mu2e {

  float CalSimParams::LRU(const CrystalId& Id) const {
    if (!Id.isValid()) {
      throw cet::exception("CALSIMPARAMS_RANGE")
          << "invalid LRU input id " << Id <<  "\n";
    }
    return _LRUs[Id.id()];
  }

  const std::vector<float>& CalSimParams::pePerMeVs(const CrystalId& Id) const {
    if (!Id.isValid()) {
      throw cet::exception("CALSIMPARAMS_RANGE")
          << "invalid pePerMev input id" << Id << "\n";
    }
    return _pePerMeVs[Id.id()];
  }

  const std::vector<float>& CalSimParams::ADCPerMeVs(const CrystalId& Id) const {
    if (!Id.isValid()) {
      throw cet::exception("CALSIMPARAMS_RANGE")
          << "invalid ADCPerMev input id" << Id << "\n";
    }
    return _ADCPerMeVs[Id.id()];
  }

  void CalSimParams::print(std::ostream& os) const {
    os << "\n";
    os << "\n" << "CalSimParams : \n" ;
    os << "\n" << " Crystalid  LRU   pePerMeV_0   pePerMeV_1  ADCPerMeV_0   ADCPerMeV_1 " ;
    for (uint16_t i = 0; i < CaloConst::_nCrystal; i++) {
      os << std::setw(6) << i << std::setw(6) << _LRUs[i] << std::setw(6) << _pePerMeVs[i].at(0)
         << std::setw(6) << _pePerMeVs[i].at(1) << std::setw(6) << _ADCPerMeVs[i].at(0)
         << std::setw(6) << _ADCPerMeVs[i].at(1)<<"\n";
    }
    os << "\n";
  }
}
