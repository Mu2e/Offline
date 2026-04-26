#include "Offline/DataProducts/inc/CaloConst.hh"
#include <ostream>

namespace mu2e {

std::string CaloConst::detTypeName(detType type) {
  switch (type) {
  case CaloConst::detType::CsI:
    return "CsI";
  case CaloConst::detType::CAPHRI:
    return "CAPHRI";
  case CaloConst::detType::PINDiode:
    return "PINDiode";
  case CaloConst::detType::Invalid:
    return "Invalid";
  default:
    return "UNKNOWN";
  }
}

std::ostream& operator<<(std::ostream& ost, const CaloConst::detType& type) {
  ost << CaloConst::detTypeName(type);
  return ost;
}

} // end namespace mu2e
