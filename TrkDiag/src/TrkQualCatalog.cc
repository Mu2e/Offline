#include "TrkDiag/inc/TrkQualCatalog.hh"

namespace mu2e {

  void TrkQualCatalog::print(std::ostream& os) const {
    os << "Entries in TrkQualCatalog:" << std::endl;
    for (const auto& i_entry : _entries) {
      os << "Train Name: " << i_entry._trainName << ", XML File: " << i_entry._xmlFileName << std::endl;
    }
    os << std::endl;
  }
}
