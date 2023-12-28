#ifndef CRVConditions_CRVScintYield_hh
#define CRVConditions_CRVScintYield_hh

//
// Holds scintillation yield spread of CRV counters
//

#include "Offline/DataProducts/inc/CRSScintillatorBarIndex.hh"
#include "Offline/Mu2eInterfaces/inc/ProditionsEntity.hh"
#include "cetlib_except/exception.h"
#include <cstdint>
#include <array>
#include <vector>

namespace mu2e {

class CRVScintYield : virtual public ProditionsEntity {
 public:
  typedef std::shared_ptr<CRVScintYield> ptr_t;
  typedef std::shared_ptr<const CRVScintYield> cptr_t;
  constexpr static const char* cxname = {"CRVScintYield"};

  typedef std::vector<float> ScintYieldMap;

  CRVScintYield(ScintYieldMap const& scintYieldMap) :
      ProditionsEntity(cxname), _scintYieldMap(scintYieldMap) {}

  // scintillation yield deviation (from nominal value) for a scintillator bar index
  float scintYieldDeviation(size_t barIndex) const {
    if (barIndex >= _scintYieldMap.size()) {
      throw cet::exception("CRVSCINTYIELD_BAD_SCINTILLATOR BAR INDEX")
          << "CRVScintYield::scinYieldDeviation bad barIndex requested: "
          << " barIndex=" << barIndex << "\n";
    }
    return _scintYieldMap.at(barIndex);
  }

 private:
  ScintYieldMap _scintYieldMap;
};

}  // namespace mu2e

#endif
