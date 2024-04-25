#ifndef CRVConditions_CRVPhotonYield_hh
#define CRVConditions_CRVPhotonYield_hh

//
// Holds photon yield deviations of CRV channels
//

#include "Offline/DataProducts/inc/CRSScintillatorBarIndex.hh"
#include "Offline/Mu2eInterfaces/inc/ProditionsEntity.hh"
#include "cetlib_except/exception.h"
#include <cstdint>
#include <array>
#include <vector>

namespace mu2e {

class CRVPhotonYield : virtual public ProditionsEntity {
 public:
  typedef std::shared_ptr<CRVPhotonYield> ptr_t;
  typedef std::shared_ptr<const CRVPhotonYield> cptr_t;
  constexpr static const char* cxname = {"CRVPhotonYield"};

  typedef std::vector<float> PhotonYieldVec;

  CRVPhotonYield(PhotonYieldVec const& svec) : ProditionsEntity(cxname), _svec(svec) {}

  // photon yield deviation (from nominal value) for a channel
  float photonYieldDeviation(std::uint16_t channel) const {
    if (channel >= _svec.size()) {
      throw cet::exception("CRVPHOTONYIELD_BAD_SCINTILLATOR BAR INDEX")
          << "CRVPhotonYield::photonYieldDeviation bad channel requested: "
          << " channel=" << channel << "\n";
    }
    return _svec.at(channel);
  }

 private:
  PhotonYieldVec _svec;
};

}  // namespace mu2e

#endif
