#ifndef CRVConditions_CRVStatus_hh
#define CRVConditions_CRVStatus_hh

//
// Prodition for CRV bad SiPM channel list
//

#include "Offline/DataProducts/inc/CRVId.hh"
#include "Offline/Mu2eInterfaces/inc/ProditionsEntity.hh"
#include <cstdint>
#include <map>

namespace mu2e {

class CRVStatus : virtual public ProditionsEntity {
 public:
  enum Flags{notConnected=0, ignoreChannel=1, noData=2, noPedestal=3, noCalibConst=4, noisy=5};

  typedef std::shared_ptr<CRVStatus> ptr_t;
  typedef std::shared_ptr<const CRVStatus> cptr_t;
  constexpr static const char* cxname = {"CRVStatus"};

  typedef std::map<std::uint16_t, int> StatusMap;

  CRVStatus(const StatusMap& smap) : ProditionsEntity(cxname), _smap(smap) {}

  // return status flag word for an offline channel
  int status(std::uint16_t channel) const {
    auto it = _smap.find(channel);
    if (it == _smap.end()) {
      return 0;
    } else {
      return it->second;
    }
  }

  const StatusMap& map() const { return _smap; }

 private:
  // map of status word for sparse channel number
  StatusMap _smap;
};

}  // namespace mu2e

#endif
