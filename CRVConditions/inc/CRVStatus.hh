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
  typedef std::shared_ptr<CRVStatus> ptr_t;
  typedef std::shared_ptr<const CRVStatus> cptr_t;
  constexpr static const char* cxname = {"CRVStatus"};

  typedef std::map<std::uint16_t, int> StatusMap;

  CRVStatus(const StatusMap& smap) : ProditionsEntity(cxname), _smap(smap) {}

  // return status flag word for an offline channel
  int status(std::uint16_t channel) {
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
