#ifndef CRVConditions_CRVOrdinal_hh
#define CRVConditions_CRVOrdinal_hh

//
// Holds conversion between CRV Offline and Online numbering
//

#include "Offline/CRVConditions/inc/CRVROC.hh"
#include "Offline/DataProducts/inc/CRVId.hh"
#include "Offline/Mu2eInterfaces/inc/ProditionsEntity.hh"
#include "cetlib_except/exception.h"
#include <cstdint>
#include <array>
#include <vector>

namespace mu2e {

class CRVOrdinal : virtual public ProditionsEntity {
 public:
  typedef std::shared_ptr<CRVOrdinal> ptr_t;
  typedef std::shared_ptr<const CRVOrdinal> cptr_t;
  constexpr static const char* cxname = {"CRVOrdinal"};

  // online numbers for each offline number
  typedef std::vector<CRVROC> OnlineMap;
  // this is a 3-dim array: offline number = x[ROC][FEB][FEBchan]
  typedef std::array<std::array<std::array<std::uint16_t, CRVId::nChanPerFEB>,
                                CRVId::nFEBPerROC+1>,   //FEB numbers start at 1
                                CRVId::nROC+1>          //ROC numbers start at 1
                                OfflineMap;

  CRVOrdinal(OnlineMap const& onMap, OfflineMap const& offMap) :
      ProditionsEntity(cxname), _onMap(onMap), _offMap(offMap) {}

  // check, if a channel exists
  bool onlineExists(std::uint16_t channel) const {
    if (channel >= _onMap.size()) return false;
    if (_onMap.at(channel).FEBchannel() >= CRVId::nChanPerFEB) return false;
    return true;
  }

  // online numbering triplet for an offline channel number
  const CRVROC& online(std::uint16_t channel) const {
    if (_onMap.at(channel).FEBchannel() >= CRVId::nChanPerFEB) {
      throw cet::exception("CRVORDINAL_BAD_OFFLINE_CHANNEL")
          << "CRVOrdinal::online bad channel requested: "
          << " channel=" << channel << "\n";
    }
    return _onMap[channel];
  }

  // offline channel number for online numbering triplet
  const std::uint16_t offline(const CRVROC& online) const {
    std::uint16_t offline =
        _offMap.at(online.ROC()).at(online.FEB()).at(online.FEBchannel());
    if (offline >= CRVId::nChannels) {
      throw cet::exception("CRVORDINAL_BAD_ONLINE_CHANNEL")
          << "CRVOrdinal::offline bad channel requested: "
          << " ROC=" << online.ROC() << " FEB=" << online.FEB()
          << " FEBchannel=" << online.FEBchannel() << "\n";
    }
    return offline;
  }

 private:
  OnlineMap _onMap;
  OfflineMap _offMap;
};

}  // namespace mu2e

#endif
