#ifndef CRVDQM_inc_CRVDQMLayout_hh
#define CRVDQM_inc_CRVDQMLayout_hh
// The layout a CRV DQM module injects into the clients, from geometry and
// conditions objects the module looked up (this header looks nothing up, so
// the clients stay service-free). Shared by the offline DqmCrv* modules and the
// otsdaq online modules. The clients' histogram set is frozen in CRVDQMRun1.hh;
// the geometry only chooses which configuration's family a job fills, and must
// match it: an unknown geometry, or a sector with no histograms, throws.
//
// Header-only: a module using it links CosmicRayShieldGeom and CRVConditions.
//
// Original Author: R. Mina

#include "Offline/CRVConditions/inc/CRVOrdinal.hh"
#include "Offline/CRVConditions/inc/CRVStatus.hh"
#include "Offline/CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "Offline/CRVDQM/inc/CRVDQMRun1.hh"
#include "Offline/CRVDQM/inc/CRVDigiDQM.hh"
#include "Offline/DataProducts/inc/CRSScintillatorBarIndex.hh"
#include "Offline/DataProducts/inc/CRVId.hh"

#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <bitset>
#include <cstdint>
#include <string>
#include <vector>

namespace mu2e {
namespace CRVDQMLayout {

// CRVDQMRun1 configuration of this geometry.
inline int configuration(const CosmicRayShield& crs) {
  const int c = CRVDQMRun1::configurationFor(crs.getName());
  if (c < 0) {
    throw cet::exception("CRVDQMLayout")
        << "CRV geometry \"" << crs.getName() << "\" is not a CRVDQMRun1 configuration. "
        << "Adding one is a binning-version change in Offline/DQMHelpers.\n";
  }
  return c;
}

// Offline channel -> index into the configuration's sector list; -1 for a
// channel flagged notConnected in the CRV status conditions (none are dropped
// when sipmStatus is null, e.g. a job without CRV status conditions).
inline std::vector<int> channelToSector(const CosmicRayShield& crs, int configuration,
                                        const CRVStatus* sipmStatus) {
  const auto& sectors = crs.getCRSScintillatorShields();
  std::vector<int> sectorIndex(sectors.size(), -1);
  for (std::size_t i = 0; i < sectors.size(); ++i) {
    sectorIndex[i] = CRVDQMRun1::sectorIndex(configuration, sectors[i].getName());
    if (sectorIndex[i] < 0) {
      throw cet::exception("CRVDQMLayout")
          << "CRV sector \"" << sectors[i].getName() << "\" of geometry \"" << crs.getName()
          << "\" has no histograms in configuration "
          << CRVDQMRun1::configurationName(configuration) << ".\n";
    }
  }
  const std::size_t nChannels = crs.getAllCRSScintillatorBars().size() * CRVId::nChanPerBar;
  std::vector<int> out(nChannels, -1);
  for (std::size_t channel = 0; channel < nChannels; ++channel) {
    if (sipmStatus != nullptr &&
        std::bitset<16>(sipmStatus->status(channel)).test(CRVStatus::Flags::notConnected)) {
      continue;
    }
    const CRSScintillatorBarIndex bar(channel / CRVId::nChanPerBar);
    out[channel] = sectorIndex.at(crs.getBar(bar).id().getShieldNumber());
  }
  return out;
}

// Which sector, module and side each FEB port reads, from the channel map, and
// the layer of each offline channel. One FEB reads two layers on one side of
// one module; a FEB that maps to more than one module side gets no partner
// timing, with one warning.
inline void febTopology(const CosmicRayShield& crs, const CRVOrdinal& channelMap,
                        std::vector<CRVDigiDQM::FebTopology>& topology,
                        std::vector<int>& channelToLayer) {
  const std::size_t nChannels = crs.getAllCRSScintillatorBars().size() * CRVId::nChanPerBar;
  topology.assign(CRVDQMRun1::kNFebPorts, CRVDigiDQM::FebTopology{});
  channelToLayer.assign(nChannels, -1);
  std::vector<bool> conflicting(CRVDQMRun1::kNFebPorts, false);
  for (std::size_t channel = 0; channel < nChannels; ++channel) {
    const CRSScintillatorBarIndex bar(channel / CRVId::nChanPerBar);
    const auto& id = crs.getBar(bar).id();
    channelToLayer[channel] = id.getLayerNumber();

    const auto offline = static_cast<std::uint16_t>(channel);
    if (!channelMap.onlineExists(offline)) continue;  // takes the offline channel
    const CRVROC& online = channelMap.online(offline);
    if (!CRVDQMRun1::onlineIdInRange(online.ROC(), online.FEB(), online.FEBchannel())) continue;
    const int port = CRVDQMRun1::febPort(online.ROC(), online.FEB());
    if (conflicting[port]) continue;

    const int side = static_cast<int>(channel % CRVId::nChanPerBar % CRVId::nSidesPerBar);
    CRVDigiDQM::FebTopology& t = topology[port];
    if (!t.valid) {
      t = {id.getShieldNumber(), id.getModuleNumber(), side, true};
    } else if (t.sector != id.getShieldNumber() || t.module != id.getModuleNumber() ||
               t.side != side) {
      mf::LogWarning("CRVDQMLayout") << "FEB port " << port << " spans more than one module side; "
                               << "partner timing for it is disabled.";
      t = CRVDigiDQM::FebTopology{};
      conflicting[port] = true;
    }
  }
}

}  // namespace CRVDQMLayout
}  // namespace mu2e

#endif /* CRVDQM_inc_CRVDQMLayout_hh */
