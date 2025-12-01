//
//
//
//#include "Offline/AnalysisConditions/inc/TrkQualCatalogCache.hh"
#include "Offline/DbService/inc/DbHandle.hh"
#include "Offline/CRVConditions/inc/CRVCalibCache.hh"
#include "Offline/CRVConditions/inc/CRVOrdinalCache.hh"
#include "Offline/CRVConditions/inc/CRVPhotonYieldCache.hh"
#include "Offline/CRVConditions/inc/CRVStatusCache.hh"
#include "Offline/CaloConditions/inc/CaloDAQMapCache.hh"
#include "Offline/DAQConditions/inc/EventTimingCache.hh"
#include "Offline/DbService/inc/DbService.hh"
#include "Offline/GeometryService/inc/GeometryService.hh"
#include "Offline/ProditionsService/inc/ProditionsService.hh"

#include "Offline/CaloConditions/inc/CaloDAQMapCache.hh"
#include "Offline/CaloConditions/inc/CalCalibCache.hh"

#include "Offline/DAQConditions/inc/EventTimingCache.hh"
#include "Offline/STMConditions/inc/STMEnergyCalibCache.hh"
#include "Offline/SimulationConditions/inc/SimBookkeeperCache.hh"

#include "Offline/TrackerConditions/inc/AlignedTrackerCache.hh"
#include "Offline/TrackerConditions/inc/FullReadoutStrawCache.hh"
#include "Offline/TrackerConditions/inc/Mu2eDetectorCache.hh"
#include "Offline/TrackerConditions/inc/Mu2eMaterialCache.hh"
#include "Offline/TrackerConditions/inc/StrawDriftCache.hh"
#include "Offline/TrackerConditions/inc/StrawElectronicsCache.hh"
#include "Offline/TrackerConditions/inc/StrawPhysicsCache.hh"
#include "Offline/TrackerConditions/inc/StrawResponseCache.hh"
#include "Offline/TrackerConditions/inc/TrackerStatusCache.hh"
#include "Offline/TrackerConditions/inc/TrackerPanelMapCache.hh"

#include "Offline/AnalysisConditions/inc/TrkQualCatalogCache.hh"
#include "Offline/SimulationConditions/inc/SimBookkeeperCache.hh"

#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"
#include <iostream>
#include <typeinfo>

using namespace std;

namespace mu2e {

ProditionsService::ProditionsService(Parameters const& sTable,
                                     art::ActivityRegistry& iRegistry) :
    _config(sTable()) {
  // create this here to force DbService to be active before Proditions
  art::ServiceHandle<DbService> d;
  // and then Geometry
  art::ServiceHandle<GeometryService> g;

  auto cor = std::make_shared<mu2e::CRVOrdinalCache>(_config.crvOrdinal());
  _caches[cor->name()] = cor;
  auto csy = std::make_shared<mu2e::CRVPhotonYieldCache>(_config.crvPhotonYield());
  _caches[csy->name()] = csy;
  auto cst = std::make_shared<mu2e::CRVStatusCache>(_config.crvStatus());
  _caches[cst->name()] = cst;
  auto cca = std::make_shared<mu2e::CRVCalibCache>(_config.crvCalib());
  _caches[cca->name()] = cca;
  auto etc = std::make_shared<mu2e::EventTimingCache>(_config.eventTiming());
  _caches[etc->name()] = etc;
  auto sep =
      std::make_shared<mu2e::STMEnergyCalibCache>(_config.stmEnergyCalib());
  _caches[sep->name()] = sep;
  auto frc =
      std::make_shared<mu2e::FullReadoutStrawCache>(_config.fullReadoutStraw());
  _caches[frc->name()] = frc;
  auto tsc = std::make_shared<mu2e::TrackerStatusCache>(_config.trackerStatus());
  _caches[tsc->name()] = tsc;
  auto tpm = std::make_shared<mu2e::TrackerPanelMapCache>(_config.trackerPanelMap());
  _caches[tpm->name()] = tpm;
  auto sdc = std::make_shared<mu2e::StrawDriftCache>(_config.strawDrift());
  _caches[sdc->name()] = sdc;
  auto spc = std::make_shared<mu2e::StrawPhysicsCache>(_config.strawPhysics());
  _caches[spc->name()] = spc;
  auto sec =
      std::make_shared<mu2e::StrawElectronicsCache>(_config.strawElectronics());
  _caches[sec->name()] = sec;
  auto src =
      std::make_shared<mu2e::StrawResponseCache>(_config.strawResponse());
  _caches[src->name()] = src;
  // tracker alignment has two templated variants, for reco and simulation
  auto atc =
      std::make_shared<mu2e::AlignedTrackerCacheReco>(_config.alignedTracker());
  _caches[atc->name()] = atc;
  auto atcs =
    std::make_shared<mu2e::AlignedTrackerCacheSim>(_config.alignedTrackerSim());
  _caches[atcs->name()+"Sim"] = atcs;
  auto mmc = std::make_shared<mu2e::Mu2eMaterialCache>(_config.mu2eMaterial());
  _caches[mmc->name()] = mmc;
  auto mdc = std::make_shared<mu2e::Mu2eDetectorCache>(_config.mu2eDetector());
  _caches[mdc->name()] = mdc;
  auto cdc =
      std::make_shared<mu2e::CaloDAQMapCache>(_config.caloDAQConditions());
  _caches[cdc->name()] = cdc;
//  auto tqc =
//      std::make_shared<mu2e::TrkQualCatalogCache>(_config.trkQualCatalog());
//  _caches[tqc->name()] = tqc;
  auto bkc =
      std::make_shared<mu2e::SimBookkeeperCache>(_config.simbookkeeper());
  _caches[bkc->name()] = bkc;
  auto cec =
      std::make_shared<mu2e::CalCalibCache>(_config.calCalib());
  _caches[cec->name()] = cec;
  if (_config.verbose() > 0) {
    cout << "Proditions built caches:" << endl;
    for (auto const& cc : _caches) {
      cout << "  " << cc.first << endl;
    }
  }
}

}  // namespace mu2e
