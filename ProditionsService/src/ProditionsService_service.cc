//
//
//
#include <iostream>
#include <typeinfo>
#include "DbService/inc/DbService.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "ProditionsService/inc/ProditionsService.hh"

#include "TrackerConditions/inc/FullReadoutStrawCache.hh"
#include "TrackerConditions/inc/DeadStrawCache.hh"
#include "TrackerConditions/inc/StrawDriftCache.hh"
#include "TrackerConditions/inc/StrawPhysicsCache.hh"
#include "TrackerConditions/inc/StrawElectronicsCache.hh"
#include "TrackerConditions/inc/StrawResponseCache.hh"
#include "TrackerConditions/inc/AlignedTrackerCache.hh"
#include "TrackerConditions/inc/Mu2eMaterialCache.hh"
#include "TrackerConditions/inc/Mu2eDetectorCache.hh"

using namespace std;

namespace mu2e {

  ProditionsService::ProditionsService(Parameters const& sTable,
             art::ActivityRegistry& iRegistry) : _config(sTable())
  {

    // create this here to force DbService to be active before Proditions
    art::ServiceHandle<DbService> d;
    // and then Geometry
    art::ServiceHandle<GeometryService> g;

    auto frc = std::make_shared<mu2e::FullReadoutStrawCache>(_config.fullReadoutStraw());
    _caches[frc->name()] = frc;
    auto dsc = std::make_shared<mu2e::DeadStrawCache>(_config.deadStraw());
    _caches[dsc->name()] = dsc;
    auto sdc = std::make_shared<mu2e::StrawDriftCache>(_config.strawDrift());
    _caches[sdc->name()] = sdc;
    auto spc = std::make_shared<mu2e::StrawPhysicsCache>(_config.strawPhysics());
    _caches[spc->name()] = spc;
    auto sec = std::make_shared<mu2e::StrawElectronicsCache>(_config.strawElectronics());
    _caches[sec->name()] = sec;
    auto src = std::make_shared<mu2e::StrawResponseCache>(_config.strawResponse());
    _caches[src->name()] = src;
    auto atc = std::make_shared<mu2e::AlignedTrackerCache>(_config.alignedTracker());
    _caches[atc->name()] = atc;
    auto mmc = std::make_shared<mu2e::Mu2eMaterialCache>(_config.mu2eMaterial());
    _caches[mmc->name()] = mmc;
    auto mdc = std::make_shared<mu2e::Mu2eDetectorCache>(_config.mu2eDetector());
    _caches[mdc->name()] = mdc;

    if( _config.verbose()>0) {
      cout << "Proditions built caches:" << endl;
      for( auto cc : _caches) {
	cout << "  " << cc.first << endl;
      }
    }

  }

}

DEFINE_ART_SERVICE(mu2e::ProditionsService);
