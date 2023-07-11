// C++ includes
#include <iostream>

// Mu2e includes
#include "Offline/TrackerConditions/inc/TrackerStatus.hh"
#include "Offline/DbTables/inc/TrkElementStatus.hh"
#include <sstream>

using namespace std;

namespace mu2e {

  void TrackerStatus::addStatus(StrawId const& sid, StrawIdMask const& mask, StrawStatus const& status) {
    StrawId mid = mask.maskStrawId(sid);
    auto ifnd = _status.find(mask);
    if(ifnd == _status.end())
      _status[mask][mid] = status;
    else {
      auto& smap = ifnd->second;
      auto jfnd = smap.find(mid);
      if(jfnd == smap.end())
        smap[mid] = status;
      else
        jfnd->second.merge(status);
    }
  }

  void TrackerStatus::print( ostream& out) const{
    for( auto ielem = _status.begin(); ielem != _status.end(); ++ielem) {
      auto const& mask = ielem->first;
      out << "Tracker level " << mask.levelName() << " has the following elements:" << std::endl;
      auto const& stati = ielem->second;
      for(auto istat = stati.begin(); istat != stati.end(); ++istat){
        out << " Element with Id " << istat->first <<
          " has status " << istat->second << std::endl;
      }
    }
  }

  StrawStatus TrackerStatus::strawStatus(StrawId const& sid) const {
    StrawStatus sstat;
    for( auto ielem = _status.begin(); ielem != _status.end(); ++ielem) {
      auto const& mask = ielem->first;
      auto const& stati = ielem->second;
      for(auto istat = stati.begin(); istat != stati.end(); ++istat){
        auto const& elem = istat->first;
        auto const& stat = istat->second;
        if(mask.equal(sid,elem))sstat.merge(stat);
      }
    }
    return sstat;
  }

  StrawStatus TrackerStatus::panelStatus(StrawId const& sid) const {
    StrawStatus sstat;
    for( auto ielem = _status.begin(); ielem != _status.end(); ++ielem) {
      auto const& mask = ielem->first;
      if(mask.level() != StrawIdMask::straw && mask.level() != StrawIdMask::uniquestraw){ // plane status also affects a panel
        auto const& stati = ielem->second;
        for(auto istat = stati.begin(); istat != stati.end(); ++istat){
          auto const& elem = istat->first;
          auto const& stat = istat->second;
          if(mask.equal(sid,elem))sstat.merge(stat);
        }
      }
    }
    return sstat;
  }

  StrawStatus TrackerStatus::planeStatus(StrawId const& sid) const {
    StrawStatus sstat;
    for( auto ielem = _status.begin(); ielem != _status.end(); ++ielem) {
      auto const& mask = ielem->first;
      if(mask.level() == StrawIdMask::plane || mask.level() == StrawIdMask::tracker){
        auto const& stati = ielem->second;
        for(auto istat = stati.begin(); istat != stati.end(); ++istat){
          auto const& elem = istat->first;
          auto const& stat = istat->second;
          if(mask.equal(sid,elem))sstat.merge(stat);
        }
      }
    }
    return sstat;
  }

  // no signal expected from this straw
  bool TrackerStatus::noSignal(StrawId const& sid) const {
    static StrawStatus mask = StrawStatus(StrawStatus::absent) |
      StrawStatus(StrawStatus::nowire) |
      StrawStatus(StrawStatus::noHV) |
      StrawStatus(StrawStatus::noLV) |
      StrawStatus(StrawStatus::nogas) |
      StrawStatus(StrawStatus::lowgasgain) |
      StrawStatus(StrawStatus::noHVPreamp) |
      StrawStatus(StrawStatus::noCalPreamp) |
      StrawStatus(StrawStatus::disabled);
    StrawStatus status = strawStatus(sid);
    return status.hasAnyProperty(mask);
  }

  bool TrackerStatus::noisy(StrawId const& sid) const {
    static StrawStatus mask = StrawStatus(StrawStatus::sparking) |
      StrawStatus(StrawStatus::noise) |
      StrawStatus(StrawStatus::pickup);
    StrawStatus status = strawStatus(sid);
    return status.hasAnyProperty(mask);
  }

  bool TrackerStatus::suppress(StrawId const& sid) const {
    static StrawStatus mask(StrawStatus::suppress);
    StrawStatus status = strawStatus(sid);
    return status.hasAnyProperty(mask);
  }
  //
  bool TrackerStatus::noMaterial(StrawId const& sid) const {
    static StrawStatus mask = StrawStatus(StrawStatus::absent);
    StrawStatus status = strawStatus(sid);
    return status.hasAnyProperty(mask);
  }
} // namespace mu2e
