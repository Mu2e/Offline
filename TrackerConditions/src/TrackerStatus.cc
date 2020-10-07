// C++ includes
#include <iostream>

// Mu2e includes
#include "TrackerConditions/inc/TrackerStatus.hh"
#include <sstream>

using namespace std;

namespace mu2e {

  void TrackerStatus::print( ostream& out) const{
    for( auto const& estat :  _estatus) {
      out << "Tracker Element with Id " << estat.sid_
	<< " level " << estat.mask_.levelName()  
	<< " status " << estat.status_ << std::endl;
    }
  }

  StrawStatus TrackerStatus::strawStatus(StrawId const& sid) const {
    StrawStatus sstat;
    for(auto const& stat : _estatus){
      if(stat.mask_.equal(sid,stat.sid_))sstat.merge(stat.status_);
    }
    return sstat;
  }

  StrawStatus TrackerStatus::panelStatus(StrawId const& sid) const {
    StrawStatus sstat;
    for(auto const& stat : _estatus){
      // don't count status specific to a straw
      if((stat.mask_.level() == StrawIdMask::panel || stat.mask_.level() == StrawIdMask::uniquepanel)
	  && stat.mask_.equal(sid,stat.sid_))sstat.merge(stat.status_);
    }
    return sstat;
  }

  StrawStatus TrackerStatus::planeStatus(StrawId const& sid) const {
    StrawStatus sstat;
    for(auto const& stat : _estatus){
      // don't count status specific to a straw or panel
      if((stat.mask_.level() == StrawIdMask::plane) && stat.mask_.equal(sid,stat.sid_))sstat.merge(stat.status_);
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
      StrawStatus(StrawStatus::noPreamp);
    StrawStatus status = strawStatus(sid);
    return status.hasAnyProperty(mask);
  }

  bool TrackerStatus::suppress(StrawId const& sid) const { 
    static StrawStatus mask = StrawStatus(StrawStatus::sparking) |
      StrawStatus(StrawStatus::noise) | 
      StrawStatus(StrawStatus::pickup) | 
      StrawStatus(StrawStatus::suppress);
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
