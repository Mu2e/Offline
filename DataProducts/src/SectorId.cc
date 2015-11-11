#include "DataProducts/inc/SectorId.hh"

namespace mu2e {

  SectorId::isep SectorId::separation(SectorId const& other) const {
    isep retval=apart;
    // same station
    if(other.getDeviceId()/2 == getDeviceId()/2){
      int plane1 = getSector()%2;
      int plane2 = other.getSector()%2;
      int dp = plane2 - plane1;
      if(other.getDeviceId() == getDeviceId()){
	if(dp == 0)
	  retval = same;
	else
	  retval = device;
      } else {
	int dd = other.getDeviceId() - getDeviceId();
	if(dp == 0)
	  retval = station2;
	else if(dd*dp>0)
	  retval = station3;
	else
	  retval = station1;
      }	
    }
    return retval;
  }
}

