#include "DataProducts/inc/PanelId.hh"

namespace mu2e {

  PanelId::isep PanelId::separation(PanelId const& other) const {
    isep retval=apart;
    // same station
    if(other.getPlaneId()/2 == getPlaneId()/2){
      int plane1 = getPanel()%2;
      int plane2 = other.getPanel()%2;
      int dp = plane2 - plane1;
      if(other.getPlaneId() == getPlaneId()){
	if(dp == 0)
	  retval = same;
	else
	  retval = plane;
      } else {
	int dd = other.getPlaneId() - getPlaneId();
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

