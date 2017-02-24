// C++ includes
#include <iostream>
#include <vector>
// art includes
#include "cetlib_except/exception.h"
// Mu2e includes
#include "GeneralUtilities/inc/splitLine.hh"
#include "DataProducts/inc/PanelId.hh"
using namespace std;

namespace mu2e {
// anonymous namespace for string functions
  namespace {
    PanelId panelIdFromString ( std::string const& s ){
      vector<string> v;
      splitLine( s, "_", v);
      if ( v.size() != 2 ){
        throw cet::exception("CONFIG")
          << "panelIdFromString: expected two parts but found: "
          << v.size()
          << "\n";
      }

      istringstream sdev(v[0]);
      istringstream ssec(v[1]);
      int dev, sec;
      sdev >> dev;
      ssec >> sec;
      return PanelId(dev,sec);
    }
  }
 
  PanelId::PanelId(std::string const& s) {
    *this = panelIdFromString(s);
  }

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

