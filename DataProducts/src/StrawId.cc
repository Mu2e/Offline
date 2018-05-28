// C++ includes
#include <iostream>
#include <vector>
// art includes
#include "cetlib_except/exception.h"
// Mu2e includes
#include "DataProducts/inc/StrawId.hh"
#include "GeneralUtilities/inc/splitLine.hh"

using namespace std;

namespace mu2e {
// statics
//
//
  StrawId::StrawId( unsigned short plane,
      unsigned short panel,
      unsigned short straw) : _sid(0) {
    setPlane(plane);
    setPanel(panel);
    setStraw(straw);
  }

  // anonymous namespace for string functions
  namespace {
    StrawId strawIdFromString ( std::string const& s ){
      vector<string> v;
      splitLine( s, "_", v);
      if ( v.size() != 3 ){
	throw cet::exception("CONFIG")
	  << "strawIdFromString: expected three parts but found: "
	  << v.size()
	  << "\n";
      }

      istringstream sdev(v[0]);
      istringstream ssec(v[1]);
      istringstream sstr(v[2]);
      unsigned short plane, panel, straw;
      sdev >> plane;
      ssec >> panel;
      sstr >> straw;
      return StrawId(plane,panel,straw);
    }
  } // end anonymous namespace

  StrawId::StrawId(std::string const& asstring) {
    *this = strawIdFromString(asstring);
  }

  void StrawId::setStraw(unsigned short istraw) {
    if(validStraw(istraw))
      _sid |= istraw;
    else
      throw cet::exception("CONFIG") << "invalid straw " << istraw << "\n";
  }

  void StrawId::setPanel(unsigned short ipanel) {
    if(validPanel(ipanel))
      _sid |= ipanel << _panelsft;
    else
      throw cet::exception("CONFIG") << "invalid panel " << ipanel << "\n";
  }

  void StrawId::setPlane(unsigned short iplane) {
    if(validPlane(iplane))
      _sid |= iplane << _planesft;
    else
      throw cet::exception("CONFIG") << "invalid plane " << iplane << "\n";
  }

  StrawId::isep StrawId::separation(StrawId const& other) const { // was PanelId
    isep retval=apart;
    // same station
    if(other.getPlane()/2 == getPlane()/2){
      int pln1 = getPanel()%2;
      int pln2 = other.getPanel()%2;
      int dp = pln2 - pln1;
      if(other.getPlane() == getPlane()){
	if(dp == 0)
	  retval = same;
	else
	  retval = plane1;
      } else {
	int dd = other.getPlane() - getPlane();
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

  std::ostream& operator<<(std::ostream& ost,
                           const StrawId& s ){
    ost << s.plane() << "_"
        << s.panel() << "_"
        << s.straw();
    return ost;
  }

} // end namespace mu2e
