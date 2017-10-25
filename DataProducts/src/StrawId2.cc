// C++ includes
#include <iostream>
#include <vector>
// art includes
#include "cetlib_except/exception.h"
// Mu2e includes
#include "DataProducts/inc/StrawId2.hh"
#include "GeneralUtilities/inc/splitLine.hh"

using namespace std;

namespace mu2e {
// statics
//
//
  StrawId2::StrawId2( unsigned short plane,
      unsigned short panel,
      unsigned short straw) : _sid(0) {
    bool goodinput =
      setPlane(plane) &&
      setPanel(panel) &&
      setStraw(straw);
    if(!goodinput) throw cet::exception("CONFIG");
  }

  // anonymous namespace for string functions
  namespace {
    StrawId2 strawIdFromString ( std::string const& s ){
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
      return StrawId2(plane,panel,straw);
    }
  }

  StrawId2::StrawId2(std::string const& asstring) {
    *this = strawIdFromString(asstring);
  }

  bool StrawId2::setStraw(unsigned short istraw) {
   bool retval = validStraw(istraw);
    if(retval)
    _sid |= istraw << _strawsft;
    return retval;
  }

  bool StrawId2::setPanel(unsigned short ipanel) {
    bool retval = validPanel(ipanel);
    if(retval)
    _sid |= ipanel << _panelsft;
    return retval;
  }

  bool StrawId2::setPlane(unsigned short iplane) {
    bool retval = validPlane(iplane);
    if(retval)
    _sid |= iplane << _planesft;
    return retval;
  }
}

