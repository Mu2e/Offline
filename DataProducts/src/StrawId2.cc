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
    setPlane(plane);
    setPanel(panel);
    setStraw(straw);
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

  void StrawId2::setStraw(unsigned short istraw) {
    if(validStraw(istraw))
      _sid |= istraw;
    else
      throw cet::exception("CONFIG") << "invalid straw " << istraw << "\n";
  }

  void StrawId2::setPanel(unsigned short ipanel) {
    if(validPanel(ipanel))
      _sid |= ipanel << _panelsft;
    else
      throw cet::exception("CONFIG") << "invalid panel " << ipanel << "\n";
  }

  void StrawId2::setPlane(unsigned short iplane) {
    if(validPlane(iplane))
      _sid |= iplane << _planesft;
    else
      throw cet::exception("CONFIG") << "invalid plane " << iplane << "\n";
  }
}

