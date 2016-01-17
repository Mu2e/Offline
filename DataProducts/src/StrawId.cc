// C++ includes
#include <iostream>
#include <vector>
// art includes
#include "cetlib/exception.h"
// Mu2e includes
#include "DataProducts/inc/StrawId.hh"
#include "GeneralUtilities/inc/splitLine.hh"

using namespace std;

namespace mu2e {
// anonymous namespace for string functions
  namespace {
    StrawId strawIdFromString ( std::string const& s ){
      vector<string> v;
      splitLine( s, "_", v);
      if ( v.size() != 4 ){
        throw cet::exception("CONFIG")
          << "strawIdFromString: expected four parts but found: "
          << v.size()
          << "\n";
      }

      istringstream sdev(v[0]);
      istringstream ssec(v[1]);
      istringstream slay(v[2]);
      istringstream sstr(v[3]);
      int dev, sec, lay, str;
      sdev >> dev;
      ssec >> sec;
      slay >> lay;
      sstr >> str;
      return StrawId(dev,sec,lay,str);
    }
  }

  StrawId::StrawId(std::string const& asstring) {
    *this = strawIdFromString(asstring);
  }
}

