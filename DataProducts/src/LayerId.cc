// C++ includes
#include <iostream>
#include <vector>
// art includes
#include "cetlib/exception.h"
// Mu2e includes
#include "GeneralUtilities/inc/splitLine.hh"
#include "DataProducts/inc/LayerId.hh"
using namespace std;

namespace mu2e {
// anonymous namespace for string functions
  namespace {
    LayerId layerIdFromString ( std::string const& s ){
      vector<string> v;
      splitLine( s, "_", v);
      if ( v.size() != 3 ){
        throw cet::exception("CONFIG")
          << "layerIdFromString: expected three parts but found: "
          << v.size()
          << "\n";
      }

      istringstream sdev(v[0]);
      istringstream ssec(v[1]);
      istringstream slay(v[2]);
      int dev, sec, lay;
      sdev >> dev;
      ssec >> sec;
      slay >> lay;
      return LayerId(dev,sec,lay);
    }
  }
  LayerId::LayerId(std::string const& s) {
    *this = layerIdFromString(s);
  }
}


