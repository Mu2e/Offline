// C++ includes
#include <iostream>
#include <sstream>
#include <vector>
// art includes
#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
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

  StrawId::StrawId( LayerId layer,
                           int n
                           ):
    _lid(layer),
    _n(n){

    // get the layer number to check if the starw number is correct
    std::ostringstream os;
    os << layer;
    vector<string> v;
    splitLine( os.str(), "_", v);
    istringstream slay(v[2]);
    int lay;
    slay >> lay;

    if ( n%2!=lay%2 ) {
      std::cerr << "CONFIG " 
        //      mf::LogWarning("CONFIG")
        //      throw cet::exception("CONFIG")
                << "StrawId(LayerId,int): incorrect straw in layer "
                << layer  << "_" << n
                << "\n";
    }

  }

}

