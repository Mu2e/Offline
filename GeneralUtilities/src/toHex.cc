//
//  Format a number in hex.
//
//
// Original author Rob Kutschke
//

#include "GeneralUtilities/inc/toHex.hh"

#include<sstream>
#include<iomanip>

using namespace std;

namespace mu2e {

  std::string toHex( int i ){
    ostringstream os;
    os<< std::showbase << std::internal << std::hex << i;
    return os.str();
  }

  std::string toHex( long i ){
    ostringstream os;
    os<< std::showbase << std::internal << std::hex << i;
    return os.str();
  }

  std::string toHex( unsigned i ){
    ostringstream os;
    os<< std::showbase << std::internal << std::hex << i;
    return os.str();
  }

  std::string toHex( unsigned long i ){
    ostringstream os;
    os<< std::showbase << std::internal << std::hex << i;
    return os.str();
  }

}
