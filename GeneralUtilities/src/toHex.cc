//
//  Format a number in hex.
//
// $Id: toHex.cc,v 1.2 2013/03/03 18:08:28 kutschke Exp $
// $Author: kutschke $
// $Date: 2013/03/03 18:08:28 $
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
