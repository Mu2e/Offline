#ifndef GeneralUtilities_toHex_hh
#define GeneralUtilities_toHex_hh

//
//  Format a number in hex.
//
// $Id: toHex.hh,v 1.1 2013/03/01 01:21:37 kutschke Exp $
// $Author: kutschke $
// $Date: 2013/03/01 01:21:37 $
//
// Original author Rob Kutschke
//

#include <string>

namespace mu2e {

  std::string toHex( int i );
  std::string toHex( unsigned i );

}
#endif /* GeneralUtilities_toHex_hh */
