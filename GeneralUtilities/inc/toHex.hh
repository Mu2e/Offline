#ifndef GeneralUtilities_toHex_hh
#define GeneralUtilities_toHex_hh

//
//  Format a number in hex.
//
// $Id: toHex.hh,v 1.2 2013/03/03 18:08:28 kutschke Exp $
// $Author: kutschke $
// $Date: 2013/03/03 18:08:28 $
//
// Original author Rob Kutschke
//

#include <string>

namespace mu2e {

  std::string toHex( int  i );
  std::string toHex( long i );
  std::string toHex( unsigned      i );
  std::string toHex( unsigned long i );

}
#endif /* GeneralUtilities_toHex_hh */
