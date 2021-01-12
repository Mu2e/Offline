#ifndef GeneralUtilities_toHex_hh
#define GeneralUtilities_toHex_hh

//
//  Format a number in hex.
//
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
