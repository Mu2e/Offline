#ifndef GeneralUtilities_trimInPlace_hh
#define GeneralUtilities_trimInPlace_hh
//
// Remove leading and trailing whitespace from a string.
// It modifies the input string.
//
// $Id: trimInPlace.hh,v 1.1 2012/07/15 22:00:36 kutschke Exp $
// $Author: kutschke $
// $Date: 2012/07/15 22:00:36 $
//
// Contact person Rob Kutschke
//

#include <string>

namespace mu2e {

void trimInPlace( std::string& s);

}
#endif /* General_trimInPlace_hh */
