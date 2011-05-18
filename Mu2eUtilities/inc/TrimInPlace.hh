#ifndef Mu2eUtilities_TrimInPlace_hh
#define Mu2eUtilities_TrimInPlace_hh

/**
 *
 * Remove leading and trailing whitespace from a string.
 * It modifies the input string.
 *
 * $Id: TrimInPlace.hh,v 1.3 2011/05/18 02:27:18 wb Exp $
 * $Author: wb $
 * $Date: 2011/05/18 02:27:18 $
 *
 * Original author Rob Kutschke
 *
 */

#include <string>

namespace mu2e {

void TrimInPlace( std::string& s);

}
#endif /* Mu2eUtilities_TrimInPlace_hh */
