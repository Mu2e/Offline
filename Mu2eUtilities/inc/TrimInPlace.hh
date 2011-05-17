#ifndef Mu2eUtilities_TrimInPlace_hh
#define Mu2eUtilities_TrimInPlace_hh

/**
 * 
 * Remove leading and trailing whitespace from a string.
 * It modifies the input string.
 *
 * $Id: TrimInPlace.hh,v 1.2 2011/05/17 15:41:36 greenc Exp $
 * $Author: greenc $ 
 * $Date: 2011/05/17 15:41:36 $
 *
 * Original author Rob Kutschke
 * 
 */

#include <string>

namespace mu2e {

void TrimInPlace( std::string& s);

}
#endif /* Mu2eUtilities_TrimInPlace_hh */
