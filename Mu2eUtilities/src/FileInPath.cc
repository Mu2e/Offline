//
// Hack a subset of the old framework's FileInPath class.
//
//  $Id: FileInPath.cc,v 1.1 2011/05/18 04:21:53 kutschke Exp $
//  $Author: kutschke $
//  $Date: 2011/05/18 04:21:53 $
//
// Original author Rob Kutschke

// Framework includes
#include "cetlib/exception.h"
#include "cetlib/search_path.h"

// Mu2e includes
#include "Mu2eUtilities/inc/FileInPath.hh"

namespace mu2e {
  
  FileInPath::FileInPath( std::string const& filename):
    filename_(filename),
    fullPath_(){

    cet::search_path sp("MU2E_SEARCH_PATH");
    if( ! sp.find_file(filename, fullPath_) )
      throw "SimpleConfig c'tor: find_file failure!";  // TODO: improve exception
    
  }

} // namespace mu2e
