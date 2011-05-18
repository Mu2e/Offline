#ifndef Mu2eUtilities_FileInPath_hh
#define Mu2eUtilities_FileInPath_hh
//
// Hack a subset of the old framework's FileInPath class.
// Only implement the subset used in the Offline code.
//
//  $Id: FileInPath.hh,v 1.1 2011/05/18 04:21:59 kutschke Exp $
//  $Author: kutschke $
//  $Date: 2011/05/18 04:21:59 $
//
// Original author Rob Kutschke
//
// Notes:
//  1) The old version of this was driven by 4 environment variables.
//     The new version is driven by a new environment variable,
//     MU2E_SEARCH_PATH, which is a colon separated list of absolute paths.
//
#include <string>

namespace mu2e {

  class FileInPath {

  public:
    FileInPath( std::string const& filename);

    std::string fullPath() const{ return fullPath_;}

  private:

    std::string filename_;
    std::string fullPath_;

  };

} // namespace mu2e

#endif /* Mu2eUtilities_FileInPath_hh */
