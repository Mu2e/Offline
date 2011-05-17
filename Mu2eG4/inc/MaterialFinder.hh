#ifndef Mu2eG4_MaterialFinder_hh
#define Mu2eG4_MaterialFinder_hh
//
// Manage lookup of G4Material from a name found in a geometry file.
// Throws if operation cannot be successfully completed.
//
// $Id: MaterialFinder.hh,v 1.3 2011/05/17 15:41:36 greenc Exp $
// $Author: greenc $
// $Date: 2011/05/17 15:41:36 $
//
// Original author Rob Kutschke

// C++ includes.
#include <string>

// Forward references.
class G4Material;

namespace mu2e {

  // Forward references.
  class SimpleConfig;

  class MaterialFinder{

  public:
    MaterialFinder( SimpleConfig const& config );

    // This call has no default; throws if key not found.
    G4Material* get(  std::string const& key );

    // This call has a default.
    G4Material* get(  std::string const& key, std::string const& defaultValue );

  private:
    
    SimpleConfig const * _config;

  };
  
}

#endif /* Mu2eG4_MaterialFinder_hh */



