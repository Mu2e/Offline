#ifndef getMaterial_HH
#define getMaterial_HH
//
// Manage lookup of G4Material from a name found in a geometry file.
// Throws if operation cannot be successfully completed.
//
// $Id: MaterialFinder.hh,v 1.2 2010/04/16 14:46:05 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/04/16 14:46:05 $
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

#endif



