#ifndef getMaterial_HH
#define getMaterial_HH
//
// Manage lookup of G4Material from a name found in a geometry file.
// Throws if operation cannot be successfully completed.
//
// $Id: MaterialFinder.hh,v 1.1 2010/04/13 23:09:54 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/04/13 23:09:54 $
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

    // key is the name of variable from a geometry file.
    G4Material* get(  std::string const& key );

  private:
    
    SimpleConfig const * _config;

  };
  
}

#endif



