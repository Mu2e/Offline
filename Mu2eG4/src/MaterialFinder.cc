//
// Manage lookup of G4Material from a name found in a geometry file.
// Throws if operation cannot be successfully completed.
//
// $Id: MaterialFinder.cc,v 1.3 2011/05/17 15:36:00 greenc Exp $
// $Author: greenc $
// $Date: 2011/05/17 15:36:00 $
//
// Original author Rob Kutschke
//

#include "Mu2eG4/inc/MaterialFinder.hh"
#include "Mu2eG4/inc/findMaterialOrThrow.hh"
#include "Mu2eUtilities/inc/SimpleConfig.hh"
#include "GeometryService/inc/GeomHandle.hh"

#include "G4Material.hh"

using namespace std;

namespace mu2e {

  MaterialFinder::MaterialFinder( SimpleConfig const& config):
    _config(&config){
  }
  
  // This call has no default; throws if key not found.
  G4Material* MaterialFinder::get( string const& key ){
    string materialName = _config->get<std::string>(key);
    return findMaterialOrThrow(materialName);
  }

  // This call has a default.
  G4Material* MaterialFinder::get( string const& key, string const& defaultValue ){
    string materialName = _config->get<std::string>(key,defaultValue);
    return findMaterialOrThrow(materialName);
  }


} // end namespace mu2e
