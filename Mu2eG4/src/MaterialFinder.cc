//
// Manage lookup of G4Material from a name found in a geometry file.
// Throws if operation cannot be successfully completed.
//
// $Id: MaterialFinder.cc,v 1.1 2010/04/13 23:09:55 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/04/13 23:09:55 $
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
  
  // key is the name of variable from a geometry file.
  G4Material* MaterialFinder::get( string const& key ){
    string materialName = _config->getString(key);
    return findMaterialOrThrow(materialName);
  }


} // end namespace mu2e
