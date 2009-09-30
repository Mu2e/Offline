#ifndef findMaterialOrThrow_H
#define findMaterialOrThrow_H 1
//
// Free function wrapper around
//   G4NistManager::FindOrBuildMaterial
// The wrapper does the job of throwing if the pointer comes
// back null.
//
// $Id: findMaterialOrThrow.hh,v 1.1 2009/09/30 22:57:47 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2009/09/30 22:57:47 $
//
// Original author Rob Kutschke
//
//

class G4Material;
class G4String;

namespace mu2e {
  
  G4Material* findMaterialOrThrow( G4String const& name);

} // end namespace mu2e
#endif
