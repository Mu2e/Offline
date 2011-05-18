#ifndef Mu2eG4_findMaterialOrThrow_hh
#define Mu2eG4_findMaterialOrThrow_hh
//
// Free function wrapper around
//   G4NistManager::FindOrBuildMaterial
// The wrapper does the job of throwing if the pointer comes
// back null.
//
// $Id: findMaterialOrThrow.hh,v 1.3 2011/05/18 02:27:17 wb Exp $
// $Author: wb $
// $Date: 2011/05/18 02:27:17 $
//
// Original author Rob Kutschke
//
//

class G4Material;
class G4String;

namespace mu2e {

  G4Material* findMaterialOrThrow( G4String const& name);

} // end namespace mu2e
#endif /* Mu2eG4_findMaterialOrThrow_hh */
