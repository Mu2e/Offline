#ifndef Mu2eG4_physicsListDecider_hh
#define Mu2eG4_physicsListDecider_hh
//
// Decide which physics list to use.
//
// $Id: physicsListDecider.hh,v 1.2 2011/05/17 15:41:36 greenc Exp $
// $Author: greenc $
// $Date: 2011/05/17 15:41:36 $
//
// Original author Rob Kutschke
//

// Forward declarations
class G4VUserPhysicsList;

namespace mu2e{

  // This only needs to be templated to share
  // code for fhicl::ParameterSet and SimpleConfig cases.

  // The returned pointer should be passed to G4, which will
  // take ownership of it.
  template<class Config> G4VUserPhysicsList* physicsListDecider (const Config& config);


}  // end namespace mu2e
#endif /* Mu2eG4_physicsListDecider_hh */
