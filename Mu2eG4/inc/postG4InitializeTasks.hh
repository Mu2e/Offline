#ifndef Mu2eG4_postG4InitializeTasks_hh
#define Mu2eG4_postG4InitializeTasks_hh
//
// Steering routine to G4 call customization routines that must be called
// only after the call to Mu2eG4RunManager::Initialize.
//
// $Id: postG4InitializeTasks.hh,v 1.1 2012/06/04 19:28:01 kutschke Exp $
// $Author: kutschke $
// $Date: 2012/06/04 19:28:01 $
//

namespace mu2e{

  // This only needs to be templated to share
  // code for fhicl::ParameterSet and SimpleConfig cases.
  template<class Config> void postG4InitializeTasks(const Config& config);

}  // end namespace mu2e

#endif /* Mu2eG4_postG4InitializeTasks_hh */


