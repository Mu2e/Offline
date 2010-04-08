//
// Decide which physics list to use.
//
// $Id: physicsListDecider.cc,v 1.2 2010/04/08 21:25:42 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/04/08 21:25:42 $
//
// Original author Rob Kutschke 
//
// Notes:
// 1) Extract the name of the requested physics list from the
//    config file, instantiate the requested physics list and
//    return a bare pointer to it.
//
// 2) The caller receives the pointer and immediately passes it 
//    to G4, which takes ownership of the physics list object.
//    The G4 interface requires a bare pointer.
//
// 3) There are two special names:
//     Minimal - the original Mu2e minimal physics list
//     N02     - the physics list copied from the G4 novice example N02.
// 
// 4) All other names are presumed to be valid names for physics lists that
//    can be created by the PhysListFactory.  At this writing ( April 2010),
//    the PhysListFactory has a default if it does not recognize name as
//    one of its known lists.  The default is QGSP_BERT 3.3.
//

// C++ includes
#include <string>

// Framework incldues
#include "FWCore/Utilities/interface/Exception.h"

// Mu2e includes
#include "Mu2eG4/inc/physicsListDecider.hh"
#include "Mu2eG4/inc/MinimalPhysicsList.hh"
#include "Mu2eG4/inc/PhysicsList.hh"
#include "Mu2eUtilities/inc/SimpleConfig.hh"

// G4 includes
#include "G4PhysListFactory.hh"
#include "G4VUserPhysicsList.hh"

using namespace std;

namespace mu2e{

  G4VUserPhysicsList* physicsListDecider ( const SimpleConfig& config ){

    G4VUserPhysicsList* physicsList(0);

    string name = config.getString("g4.physicsListName","N02");

    // Two special cases
    if ( name  == "Minimal" ) {
      physicsList = dynamic_cast<G4VUserPhysicsList*>(new MinimalPhysicsList );
    }

    else if ( name == "N02" ){
      physicsList = dynamic_cast<G4VUserPhysicsList*>(new PhysicsList(config) );
    }

    // General case
    else {
      G4PhysListFactory physListFactory;
      physicsList = physListFactory.GetReferencePhysList(name);
    }

    if ( !physicsList ){
      throw cms::Exception("G4CONTROL")
        << "Unable to load physics list named: "
        << name
        << "\n";
    }

    return physicsList;

  }

} // end namespace mu2e
  
