//
// Decide which physics list to use.
//
// $Id: physicsListDecider.cc,v 1.11 2012/03/14 21:46:49 genser Exp $
// $Author: genser $
// $Date: 2012/03/14 21:46:49 $
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

// Framework includes
#include "cetlib/exception.h"

// Mu2e includes
#include "Mu2eG4/inc/physicsListDecider.hh"
#include "Mu2eG4/inc/MinimalPhysicsList.hh"
#include "Mu2eG4/inc/QGSP_BERT_MU2E00.hh"
#include "Mu2eG4/inc/PhysicsList.hh"
#include "Mu2eG4/inc/StepLimiterPhysConstructor.hh"
#include "Mu2eUtilities/inc/SimpleConfig.hh"

//tmp arrangement
#include "Mu2eG4/inc/QGSP_BERT_HP_MU2E00.hh"

// G4 includes
#include "G4PhysListFactory.hh"
#include "G4VUserPhysicsList.hh"
#include "QGSP.hh"

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

    else if ( name == "QGSP" ){
      G4VModularPhysicsList* tmp = new QGSP();
      tmp->RegisterPhysics( new StepLimiterPhysConstructor() );
      physicsList = tmp;
    }

    else if ( name == "QGSP_BERT_MU2E00" ){
      G4VModularPhysicsList* tmp = new QGSP_BERT_MU2E00();
      tmp->RegisterPhysics( new StepLimiterPhysConstructor() );
      physicsList = tmp;
    }

    else if ( name == "QGSP_BERT_HP_MU2E00" ){
      G4VModularPhysicsList* tmp = new QGSP_BERT_HP_MU2E00();
      tmp->RegisterPhysics( new StepLimiterPhysConstructor() );
      physicsList = tmp;
    }

    // General case
    else {
      G4PhysListFactory physListFactory;
      G4VModularPhysicsList* tmp = physListFactory.GetReferencePhysList(name);

      // The modular physics list takes ownership of the StepLimiterPhysConstructor.
      tmp->RegisterPhysics( new StepLimiterPhysConstructor() );

      physicsList = tmp;
    }

    if ( !physicsList ){
      throw cet::exception("G4CONTROL")
        << "Unable to load physics list named: "
        << name
        << "\n";
    }

    return physicsList;

  }

} // end namespace mu2e

