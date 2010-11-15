//
// G4Helper plugin.
//
// $Id: G4Helper_plugin.cc,v 1.1 2010/11/15 23:20:06 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2010/11/15 23:20:06 $
//
// Original author Rob Kutschke
//

// Framework includes
#include "FWCore/ServiceRegistry/interface/ServiceMaker.h"

// Mu2e includes
#include "Mu2eG4/inc/G4Helper.hh"

using namespace std;

namespace mu2e {

  G4Helper::G4Helper(edm::ParameterSet const& iPS, 
                     edm::ActivityRegistry&iRegistry){
  }

  G4Helper::~G4Helper(){
  }


} // end namespace mu2e

using mu2e::G4Helper;
DEFINE_FWK_SERVICE(G4Helper);
