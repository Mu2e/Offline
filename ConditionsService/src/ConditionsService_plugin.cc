//
// Plugin for the Conditions Service.
//
// $Id: ConditionsService_plugin.cc,v 1.1 2009/11/12 00:51:08 kutschke Exp $
// $Author: kutschke $
// $Date: 2009/11/12 00:51:08 $
//
// Original author Rob Kutschke
//

#include "FWCore/ServiceRegistry/interface/ServiceMaker.h"
#include "ConditionsService/inc/ConditionsService.hh"

using mu2e::ConditionsService;
DEFINE_FWK_SERVICE(ConditionsService);
